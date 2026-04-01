
function [passTest, confidence, contamination, timeOfLowestCont, ...
    nViolShort, confMatrix, cont, rp, nACG] = slidingRP(spikeTimes, varargin)
% SLIDINGRP  Compute the Sliding Refractory Period quality metric for one cluster.
%
%   [passTest, confidence, contamination, timeOfLowestCont, nViolShort, ...
%    confMatrix, cont, rp, nACG] = slidingRP(spikeTimes)
%   [...] = slidingRP(spikeTimes, params)
%
%   Tests whether a spike train has contaminated refractory periods without
%   assuming a fixed RP duration. For each candidate RP duration tau_r,
%   the method computes the statistical confidence that contamination is
%   below a user-defined threshold. A unit passes if that confidence
%   exceeds the required level for at least one tau_r > tau_min.
%
%   See Roth, Chapuis, Winter et al. (2026) for full algorithm description.
%   Requires histdiff from cortex-lab/spikes on the MATLAB path.
%
%   INPUTS
%     spikeTimes - [N x 1] vector of spike times in seconds.
%     params     - (optional) struct with any of the following fields:
%
%       .contaminationThresh - Maximum acceptable contamination (%).
%                              Default: 10.
%       .confidenceThresh    - Minimum confidence (%) required to pass.
%                              Default: 90.
%       .recDur              - Recording duration (seconds). Defaults to
%                              max(spikeTimes) if not provided. Recommended
%                              to specify explicitly.
%       .rpReject            - Minimum tau_r considered for pass/fail (s).
%                              Bins below this are excluded to avoid noise
%                              artefacts. Default: 0.0005 (0.5 ms).
%       .cont                - Vector of contamination levels (%) to test.
%                              Passing a single value equal to
%                              contaminationThresh speeds up computation
%                              but disables contamination estimation.
%                              Default: 0.5:0.5:35.
%       .nViolShortThresh    - Threshold for nViolShort output (s).
%                              Default: 0.002 (2 ms).
%       Additional fields are passed through to computeMatrix().
%
%   OUTPUTS
%     passTest         - Logical. True if confidence >= confidenceThresh for
%                        any tau_r > rpReject at contaminationThresh.
%     confidence       - Maximum confidence (%) that contamination is below
%                        contaminationThresh, across all tested tau_r.
%     contamination    - Minimum contamination level (%) confirmable at
%                        confidenceThresh. NaN if no level within the tested
%                        range (0.5–35%) reaches the required confidence
%                        (implying contamination is likely > 35%).
%     timeOfLowestCont - tau_r (seconds) at which minimum contamination is
%                        achieved. This is an estimate of the unit's true RP
%                        duration. NaN when contamination is NaN.
%     nViolShort       - Number of ACG counts with ISI < nViolShortThresh.
%                        Useful for identifying low-firing-rate units that
%                        failed only because of insufficient statistical
%                        power (nViolShort == 0 means no observed violations).
%     confMatrix       - [nCont x nRP] confidence matrix. See computeMatrix.
%     cont             - Vector of contamination levels tested (%).
%     rp               - Vector of tau_r values tested (seconds).
%     nACG             - ACG counts at each tau_r (counts per bin).

if nargin > 1
    params = varargin{1};
else
    params = struct();
end

if isfield(params, 'contaminationThresh')
    contThresh = params.contaminationThresh;
else
    contThresh = 10;
end
if isfield(params, 'confidenceThresh')
    confThresh = params.confidenceThresh;
else
    confThresh = 90;
end
if isfield(params, 'nViolShortThresh')
    nViolShortThresh = params.nViolShortThresh;
else
    nViolShortThresh = 0.002;
end
if isfield(params, 'rpReject')
    rpReject = params.rpReject;
else
    rpReject = 0.0005;
end

[confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params);
% confMatrix is [nCont x nRP]

% Exclude tau_r bins below rpReject from pass/fail decisions
testTimes = rp > rpReject;

% Maximum confidence at the user-defined contamination threshold
confidence = max(confMatrix(find(cont >= contThresh, 1), testTimes));

% Minimum contamination that can be confirmed at confThresh
[ii, ~] = find(confMatrix(:, testTimes) > confThresh);
[minI, ~] = min(ii);
contamination = cont(minI);
if isempty(contamination); contamination = NaN; end

% tau_r at which minimum contamination is achieved (RP duration estimate)
[~, minRP] = max(confMatrix(minI, testTimes));
timeOfLowestCont = rp(minRP + find(testTimes, 1) - 1);
if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

% Count short-ISI spikes (legacy diagnostic for low-rate units)
nViolShort = sum(nACG(1:find(rp > nViolShortThresh, 1)));

passTest = confidence > confThresh;
