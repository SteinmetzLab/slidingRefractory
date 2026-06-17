
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
%       .correction          - Logical. If true, apply the family-wise
%                              multiple-comparisons correction across tau_r
%                              (exact Poisson first-passage; Fig. S3). Slow and
%                              over-conservative for short-RP units; intended for
%                              Fig. S3 only. Default: false (standard metric).
%       Additional fields are passed through to computeMatrix().
%
%   OUTPUTS
%     passTest         - Logical. True if confidence >= confidenceThresh for
%                        any tau_r > rpReject at contaminationThresh.
%     confidence       - Maximum confidence (%) that contamination is below
%                        contaminationThresh, across all tested tau_r.
%     contamination    - Minimum contamination level (%) confirmable at
%                        confidenceThresh, computed analytically (continuous,
%                        not grid-rounded). NaN if it exceeds 35% (implying
%                        contamination is likely > 35%).
%     timeOfLowestCont - tau_r (seconds) at which minimum contamination is
%                        achieved. This is an estimate of the unit's true RP
%                        duration. NaN when contamination is NaN.
%     nViolShort       - Number of ACG counts with ISI < nViolShortThresh.
%                        Useful for identifying low-firing-rate units that
%                        failed only because of insufficient statistical
%                        power (nViolShort == 0 means no observed violations).
%     confMatrix       - [nCont x nRP] confidence matrix, built only when this
%                        output is requested (nargout >= 6); [] otherwise. The
%                        scalar metrics above never require it. See computeMatrix.
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

useCorrection = isfield(params, 'correction') && params.correction;

if useCorrection
    % ---- FWER-corrected variant (Fig. S3) -------------------------------
    % Build the family-wise-corrected confidence matrix and derive the scalar
    % metrics from it (grid-based; there is no analytical shortcut for the
    % corrected confidence). This is the expensive, non-default path.
    [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params);
    testTimes = rp > rpReject;

    confidence = max(confMatrix(find(cont >= contThresh, 1), testTimes));

    [ii, ~] = find(confMatrix(:, testTimes) > confThresh);
    [minI, ~] = min(ii);
    contamination = cont(minI);
    if isempty(contamination); contamination = NaN; end

    [~, minRP] = max(confMatrix(minI, testTimes));
    timeOfLowestCont = rp(minRP + find(testTimes, 1) - 1);
    if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

    nViolShort = sum(nACG(1:find(rp > nViolShortThresh, 1)));
    passTest = confidence > confThresh;
    return;
end

% ---- Standard (analytical fast) path ------------------------------------
if isfield(params, 'recDur');     recDur     = params.recDur;     else; recDur     = max(spikeTimes); end
if isfield(params, 'acgBinSize'); acgBinSize = params.acgBinSize; else; acgBinSize = 1/30000;          end
if isfield(params, 'testWindow'); testWindow = params.testWindow; else; testWindow = 0.01;             end

spikeCount = numel(spikeTimes);
rpEdges    = 0:acgBinSize:testWindow;
[nACG, rp] = histdiff(spikeTimes, spikeTimes, rpEdges);
obsViol    = cumsum(nACG);          % cumulative observed violations up to each tau_r
refDur     = rp + acgBinSize/2;     % right bin edge = tested tau_r

% Exclude tau_r bins below rpReject from pass/fail decisions
testTimes = rp > rpReject;

% Maximum confidence at the user-defined contamination threshold. This is the
% same quantity computeMatrix would give at the contThresh row, but computed as
% a single vector over tau_r (no full contamination grid needed).
confAtThresh = 100 * computeViol(obsViol, [], spikeCount, refDur, contThresh/100, recDur);
confidence = max(confAtThresh(testTimes));

% Minimum contamination confirmable at confThresh, computed analytically
% (continuous; equivalent to the grid search to within grid resolution).
[contamination, timeOfLowestCont] = computeMinContamination(...
    obsViol, spikeCount, refDur, rp, recDur, confThresh, rpReject);

% Count short-ISI spikes (diagnostic for low-rate units)
nViolShort = sum(nACG(1:find(rp > nViolShortThresh, 1)));

passTest = confidence > confThresh;

% Full confidence matrix is only built when the caller requests it (outputs
% 6-7). The scalar metrics above never need it.
confMatrix = [];
cont = [];
if nargout >= 6
    [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params);
end
