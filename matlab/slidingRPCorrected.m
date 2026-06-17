

function [passTest, confidence, contamination, timeOfLowestCont,...
    nACGBelow2, confMatrix, cont, rp, nACG] ...
    = slidingRPCorrected(spikeTimes, varargin)
% SLIDINGRPCORRECTED  Sliding RP metric with a multiple-comparisons correction.
%
%   [passTest, confidence, contamination, timeOfLowestCont, nACGBelow2, ...
%    confMatrix, cont, rp, nACG] = slidingRPCorrected(spikeTimes)
%   [...] = slidingRPCorrected(spikeTimes, params)
%
%   Identical in interface to slidingRP(), but the confidence matrix is
%   computed by computeMatrixCorrected() instead of computeMatrix(). That
%   function applies a family-wise-error-rate correction across the swept
%   tau_r values (an exact Poisson first-passage / Markov-chain calculation;
%   see Methods, "False acceptance rate in the Sliding RP metric", Fig. S3),
%   so the false-passing rate is held to (1 - confidenceThresh) when the true
%   RP duration matches the full tested window.
%
%   This variant is provided for the manuscript's Fig. S3 only. It is NOT the
%   default metric: the correction is computationally expensive and, because
%   it cannot be made exact without knowing each unit's true RP duration, it
%   over-penalises units whose true RP is shorter than the tested window. Use
%   slidingRP() for standard analyses.
%
%   INPUTS
%     spikeTimes - [N x 1] vector of spike times in seconds.
%     params     - (optional) struct, same fields as slidingRP():
%       .contaminationThresh - Max acceptable contamination (%). Default 10.
%       .confidenceThresh    - Min confidence (%) to pass. Default 90.
%       .recDur              - Recording duration (s). Defaults to
%                              max(spikeTimes); strongly recommended to set.
%       .cont                - Vector of contamination levels (%) to test.
%                              A single value (= contaminationThresh) speeds
%                              up the computation but then only 'confidence'
%                              is meaningful, not 'contamination'.
%       .rpReject            - Min tau_r considered for pass/fail (s).
%                              Default 0.0005 (0.5 ms).
%
%   OUTPUTS
%     passTest         - Logical. confidence > confidenceThresh.
%     confidence       - Max corrected confidence (%) at contaminationThresh.
%     contamination    - Min contamination (%) confirmable at confidenceThresh;
%                        NaN if none within the tested range.
%     timeOfLowestCont - tau_r (s) at which that minimum contamination occurs.
%     nACGBelow2       - ACG count with ISI < 2 ms (diagnostic: 0 flags a
%                        low-rate unit rejected only for lack of power).
%     confMatrix, cont, rp, nACG - see computeMatrixCorrected().
%
%   SEE ALSO
%     slidingRP, computeMatrixCorrected

if nargin>1
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

if isfield(params, 'rpReject')
    rpReject = params.rpReject;
else
    rpReject = 0.0005;
end

[confMatrix, cont, rp, nACG] = computeMatrixCorrected(spikeTimes, params);
% confMatrix is [nCont x nRP]

testTimes = rp>rpReject; % exclude tau_r below rpReject from the pass/fail decision

% Maximum confidence at the contamination threshold row, across tested tau_r
confidence = max(confMatrix(find(cont>=contThresh,1), testTimes));

% Minimum contamination level reaching confThresh anywhere in the tested tau_r
[ii,~] = find(confMatrix(:,testTimes)>confThresh);
[minI, ~] = min(ii);
contamination = cont(minI);
if isempty(contamination); contamination = NaN; end

% tau_r at which that minimum contamination is most confidently confirmed
% (an estimate of the unit's RP duration); offset maps back to the full rp axis
[~,minRP] = max(confMatrix(minI,testTimes));
timeOfLowestCont = rp(minRP+find(testTimes,1)-1);
if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

% Diagnostic: ACG count with ISI < 2 ms (0 => low-rate unit, no observed violations)
nACGBelow2 = sum(nACG(1:find(rp>0.002,1)));

passTest = confidence>confThresh;