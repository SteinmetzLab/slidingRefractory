

function [passTest, confidence, contamination, timeOfLowestCont,...
    nACGBelow2, confMatrix, cont, rp, nACG] ...
    = slidingRP(spikeTimes, varargin)
% compute the metric for a single cluster in a recording
%
% Inputs:
% - spikeTimes, a vector of times in seconds
% - params, a struct which may contain: 
%   - contaminationThresh, the amount of contamination that will be
%   acceptable to pass the neuron, default 10 (%)
%   - confidenceThresh, the amount of confidence that will be required to
%   decide that a given contamination level is not breached, default 90 (%)
%   - recDur, recording duration in seconds, if not specified will be taken
%   as the max spike time - highly recommended to specify
%   - cont, a vector of contamination levels at which to test the neuron.
%   The computation can be accelerated by setting this to only a single
%   value (which would need to be contaminationThresh), but this will not
%   allow a useful estimate of the contamination level, i.e. the returned
%   variable 'confidence' will be correct but 'contamination' will not. 
%   - rpReject
%
%
% Outputs:
% - passTest, binary indicating whether the neuron's contamination level is
% below contaminationThresh with at least confidenceThresh level of
% confidence
% - confidence, the confidence that you have less than the threshold level
% of contamination
% - contamination, the minimum contamination for which you have at least
% the threshold level of confidence. A value of NaN indicates high
% contamination
% - timeOfLowestCont, the time at which best score happens
% - nACGBelow2, the number of spikes in the ACG below 2 ms, so you can
% check whether this is 0 for low firing rate neurons that have been
% rejected by the normal metric
% - for other outputs, see documentation of computeMatrix, as these are
% passed through from there. 

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

[confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params); 
% matrix is [nCont x nRP]

testTimes = rp>rpReject; % don't consider any times below this for analysis

confidence = max(confMatrix(find(cont>=contThresh,1), testTimes)); 

[ii,jj] = find(confMatrix(:,testTimes)>confThresh); 
[minI, idx] = min(ii); 
contamination = cont(minI);  
if isempty(contamination); contamination = NaN; end

% TODO comment 
[~,minRP] = max(confMatrix(minI,testTimes)); 
timeOfLowestCont = rp(minRP+find(testTimes,1)-1);
if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

nACGBelow2 = sum(nACG(1:find(rp>0.002,1)));

passTest = confidence>confThresh;