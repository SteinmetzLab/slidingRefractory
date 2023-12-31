
function [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params)
% function [confMatrix, cont, rp, nACG] = computeMatrix(spikeTimes, params)
%
% Core function for Sliding RP metric. 
% 
% Inputs: 
% - spikeTimes, a vector of times in seconds
% - params, a struct which may have: 
%    - recDur, recording duration in seconds, if not specified will be taken
%   as the max spike time - highly recommended to specify
%   - cont, a vector of contamination levels at which to test the neuron.
%   The computation can be accelerated by setting this to only a single
%   value (which would need to be contaminationThresh), but this will not
%   allow a useful estimate of the contamination level, i.e. the returned
%   variable 'confidence' will be correct but 'contamination' will not. 
%
% Outputs: 
% - confMatrix - result of the test, dimensions [contamination levels, RP duration]
%    each entry reflects the confidence that the contamination is less than
%    that level when testing at that RP duration. 
% - cont - scalar or vector of the contamination levels tested
% - rp - vector of the RP durations tested
% - nACG - the ACG function of the neuron at the values in rp. Units are
%   counts


if nargin>1 && isfield(params, 'cont')
    cont = params.cont;
else
    cont = 0.5:0.5:35; % contamination levels at which to assess confidence
end
if nargin>1 && isfield(params, 'recDur')
    recDur = params.recDur;
else
    recDur = max(spikeTimes); 
end
% rpEdges = (0:10*30)/30000; % 10 ms at 30kHz
rpBinSize = 1/30000;
rpEdges = 0:rpBinSize:10/1000;

% compute firing rate and spike count
spikeCount = numel(spikeTimes); 
% [nACG,~] = histdiff(spikeTimes, spikeTimes, [0 1 2]);
% firingRate = nACG(1,2)/spikeCount;
firingRate = []; % note - firing rate no longer required with updated calculation

[nACG,rp] = histdiff(spikeTimes, spikeTimes, rpEdges);

confMatrix = nan(numel(cont), numel(rp));
for rpIdx = 1:numel(rp)
    
    % compute observed violations
    obsViol = sum(nACG(1,1:rpIdx));
    
    for cidx = 1:numel(cont)
        
        confMatrix(cidx, rpIdx) = 100*computeViol(...
            obsViol, firingRate, spikeCount, rp(rpIdx)+rpBinSize/2, cont(cidx)/100, recDur);
    end
end
