
function [confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params)

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
[nACG,~] = histdiff(spikeTimes, spikeTimes, [0 1 2]);
firingRate = nACG(1,2)/spikeCount;

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
