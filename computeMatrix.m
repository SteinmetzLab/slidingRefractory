
function [confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params)

cont = 0.5:0.5:35; 
rpEdges = (0:10*30)/30000; % 10 ms at 30kHz

% compute firing rate and spike count
spikeCount = numel(spikeTimes); 
[nACG,~] = histdiff(spikeTimes, spikeTimes, [0 1 2]);
firingRate = nACG(1,2)/spikeCount;


[nACG,rp] = histdiff(spikeTimes, spikeTimes, rpEdges);

confMatrix = nan(numel(cont), numel(rp));
for rpIdx = 1:numel(rp)
    
    % compute observed violations
    obsViol = sum(nACG(1,1:rpIdx)); % to check: rpIdx or rpIdx-1?? 
    
    for cidx = 1:numel(cont)
        
        confMatrix(cidx, rpIdx) = 100*computeViol(...
            obsViol, firingRate, spikeCount, rp(rpIdx), cont(cidx)/100);
    end
end
