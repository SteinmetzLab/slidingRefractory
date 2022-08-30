
function [rpMetrics, cont, rp] = slidingRP_all(spikeTimes, spikeClusters, params)
% compute the metric for each cluster in a recording

if nargin<3
    params = struct();
end

if ~isempty(params) && isfield(params, 'returnMatrix')
    returnMatrix = params.returnMatrix; 
else
    returnMatrix = false;
end

if ~isempty(params) && isfield(params, 'verbose')
    verbose = params.verbose; 
else
    verbose = true;
end

cids = unique(spikeClusters); 

rpMetrics = struct();

if verbose
    fprintf(1, 'Computing metrics for %d clusters\n', numel(cids)); 
end
for cidx = 1:numel(cids)
    st = spikeTimes(spikeClusters==cids(cidx)); 
    
    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
        nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate] ...
        = slidingRP(st, params);

    rpMetrics(cidx).cid = cids(cidx); 
    rpMetrics(cidx).maxConfidenceAt10Cont = maxConfidenceAt10Cont;
    rpMetrics(cidx).minContWith90Confidence = minContWith90Confidence;
    rpMetrics(cidx).timeOfLowestCont = timeOfLowestCont;
    rpMetrics(cidx).nSpikesBelow2 = nSpikesBelow2;
    if returnMatrix
        rpMetrics(cidx).confMatrix = confMatrix;
    end
    if verbose
        if minContWith90Confidence<=10; pfstring = 'PASS'; else pfstring = 'FAIL'; end;
        fprintf(1, '  %d: contamination = %.1f%%, %s max conf = %.2f%%, time = %.2f ms, n below 2 ms = %d\n', ...
            cids(cidx), pfstring, minContWith90Confidence, maxConfidenceAt10Cont, ...
             timeOfLowestCont*1000, nSpikesBelow2);
    end
end