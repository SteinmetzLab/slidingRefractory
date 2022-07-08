
function [rpMetrics, cont, rp] = slidingRP_all(spikeTimes, spikeClusters, params)
% compute the metric for each cluster in a recording

if ~isempty(params) && isfield(params, 'returnMatrix')
    returnMatrix = params.returnMatrix; 
else
    returnMatrix = false;
end

if ~isempty(params) && isfield(params, 'verbose')
    verbose = params.verbose; 
else
    verbose = false;
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
        if maxConfidenceAt10Cont>=90; pfstring = 'PASS'; else pfstring = 'FAIL'; end;
        fprintf(1, '  %d: %s max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms, n below 2 ms = %d\n', ...
            cids(cidx), pfstring, maxConfidenceAt10Cont, ...
            minContWith90Confidence, timeOfLowestCont*1000, nSpikesBelow2);
    end
end