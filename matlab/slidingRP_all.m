
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

if ~isempty(params) && isfield(params, '2msNoSpikesCondition')
    2msNoSpikesCondition = params.2msNoSpikesCondition;
    if ~isempty(params) && isfield(params, 'FRthresh')
        FRthresh = params.FRthresh;
    else
        FRthresh = 0.5; %default FR thresh 
    end   
else
    2msNoSpikesCondition = false;
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
    


    %Add returned value of metric (pass or fail)
    if params['2msNoSpikesCondition']
        %In this case, reject neurons below FRthresh and accept if no
        %spikes below 2 (otherwise follow the regular behavior of the
        %metric)
        
            if firingRate < FRthresh
                rpMetrics(cidx).value = 0;
            else
                if nSpikesBelow2 == 0
                    rpMetrics(cidx).value = 1;
                else
                    if minContWith90Confidence <= 10
                        rpMetrics(cidx).value = 1;
                    else
                         rpMetrics(cidx).value = 0;
                    end
                end
            end
    else 
        %regular behavior of the metric, disregarding whether neurons have
        %spikes below 2ms  
        if minContWith90Confidence <=10
            rpMetrics(cidx).value = 1;
        else
            rpMetrics(cidx).value = 0;
        end
    end

    
    
    
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