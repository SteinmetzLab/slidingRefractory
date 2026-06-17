
function rpMetrics = slidingRP_all(spikeTimes, spikeClusters, params)
% SLIDINGRP_ALL  Compute the Sliding Refractory Period quality metric for
%   every cluster in a recording.
%
%   rpMetrics = slidingRP_all(spikeTimes, spikeClusters)
%   rpMetrics = slidingRP_all(spikeTimes, spikeClusters, params)
%
%   For each unique cluster ID in spikeClusters, this function extracts the
%   corresponding spike times and passes them to slidingRP(), which
%   estimates contamination without assuming a fixed refractory period
%   duration. See Roth, Chapuis, Winter et al. (2026) for details.
%
%   INPUTS
%     spikeTimes    - [N x 1] vector of spike times in seconds for all
%                     clusters, where N is the total number of spikes.
%     spikeClusters - [N x 1] vector of cluster IDs, same length as
%                     spikeTimes. Each element identifies which cluster
%                     the corresponding spike belongs to.
%     params        - (optional) scalar struct of parameters. Any fields
%                     not specified will use their default values. Fields:
%
%       .contaminationThresh - Maximum acceptable contamination proportion.
%                              Passed to slidingRP(). Default: 0.1 (10%).
%       .confidenceThresh    - Minimum confidence required to accept a
%                              unit. Passed to slidingRP(). Default: 90.
%       .tauMin              - Minimum RP duration to consider (seconds).
%                              Passed to slidingRP(). Default: 0.0005.
%       .returnMatrix        - Logical. If true, the full confidence matrix
%                              from slidingRP() is stored in each element
%                              of rpMetrics. Default: false.
%       .verbose             - Logical. If true, print progress to the
%                              command window. Default: true.
%       .useParallel         - Logical. If true, use parfor when a parallel
%                              pool is already running. If no pool exists,
%                              the loop runs serially regardless. Set to
%                              false to force serial execution even when a
%                              pool is available. Default: true.
%
%   OUTPUTS
%     rpMetrics - [1 x nClusters] struct array with fields:
%       .cid              - Cluster ID.
%       .confidence       - Maximum confidence (%) that contamination is
%                           below the threshold, across all tested RP
%                           durations.
%       .contamination    - Minimum contamination level (%) that can be
%                           confirmed at the specified confidence level.
%       .timeOfLowestCont - RP duration (seconds) at which the minimum
%                           confirmed contamination was achieved.
%       .confMatrix       - Full confidence matrix (only populated when
%                           params.returnMatrix is true; empty otherwise).
%
%   EXAMPLES
%     % Basic usage with default parameters
%     rpMetrics = slidingRP_all(spikeTimes, spikeClusters);
%
%     % Custom thresholds, serial execution, with confidence matrices
%     params.contaminationThresh = 0.05;
%     params.confidenceThresh    = 95;
%     params.returnMatrix        = true;
%     params.useParallel         = false;
%     rpMetrics = slidingRP_all(spikeTimes, spikeClusters, params);

% ---- Parse optional parameters with defaults ----------------------------
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

if ~isfield(params, 'useParallel')
    params.useParallel = true;  % default: allow parallel
end

% ---- Determine parallel execution mode ----------------------------------
%   parfor with nWorkers=0 runs as a plain for loop, so we only need one
%   copy of the loop body.
if params.useParallel
    p = gcp('nocreate');  % check for existing pool, don't start one
    if isempty(p)
        nWorkers = 0;     % no pool running → run serial
    else
        nWorkers = p.NumWorkers;
    end
else
    nWorkers = 0;
end

% ---- Identify clusters and preallocate output ----------------------------
cids = unique(spikeClusters); 

rpMetrics(1:numel(cids)) = struct('cid',[], 'confidence',[], ...
    'contamination',[], 'timeOfLowestCont',[], 'confMatrix',[]);

if verbose
    fprintf(1, 'Computing metrics for %d clusters\n', numel(cids)); 
end

% ---- Pre-slice spike times by cluster ------------------------------------
%   Indexing into the full spikeTimes/spikeClusters arrays inside a parfor
%   would broadcast the entire arrays to every worker. Pre-slicing into a
%   cell array makes each element a proper sliced variable.
stCell = cell(numel(cids), 1);
for cidx = 1:numel(cids)
    stCell{cidx} = spikeTimes(spikeClusters == cids(cidx));
end

% ---- Main loop -----------------------------------------------------------
parfor (cidx = 1:numel(cids), nWorkers)
    st = stCell{cidx};   
    
    [passTest, confidence, contamination, timeOfLowestCont,...
        nViolShort, confMatrix, cont, rp, nACG] ...
        = slidingRP(st, params);
    
    s = struct();
    s.cid = cids(cidx);
    s.confidence = confidence;
    s.contamination = contamination;
    s.timeOfLowestCont = timeOfLowestCont;

    s.confMatrix = [];
    if returnMatrix
        s.confMatrix = confMatrix;
    end
    
    rpMetrics(cidx) = s;   
    
    if verbose
        if passTest; pfstring = 'PASS'; else pfstring = 'FAIL'; end;
        fprintf(1, '  %d: %s - contamination = %.1f%%, max conf = %.2f%%, time = %.2f ms\n', ...
            cids(cidx), pfstring, contamination, confidence, ...
             timeOfLowestCont*1000);
    end
end