

% script for testing sliding refractory period code

%% paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory')))

%% load

allst = readNPY('D:\Hopkins\spike_times.npy'); 
clu = readNPY('D:\Hopkins\spike_clusters.npy'); 

allst = double(allst)/30000;
cids = unique(clu);

%%
cidx = 0;
st = allst(clu==cidx);

[confMatrix, cont, rp] = computeMatrix(spikeTimes, []);
figure; imagesc(rp, cont, confMatrix)

[maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont] ...
    = slidingRP(spikeTimes, [])

%% test plotting

cidx =26;
st = allst(clu==cids(cidx));
params = struct();
params.cidx = cids(cidx);

plotSlidingRP(st, params)

%% test running it on all
params = [];
params.verbose = true; 
[rpMetrics, cont, rp] = slidingRP_all(allst, clu, params);