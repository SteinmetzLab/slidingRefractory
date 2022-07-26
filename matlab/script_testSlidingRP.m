

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

[confMatrix, cont, rp] = computeMatrix(st, []);
figure; imagesc(rp, cont, confMatrix)

[maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont] ...
    = slidingRP(st, [])

%% test plotting

% cidx =1;
cidx = find(cids==0);
st = allst(clu==cids(cidx));
params = struct();
params.cidx = cids(cidx);

plotSlidingRP(st, params)

%% test running it on all
params = [];
params.verbose = true; 
[rpMetrics, cont, rp] = slidingRP_all(allst, clu, params);

%% quick test of firing rate estimate from ACG versus total spike count
fr = zeros(numel(cids), 2); 
recDur = allst(end);

for cidx = 1:numel(cids)
    st = allst(clu==cids(cidx));
    spikeCount = numel(st); 
    recDur = st(end); 
    fr(cidx,1) = spikeCount/recDur; 
    [nACG,~] = histdiff(st, st, [0 1 2]);
    fr(cidx,2) = nACG(1,2)/spikeCount;
end

figure; plot(fr(:,1), fr(:,2),'.')
xlabel('from spikecount/recdur'); 
ylabel('from acg'); 
addXeqY()
% axis square