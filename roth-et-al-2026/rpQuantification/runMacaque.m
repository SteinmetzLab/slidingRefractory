

% script to analyze Horwitz macaque data

%% paths

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))

rootd = 'D:\Horwitz'; % where you put the data from https://doi.org/10.6084/m9.figshare.31587391.v1

%% iterate over all macaque datasets to compute RPs

names = {'LGN_1', 'LGN_2', 'LGN_3', 'V1_1', 'V1_2', 'V1_3', 'V1_4'};

for d = 1:numel(names)
    dirs{d} = fullfile(rootd, names{d});
end


Fs = 30000;

rpBinSize = 1/Fs;
rpEdges = 0:rpBinSize:10/1000;

f = figure; 
for d = 1:numel(dirs)

    fprintf(1, 'dataset %s\n', dirs{d}); 

    % load

    clu = readNPY(fullfile(dirs{d}, 'spike_clusters.npy'));
    st = readNPY(fullfile(dirs{d},'spike_times.npy'));    

    st = double(st)/Fs;
    ucid = unique(clu);

    params = struct(); params.cont = 10; params.recDur = max(st); % for slidingRP below

    inclACG = false(numel(ucid),1);
    for u = 1:numel(ucid)

        % compute ACG
        theseST = st(clu==ucid(u));
        [nACG,rp] = histdiff(theseST, theseST, rpEdges);

        if u==1
            allACG = zeros(numel(ucid), numel(nACG));
            normACG = zeros(numel(ucid), numel(nACG));
            allRP = zeros(numel(ucid),1);
            allFits = {};
            allSlidingRPpass = zeros(numel(ucid),1);

        end

        % estimate RP
        [RP_estimated, sigmoid_params, fit_quality] = estimate_refractory_period(nACG, rp, rpBinSize);

        % calculate possible inclusion determination

        inclACG(u) = sum(nACG(rp<0.00075),2)==0 & sum(nACG)>100;

        allACG(u,:) = nACG;
        normACG(u,:) = nACG./prctile(nACG, 90);
        allRP(u) = RP_estimated; 
        allFits{u} = {sigmoid_params, fit_quality};

        % also calculate the sliding RP metric        
        [passTest, confidence, contamination, timeOfLowestCont,...
            nACGBelow2, confMatrix, cont, rp, nACG] ...
            = slidingRP(theseST, params);

        allSlidingRPpass(u) = passTest;

    end

    % save it out

   save(names{d}, 'allACG', 'normACG', 'allRP', 'allFits', 'allSlidingRPpass');

   subplot(numel(dirs),2, (d-1)*2+1);
   hist(allRP, 0:0.2:6)
   title(names{d}); 

   subplot(numel(dirs),2, (d-1)*2+2);
   hist(allRP(inclACG), 0:0.2:6)

   drawnow;
end



