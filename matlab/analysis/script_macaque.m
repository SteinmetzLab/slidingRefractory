

% script to analyze Horwitz macaque data

%% paths

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))

%% load
clear all; close all; 

% cd D:\Dropbox\code\slidingRefractory\python\slidingRP\Horowitz_data
% cd D:\Horwitz\LGN_Kilosort_Output\2024-04-24_Osiris\sorter_output
% cd D:\Horwitz\LGN_Kilosort_Output\2024-05-07_Osiris\sorter_output
% cd D:\Horwitz\LGN_Kilosort_Output\2024-05-14_Osiris\sorter_output
% cd('D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R01\kilosort3')
% cd('D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R02\kilosort3')
% cd('D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R03\kilosort3')
% cd('D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R04\kilosort3')
% cd D:\Horwitz\V1_NPX_ForNick\Q101525
% cd D:\Horwitz\V1_NPX_ForNick\Q122625
% cd D:\Horwitz\V1_NPX_ForNick\W040925
cd D:\Horwitz\V1_NPX_ForNick\W123024


clu = readNPY('spike_clusters.npy');
st = readNPY('spike_times.npy');

Fs = 30000; 

st = double(st)/Fs;
ucid = unique(clu);
rpBinSize = 2/Fs;
rpEdges = 0:rpBinSize:10/1000;


%% test ACGs

 


figure; pidx = 1; nstart = 370; 
for u = nstart:nstart+8
    subplot(3,3,pidx);

    theseST = st(clu==ucid(u));
    [nACG,rp] = histdiff(theseST, theseST, rpEdges);

    [xs, ys] = stairs(rp*1000, nACG);
    fill([xs; xs(end); xs(1)], [ys; 0; 0], 'k', 'EdgeColor', 'none');

    title(ucid(u))
    h = addX(2); h.Color = 'r';

    pidx = pidx+1;
end

%% image of them
inclACG = false(numel(ucid),1);
for u = 1:numel(ucid)

    theseST = st(clu==ucid(u));
    [nACG,rp] = histdiff(theseST, theseST, rpEdges);

    if u==1
        allACG = zeros(numel(ucid), numel(nACG)); 
        normACG = zeros(numel(ucid), numel(nACG)); 
    end

    inclACG(u) = sum(nACG(rp<0.00075),2)==0 & sum(nACG)>100;

    allACG(u,:) = nACG; 
    normACG(u,:) = nACG./prctile(nACG, 90);

end



figure; 
% imagesc(rp*1000, [], allACG(inclACG,:))
imagesc(rp*1000, [], normACG(inclACG,:))
h = addX(2); h.Color = 'r';

%% test the refractory period fitting

bin_size = rpBinSize;
bin_centers = rp;

inclIdx = find(inclACG); 
for nn = 1:numel(inclIdx)
    thisCID = inclIdx(nn); 
    acg = allACG(thisCID,:);


    [RP_estimated, sigmoid_params, fit_quality] = estimate_refractory_period(acg, bin_centers, bin_size);
    
    if fit_quality.success
        hold on; 
        plot(RP_estimated,nn, 'rx'); 
    end
end
%%

% Plot results
figure;
subplot(2,1,1);
plot(bin_centers*1000, acg, 'b.-', 'DisplayName', 'Original ACG');
hold on;

% Plot filtered ACG
filter_window_bins = round(0.83e-3 / bin_size);
if mod(filter_window_bins, 2) == 0; filter_window_bins = filter_window_bins + 1; end
acg_filtered = medfilt1(acg, filter_window_bins);
plot(bin_centers*1000, acg_filtered, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered ACG');

xlabel('Time (ms)');
ylabel('Firing rate');
title('Single-Sided Autocorrelogram');
legend('Location', 'best');
grid on;
xlim([0 10]);

subplot(2,1,2);
plot(bin_centers*1000, acg_filtered, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered ACG');
hold on;

if fit_quality.success
    params = sigmoid_params;
    RP_est = RP_estimated;
    sigmoid_func = @(x) params(1) + (params(2) - params(1)) ./ ...
        (1 + exp(-params(3) * (x - params(4))));
    x_plot = linspace(0, max(bin_centers), 500);
    plot(x_plot*1000, sigmoid_func(x_plot), 'g-', 'LineWidth', 2, 'DisplayName', 'Sigmoid fit');
    
    % Mark estimated RP
    RP_value = sigmoid_func(RP_est/1000);
    plot(RP_est, RP_value, 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'y', 'DisplayName', sprintf('RP_{est} = %.2f ms', RP_est));
    
    % Mark true RP for comparison
    % plot(true_RP*1000, sigmoid_func(true_RP), 'ro', 'MarkerSize', 10, ...
    %     'MarkerFaceColor', 'r', 'DisplayName', sprintf('RP_{true} = %.2f ms', true_RP*1000));
end

xlabel('Time (ms)');
ylabel('Firing rate');
title('Sigmoid Fit and RP Estimation');
legend('Location', 'southeast');
grid on;
xlim([0 5]); % Zoom in on recovery region


%% ok, now iterate over all macaque datasets to compute RPs

d=1;
% dirs{d} = 'D:\Dropbox\code\slidingRefractory\python\slidingRP\Horowitz_data'; names{d} = 'HorwitzV1'; d=d+1;
dirs{d} = 'D:\Horwitz\LGN_Kilosort_Output\2024-04-24_Osiris\sorter_output'; names{d} = 'HorwitzLGN_1'; d=d+1;
dirs{d} = 'D:\Horwitz\LGN_Kilosort_Output\2024-05-07_Osiris\sorter_output'; names{d} = 'HorwitzLGN_2'; d=d+1;
dirs{d} = 'D:\Horwitz\LGN_Kilosort_Output\2024-05-14_Osiris\sorter_output'; names{d} = 'HorwitzLGN_3'; d=d+1;
% dirs{d} = 'D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R01\kilosort3'; names{d} = 'Shude_1'; d=d+1;
% dirs{d} = 'D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R02\kilosort3'; names{d} = 'Shude_2'; d=d+1;
% dirs{d} = 'D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R03\kilosort3'; names{d} = 'Shude_3'; d=d+1;
% dirs{d} = 'D:\Dropbox\papers\2026_SlidingRP\Shude_data\Shude_Macaque_V1_ks3 Shude Zhu\R04\kilosort3'; names{d} = 'Shude_4'; d=d+1;
dirs{d} = 'D:\Horwitz\V1_NPX_ForNick\Q101525'; names{d} = 'HorwitzV1_Q10'; d=d+1;
dirs{d} = 'D:\Horwitz\V1_NPX_ForNick\Q122625'; names{d} = 'HorwitzV1_Q12'; d=d+1;
dirs{d} = 'D:\Horwitz\V1_NPX_ForNick\W040925'; names{d} = 'HorwitzV1_W04'; d=d+1;
dirs{d} = 'D:\Horwitz\V1_NPX_ForNick\W123024'; names{d} = 'HorwitzV1_W12'; d=d+1;

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

% print(f, 'RPdists', '-dpdf', '-vector');


