

% Fig 1 plots of comparison of estmiated RPs across datasets,
% using Noam's table of the fit values. 

%%

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'roth-et-al-2026')))
addpath(genpath(fullfile(githubDir, 'spikes'))) % github/cortex-lab/spikes


t = readtable(fullfile(githubDir, 'slidingRefractory', 'roth-et-al-2026','rpQuantification', ...
    'df_units__steinmetz_allen_ibl__Isocortex_TH_HPF__fr2.csv'));


%% plot

dsets = unique(t.dataset);
regs = unique(t.region);

plotD = 2; % plot IBL dataset

f = figure; f.Color = 'w';
colors = colororder();

subplot(1,2,1); 

for r = 1:3

theseRP = t.estimatedRP(strcmp(t.dataset, dsets{plotD})&strcmp(t.region, regs{r}));
theseRP = theseRP(theseRP>1); 

[n,x] = hist(theseRP, 0:0.2:5);
nNorm = n./sum(n);
H = stairs(x,nNorm);
H.Color = colors(r,:); H.LineWidth = 2;
legH{r} = H;

mns(r) = median(theseRP); 
hold on; 

end

for r = 1:3
    h = addX(mns(r)); h.Color = colors(r,:); h.LineWidth = 2;
end
box off;
xlabel('Estimated RP duration (ms)')
ylabel('Proportion of neurons');
title(sprintf('%s dataset', upper(dsets{plotD})));
leg = legend([legH{:}], regs); leg.Box = 'off';
axis square;


subplot(1,2,2); hold on;

for r = 1:numel(regs)
    for d = 1:numel(dsets)
        
        theseRP = t.estimatedRP(strcmp(t.dataset, dsets{d})&strcmp(t.region, regs{r}));
        theseRP = theseRP(theseRP>1); 
        
        fprintf(1, '%s, %s, n=%d\n', regs{r}, dsets{d}, numel(theseRP));
        mn = median(theseRP); 
        %stderr = std(theseRP)./sqrt(numel(theseRP));
        ci = bootci(10000, @median, theseRP);

        yplot = sub2ind([3 3], r, d);
        plot(mn, yplot, 'o', 'Color', colors(r,:)); 
        %plot(mn+stderr*[-1 1], yplot*[1 1],'LineWidth', 2.0, 'Color', colors(r,:));
        plot(ci, yplot*[1 1],'LineWidth', 2.0, 'Color', colors(r,:));
        labels{yplot} = sprintf('%s, %s', dsets{d}, regs{r});
    end
end
set(gca, 'YTick', 1:9, 'YTickLabel', labels);
box off; 
xlabel('Estimated RP duration (ms)')
xlim([1 4])
axis square;

%% stats

% [P,T,STATS,TERMS]=anovan(t.estimatedRP, {t.dataset t.region}, ...
%     'model', 'interaction', 'varnames', {'dataset', 'region'})

% definitely significant: 
%   Source           Sum Sq.   d.f.    Mean Sq.     F        Prob>F   
% --------------------------------------------------------------------
%   dataset            66.52       2    33.259     48.52   1.04043e-21
%   region            280.04       2   140.022    204.29    7.3425e-88
%   dataset*region     50.49       4    12.623     18.42   4.27811e-15
%   Error            7636.12   11141     0.685                        
%   Total            8504.21   11149                                  

[P,anovatab,STATS]=anova1(t.estimatedRP, t.region)

% Source     SS       df       MS        F         Prob>F   
% ----------------------------------------------------------
% Groups    732.85       2   366.425   525.59   7.41547e-219
% Error    7771.36   11147     0.697                        
% Total    8504.21   11149                       

%% make a similar plot for the macaque data

% Need to run "runMacaque.m" first, then make rootd here the same as what
% you used there
rootd = 'D:\Horwitz';

names = {'LGN_1', 'LGN_2', 'LGN_3', 'V1_1', 'V1_2', 'V1_3', 'V1_4'};

for d = 1:numel(names)
    dirs{d} = fullfile(rootd, names{d});
end

f = figure; f.Color = 'w';
clear nameLabels
for n = 1:numel(names)
    load(names{n})
    theseRP = allRP;
    
    inclUnits = false(size(allRP));
    for u = 1:size(allACG,1)
        if ~isnan(allRP(u))
            preRate = mean(allACG(u,rp*1000<allRP(u)));
            postRate = mean(allACG(u,rp*1000>allRP(u)));
            % inclUnits(u) =  postRate>5*preRate & allRP(u)>1;            
            % inclUnits(u) =  allRP(u)>1;  
            %inclUnits(u) =  allRP(u)>1 & allSlidingRPpass(u); 
            inclUnits(u) =  postRate>5*preRate & allRP(u)>1 & allSlidingRPpass(u);
        end
    end

    fprintf(1, '%s, n=%d\n', names{n}, sum(inclUnits));
    nameLabels{n} = sprintf('%s n=%d', names{n}, sum(inclUnits));
    
    theseRP = theseRP(inclUnits); 
    
    mn = nanmedian(theseRP); 
    %stderr = nanstd(theseRP)./sqrt(sum(~isnan(theseRP)));
    ci = bootci(10000, @median, theseRP);
    
    subplot(1,2,1);
    yplot = n;
    plot(mn, yplot, 'ko'); hold on;
    %plot(mn+stderr*[-1 1], yplot*[1 1], 'LineWidth', 2.0);
    plot(ci, yplot*[1 1], 'k', 'LineWidth', 2.0);
    
    subplot(1,2,2); 
    [n,x] = hist(theseRP, 0:0.2:5);
    nNorm = n./sum(n);
    H = stairs(x,nNorm*3+yplot); hold on;
    H.LineWidth = 2;
    
    drawnow;
end

subplot(1,2,1);
box off; 
xlim([1 4]); ylim([0 numel(names)+1]);
xlabel('Estimated RP duration (ms)')
set(gca, 'YTick', 1:numel(names), 'YTickLabel', nameLabels);
axis square;

subplot(1,2,2);
box off; 
ylim([0 numel(names)+1]);
set(gca, 'YTick', 1:numel(names), 'YTickLabel', nameLabels);
xlabel('Estimated RP duration (ms)')
axis square;


%% more detailed macaque plots

Fs = 30000;

rpBinSize = 1/Fs;
rpEdges = 0:rpBinSize:10/1000;
rp = rpEdges(1:end-1)+rpBinSize/2;

f = figure; f.Color = 'w';
f2 = figure; f2.Color = 'w';

theseIdx = 1:numel(names);
for n = theseIdx%1:numel(names)
    load(names{n});
    
    %inclUnits = sum(allACG(:,rp<0.00075),2)==0 & sum(allACG,2)>100;
    inclUnits = false(size(allRP));
    for u = 1:size(allACG,1)
        if ~isnan(allRP(u))
            preRate = mean(allACG(u,rp*1000<allRP(u)));
            postRate = mean(allACG(u,rp*1000>allRP(u)));
            inclUnits(u) =  postRate>5*preRate & allRP(u)>1;
        end
    end
    
    figure(f)
    subplot(1,numel(names),n); hold on;
    
    inclACG = normACG(inclUnits,:); 
    inclRP = allRP(inclUnits); 
    [sRP,ii] = sort(inclRP); 
    sortACG = inclACG(ii,:);
    imagesc(rp*1000, [], sortACG); hold on;
    plot(sRP, 1:numel(sRP), 'rx');
    h =addX(2); h.Color = 'g';h.LineWidth = 2;
    caxis([0 5]);
    xlim([0 6]); 
    title([names{n} ' ' num2str(sum(inclUnits))])
    

    figure(f2); 
    plot(sRP, (1:numel(sRP))/numel(sRP), 'o-', 'LineWidth', 2.0); hold on; 

    nU = sum(inclUnits); 
    if contains(names{n}, 'V1'); regAcr = 'Isocortex'; else; regAcr = 'TH'; end;
    newRows = table(ones(nU,1), ones(nU,1), repmat({names{n}},nU,1), repmat({regAcr}, nU,1), ...
        ones(nU,1), inclRP, repmat({'monkey'}, nU,1),...
        'VariableNames', {'number','index','dataset','region','mean_fr','estimatedRP','species'});
    t = [t; newRows];
end
legend(names(theseIdx))
colororder('gem12')
box off; 
xlabel('Estimated RP duration (ms)')
ylabel('Neurons'); 

%% mouse vs monkey stats

% need to run the cell at the top (loading table t of mouse data) and cell
% right above this one to put macaque data in it. 

% this one doesn't work presumably because not all regions are in all
% species? and not all regions are in all datasets?
% [P,T,STATS,TERMS]=anovan(t.estimatedRP, {t.dataset t.region, t.species}, ...
%     'model', 'interaction', 'varnames', {'dataset', 'region', 'species'})

tSub = t(strcmp(t.region, 'TH') | strcmp(t.region, 'Isocortex'), :);
[P,T,STATS,TERMS]=anovan(tSub.estimatedRP, { tSub.region, tSub.species}, ...
    'model', 'interaction', 'varnames', {'region', 'species'})

%   Source           Sum Sq.   d.f.    Mean Sq.     F        Prob>F   
% --------------------------------------------------------------------
%   region             27.68       1   27.6809     42.96   5.87394e-11
%   species            30          1   30.0015     46.56   9.41204e-12
%   region:species     81.81       1   81.8095    126.95   2.84329e-29
%   Error            6516.95   10113    0.6444                        
%   Total            6895.86   10116   

% this ignores dataset but finds a significant effect of region, of
% species, and interaction. 

%% also plot example ACGs 

names = {'LGN_1', 'LGN_2', 'LGN_3', 'V1_1', 'V1_2', 'V1_3', 'V1_4'};

for d = 1:numel(names)
    dirs{d} = fullfile(rootd, names{d});
end


rebin = 5;
startCID = 28;
theseIdx = 1:numel(names);
for n = theseIdx%1:numel(names)
    load(dirs{n});
    
    %inclUnits = sum(allACG(:,rp<0.00075),2)==0 & sum(allACG,2)>100;
    inclUnits = false(size(allRP));
    for u = 1:size(allACG,1)
        if ~isnan(allRP(u))
            preRate = mean(allACG(u,rp*1000<allRP(u)));
            postRate = mean(allACG(u,rp*1000>allRP(u)));
            inclUnits(u) =  postRate>5*preRate & allRP(u)>1;
        end
    end
    
    inclACG = normACG(inclUnits,:); 
    inclRP = allRP(inclUnits);
    
    f = figure; f.Color = 'w';
    
    for q = 1:9
        subplot(3,3,q);
        
        nACG = inclACG(startCID-1+q,:);
        rebinACG = mean(reshape(nACG, rebin, []));
        rebinRP = mean(reshape(rp, rebin, []));
        
        [xs, ys] = stairs(rebinRP*1000, rebinACG);
        fill([xs; xs(end); xs(1)], [ys; 0; 0], 'k', 'EdgeColor', 'none'); hold on;
        
        title(sprintf('%s, %d, %.2f ms', names{n}, startCID-1+q, inclRP(startCID-1+q)))
        xlim([0 6]); 
        for gridL = 0:6; h = addX(gridL); h.Color = 0.7*[1 1 1]; h.LineStyle = '-'; end;
        
        h = addX(inclRP(startCID-1+q)); h.Color = 'r';
        
        box off; axis off;
        
        pidx = pidx+1;
    end
    drawnow;
end

%% example ACGs for mouse data
% 

addpath(genpath(fullfile(githubDir, 'steinmetz-et-al-2019'))) % nsteinme/steinmetz-et-al-2019/
addpath(genpath(fullfile(githubDir, 'allenCCF'))) % github/cortex-lab/allenCCF
addpath(genpath(fullfile(githubDir, 'spikes')))% github/cortex-lab/spikes

clear all; close all;

st = loadStructureTree();

recName = 'Cori_2016-12-18';

% specify location of data from https://doi.org/10.6084/m9.figshare.9598406
rootD = fullfile('x:/dataLocation', recName);

s = loadSession(rootD);

%%

clusterCCF = s.channels.brainLocation.allen_ontology(s.clusters.peakChannel,:); 

% identify cortex, hippocampus, and thalamus clusters
clear isCTX isHPF isTH fr
for q = 1:size(clusterCCF,1)
    c = clusterCCF(q,:);
    spidx = find(c==' ',1); if ~isempty(spidx); c = c(1:spidx-1); end;
    idx = find(contains(st.acronym, c));
    isCTX(q) = contains(st.structure_id_path(idx(1)), '/315/');
    isHPF(q) = contains(st.structure_id_path(idx(1)), '/1089/');
    isTH(q) = contains(st.structure_id_path(idx(1)), '/549/');
    theseST = s.spikes.times(s.spikes.clusters==q);
    fr(q) = numel(theseST)./recDur;
end

fprintf(1, 'CTX: %d, HPF: %d, TH: %d\n', sum(isCTX&fr<2), sum(isHPF&fr<2), sum(isTH&fr<2));

%% plot example ACGs from each region
startCID = 10;
ctxID = find(isCTX); 
names = {'CTX', 'HPF', 'TH'};
ids = {find(isCTX&fr>2), find(isHPF&fr>2), find(isTH&fr>2)};

rebin = 5;
Fs = 30000;
rpBinSize = 1/Fs;
rpEdges = 0:rpBinSize:10/1000;
recDur = s.spikes.times(end);
theseIdx = 1:numel(names);
for n = theseIdx%1:numel(names)
    
    f = figure; f.Color = 'w'; f.Name = names{n};
    
    for q = 1:9
        
        if q+startCID-1<=numel(ids{n})
            cid = ids{n}(q+startCID-1);
            theseST = s.spikes.times(s.spikes.clusters==cid);
            
            % compute ACG and FR and RPdur
            [nACG,rp] = histdiff(theseST, theseST, rpEdges);            
            [RPest, sigmoid_params, fit_quality] = estimate_refractory_period(nACG, rp, rpBinSize, 0.1);
            
            subplot(3,3,q);
            
            rebinACG = mean(reshape(nACG, rebin, []));
            rebinRP = mean(reshape(rp, rebin, []));
            
            [xs, ys] = stairs(rebinRP*1000, rebinACG);
            fill([xs; xs(end); xs(1)], [ys; 0; 0], 'k', 'EdgeColor', 'none'); hold on;
            
            title(sprintf('%s, %d, %.2f ms, %.2f sp/s', names{n}, cid, RPest, fr(cid)))
            xlim([0 6]);
            for gridL = 0:6; h = addX(gridL); h.Color = 0.7*[1 1 1]; h.LineStyle = '-'; end;
            
            h = addX(RPest); h.Color = 'r';
            
            box off; axis off;
        end
    end
    drawnow;
end
