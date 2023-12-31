
% notes: 
% - can put nsim loop outside to save on generating spike trains
% - can build a big table with all parameter values - this would allow
% quickly re-plotting

%% paths
addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'nickBox')))

figSaveDir = fullfile(githubDir, 'slidingRefractory', 'matlab', 'simulations');

%% Fig 2 example neuron test



%% Fig 3a: performance across firing rates

saveFig = true;
nSim = 1000; % number of simulations
RPdur = 0.0025; % true RP duration, s
recDur = 3600*2; % recording duration, s
contProp = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20]/100; % simulated proportion contamination
baseRates = [0.5 1 2 5 10 20]; % rate of the true neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail
confThresh = 90; % confidence we need to accept a neuron

params = struct(); params.cont = contThresh;
params.contaminationThresh = contThresh;
params.confidenceThresh = confThresh;

f = figure; f.Color = 'w';
colors = myCopper(0.6, numel(baseRates)+1);
colors = colors(2:end,:);

passPct = zeros(numel(baseRates), numel(contProp));
passErr = zeros(numel(baseRates), numel(contProp),2);

for bidx = 1:numel(baseRates)
    totalRate = baseRates(bidx)
    for c = 1:numel(contProp)
        
        % this calculation ensures that the contaminating spikes generated
        % at this rate do form the correct proportion of the total
        % contRate = contProp(c)*baseRate/(1-contProp(c));
        baseRate = (1-contProp(c))*totalRate;
        contRate = contProp(c)*totalRate;
        
        params.recDur = recDur;
        
        simRes = zeros(nSim, 1);
        for n = 1:nSim
            
            st = genST(baseRate, recDur, RPdur); % true spikes
            contST = genST(contRate,recDur, 0); % contaminating spikes
            combST = sort([st; contST]); % combined spike train
 
%             [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
%             
%             simRes(n) = max(confMatrix(rpTestVals>0.0005))>conf;
            
            [passTest, confidence, contamination, timeOfLowestCont,...
                nSpikesBelow2, confMatrix, cont, rp, nACG] ...
                = slidingRP(combST, params);
            
            simRes(n) = passTest;

        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        [~, pci] = binofit(sum(simRes),nSim,0.05); % 95% confidence interval
        passErr(bidx,c,:) = pci*100;
        
    end
    % legH(bidx) = plot(contProp*100, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
    legH(bidx) = plotWithErrUL(contProp*100, passPct(bidx,:), squeeze(passErr(bidx,:,:)), colors(bidx,:)); 
    legH(bidx).Marker = 'o'; legH(bidx).MarkerFaceColor = colors(bidx,:);
        
    hold on; drawnow;
end

xlabel('Contamination (%)'); 
ylabel('Percent pass (%)');
addX(10); 
leg = legend(legH, array2stringCell(baseRates)); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off

if saveFig
    print(f, fullfile(figSaveDir, 'Fig3a.pdf'), '-dpdf');
end

%% Fig 3c: performance by recording duration
% TODO: Update the base rate calculation method

nSim = 1000; % number of simulations
RPdur = 0.0025; % true RP duration, s
recDurs = [0.5 1 2 3]*3600; % recording duration, s
contProp = [0:22]/100; % simulated proportion contamination
baseRate = 5; % rate of the true neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail
conf = 90; % confidence we need to accept a neuron

params = struct(); params.cont = contThresh;

testSliding = true; % alternative is to test at a single cutoff value

f = figure; f.Color = 'w';
colors = myCopper(0.3, numel(recDurs)+1);
colors = colors(2:end,:);

passPct = zeros(numel(recDurs), numel(contProp));
passErr = zeros(numel(recDurs), numel(contProp),2);

for bidx = 1:numel(recDurs)
    recDur = recDurs(bidx)
    for c = 1:numel(contProp)
        
        % this calculation ensures that the contaminating spikes generated
        % at this rate do form the correct proportion of the total
        contRate = contProp(c)*baseRate/(1-contProp(c));
        
        simRes = zeros(nSim, 1);
        for n = 1:nSim
            
            st = genST(baseRate, recDur, RPdur); 
            contST = genST(contRate,recDur, 0);
            combST = sort([st; contST]);
 
            [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
            
            if testSliding
                % normal sliding rp test result
                simRes(n) = max(confMatrix(rpTestVals>0.0005))>conf;
            else
                % testing only at a single rp duration
                simRes(n) = confMatrix(find(rpTestVals>0.002,1))>conf;
            end
        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        [~, pci] = binofit(sum(simRes),nSim,0.05); % 95% confidence interval
        passErr(bidx,c,:) = pci*100;
        
    end
    % legH(bidx) = plot(contProp*100, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
    legH(bidx) = plotWithErrUL(contProp*100, passPct(bidx,:), squeeze(passErr(bidx,:,:)), colors(bidx,:)); 
    legH(bidx).Marker = 'o'; legH(bidx).MarkerFaceColor = colors(bidx,:);
        
    hold on; drawnow;
end

xlabel('Contamination (%)'); 
ylabel('Percent pass (%)');
addX(10); 
leg = legend(legH, array2stringCell(recDurs/3600)); 
leg.Title.String = 'Recording Duration (h)';
legend boxoff
box off

print(f, fullfile(figSaveDir, 'Fig3c.pdf'), '-dpdf');

%% Fig 3d: performance by RP duration
% TODO: update base rate calculation

nSim = 1000; % number of simulations
RPdurs = [1.5 2 3 4 5 6]/1000; % true RP duration, s
recDur = 2*3600; % recording duration, s
contProp = [0:22]/100; % simulated proportion contamination
baseRate = 5; % rate of the true neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail
conf = 90; % confidence we need to accept a neuron

params = struct(); params.cont = contThresh;

testSliding = true; % alternative is to test at a single cutoff value

f = figure; f.Color = 'w';
colors = myCopper(0.8, numel(RPdurs)+1);
colors = colors(2:end,:);

passPct = zeros(numel(RPdurs), numel(contProp));
passErr = zeros(numel(RPdurs), numel(contProp),2);

for bidx = 1:numel(RPdurs)
    RPdur = RPdurs(bidx)
    for c = 1:numel(contProp)
        
        % this calculation ensures that the contaminating spikes generated
        % at this rate do form the correct proportion of the total
        contRate = contProp(c)*baseRate/(1-contProp(c));
        
        simRes = zeros(nSim, 1);
        for n = 1:nSim
            
            st = genST(baseRate, recDur, RPdur); 
            contST = genST(contRate,recDur, 0);
            combST = sort([st; contST]);
 
            [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
            
            if testSliding
                % normal sliding rp test result
                simRes(n) = max(confMatrix(rpTestVals>0.0005))>conf;
            else
                % testing only at a single rp duration
                simRes(n) = confMatrix(find(rpTestVals>0.002,1))>conf;
            end
        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        [~, pci] = binofit(sum(simRes),nSim,0.05); % 95% confidence interval
        passErr(bidx,c,:) = pci*100;
        
    end
    % legH(bidx) = plot(contProp*100, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
    legH(bidx) = plotWithErrUL(contProp*100, passPct(bidx,:), squeeze(passErr(bidx,:,:)), colors(bidx,:)); 
    legH(bidx).Marker = 'o'; legH(bidx).MarkerFaceColor = colors(bidx,:);
        
    hold on; drawnow;
end

xlabel('Contamination (%)'); 
ylabel('Percent pass (%)');
addX(10); 
leg = legend(legH, array2stringCell(RPdurs*1000)); 
leg.Title.String = 'RP Duration (ms)';
legend boxoff
box off

print(f, fullfile(figSaveDir, 'Fig3d.pdf'), '-dpdf');

%% Run all simulations

tic
nSim = 10; % number of simulations
% RPdurs = [1.5 2 3 4 5 6]/1000; % true RP duration, s
% recDurs = [0.5 1 2 3]*3600; % recording duration, s
% contProp = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20]/100; % simulated proportion contamination
% baseRates = [0.5 1 2 5 10 20]; % rate of the true neuron
% confThreshes = [50:10:90]; % confidence we need to accept a neuron

RPdurs = [1.5 3 6]/1000; % true RP duration, s
recDurs = [0.5  2 ]*3600; % recording duration, s
contProp = [0  4   8  10  12   16  20]/100; % simulated proportion contamination
baseRates = [0.5 1  5 ]; % rate of the true neuron
confThreshes = [75 90]; % confidence we need to accept a neuron

contThresh = 10; % acceptable percentage of contamination when determining pass/fail

params = struct(); params.cont = contThresh;
params.contaminationThresh = contThresh;

paramsCompare = struct(); 
paramsCompare.contaminationThresh = contThresh;

totalidx = 1; 
totaln = numel(baseRates)*numel(contProp)*numel(RPdurs)*numel(recDurs)*numel(confThreshes);
passPct = nan(totaln,1); 
passPctLlobet1_5 = nan(totaln,1); 
passPctLlobet2 = nan(totaln,1); 
passPctLlobet3 = nan(totaln,1); 
passPctHill1_5 = nan(totaln,1); 
passPctHill2 = nan(totaln,1); 
passPctHill3 = nan(totaln,1); 
base_rate = nan(totaln,1); cont_prop = nan(totaln,1); RP_dur = nan(totaln,1); 
rec_dur = nan(totaln,1); conf_level = nan(totaln,1); 

% - can put recDur on the inside and just select subsets of spikes.
% - contaminating spikes don't depend on RPdur- could in principle generate
% these just once
% - could do some parfor

for bidx = 1:numel(baseRates)
    totalRate = baseRates(bidx);
    for cidx = 1:numel(contProp)
        
        % this calculation ensures that the contaminating spikes generated
        % at this rate do form the correct proportion of the total
        % contRate = contProp(cidx)*baseRate/(1-contProp(cidx));
        baseRate = (1-contProp(cidx))*totalRate;
        contRate = contProp(cidx)*totalRate;
        
        for RPidx = 1:numel(RPdurs)
            RPdur = RPdurs(RPidx);
            fprintf(1, 'br %d/%d cp %d/%d rp %d/%d\n', bidx, numel(baseRates), ...
                cidx, numel(contProp), RPidx, numel(RPdurs)); 
          
            
            for ridx = 1:numel(recDurs)
                recDur = recDurs(ridx); 
                params.recDur = recDur;
                paramsCompare.recDur = recDur;

                for confIdx = 1:numel(confThreshes)
                    confThresh = confThreshes(confIdx);
                    params.confidenceThresh = confThresh;
        
                    simRes = zeros(nSim, 7);
                    for n = 1:nSim

                        st = genST(baseRate, recDur, RPdur); % true spikes
                        contST = genST(contRate,recDur, 0); % contaminating spikes
                        combST = sort([st; contST]); % combined spike train

                        [passTest, confidence, contamination, timeOfLowestCont,...
                            nSpikesBelow2, confMatrix, cont, rp, nACG] ...
                            = slidingRP(combST, params);

                        simRes(n,1) = passTest;
                        
                        paramsCompare.rp = rp; paramsCompare.nACG = nACG; paramsCompare.spikeCount = numel(combST);
                        
                        paramsCompare.metricType = 'Llobet';
                        paramsCompare.RPdur = 0.0015;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,2) = passTest;

                        paramsCompare.RPdur = 0.002;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,3) = passTest;

                        paramsCompare.RPdur = 0.003;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,4) = passTest;

                        paramsCompare.metricType = 'Hill';
                        paramsCompare.RPdur = 0.0015;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,5) = passTest;

                        paramsCompare.RPdur = 0.002;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,6) = passTest;

                        paramsCompare.RPdur = 0.003;
                        [passTest, estContam] = RPmetric_Classic([], paramsCompare);
                        simRes(n,7) = passTest;


                    end

                    passPct(totalidx) = sum(simRes(:,1))/nSim*100;
                    passPctLlobet1_5(totalidx) = sum(simRes(:,2))/nSim*100;
                    passPctLlobet2(totalidx) = sum(simRes(:,3))/nSim*100;
                    passPctLlobet3(totalidx) = sum(simRes(:,4))/nSim*100;
                    passPctHill1_5(totalidx) = sum(simRes(:,5))/nSim*100;
                    passPctHill2(totalidx) = sum(simRes(:,6))/nSim*100;
                    passPctHill3(totalidx) = sum(simRes(:,7))/nSim*100;
                    base_rate(totalidx) = baseRate; 
                    cont_prop(totalidx) = contProp(cidx); 
                    RP_dur(totalidx) = RPdur; 
                    rec_dur(totalidx) = recDur; 
                    conf_level(totalidx) = confThresh; 
                    totalidx = totalidx+1; 
                    %[~, pci] = binofit(sum(simRes),nSim,0.05); % 95% confidence interval
                    %passErr(bidx,c,:) = pci*100;
                    
                end
            end
        end
    end
end

simDat = table(base_rate, cont_prop, RP_dur, rec_dur, conf_level, passPct, ...
    passPctLlobet1_5, passPctLlobet2, passPctLlobet3, passPctHill1_5, ...
    passPctHill2, passPctHill3);
% save simDat.mat simDat nSim
toc
%% example plot from simDat -- SEE simDatFigure.m for updated/better version

figure; 

br = unique(simDat.base_rate);
cont = unique(simDat.cont_prop); 
recDurs = unique(simDat.rec_dur); 

colors = myCopper(0.6, numel(baseRates)+1);
colors = colors(2:end,:);

clear legH
for bidx = 1:numel(br)
    incl = simDat.base_rate==br(bidx) & simDat.RP_dur == 0.003 & simDat.rec_dur==7200 & simDat.conf_level==90;
    legH(bidx) = plot(cont*100, simDat.passPct(incl), '.-', 'Color', colors(bidx,:)); hold on;
end

xlabel('Contamination (%)'); 
ylabel('Percent pass (%)');
addX(10); 
leg = legend(legH, array2stringCell(br)); 
leg.Title.String = 'Base rate (sp/s)';
legend boxoff
box off


%% Comparing Sliding RP to Llobet and Hill

% TODO: Is estContam returning the right thing?? 

clearvars -except figSaveDir
saveFig = false;
nSim = 500; % number of simulations
RPdur = 0.0025; % true RP duration, s
recDur = 3600*2; % recording duration, s
% contProp = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20 25 30]/100; % simulated proportion contamination
contProp = [0 4 10 14 16 18 19 20 21 22 25 30]/100; % simulated proportion contamination
%baseRates = [0.5 1 2 5 10 20]; % rate of the true neuron
baseRates = [0.5 1 10]; % rate of the true neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail
confThresh = 75; % confidence we need to accept a neuron

params = struct(); params.cont = contThresh;
params.contaminationThresh = contThresh;
params.confidenceThresh = confThresh;

paramsCompare = struct(); 
paramsCompare.metricType = 'Llobet';
paramsCompare.contaminationThresh = contThresh;
paramsCompare.recDur = recDur; 
paramsCompare.RPdur = 0.002;

plotTitles = {}; plotTitles{1} = 'Sliding RP';
plotTitles{2} = sprintf('%s %dms', paramsCompare.metricType, 1000*paramsCompare.RPdur);

f = figure; f.Color = 'w';
colors = myCopper(0.6, numel(baseRates)+1);
colors = colors(2:end,:);

passPct = zeros(numel(baseRates), numel(contProp));
passErr = zeros(numel(baseRates), numel(contProp),2);
passPctCompare = zeros(numel(baseRates), numel(contProp));
passErrCompare = zeros(numel(baseRates), numel(contProp),2);

for bidx = 1:numel(baseRates)
    %baseRate = baseRates(bidx)
    totalRate = baseRates(bidx)
    for c = 1:numel(contProp)
        
        % this calculation ensures that the contaminating spikes generated
        % at this rate do form the correct proportion of the total
        % contRate = contProp(c)*baseRate/(1-contProp(c));
        contRate = contProp(c)*totalRate;
        baseRate = (1-contProp(c))*totalRate;
        
        params.recDur = recDur;
        
        simRes = zeros(nSim, 1);
        simResCompare = zeros(nSim, 1);
        for n = 1:nSim
            
            st = genST(baseRate, recDur, RPdur); % true spikes
            contST = genST(contRate,recDur, 0); % contaminating spikes
            combST = sort([st; contST]); % combined spike train
 
%             [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
%             
%             simRes(n) = max(confMatrix(rpTestVals>0.0005))>conf;
            
            [passTest, confidence, contamination, timeOfLowestCont,...
                nSpikesBelow2, confMatrix, cont, rp, nACG] ...
                = slidingRP(combST, params);
            
            simRes(n) = passTest;

            paramsCompare.rp = rp; paramsCompare.nACG = nACG; paramsCompare.spikeCount = numel(combST);
            [passTest, estContam] = RPmetric_Classic([], paramsCompare);
            simResCompare(n) = passTest;
            %simResCompare(n) = estContam<=contThresh;

        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        [~, pci] = binofit(sum(simRes),nSim,0.05); % 95% confidence interval
        passErr(bidx,c,:) = pci*100;

        passPctCompare(bidx,c) = sum(simResCompare)/nSim*100;
        [~, pci] = binofit(sum(simResCompare),nSim,0.05); % 95% confidence interval
        passErrCompare(bidx,c,:) = pci*100;
        
    end

    subplot(1,2,1); 
    legH(bidx) = plotWithErrUL(contProp*100, passPct(bidx,:), squeeze(passErr(bidx,:,:)), colors(bidx,:)); 
    legH(bidx).Marker = 'o'; legH(bidx).MarkerFaceColor = colors(bidx,:);
    hold on; 
    
    subplot(1,2,2); 
    legH(bidx) = plotWithErrUL(contProp*100, passPctCompare(bidx,:), squeeze(passErrCompare(bidx,:,:)), colors(bidx,:)); 
    legH(bidx).Marker = 'o'; legH(bidx).MarkerFaceColor = colors(bidx,:);
    hold on; 
    
    
    drawnow;
end

for sp = 1:2
    subplot(1,2,sp)
    xlabel('Contamination (%)'); 
    ylabel('Percent pass (%)');
    ylim([0 100]); 
    addX(10); 
    box off
    title(plotTitles{sp})
    
end

leg = legend(legH, array2stringCell(baseRates)); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff

if saveFig
    print(f, fullfile(figSaveDir, 'FigLlobetCompare.pdf'), '-dpdf');
end

