
% TODO: define corrected metric as a parameter of slidingRP function 

%% paths

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'roth-et-al-2026')))

addpath(genpath(fullfile(githubDir, 'spikes')))

%%


tic
nSim = 1000; % number of simulations

RPdurs = [1.5 10]/1000; % true RP duration, s
recDurs = [2 ]*3600; % recording duration, s
contProp = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20]/100; % simulated proportion contamination
baseRates = [5]; % rate of the true neuron
confThreshes = [90]; % confidence we need to accept a neuron

contThresh = 10; % acceptable percentage of contamination when determining pass/fail

params = struct(); params.cont = contThresh;
params.contaminationThresh = contThresh;

paramsCompare = struct(); 
paramsCompare.contaminationThresh = contThresh;

totalidx = 1; 
totaln = numel(baseRates)*numel(contProp)*numel(RPdurs)*numel(recDurs)*numel(confThreshes);
passPct = nan(totaln,1); 
passPctCorrected = nan(totaln,1); 
 
total_rate = nan(totaln,1); cont_prop = nan(totaln,1); RP_dur = nan(totaln,1); 
rec_dur = nan(totaln,1); conf_level = nan(totaln,1); 

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
        


                    simRes1 = zeros(nSim, 1);
                    simRes2 = zeros(nSim, 1);
                    
                    % >>> PARFOR LOOP STARTS HERE <<<
                    parfor n = 1:nSim
                        % Generate spikes
                        st = genST(baseRate, recDur, RPdur); 
                        contST = genST(contRate, recDur, 0); 
                        combST = sort([st; contST]); 

                        % Only request the first output (passTest) to save IPC overhead
                        paramsCorr = params; paramsCorr.correction = true;
                        pass1 = slidingRP(combST, params);            % standard
                        pass2 = slidingRP(combST, paramsCorr);        % FWER-corrected

                        simRes1(n) = pass1;
                        simRes2(n) = pass2;
                    end

                    passPct(totalidx) = sum(simRes1)/nSim*100;
                    passPctCorrected(totalidx) = sum(simRes2)/nSim*100;
                    
                    total_rate(totalidx) = totalRate; 
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

simDat = table(total_rate, cont_prop, RP_dur, rec_dur, conf_level, passPct, ...
    passPctCorrected);

toc

%% plot

f = figure; f.Color = 'w';

subplot(1,2,1); hold on;

thisRecDur = 7200; 
thisRPdur = 0.01; 
thisConfLevel = 90; 

theseFR = unique(simDat.total_rate);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.6, numel(theseFR)+1);
colors = colors(2:end,:);
colors2 = myCopper(0.3, numel(theseFR)+1);
colors2 = colors2(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
pidx = 1;
for rr = 1:numel(theseFR)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(pidx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(pidx).Marker = 'o'; legH(pidx).MarkerFaceColor = colors(rr,:); legH(pidx).MarkerSize = 3;
    pidx = pidx+1;
    
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctCorrected(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(pidx) = plotWithErrUL(theseCont*100, passPct, passCI, colors2(rr,:)); 
    legH(pidx).Marker = 'o'; legH(pidx).MarkerFaceColor = colors2(rr,:); legH(pidx).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Conf. = %d', thisRecDur/3600, thisRPdur*1000, thisConfLevel)); 
addX(10); addY(10);
% leg = legend(legH, array2stringCell(theseFR)); 
leg = legend(legH, {'Sliding RP', 'Mult Compare Corr'});
% leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off
axis square;


subplot(1,2,2); hold on;

thisRecDur = 7200; 
thisRPdur = 0.0015; 
thisConfLevel = 90; 

theseFR = unique(simDat.total_rate);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.6, numel(theseFR)+1);
colors = colors(2:end,:);
colors2 = myCopper(0.3, numel(theseFR)+1);
colors2 = colors2(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
pidx = 1;
for rr = 1:numel(theseFR)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(pidx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(pidx).Marker = 'o'; legH(pidx).MarkerFaceColor = colors(rr,:); legH(pidx).MarkerSize = 3;
    pidx = pidx+1;
    
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctCorrected(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(pidx) = plotWithErrUL(theseCont*100, passPct, passCI, colors2(rr,:)); 
    legH(pidx).Marker = 'o'; legH(pidx).MarkerFaceColor = colors2(rr,:); legH(pidx).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Conf. = %d', thisRecDur/3600, thisRPdur*1000, thisConfLevel)); 
addX(10); addY(10);
% leg = legend(legH, array2stringCell(theseFR)); 
leg = legend(legH, {'Sliding RP', 'Mult Compare Corr'});
% leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off

axis square;