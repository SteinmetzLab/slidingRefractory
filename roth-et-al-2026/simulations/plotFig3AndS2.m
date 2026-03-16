

% Fig 3 


%% load

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'roth-et-al-2026')))
addpath(genpath(fullfile(githubDir, 'spikes'))) % github/cortex-lab/spikes

% Using simDat from runSimulations
load simDat.mat

%%


% panels as follows, in a 3x3 grid: 
% 3a, perf by firing rate
% S2a, perf by rec dur
% 3b, perf by refractory period
% 3c, Llobet by firing rate
% S2b, Llobet by rec dur
% 3d, Llobet by refractory period
% 3e, overlay Sliding/Llobet by firing rate
% S2c, overlay by rec dur
% 3f, overlay by RP dur

f = figure; f.Color = 'w'; 

%% sliding RP by FR

subplot(3,3,1); hold on; 

thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 

theseFR = unique(simDat.total_rate);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.6, numel(theseFR)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
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

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Conf. = %d', thisRecDur/3600, thisRPdur*1000, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell(theseFR)); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off

%% sliding RP by rec dur

subplot(3,3,2); hold on; 

% thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRecDur = unique(simDat.rec_dur);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.3, numel(theseRecDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseRecDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==theseRecDur(rr) & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; RP dur. = %.1f ms; Conf. = %d', thisFR, thisRPdur*1000, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell(theseRecDur/3600)); 
leg.Title.String = 'Rec. dur. (hr)';
legend boxoff
box off

%% sliding RP by RP dur

subplot(3,3,3); hold off; 

thisRecDur = 7200; 
% thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRPDur = unique(simDat.RP_dur);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.8, numel(theseRPDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseRPDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==theseRPDur(rr) & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); hold on;
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; Rec. dur. = %.1f hr; Conf. = %d', thisFR, thisRecDur/3600, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell(theseRPDur*1000)); 
leg.Title.String = 'RP dur. (ms)';
legend boxoff
box off

%% llobet by FR

subplot(3,3,4); hold on; 

thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 

theseFR = unique(simDat.total_rate);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.1, numel(theseFR)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseFR)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms', thisRecDur/3600, thisRPdur*1000)); 
addX(10); 
leg = legend(legH, array2stringCell(theseFR)); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off

%% llobet by rec dur

subplot(3,3,5); hold on; 

% thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRecDur = unique(simDat.rec_dur);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.9, numel(theseRecDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseRecDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==theseRecDur(rr) & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; RP dur. = %.1f ms', thisFR, thisRPdur*1000)); 
addX(10); 
leg = legend(legH, array2stringCell(theseRecDur/3600)); 
leg.Title.String = 'Rec. dur. (hr)';
legend boxoff
box off

%% llobet by RP dur


subplot(3,3,6); hold on; 

thisRecDur = 7200; 
% thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRPDur = unique(simDat.RP_dur);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.4, numel(theseRPDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseRPDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==theseRPDur(rr) & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; Rec. dur. = %.1f hr', thisFR, thisRecDur/3600)); 
addX(10); 
leg = legend(legH, array2stringCell(theseRPDur*1000)); 
leg.Title.String = 'RP dur. (ms)';
legend boxoff
box off

%% overlay by FR
% include 1 and 10 sp/s

subplot(3,3,7); hold off; 

thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 

theseFR = unique(simDat.total_rate);
plotFR = [1 10];
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.6, numel(theseFR)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH; plotIdx = 1;
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

    if ismember(theseFR(rr),plotFR)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); hold on;
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

colors = myCopper(0.1, numel(allFR)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
% clear legH
for rr = 1:numel(theseFR)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==theseFR(rr) & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    if ismember(theseFR(rr),plotFR)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Conf. = %d', thisRecDur/3600, thisRPdur*1000, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell([plotFR plotFR])); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off

%% overlay by rec dur

subplot(3,3,8); hold off; 

% thisRecDur = 7200; 
thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRecDur = unique(simDat.rec_dur);
plotRecDur = 3600*[0.5 3];
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.3, numel(theseRecDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH; plotIdx = 1;
for rr = 1:numel(theseRecDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==theseRecDur(rr) & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    if ismember(theseRecDur(rr), plotRecDur)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); hold on; 
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

colors = myCopper(0.9, numel(theseRecDur)+1);
colors = colors(2:end,:);

for rr = 1:numel(theseRecDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==theseRecDur(rr) & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    if ismember(theseRecDur(rr), plotRecDur)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:));
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; RP dur. = %.1f ms; Conf. = %d', thisFR, thisRPdur*1000, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell([plotRecDur plotRecDur]/3600)); 
leg.Title.String = 'Rec. dur. (hr)';
legend boxoff
box off

%% overlay by RP dur

subplot(3,3,9); hold off; 

thisRecDur = 7200; 
% thisRPdur = 0.003; 
thisConfLevel = 90; 
thisFR = 5; 

theseRPDur = unique(simDat.RP_dur);
plotRPDur = [0.002 0.005];
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.8, numel(theseRPDur)+1);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH; plotIdx = 1;
for rr = 1:numel(theseRPDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==theseRPDur(rr) & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    if ismember(theseRPDur(rr), plotRPDur)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); hold on;
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

colors = myCopper(0.4, numel(theseRPDur)+1);
colors = colors(2:end,:);

for rr = 1:numel(theseRPDur)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==theseRPDur(rr) & ...
            simDat.conf_level==thisConfLevel & ...
            simDat.total_rate==thisFR & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPctLlobet3(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    if ismember(theseRPDur(rr), plotRPDur)
        legH(plotIdx) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); hold on;
        legH(plotIdx).Marker = 'o'; legH(plotIdx).MarkerFaceColor = colors(rr,:); legH(plotIdx).MarkerSize = 3;
        plotIdx = plotIdx+1;
    end
end

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('FR = %.1f sp/s; Rec. dur. = %.1f hr; Conf. = %d', thisFR, thisRecDur/3600, thisConfLevel)); 
addX(10); 
leg = legend(legH, array2stringCell([plotRPDur plotRPDur]*1000)); 
leg.Title.String = 'RP dur. (ms)';
legend boxoff
box off

