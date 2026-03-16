%% Fig 4 (confidence levels)

% Using simDat from runSimulations
load simDat.mat

f = figure; f.Color = 'w'; 

%% sliding RP by Confidence

subplot(2,3,1); hold on; 

thisRecDur = 3600; 
thisRPdur = 0.003; 
% thisConfLevel = 90; 
thisRate = 5.0;

theseConf = unique(simDat.conf_level);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.5, numel(theseConf)+2);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseConf)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==theseConf(rr) & ...
            simDat.total_rate==thisRate & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

% add Llobet 
for cc = 1:numel(theseCont)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==theseCont(cc);
    thisPP = simDat.passPctLlobet2(thisRow);
    passPct(cc)=thisPP;
    [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
    passCI(cc,:) = ci*100;
end

rr = rr+1;
legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, 'r');
legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Rate = %.2f sp/s', thisRecDur/3600, thisRPdur*1000, thisRate)); 
addX(10); 
leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}]); 
leg.Title.String = 'Confidence';
legend boxoff
box off
axis square

%% ROC version

% picking false positive rate and true positive rate at 8 and 12%
% contamination, for different confidence values

subplot(2,3,4); hold on;

% implicitly using parameters from the cell above so that it matches
testContLower = 0.08; testContUpper = 0.12;

clear legH
for rr = 1:numel(theseConf)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContLower;
    lowerPP = simDat.passPct(thisRow);
    lowerLlobet = simDat.passPctLlobet2(thisRow);
    
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContUpper;
    upperPP = simDat.passPct(thisRow);
    upperLlobet = simDat.passPctLlobet2(thisRow);
    
    legH(rr) = plot(upperPP, lowerPP, 'o', 'MarkerFaceColor', colors(rr,:), 'MarkerEdgeColor', [1 1 1]);        
    
end

rr = rr+1;
legH(rr) = plot(upperLlobet, lowerLlobet, 'rx');
addXeqY();

leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}],'Location','southeast'); 
leg.Title.String = 'Confidence';
legend boxoff

xlabel('False passing rate'); 
ylabel('True passing rate'); 
xlim([0 100]); 
ylim([0 100]);
box off; 
axis square;


%% by confidence, for different parameters

subplot(2,3,2); hold on; 

thisRecDur = 3600; 
thisRPdur = 0.003; 
% thisConfLevel = 90; 
thisRate = 2;

theseConf = unique(simDat.conf_level);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.5, numel(theseConf)+2);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseConf)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==theseConf(rr) & ...
            simDat.total_rate==thisRate & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

% add Llobet 
for cc = 1:numel(theseCont)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==theseCont(cc);
    thisPP = simDat.passPctLlobet2(thisRow);
    passPct(cc)=thisPP;
    [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
    passCI(cc,:) = ci*100;
end

rr = rr+1;
legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, 'r');
legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Rate = %.2f sp/s', thisRecDur/3600, thisRPdur*1000, thisRate)); 
addX(10); 
leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}]); 
leg.Title.String = 'Confidence';
legend boxoff
box off
axis square

%%
subplot(2,3,5); hold on;

% implicitly using parameters from the cell above so that it matches
testContLower = 0.08; testContUpper = 0.12;

clear legH
for rr = 1:numel(theseConf)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContLower;
    lowerPP = simDat.passPct(thisRow);
    lowerLlobet = simDat.passPctLlobet2(thisRow);
    
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContUpper;
    upperPP = simDat.passPct(thisRow);
    upperLlobet = simDat.passPctLlobet2(thisRow);
    
    legH(rr) = plot(upperPP, lowerPP, 'o', 'MarkerFaceColor', colors(rr,:), 'MarkerEdgeColor', [1 1 1]);        
    
end


rr = rr+1;
legH(rr) = plot(upperLlobet, lowerLlobet, 'rx');
addXeqY();

leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}],'Location','southeast'); 
leg.Title.String = 'Confidence';
legend boxoff

xlabel('False passing rate'); 
ylabel('True passing rate'); 
xlim([0 100]); 
ylim([0 100]);
box off; 
axis square;

%%
subplot(2,3,3); hold on; 

thisRecDur = 3600; 
thisRPdur = 0.003; 
% thisConfLevel = 90; 
thisRate = 1;

theseConf = unique(simDat.conf_level);
theseCont = unique(simDat.cont_prop);

colors = myCopper(0.5, numel(theseConf)+2);
colors = colors(2:end,:);

passPct = zeros(numel(theseCont),1); passCI = zeros(numel(theseCont),2); 
clear legH
for rr = 1:numel(theseConf)
    for cc = 1:numel(theseCont)
        thisRow = simDat.rec_dur==thisRecDur & ...
            simDat.RP_dur==thisRPdur & ...
            simDat.conf_level==theseConf(rr) & ...
            simDat.total_rate==thisRate & ...
            simDat.cont_prop==theseCont(cc);
        thisPP = simDat.passPct(thisRow);
        passPct(cc)=thisPP; 
        [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
        passCI(cc,:) = ci*100; 
    end

    legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, colors(rr,:)); 
    legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;
end

% add Llobet 
for cc = 1:numel(theseCont)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==theseCont(cc);
    thisPP = simDat.passPctLlobet2(thisRow);
    passPct(cc)=thisPP;
    [~,ci] = binofit(round(thisPP*nSim/100), nSim, 0.05);
    passCI(cc,:) = ci*100;
end

rr = rr+1;
legH(rr) = plotWithErrUL(theseCont*100, passPct, passCI, 'r');
legH(rr).Marker = 'o'; legH(rr).MarkerFaceColor = colors(rr,:); legH(rr).MarkerSize = 3;

xlabel('True simulated contamination (%)'); 
ylabel('Percent pass (%)');
title(sprintf('Rec. dur. = %.1f hr; RP dur. = %.1f ms; Rate = %.2f sp/s', thisRecDur/3600, thisRPdur*1000, thisRate)); 
addX(10); 
leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}]); 
leg.Title.String = 'Confidence';
legend boxoff
box off
axis square

%%

subplot(2,3,6); hold on;

% implicitly using parameters from the cell above so that it matches
testContLower = 0.08; testContUpper = 0.12;

clear legH
for rr = 1:numel(theseConf)
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContLower;
    lowerPP = simDat.passPct(thisRow);
    lowerLlobet = simDat.passPctLlobet2(thisRow);
    
    thisRow = simDat.rec_dur==thisRecDur & ...
        simDat.RP_dur==thisRPdur & ...
        simDat.conf_level==theseConf(rr) & ...
        simDat.total_rate==thisRate & ...
        simDat.cont_prop==testContUpper;
    upperPP = simDat.passPct(thisRow);
    upperLlobet = simDat.passPctLlobet2(thisRow);
    
    legH(rr) = plot(upperPP, lowerPP, 'o', 'MarkerFaceColor', colors(rr,:), 'MarkerEdgeColor', [1 1 1]);        
    
end

rr = rr+1;
legH(rr) = plot(upperLlobet, lowerLlobet, 'rx');
addXeqY();

leg = legend(legH, [array2stringCell(theseConf); {'Llobet'}],'Location','southeast'); 
leg.Title.String = 'Confidence';
legend boxoff

xlabel('False passing rate'); 
ylabel('True passing rate'); 
xlim([0 100]); 
ylim([0 100]);
box off; 
axis square;
