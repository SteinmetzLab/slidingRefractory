

function plotSimDat(hCalling,keyDat)


f = hCalling.Parent;
ud = f.UserData;
simDat = ud.simDat;
bg = ud.bg;
ax = ud.ax;


br = unique(simDat.base_rate);
cont = unique(simDat.cont_prop); 
rpd = unique(simDat.RP_dur); 
recDurs = unique(simDat.rec_dur); 
conf = unique(simDat.conf_level);

% round sliders and get parameter selections
brSlider = ud.brSlider;
brSlider.Value = round(brSlider.Value);
defBR = br(brSlider.Value);

rpSlider = ud.rpSlider;
rpSlider.Value = round(rpSlider.Value);
defRPd = rpd(rpSlider.Value);

rdSlider = ud.rdSlider;
rdSlider.Value = round(rdSlider.Value);
defRD = recDurs(rdSlider.Value);

confSlider = ud.confSlider;
confSlider.Value = round(confSlider.Value);
defConf = conf(confSlider.Value);


switch bg.SelectedObject.Text
    case 'Base rate'
        hue = 0.6; 
        splitVar = br;    
        legVar = br; 
        titleStr = 'Base rate (sp/s)';
    case 'RP duration'
        hue = 0.8; 
        splitVar = rpd; 
        legVar = rpd*1000; 
        titleStr = 'RP duration (ms)';
    case 'Recording duration'
        hue = 0.3; 
        splitVar = recDurs; 
        legVar = recDurs/3600;
        titleStr = 'Recording duration (h)';
    case 'Confidence thresh'
        hue = 0.1; 
        splitVar = conf;
        legVar = conf;
        titleStr = 'Confidence level (%)';
end

colors = myCopper(hue, numel(splitVar)+1);
colors = colors(2:end,:);
hold(ax, 'off');

clear legH
for bidx = 1:numel(splitVar)

    switch bg.SelectedObject.Text
        case 'Base rate'
            incl = simDat.base_rate==splitVar(bidx) & simDat.RP_dur == defRPd & simDat.rec_dur==defRD & simDat.conf_level==defConf;
        case 'RP duration'
            incl = simDat.base_rate==defBR & simDat.RP_dur == splitVar(bidx) & simDat.rec_dur==defRD & simDat.conf_level==defConf;
        case 'Recording duration'
            incl = simDat.base_rate==defBR & simDat.RP_dur == defRPd & simDat.rec_dur==splitVar(bidx) & simDat.conf_level==defConf;
        case 'Confidence thresh'
            incl = simDat.base_rate==defBR & simDat.RP_dur == defRPd & simDat.rec_dur==defRD & simDat.conf_level==splitVar(bidx);
    end
    legH(bidx) = plot(ax,cont*100, simDat.passPct(incl), 'o-', 'LineWidth', 2.0, 'Color', colors(bidx,:));
    hold(ax,'on');
end

xlabel(ax,'Contamination (%)'); 
ylabel(ax,'Percent pass (%)');
plot(ax, [10 10], [0 100], 'k--');
leg = legend(legH, array2stringCell(legVar)); 
leg.Title.String = titleStr;
legend(ax,'boxoff');
box(ax,'off');