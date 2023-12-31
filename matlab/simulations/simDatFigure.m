
%%

load simDat.mat

% Create UIFigure and hide until all components are created
f = uifigure('Visible', 'off');
f.Position = [100 100 628 552];
f.Name = 'Sliding RP Simulations';


% Create UIAxes
ax = uiaxes(f);
title(ax, 'Metric performance')
xlabel(ax, 'Contamination level (%)')
ylabel(ax, 'Percent passing (%)')
ax.Position = [72 256 496 282];

% Create SplitbyButtonGroup
bg = uibuttongroup(f);
bg.Title = 'Split by';
bg.Position = [353 18 145 208];
bg.SelectionChangedFcn = @(f,k)plotSimDat(f,k);

% Create BaserateButton
BaserateButton = uiradiobutton(bg);
BaserateButton.Text = 'Base rate';
BaserateButton.Position = [11 159 73 22];
BaserateButton.Value = true;

% Create RPdurationButton
RPdurationButton = uiradiobutton(bg);
RPdurationButton.Text = 'RP duration';
RPdurationButton.Position = [11 114 86 22];

% Create RecordingdurationButton
RecordingdurationButton = uiradiobutton(bg);
RecordingdurationButton.Text = 'Recording duration';
RecordingdurationButton.Position = [11 67 125 22];

% Create ConfidencethreshButton
ConfidencethreshButton = uiradiobutton(bg);
ConfidencethreshButton.Text = 'Confidence thresh';
ConfidencethreshButton.Position = [11 23 120 22];

% Create Sliders
br = unique(simDat.base_rate);
rpd = unique(simDat.RP_dur); 
recDurs = unique(simDat.rec_dur); 
conf = unique(simDat.conf_level);

BaserateSlider = uislider(f);
BaserateSlider.Position = [181 187 150 3];
BaserateSlider.ValueChangedFcn = @(f,k)plotSimDat(f,k);
BaserateSlider.MajorTicks = 1:numel(br);
BaserateSlider.MajorTickLabels = array2stringCell(br);
BaserateSlider.Limits = [1 numel(br)];
BaserateSlider.MinorTicks = [];

RPdurationSlider = uislider(f);
RPdurationSlider.Position = [181 142 150 3];
RPdurationSlider.ValueChangedFcn = @(f,k)plotSimDat(f,k);
RPdurationSlider.MajorTicks = 1:numel(rpd);
RPdurationSlider.MajorTickLabels = array2stringCell(rpd*1000);
RPdurationSlider.Limits = [1 numel(rpd)];
RPdurationSlider.MinorTicks = [];

RecordingdurationSlider = uislider(f);
RecordingdurationSlider.Position = [181 95 150 3];
RecordingdurationSlider.ValueChangedFcn = @(f,k)plotSimDat(f,k);
RecordingdurationSlider.MajorTicks = 1:numel(recDurs);
RecordingdurationSlider.MajorTickLabels = array2stringCell(recDurs/3600);
RecordingdurationSlider.Limits = [1 numel(recDurs)];
RecordingdurationSlider.MinorTicks = [];

ConfidencethreshSlider = uislider(f);
ConfidencethreshSlider.Position = [181 51 150 3];
ConfidencethreshSlider.ValueChangedFcn = @(f,k)plotSimDat(f,k);
ConfidencethreshSlider.MajorTicks = 1:numel(conf);
ConfidencethreshSlider.MajorTickLabels = array2stringCell(conf);
ConfidencethreshSlider.Limits = [1 numel(conf)];
ConfidencethreshSlider.MinorTicks = [];

% Show the figure after all components are created
f.Visible = 'on';

ud = struct();
ud.simDat = simDat;
ud.bg = bg;
ud.brSlider = BaserateSlider;
ud.rpSlider = RPdurationSlider;
ud.rdSlider = RecordingdurationSlider;
ud.confSlider = ConfidencethreshSlider;
ud.ax = ax;
f.UserData = ud;