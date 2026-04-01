%% script_runAndPlotSimulations
% Run the sliding-RP simulation battery and generate figures.
%
% Steps:
%   1. Run simulations  → saves simDat.mat in this directory
%   2. Plot Fig 3 / S2  (pass-rate vs contamination, multiple panels)
%   3. Plot Fig 4       (confidence threshold / ROC panels)
%
% For a quick functional check use nSim = 10.
% For publication-quality results use nSim = 1000.

%% Paths
thisDir   = fileparts(mfilename('fullpath'));
matlabDir = fullfile(thisDir, '..');           % matlab/
repoRoot  = fullfile(matlabDir, '..');         % repo root
simFigDir = fullfile(repoRoot, 'roth-et-al-2026', 'simulations');

addpath(genpath(matlabDir));
addpath(genpath(simFigDir));
addpath(genpath(fullfile(githubDir, 'spikes'))); % cortex-lab/spikes (histdiff)

%% Parameters
nSim     = 10;    % set to 1000 for publication results
savePath = fullfile(thisDir, 'simDat.mat');

%% 1. Run simulations
simDat = runSimulations(nSim, savePath); %#ok<NASGU>

%% 2 & 3. Plot figures (scripts use load simDat.mat from the working dir)
prevDir = pwd;
cd(thisDir);

run(fullfile(simFigDir, 'plotFig3AndS2.m'));
run(fullfile(simFigDir, 'plotFig4.m'));

cd(prevDir);
