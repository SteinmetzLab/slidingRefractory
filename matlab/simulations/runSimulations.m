
function simDat = runSimulations(nSim, savePath)
% RUNSIMULATIONS  Run the full sliding-RP simulation battery and save results.
%
%   simDat = runSimulations(nSim)
%   simDat = runSimulations(nSim, savePath)
%
%   Sweeps over all combinations of firing rate, contamination proportion,
%   RP duration, recording duration, and confidence threshold.  For each
%   parameter combination, nSim synthetic spike trains are generated and
%   evaluated by slidingRP() and by RPmetric_Classic() at three fixed RP
%   durations (1.5, 2, 3 ms) for both the Llobet and Hill formulations.
%
%   A single parfor block over all parameter combinations is used (rather
%   than parfor on the inner nSim loop) to avoid repeated parfor startup
%   overhead and to maximise worker utilisation.
%
%   INPUTS
%     nSim      - Number of simulations per parameter combination.
%                 Use nSim=10 for a quick functional check; nSim=1000 for
%                 publication-quality results.
%     savePath  - (optional) Full path for the output .mat file.
%                 Default: <this file's directory>/simDat.mat
%
%   OUTPUT
%     simDat - Table with one row per parameter combination and columns:
%                total_rate, cont_prop, RP_dur, rec_dur, conf_level,
%                passPct, passPctLlobet1_5, passPctLlobet2, passPctLlobet3,
%                passPctHill1_5, passPctHill2, passPctHill3
%              Also saved to savePath as 'simDat' and 'nSim'.
%
%   REQUIRES
%     slidingRP, RPmetric_Classic, genST  (matlab/ directory)
%     histdiff from cortex-lab/spikes on the MATLAB path

if nargin < 2 || isempty(savePath)
    savePath = fullfile(fileparts(mfilename('fullpath')), 'simDat.mat');
end

%% Simulation parameters
RPdurs       = [1.5 2 3 4 5 6] / 1000;   % true RP duration (s)
recDurs      = [0.5 1 2 3] * 3600;        % recording duration (s)
contProp     = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20] / 100;
baseRates    = [0.5 1 2 5 10 20];         % total spike rate (sp/s)
confThreshes = 50:10:90;                  % confidence threshold (%)
contThresh   = 10;                        % contamination threshold (%)

%% Flatten all parameter combinations into index vectors
[bidxV, cidxV, RPidxV, ridxV, confIdxV] = ndgrid( ...
    1:numel(baseRates), 1:numel(contProp), 1:numel(RPdurs), ...
    1:numel(recDurs),   1:numel(confThreshes));
bidxV    = bidxV(:);
cidxV    = cidxV(:);
RPidxV   = RPidxV(:);
ridxV    = ridxV(:);
confIdxV = confIdxV(:);
totaln   = numel(bidxV);

fprintf('Running %d parameter combinations × %d sims = %d total simulations\n', ...
    totaln, nSim, totaln * nSim);

%% Base parameter structs
params0 = struct('cont', contThresh, 'contaminationThresh', contThresh);
paramsCompare0 = struct('contaminationThresh', contThresh);

%% Pre-allocate result arrays (rows = parameter combos)
pctAll   = zeros(totaln, 7);  % [passPct, Llobet1.5/2/3, Hill1.5/2/3]
paramOut = zeros(totaln, 5);  % [total_rate, cont_prop, RP_dur, rec_dur, conf_level]

%% Progress tracking via DataQueue
t0 = tic;
initProgress(totaln, t0);
q = parallel.pool.DataQueue;
afterEach(q, @printProgress);

%% Single parfor over all parameter combinations
parfor idx = 1:totaln
    bI  = bidxV(idx);    cI  = cidxV(idx);
    rpI = RPidxV(idx);   rI  = ridxV(idx);   cfI = confIdxV(idx);

    totalRate  = baseRates(bI);
    baseRate   = (1 - contProp(cI)) * totalRate;
    contRate   =      contProp(cI)  * totalRate;
    RPdur      = RPdurs(rpI);
    recDur     = recDurs(rI);
    confThresh = confThreshes(cfI);

    p = params0;
    p.recDur           = recDur;
    p.confidenceThresh = confThresh;

    simRes = zeros(nSim, 7);
    for n = 1:nSim
        st     = genST(baseRate, recDur, RPdur);
        contST = genST(contRate,  recDur, 0);
        combST = sort([st; contST]);

        [passTest, ~, ~, ~, ~, ~, ~, rp, nACG] = slidingRP(combST, p);

        pc            = paramsCompare0;
        pc.recDur     = recDur;
        pc.rp         = rp;
        pc.nACG       = nACG;
        pc.spikeCount = numel(combST);

        pc.metricType = 'Llobet';
        pc.RPdur = 0.0015; r2 = RPmetric_Classic([], pc);
        pc.RPdur = 0.002;  r3 = RPmetric_Classic([], pc);
        pc.RPdur = 0.003;  r4 = RPmetric_Classic([], pc);

        pc.metricType = 'Hill';
        pc.RPdur = 0.0015; r5 = RPmetric_Classic([], pc);
        pc.RPdur = 0.002;  r6 = RPmetric_Classic([], pc);
        pc.RPdur = 0.003;  r7 = RPmetric_Classic([], pc);

        simRes(n, :) = [passTest, r2, r3, r4, r5, r6, r7];
    end

    pctAll(idx, :)   = sum(simRes) / nSim * 100;
    paramOut(idx, :) = [totalRate, contProp(cI), RPdur, recDur, confThresh];

    send(q, idx);
end

elapsed = toc(t0);
fprintf('Done. Total time: %.0f s (%.1f min)\n', elapsed, elapsed/60);

%% Assemble output table
total_rate       = paramOut(:,1);
cont_prop        = paramOut(:,2);
RP_dur           = paramOut(:,3);
rec_dur          = paramOut(:,4);
conf_level       = paramOut(:,5);
passPct          = pctAll(:,1);
passPctLlobet1_5 = pctAll(:,2);
passPctLlobet2   = pctAll(:,3);
passPctLlobet3   = pctAll(:,4);
passPctHill1_5   = pctAll(:,5);
passPctHill2     = pctAll(:,6);
passPctHill3     = pctAll(:,7);

simDat = table(total_rate, cont_prop, RP_dur, rec_dur, conf_level, passPct, ...
    passPctLlobet1_5, passPctLlobet2, passPctLlobet3, ...
    passPctHill1_5, passPctHill2, passPctHill3);

save(savePath, 'simDat', 'nSim');
fprintf('Saved results to %s\n', savePath);

end % main function


%% --- Local helper functions ---

function initProgress(totaln, t0)
% Initialise persistent state for printProgress.
persistent pTotaln pT0 pCount pReportEvery
pTotaln     = totaln;
pT0         = t0;
pCount      = 0;
pReportEvery = max(1, floor(totaln / 50));  % report ~every 2%
end

function printProgress(~)
% Called by DataQueue afterEach; uses persistent state set by initProgress.
persistent pTotaln pT0 pCount pReportEvery
pCount = pCount + 1;
if mod(pCount, pReportEvery) == 0 || pCount == pTotaln
    elapsed = toc(pT0);
    eta     = elapsed / pCount * (pTotaln - pCount);
    fprintf('[%5.1f%%] %d/%d combos | %4.0f s elapsed | ETA ~%4.0f s (~%.1f min)\n', ...
        pCount/pTotaln*100, pCount, pTotaln, elapsed, eta, eta/60);
end
end
