
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
%   The inner nSim loop is parallelised with parfor (requires Parallel
%   Computing Toolbox; falls back to a serial for-loop otherwise).
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
RPdurs      = [1.5 2 3 4 5 6] / 1000;   % true RP duration (s)
recDurs     = [0.5 1 2 3] * 3600;        % recording duration (s)
contProp    = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20] / 100;
baseRates   = [0.5 1 2 5 10 20];         % total spike rate (sp/s)
confThreshes = 50:10:90;                 % confidence threshold (%)
contThresh  = 10;                        % contamination threshold (%)

totaln = numel(baseRates) * numel(contProp) * numel(RPdurs) * ...
         numel(recDurs) * numel(confThreshes);

% Pre-allocate result vectors
passPct         = nan(totaln, 1);
passPctLlobet1_5 = nan(totaln, 1);
passPctLlobet2  = nan(totaln, 1);
passPctLlobet3  = nan(totaln, 1);
passPctHill1_5  = nan(totaln, 1);
passPctHill2    = nan(totaln, 1);
passPctHill3    = nan(totaln, 1);
total_rate = nan(totaln, 1);
cont_prop  = nan(totaln, 1);
RP_dur     = nan(totaln, 1);
rec_dur    = nan(totaln, 1);
conf_level = nan(totaln, 1);

%% Parameter base structs (fields set inside loops)
params0 = struct();
params0.cont                = contThresh;
params0.contaminationThresh = contThresh;

paramsCompare0 = struct();
paramsCompare0.contaminationThresh = contThresh;

%% Sweep
totalidx = 1;
tic

for bidx = 1:numel(baseRates)
    totalRate = baseRates(bidx);

    for cidx = 1:numel(contProp)
        baseRate = (1 - contProp(cidx)) * totalRate;
        contRate =      contProp(cidx)  * totalRate;

        for RPidx = 1:numel(RPdurs)
            RPdur = RPdurs(RPidx);
            fprintf('br %d/%d  cp %d/%d  rp %d/%d\n', ...
                bidx, numel(baseRates), cidx, numel(contProp), RPidx, numel(RPdurs));

            for ridx = 1:numel(recDurs)
                recDur = recDurs(ridx);

                for confIdx = 1:numel(confThreshes)
                    confThresh = confThreshes(confIdx);

                    params = params0;
                    params.recDur           = recDur;
                    params.confidenceThresh = confThresh;

                    simRes = zeros(nSim, 7);

                    % --- parallelised inner loop ---
                    parfor n = 1:nSim
                        st      = genST(baseRate, recDur, RPdur);
                        contST  = genST(contRate,  recDur, 0);
                        combST  = sort([st; contST]);

                        [passTest, ~, ~, ~, ~, ~, ~, rp, nACG] = slidingRP(combST, params);

                        % Local copy avoids parfor broadcast-modification issues
                        pc = paramsCompare0;
                        pc.recDur      = recDur;
                        pc.rp          = rp;
                        pc.nACG        = nACG;
                        pc.spikeCount  = numel(combST);

                        pc.metricType = 'Llobet';
                        pc.RPdur = 0.0015;
                        r2 = RPmetric_Classic([], pc);

                        pc.RPdur = 0.002;
                        r3 = RPmetric_Classic([], pc);

                        pc.RPdur = 0.003;
                        r4 = RPmetric_Classic([], pc);

                        pc.metricType = 'Hill';
                        pc.RPdur = 0.0015;
                        r5 = RPmetric_Classic([], pc);

                        pc.RPdur = 0.002;
                        r6 = RPmetric_Classic([], pc);

                        pc.RPdur = 0.003;
                        r7 = RPmetric_Classic([], pc);

                        % Assign whole row so parfor sees a single subscript pattern
                        simRes(n, :) = [passTest, r2, r3, r4, r5, r6, r7];
                    end

                    passPct(totalidx)          = sum(simRes(:,1)) / nSim * 100;
                    passPctLlobet1_5(totalidx) = sum(simRes(:,2)) / nSim * 100;
                    passPctLlobet2(totalidx)   = sum(simRes(:,3)) / nSim * 100;
                    passPctLlobet3(totalidx)   = sum(simRes(:,4)) / nSim * 100;
                    passPctHill1_5(totalidx)   = sum(simRes(:,5)) / nSim * 100;
                    passPctHill2(totalidx)     = sum(simRes(:,6)) / nSim * 100;
                    passPctHill3(totalidx)     = sum(simRes(:,7)) / nSim * 100;
                    total_rate(totalidx) = totalRate;
                    cont_prop(totalidx)  = contProp(cidx);
                    RP_dur(totalidx)     = RPdur;
                    rec_dur(totalidx)    = recDur;
                    conf_level(totalidx) = confThresh;
                    totalidx = totalidx + 1;
                end
            end
        end
    end
end

toc

simDat = table(total_rate, cont_prop, RP_dur, rec_dur, conf_level, passPct, ...
    passPctLlobet1_5, passPctLlobet2, passPctLlobet3, ...
    passPctHill1_5, passPctHill2, passPctHill3);

save(savePath, 'simDat', 'nSim');
fprintf('Saved results to %s\n', savePath);
