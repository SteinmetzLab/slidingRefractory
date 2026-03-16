
% TODO parfor at nSim loop

%% paths

addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'matlab')))
addpath(genpath(fullfile(githubDir, 'slidingRefractory', 'roth-et-al-2026')))
addpath(genpath(fullfile(githubDir, 'spikes'))) % github/cortex-lab/spikes

%% Run a single simulation with defined parameters and plot it

% Define parameters of simulation
RPdur = 2/1000; % true RP duration, s
recDur = [2]*3600; % recording duration, s
contProp = 95/100; % simulated proportion contamination
totalRate = 5; % rate of the true neuron, sp/s
confThresh = 90; % confidence we need to accept a neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail

params = struct(); 
params.contaminationThresh = contThresh;
params.confidenceThresh = confThresh;
params.recDur = recDur;


% simulate the spike train
baseRate = (1-contProp)*totalRate;
contRate = contProp*totalRate;

st = genST(baseRate, recDur, RPdur); % true spikes
contST = genST(contRate,recDur, 0); % contaminating spikes
combST = sort([st; contST]); % combined spike train

% run the test
[passTest, confidence, contamination, timeOfLowestCont,...
    nSpikesBelow2, confMatrix, cont, rp, nACG] ...
    = slidingRP(combST, params);

plotSlidingRP(combST, params);


%% Run all simulations

tic
nSim = 1000; % number of simulations
RPdurs = [1.5 2 3 4 5 6]/1000; % true RP duration, s
recDurs = [0.5 1 2 3]*3600; % recording duration, s
contProp = [0 2 4 6 7 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 16 18 20]/100; % simulated proportion contamination
baseRates = [0.5 1 2 5 10 20]; % rate of the true neuron
confThreshes = [50:10:90]; % confidence we need to accept a neuron

% RPdurs = [1.5 3 6]/1000; % true RP duration, s
% recDurs = [0.5  2 ]*3600; % recording duration, s
% contProp = [0  4   8  10  12   16  20]/100; % simulated proportion contamination
% baseRates = [0.5 1  5 ]; % rate of the true neuron
% confThreshes = [75 90]; % confidence we need to accept a neuron

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
total_rate = nan(totaln,1); cont_prop = nan(totaln,1); RP_dur = nan(totaln,1); 
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
    passPctLlobet1_5, passPctLlobet2, passPctLlobet3, passPctHill1_5, ...
    passPctHill2, passPctHill3);
save simDat.mat simDat nSim
toc