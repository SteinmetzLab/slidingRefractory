
%% paths

addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'nickBox')))

%% Fig 2 example neuron test



%% Fig 3a: performance across firing rates

nSim = 100; % number of simulations
RPdur = 0.0025; % true RP duration, s
recDur = 3600*2; % recording duration, s
contProp = [0:22]/100; % simulated proportion contamination
baseRates = [0.5 1 5 10]; % rate of the true neuron
contThresh = 10; % acceptable percentage of contamination when determining pass/fail
conf = 90; % confidence we need to accept a neuron

params = struct(); params.cont = contThresh;

testSliding = true; % alternative is to test at a single cutoff value

f = figure; f.Color = 'w';
colors = myCopper(0.6, numel(baseRates)+1);
colors = colors(2:end,:);
for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx)
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
        
        
    end
    legH(bidx) = plot(contProp*100, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
        
    hold on; drawnow;
end

xlabel('Contamination (%)'); 
ylabel('Percent pass (%)');
addX(10); 
leg = legend(legH, array2stringCell(baseRates)); 
leg.Title.String = 'Firing rate (sp/s)';
legend boxoff
box off
