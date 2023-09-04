

%% simulation test


clear all; 

nSim = 100; 

binSize = 1/30000;
b = (1:300)/30000;

baseRates = logspace(-0.3,1.3,8); baseRates = baseRates(3:end); 
% contPct = 0:0.02:0.2;
contPct = 0.095:0.0025:0.115;
% baseRates = logspace(-0.3,1.3,4); 
% contPct = 0:0.05:0.2;
recDur = 3600; 
rp = 0.002;

params = struct(); params.cont = 10;

figure; 
passPct = zeros(numel(baseRates), numel(contPct));
colors = hsv(numel(baseRates)); 
for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx)
    for c = 1:numel(contPct)
        
        %contRate = baseRate*contPct(c)
        contRate = contPct(c)*baseRate/(1-contPct(c))
        
        simRes = zeros(nSim, 1);
        for n = 1:nSim
            
%             st = genST(baseRate,recDur);
%             isi = diff([0; st]); isi(isi<rp) = [];
%             st = cumsum(isi);
%             
%             contST = genST(contRate,recDur);
%             
%             combST = sort([st; contST]);

            st = genST(baseRate, recDur, rp); 
            contST = genST(contRate,recDur, 0);
            combST = sort([st; contST]);
 
            [confMatrix, cont, rpTestVals, ~, ~] = computeMatrix(combST, params);
            %simRes(n) = max(confMatrix(rpTestVals>0.0005))>90;
            simRes(n) = confMatrix(find(rpTestVals>0.0015,1))>90;
        end
        
        passPct(bidx,c) = sum(simRes)/nSim*100;
        
        
    end
    legH(bidx) = plot(contPct, passPct(bidx,:), 'o-', 'Color', colors(bidx,:), 'MarkerFaceColor', colors(bidx,:)); 
        
    hold on; drawnow;
end
addX(0.1);
legend(legH, array2stringCell(baseRates));
set(gcf, 'Color', 'w'); 
box off; 
xlabel('True contamination proportion'); 
ylabel('Percentage of simulations that pass');  


%% test the 'expected count'

baseRates = logspace(-0.6,1.5,80); 
recDur = 7200; 
rp = 0.002;
contProp = 0.1;
params = struct(); params.cont = 10;

for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx);
%     contRate = baseRate*contProp;
    contRate = contProp*baseRate/(1-contProp);
    
%     st = genST(baseRate,recDur);
%     isi = diff([0; st]); isi(isi<rp) = [];
%     st = cumsum(isi);
%     
%     contST = genST(contRate,recDur);
%     
%     combST = sort([st; contST]);

    st = genST(baseRate, recDur, rp);
    contST = genST(contRate,recDur, 0);
    combST = sort([st; contST]);
    
    [confMatrix, cont, rpTestVals, nACG, firingRate] = computeMatrix(combST, params);
    
    %contaminationRate = firingRate*contProp; 
    %expectedViol = contaminationRate*rp*2*numel(combST);
    
    expectedViol = contRate*rp*2*(baseRate+contRate)*recDur;
    
    obsExp(bidx,1) = sum(nACG(rpTestVals<=rp)); 
    obsExp(bidx,2) = expectedViol;
    obsExp(bidx,3) = firingRate; 
    obsExp(bidx, 4) = 2*rp*contRate * (baseRate*recDur + (contRate*recDur-1)/2);
end

f = figure; f.Color = 'w';

subplot(1,2,1); 
plot(baseRates, obsExp(:,1), 'r.-'); hold on;
plot(baseRates, obsExp(:,2), 'k.-');
plot(baseRates, obsExp(:,4), 'b.-');

xlabel('firing rate (sp/s)')
ylabel('count'); 
legend({'observed', 'expected', 'new expected'});

subplot(1,2,2); 

contRates = contProp.*baseRates./(1-contProp);
plot(baseRates, obsExp(:,3), 'k.-'); hold on;
plot(baseRates, obsExp(:,3)-contRates', 'r.-'); 
xlabel('intended firing rate (sp/s)')
ylabel('actual firing rate (sp/s)'); 
addXeqY


%% 

% test simulation: generate uncontaminated spike train, add 1 contaminating
% spike, 2 spikes, etc - where does prediction diverge from simulation?


%% some math

% hill formula is this: 
% expectedViol = contRate*rp*2*(baseRate+contRate)*recDur;
% = R_c * P * 2 * N_t
% = 2P/D * N_c * N_t
% i.e. rate of contaminating spikes times total time available to them
% two problems with this:
% 1) can't some contaminating spikes add a count > 1? if a contaminating
% spike lands between two spikes that are alread separated by < P, then you
% go from having one contaminating event to having 3. (if counting by ACG)
% 2) the effect of each contaminating spike on the total number of
% violations is dependent on how many contaminating events are there so far

% a different calculation: 
% - when you add the first contaminating spike, you can add up to 1
% violation, and the probability that you add a violation is given by: 
% P(v) = N_b * P * 2 / D = R_b * P * 2
% so the expected value of the number of violations is now: 
% E(v) = 1*P(v) + 0*(1-P(v)) = P(v)
% when you add the next contaminating spike, you now have this expression: 
% P(v) = (N_b+1) * P * 2 / D
% since the new spike can collide with any of the base spikes as well as
% the new contaminating spike
%
% **** what is the correction factor now for a double collision? i.e. if
% the last spike did land near another spike, then the total time
% available for collisions with this spike is smaller than what's written
% there. proceeding by ignoring this for now
%
% And as above, E(v) = P(v). likewise for the next spike, 
% E(v) = (N_b+2)*P*2/D, etc
% So to get the total expected value of adding N_c contaminating spikes, we
% sum all these E(v): 
% E_total(v) = sum over i where i is (0->N_c-1) of (N_b+i)*P*2/D
% = 2P/D * sum(N_b+i)
% = 2P/D * (N_b*N_c + sum(i))
% = 2P/D * (N_b*N_c + N_c*(N_c-1)/2)
% = 2P/D * N_c * (N_b + (N_c-1)/2)
% so this is definitely a smaller expected total than we had before

% to find the difference we can rewrite hill as: 
% 2P/D * N_c * (N_c + N_b)
% = 2P/D * N_c * N_b + 2P/D * N_c^2
% and our calculation as 
% = 2P/D * N_c * N_b + 2P/D * N_c * (N_c-1)/2)
% the difference is: 
% delta = 2P/D * N_c * (N_c - (N_c-1)/2)
% = 2P/D * N_c * ( (2*N_c-N_c+1)/2 ) 
% = 2P/D * N_c * (N_c+1)/2

% let's just plot what these look like

D = 3600; 
P = 0.002; 
params = struct(); params.cont = 10; 

R_c = 0:0.05:2; 
R_b = 5; 
N_c = R_c*D;
N_b = R_b*D;

hill = 2*P/D * N_c .* (N_c + N_b);
newcalc = 2*P/D * N_c .* (N_b + (N_c-1)/2);
delta = 2*P/D * N_c .* (N_c+1)/2;
empdelta = hill - newcalc; 

figure; co = get(gca, 'ColorOrder'); 
plot(R_c, hill, 'k.-'); hold on; 
plot(R_c, newcalc, '.-', 'Color', co(1,:)); 
plot(R_c, delta, '.-', 'Color', co(2,:)); 
plot(R_c, empdelta, 'o-', 'Color', co(2,:)); 

% add a simulation
n_reps = 100; 
actViol = zeros(n_reps, numel(R_c));
actViolACG = zeros(n_reps, numel(R_c));
for q = 1:n_reps
for idx = 1:numel(R_c)
    st = genST(R_b, D, P);
    contST = genST(R_c(idx),D, 0);
    combST = sort([st; contST]);
    
    actViol(q,idx) = sum(diff(combST)<P);
    
    [confMatrix, cont, rpTestVals, nACG, firingRate] = computeMatrix(combST, params);
    actViolACG(q,idx) = sum(nACG(rpTestVals<P));
end
end

plot(R_c, mean(actViol,1), 'o-', 'Color', co(4,:)); 
plot(R_c, mean(actViolACG,1), 'o-', 'Color', co(5,:)); 


legend({'Hill', 'New calc.', 'Difference', 'Empirical Diff', 'Simulation: ISI', 'Simulation: ACG'});
xlabel('Contamination rate (sp/s)'); 
ylabel('N ISI violations'); 