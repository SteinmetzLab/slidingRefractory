

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
rp = 0.003;

params = struct(); params.cont = 10;

figure; 
passPct = zeros(numel(baseRates), numel(contPct));
colors = hsv(numel(baseRates)); 
for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx)
    for c = 1:numel(contPct)
        
        %contRate = baseRate*contPct(c)
        contRate = contPct(c)*baseRate/(1-contPct(c));
        
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
            
            % normal sliding rp test result
%             simRes(n) = max(confMatrix(rpTestVals>0.0005))>90;
            
            % testing only at a single rp duration
            simRes(n) = confMatrix(find(rpTestVals>0.002,1))>90;
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

D = 7200; 
P = 0.004; 
params = struct(); params.cont = 10; 

R_c = 0:0.2:8; 
R_b = 20; 
N_c = R_c*D;
N_b = R_b*D;
Nt = N_c+N_b;

hill = 2*P/D * N_c .* (N_c + N_b);
newcalc = 2*P/D * N_c .* (N_b + (N_c-1)/2);
delta = 2*P/D * N_c .* (N_c+1)/2;
empdelta = hill - newcalc; 

realhill = 2*P/D * N_c .* N_b;


figure; co = get(gca, 'ColorOrder'); 
plot(R_c, hill, 'k.-'); hold on; 
plot(R_c, newcalc, '.-', 'Color', co(1,:)); 
plot(R_c, realhill, '.-', 'Color', co(2,:)); 
% plot(R_c, delta, '.-', 'Color', co(2,:)); 
% plot(R_c, empdelta, 'o-', 'Color', co(2,:)); 

% add a simulation
n_reps = 10; 
actViol = zeros(n_reps, numel(R_c));
actViolACG = zeros(n_reps, numel(R_c));
for q = 1:n_reps
for idx = 1:numel(R_c)
    st = genST(R_b, D, P);
    contST = genST(R_c(idx),D, 0);
    combST = sort([st; contST]);
    
    actViol(q,idx) = sum(diff(combST)<P);
    
    [confMatrix, cont, rpTestVals, nACG, firingRate] = computeMatrix(combST, params);
    Nv =  sum(nACG(rpTestVals<P));
    actViolACG(q,idx) = Nv;
    
%     estContHill(q,idx) = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P)); 
%     estContNew(q,idx) = 1 - sqrt(P*(P*Nt^2-D*Nv))/(P*Nt); 
%     estContOld(q,idx) = Nv*D/(2*P*Nt^2);
    
end
end

plot(R_c, mean(actViol,1), 'o-', 'Color', co(4,:)); 
plot(R_c, mean(actViolACG,1), 'o-', 'Color', co(5,:)); 


% legend({'Hill', 'New calc.', 'Difference', 'Empirical Diff', 'Simulation: ISI', 'Simulation: ACG'});
legend({'Hill', 'New calc.', 'Real Hill', 'Simulation: ISI', 'Simulation: ACG'});

xlabel('Contamination rate (sp/s)'); 
ylabel('N ISI violations'); 

%% checking updated calculation

% if 
% Nv = 2*P/D * Nc * (Nb + Nc/2)
%    = 2*P/D * Nc * (Nt - Nc/2) 
% then proportion of contamination is: 
% Nc/Nt = 1 +/- sqrt(P*(P*Nt^2-D*Nv))/(P*Nt)
% (from wolfram alpha's solve)
%
% where Llobet gives: 
% Nc/Nt = 1 - sqrt(1 - Nv*D/(Nt^2*P))
% 
% and the old version from cortexlab/ecephys was: 
% Nc/Nt = Nv*D/(2*P*Nt^2)
% 
% and hill had: 
% Nv = 2*P/D * Nc * Nb
% which gives: 
% Nc/Nt = (P*Nt - sqrt(P(P*Nt^2-2*D*Nv)))/(2*P*Nt)

% so some example just to check

% D = 3600; 
% P = 0.002; 
% Nb = 4*D; 
% Nc = 1.65*D; 

D = 7200; 
P = 0.003; 
Nb = 2*D; 
Nc = 0.5*D; 

Nt = Nb+Nc; 

Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2)

C_true = Nc/Nt 

C_old = Nv*D/(2*P*Nt^2)

C_llobet =  1 - sqrt(1 - Nv*D/(Nt^2*P))

C_new = 1 - sqrt(P*(P*Nt^2-D*Nv))/(P*Nt)

C_hill = (P*Nt - sqrt(P*(P*Nt^2-2*D*Nv)))/(2*P*Nt)

C_llobet2 = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P))

% looks good. can use Llobet formula. 

%% testing the comparison between hill / old / corrected

D = 7200; 
P = 0.003; 
Nb = 2*D; 
Nc = 0.5*D; 
Nt = Nb+Nc; 
Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);

figure; set(gcf, 'Color', 'w')

subplot(2,2,1); 
Dvals = 300:300:7200;
clear C_old C_llobet C_hill
for idx = 1:numel(Dvals)
    D = Dvals(idx);
    Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);
    C_old(idx) = Nv*D/(2*P*Nt^2);
    C_llobet(idx) =  1 - sqrt(1 - Nv*D/(Nt^2*P));
    C_hill(idx) = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P));
end
co = get(gca, 'ColorOrder'); hold on;
plot(Dvals, C_llobet, '.-', 'Color', co(1,:));
plot(Dvals, C_old, '.-', 'Color', co(4,:)); 
plot(Dvals, C_hill, '.-', 'Color', co(2,:)); 
xlabel('Recording dur (s)'); ylabel('Reported contamination'); 
legend({'Llobet (true)', 'Old method', 'Hill'});

D = 7200; 
P = 0.003; 
Nb = 2*D; 
Nc = 0.5*D; 
Nt = Nb+Nc; 
Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);

subplot(2,2,2); 
Pvals = 0.001:0.001:0.008;
clear C_old C_llobet C_hill
for idx = 1:numel(Pvals)
    P = Pvals(idx); 
    Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);
    C_old(idx) = Nv*D/(2*P*Nt^2);
    C_llobet(idx) =  1 - sqrt(1 - Nv*D/(Nt^2*P));
    C_hill(idx) = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P));
end
co = get(gca, 'ColorOrder'); hold on;
plot(Pvals, C_llobet, '.-', 'Color', co(1,:));
plot(Pvals, C_old, '.-', 'Color', co(4,:)); 
plot(Pvals, C_hill, '.-', 'Color', co(2,:)); 
xlabel('Ref period duration (s)'); ylabel('Reported contamination'); 
legend({'Llobet (true)', 'Old method', 'Hill'});



D = 7200; 
P = 0.003; 
Nt = 3*D; 

subplot(2,2,3); 
contPropVals = 0:0.02:0.5;
clear C_old C_llobet C_hill isim
for idx = 1:numel(contPropVals)
    CP = contPropVals(idx); 
    Nc = CP*Nt;
    Nb = (1-CP)*Nt;    
    Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);
    C_old(idx) = Nv*D/(2*P*Nt^2);
    C_llobet(idx) =  1 - sqrt(1 - Nv*D/(Nt^2*P));
    ch = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P));
    C_hill(idx) = real(ch); 
    if ~isreal(ch); isim(idx) = true; else; isim(idx) = false; end;
end
co = get(gca, 'ColorOrder'); hold on;
plot(contPropVals, C_llobet, '.-', 'Color', co(1,:));
plot(contPropVals, C_old, '.-', 'Color', co(4,:)); 
plot(contPropVals, C_hill, '.-', 'Color', co(2,:)); 
plot(contPropVals(isim), C_hill(isim), 'ko'); 
xlabel('Proportion contamination'); ylabel('Reported contamination'); 
addXeqY
legend({'Llobet (true)', 'Old method', 'Hill', '(imaginary)'});


D = 7200; 
P = 0.003; 
CP = 0.2;

subplot(2,2,4); 
totalRateVals = [0.2:0.2:10]*D;
clear C_old C_llobet C_hill isim
for idx = 1:numel(totalRateVals)
    Nt = totalRateVals(idx);      
    Nc = CP*Nt;
    Nb = (1-CP)*Nt;    
    Nv = 2*P/D * Nc .* (Nb + (Nc-1)/2);
    C_old(idx) = Nv*D/(2*P*Nt^2);
    C_llobet(idx) =  1 - sqrt(1 - Nv*D/(Nt^2*P));
    ch = 1/2 * (1 - sqrt(1 - 2*Nv*D/Nt^2/P));
    C_hill(idx) = real(ch); 
    if ~isreal(ch); isim(idx) = true; else; isim(idx) = false; end;
end
co = get(gca, 'ColorOrder'); hold on;
plot(totalRateVals, C_llobet, '.-', 'Color', co(1,:));
plot(totalRateVals, C_old, '.-', 'Color', co(4,:)); 
plot(totalRateVals, C_hill, '.-', 'Color', co(2,:)); 
plot(totalRateVals(isim), C_hill(isim), 'ko'); 
xlabel('Total spike count'); ylabel('Reported contamination'); 
legend({'Llobet (true)', 'Old method', 'Hill', '(imaginary)'});