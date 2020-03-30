

% goal is to determine what firing rate and measured false positive rate is
% safe to include. 

% The point is that some neurons have truly really low refractory periods
% so you need to test those short windows. But if you test a short window
% the probability of getting a low score by chance is higher, so you will
% need a higher overall firing rate to accept that neuron. 

% By working out this relationship, we can fix our "false positive
% estimate" measure to not depend on chosen refractory period duration, and
% set a rational firing rate inclusion threshold, at the same time. 


% To do this, determine for each combination of firing rate and refractory
% period, what is the chance of failing to see violations by chance, for
% the upper limit of acceptable contaminating spikes. As refractory period
% is varied, each neuron will trace out a curve of observed violations. If
% for some choice of refractory period, we can conclude that we would not
% have obtained that low of observed violations by chance, then we can
% include the neuron. 


%% an example calculation

firingRate = 4; % sp/s
recDur = 3600; % s, recording duration
acceptableCont = 0.1*firingRate; % we can accept that 10% of spikes are contaminants
thresh = 0.1; % we can accept it if the true contamination rate is above the acceptable with this probability

% now let's assume for the moment that this neuron has a larger true refractory
% period than we will examine. So all spikes within the refractory period
% are contaminants. What we need to calculate is: what is the greatest
% number of violations we can see to have a sufficiently low probability
% that the level of contamination is greater than what's acceptable? 
% We want p( (trueContRate > acceptableCont) | obsContCount) < thresh. So
% we need to calculate the probability distribution of the true contaminant
% rate, given the observed contaminant count. 

% Under Bayes's Theorem and given no prior on the probability of true
% contaminant rates or on the probability of observed contamination counts,
% we just have p(trueContRate | obsContCount) = p(obsContCount |
% trueContRate). 

trueContRate = 0:0.01:firingRate; 
obsContCount = 0:100; 

refDur = 0.002; % s, refractory period duration
timeForViol = refDur * 2 * (firingRate-acceptableCont) * recDur; % total time available for violations 
probContGivenRate = zeros(numel(trueContRate), numel(obsContCount)); 
for t = 1:numel(trueContRate)
     probContGivenRate(t,:) = poisspdf(obsContCount,timeForViol*trueContRate(t));
end

figure; nsp = 3; subplot(1,nsp,1); 
imagesc(trueContRate, obsContCount, probContGivenRate'); 
h = addX(acceptableCont); h.Color = 'w';
xlabel('true cont rate'); ylabel('observed cont count');
colorbar;
title('p(count|rate)'); 
set(gca, 'YDir', 'normal')
axis square;
    
% a horizontal slice of this array gives us the probability of the true
% rate 
exampleObsContCount = 10; 
pObs = poisspdf(exampleObsContCount, trueContRate*timeForViol);
% pObs = pObs./sum(pObs)*100; 
subplot(1,nsp,2); plot(trueContRate, pObs, '.-')
addX(acceptableCont);
xlabel('true cont rate'); ylabel(sprintf('probability of seeing count=%d (%%)',exampleObsContCount)); 

% taking a uniform prior on true contaminant rate, we take this as a
% probability distribution, and can infer the probability that, for each
% count, the true rate was above acceptable. 
probUnaccept = zeros(numel(obsContCount),1);
for idx = 1:numel(obsContCount)
    pObs = poisspdf(obsContCount(idx), trueContRate*timeForViol);
    pObs = pObs./sum(pObs)*100;
    probUnaccept(idx) = sum(pObs(trueContRate>acceptableCont)); 
end
maxAcceptableCount = obsContCount(find(probUnaccept/100<thresh,1,'last')); 

subplot(1,nsp,3); 
plot(obsContCount, probUnaccept,'.-'); 

if isempty(maxAcceptableCount)
    maxAcceptableCount = -1; 
else
    addX(maxAcceptableCount); 
end
xlabel('obs cont count'); ylabel('probability that true rate is above acceptable'); 
title(sprintf('max acceptable count = %d', maxAcceptableCount));

% Now we can turn that into a function, to get the max acceptable count for
% a given set of parameters

m = maxAcceptableISIviol(firingRate, refDur, recDur, acceptableCont, thresh)

%% compute the max acceptable violations for a range of parameters

recDur = 3600; % s, recording duration
thresh = 0.1; % we can accept it if the true contamination rate is above the acceptable with this probability

fr = 0.3:0.3:6;
refDur = 0.0002:0.0003:0.006;

maxCounts = zeros(numel(fr), numel(refDur)); 
for fidx = 1:numel(fr)
    acceptableCont = 0.1*fr(fidx); % we can accept that 10% of spikes are contaminants
    for ridx = 1:numel(refDur)
        maxCounts(fidx,ridx) = maxAcceptableISIviol(fr(fidx), refDur(ridx), recDur, acceptableCont, thresh);
    end
end

figure; imagesc(refDur, fr, maxCounts); 
ylabel('firing rate'); xlabel('refractory period duration'); 
title('max acceptable count of violations')
set(gca, 'YDir', 'normal')
colorbar;
caxis([-1 10])
% hold on; 
% contour(refDur, fr, maxCounts, [0 0]); 
% contour(refDur, fr, maxCounts, [5 5]); 

% Finally, we can take a row out of that matrix corresponding to our actual
% firing rate. We can plot the calculated max acceptable violations as a
% fucntion of refractory duration, relative to the measured violations at
% each. If the true counts are below that max limit at any point, then
% we're good. The true counts will by definition start out above the line,
% because the line will say that even zero isn't good enough for a short
% enough refractory period. And they will end above the line because
% eventually you will have what looks like 100% contamination. The question
% is whether it is below the line for some range in between. 


%% Test neuron
% addpath(genpath(fullfile(githubDir, 'spikes')))
% addpath(genpath(fullfile(githubDir, 'npy-matlab')))
% sp = loadKSdir('D:\Hopkins\test\data.cortexlab.net\singlePhase3\data\Hopkins2016-07-22_metrics\ks2_out');

cid = 8;
st = sp.st(sp.clu==cid);

recDur = max(sp.st)-min(sp.st); 

fr = numel(st)/recDur;

acceptableCont = 0.1*fr; % we can accept that 10% of spikes are contaminants

refDur = 0.00005:0.0001:0.0022;

maxCounts = zeros(1, numel(refDur)); 
actualCounts = zeros(1, numel(refDur)); 
isi = sort(diff(st)); 
for ridx = 1:numel(refDur)
    maxCounts(1,ridx) = maxAcceptableISIviol(fr, refDur(ridx), recDur, acceptableCont, thresh);
    actualCounts(1,ridx) = sum(isi<=refDur(ridx)); 
end

figure; plot(refDur, maxCounts, 'o-', 'LineWidth',2.0);
hold on; plot(refDur, actualCounts, 'x-', 'LineWidth',2.0);
legend({'max acceptable count', 'actual count'});
xlabel('refractory period (s)'); 
ylabel('violation count');
set(gcf, 'Color', 'w'); box off; 

%% quick test of acg function

rate = 5; n = 10000;
st = genST(rate,n);

binSize = 0.0005;
b = binSize:binSize:0.005; 

[nACG,xACG] = histdiff(st, st, b);
% nACG = nACG./binSize; 

figure; 
stairs(xACG-binSize/2, nACG, 'LineWidth', 2.0);
box off; 
ylim([0 max(ylim())])

%% simulations


binSize = 0.0005;
b = 1e-6:binSize:0.0055;

baseRate = 15; 
recDur = 3600; 
rp = 0.004;

st = genST(rate,recDur);
isi = diff([0; st]); isi(isi<rp) = []; 
st = cumsum(isi); 

contRate = baseRate*0.02; 

contST = genST(contRate,recDur);

combST = sort([st; contST]); 

[nBase,xACG] = histdiff(st, st, b);
[nCont,xACG] = histdiff(contST, contST, b);
[nComb,xACG] = histdiff(combST, combST, b);

figure; subplot(1,2,1);hold on; 

stairs(xACG-binSize/2, nBase, 'LineWidth', 2.0);
stairs(xACG-binSize/2, nCont, 'LineWidth', 2.0);
stairs(xACG-binSize/2, nComb, 'LineWidth', 2.0);
box off; 
ylim([0 max(ylim())])
ylabel('acg counts'); 

subplot(1,2,2); hold on; 

stairs(xACG-binSize/2, cumsum(nBase), 'LineWidth', 2.0);
stairs(xACG-binSize/2, cumsum(nCont), 'LineWidth', 2.0);
stairs(xACG-binSize/2, cumsum(nComb), 'LineWidth', 2.0);

m = arrayfun(@(xx)maxAcceptableISIviol(rate, xx, recDur, rate/10, 0.1), b(2:end-1));
hold on; 
plot(b(2:end-1), m, 'o-', 'LineWidth', 2.0);

%% repeat the simulation above, but for many reps

nSim = 1000; 

usePct = [10 25 50 75 90]; 

binSize = 0.0005;
b = 1e-6:binSize:0.0055;

baseRate = 2; 
recDur = 3600; 
rp = 0.004;
contPct = 0.05;
contRate = baseRate*contPct; 


m = arrayfun(@(xx)maxAcceptableISIviol(baseRate, xx, recDur, baseRate/10, 0.1), b(2:end-1));



figure; hold on; box off;
plot(b(2:end-1), m, 'ko-', 'LineWidth', 2.0);

simRes = zeros(nSim, numel(b)-1); 
for n = 1:nSim

    st = genST(rate,recDur);
    isi = diff([0; st]); isi(isi<rp) = [];
    st = cumsum(isi);
    
    contST = genST(contRate,recDur);

    combST = sort([st; contST]); 

    [nComb,xACG] = histdiff(combST, combST, b);
    
    simRes(n,:) = cumsum(nComb); 
end

for bidx = 1:numel(b)-1
    boxp = prctile(simRes(:,bidx), usePct); 
    plot(b(bidx+1), boxp, 'bo-');
end

passPct = sum(any(simRes(:,1:end-1)<=repmat(m,nSim,1),2))/nSim*100;

title(sprintf('%.2f%% passed at contamination %.2f%%', passPct, contPct*100)); 

% I think the point here is that there is one thing that we've calculated,
% and another thing we're assessing. The thing we calculated is "what's the
% probability that the true mean was less than my acceptable contamination
% rate?" and the thing we're assessing is "what's the probability that
% contamination at the just-unacceptable rate gets rejected?" I think the
% difficulty is that the former question, while it's what we'd like to ask,
% requires some assumption about the prior probability of different
% contamination rates. So perhaps instead we should simply define the
% calculation in terms of the other question. 

%% now with another version of computing the max acceptable
clear all; 

nSim = 1000; 

usePct = [10 25 50 75 90]; 

binSize = 0.0005;
b = 1e-6:binSize:0.0055;

baseRate = 2; 
recDur = 3600; 
rp = 0.004;
contPct = 0.1;
contRate = baseRate*contPct; 

m = arrayfun(@(xx)maxAcceptableISIviol2(baseRate, xx, recDur, baseRate/10, 0.1), b(2:end-1));


figure; hold on; box off;
plot(b(2:end-1), m, 'ko-', 'LineWidth', 2.0);

simRes = zeros(nSim, numel(b)-1); 
for n = 1:nSim

    st = genST(baseRate,recDur);
    isi = diff([0; st]); isi(isi<rp) = [];
    st = cumsum(isi);
    
    contST = genST(contRate,recDur);

    combST = sort([st; contST]); 

    [nComb,xACG] = histdiff(combST, combST, b);
    
    simRes(n,:) = cumsum(nComb); 
end

for bidx = 1:numel(b)-1
    boxp = prctile(simRes(:,bidx), usePct); 
    plot(b(bidx+1), boxp, 'bo-');
end

passPct = sum(any(simRes(:,1:end-1)<=repmat(m,nSim,1),2))/nSim*100;

title(sprintf('%.2f%% passed at contamination %.2f%%', passPct, contPct*100)); 

% So, the method is correctly identifying the 10% point of the poisson
% distribution. That's good. Now, the threshold isn't working to reject 90%
% of neurons with borderline contamination, instead only rejecting about 
% 70% probably because different simulations pass below criterion at
% different bins. 
%
% One idea to sharpen this up is to require that it be below the curve for
% at least a few bins in a row, for instance. We can probably skip the
% 0.5ms bin, too. 
%
% But first let's just see what the empirical reject rate is for different
% levels of contamination

%% test reject rate for different levels of contamination

clear all; 

nSim = 1000; 

binSize = 0.0005;
b = 1e-6:binSize:0.0055;

baseRates = logspace(-0.3,1.3,8); 
recDur = 3600; 
rp = 0.004;
contPct = 0:0.02:0.2;
thresh = 0.2;

figure; 
passPct = zeros(numel(baseRates), numel(contPct));
colors = hsv(numel(baseRates)); 
for bidx = 1:numel(baseRates)
    baseRate = baseRates(bidx);
    for c = 1:numel(contPct)
        
        contRate = baseRate*contPct(c);
        
        m = arrayfun(@(xx)maxAcceptableISIviol2(baseRate, xx, recDur, baseRate/10, thresh), b(2:end-1));
        
        simRes = zeros(nSim, numel(b)-1);
        for n = 1:nSim
            
            st = genST(baseRate,recDur);
            isi = diff([0; st]); isi(isi<rp) = [];
            st = cumsum(isi);
            
            contST = genST(contRate,recDur);
            
            combST = sort([st; contST]);
            
            [nComb,xACG] = histdiff(combST, combST, b);
            
            simRes(n,:) = cumsum(nComb);
        end
        
        passPct(bidx,c) = sum(any(simRes(:,1:end-1)<=repmat(m,nSim,1),2))/nSim*100;
        
        
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
	