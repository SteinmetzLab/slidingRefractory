

function m = maxAcceptableISIviol(firingRate, refDur, recDur, acceptableCont, thresh)


% trueContRate = 0:0.01:firingRate; % sp/s
trueContRate = linspace(0,firingRate, 100); 

% factor of 2 necessary here??
timeForViol = refDur * 2 * (firingRate-acceptableCont) * recDur; % total time available for violations

expectedCountAtFullCont = firingRate*timeForViol; 

obsContCount = 0:expectedCountAtFullCont; 

% taking a uniform prior on true contaminant rate, we take this as a
% probability distribution, and can infer the probability that, for each
% count, the true rate was above acceptable. 
probUnaccept = zeros(numel(obsContCount),1);
for idx = 1:numel(obsContCount)
    pObs = poisspdf(obsContCount(idx), trueContRate*timeForViol);
    pObs = pObs./sum(pObs)*100;
    probUnaccept(idx) = sum(pObs(trueContRate>acceptableCont)); 
end
m = obsContCount(find(probUnaccept/100<thresh,1,'last')); 

if isempty(m) % no count is below threshold
    m = -1; % meaning, even a count of zero is not sufficient evidence here
end
