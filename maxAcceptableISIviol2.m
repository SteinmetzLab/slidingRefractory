

function m = maxAcceptableISIviol2(firingRate, refDur, recDur, acceptableCont, thresh)


% factor of 2 necessary here??
% timeForViol = refDur * 2 * (firingRate-acceptableCont) * recDur; % total time available for violations

timeForViol = refDur * 2 * (firingRate) * recDur; % total time available for violations


expectedCountForAcceptableLimit = timeForViol*acceptableCont;

m = poissinv(thresh, expectedCountForAcceptableLimit);

if m==0 && poisspdf(0, expectedCountForAcceptableLimit)>thresh
    m=-1; % meaning, even a count of zero is not sufficient evidence here
end

