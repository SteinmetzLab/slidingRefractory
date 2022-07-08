

function confidenceScore = computeViol(...
    obsViol, firingRate, spikeCount, refDur, contaminationProp)

% refDur is the refractory period duration in seconds
% contaminationProp is the query contamination level as a proportion,
% i.e. 0.1 for 10% contamination

contaminationRate = firingRate*contaminationProp; 
expectedViol = contaminationRate*refDur*2*spikeCount; 

confidenceScore = 1-poisscdf(obsViol, expectedViol); 

