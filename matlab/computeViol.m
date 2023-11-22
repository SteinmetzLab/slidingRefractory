

function [confidenceScore, expectedViol] = computeViol(...
    obsViol, firingRate, spikeCount, refDur, contaminationProp, recDur)

% refDur is the refractory period duration in seconds
% contaminationProp is the query contamination level as a proportion,
% i.e. 0.1 for 10% contamination

% contaminationRate = firingRate*contaminationProp; 
% expectedViol = contaminationRate*refDur*2*spikeCount; 

% Updated calculation from Llobet et al 2022
Nc = spikeCount * contaminationProp;
Nb = spikeCount * (1-contaminationProp); 
expectedViol = 2*refDur/recDur * Nc .* (Nb + (Nc-1)/2);

confidenceScore = 1-poisscdf(obsViol, expectedViol); 

