# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 11:34:59 2022

@author: Noam Roth

compute the metric for a single cluster (neuron) in a recording
"""

def slidingRP(spikeTimes, params)

''' 
inputs:
    spikeTimes: spike times vector for a single unit (in ms) 
    params: TODO
    
returns: 

maxConfidenceAt10Cont:   Max confidence that you have <= 10% contamination
minContWith90Confidence: Minimum contamination for which you have >=90% confidence
timeOfLowestCont:        Time at which best score happens
nSpikesBelow2:           Number of observed spikes that occur before 2 ms
confMatrix:              Full confidence matrix of size nCont x nRP
cont:Vector of contamination values tested 
rp: Vector of refractory period durations tested  
nACG: the autocorrelogram of the neuron
firingRate: firing rate of the cluster, computed as the average acg value from 1-2 seconds
'''

[confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params)
# matrix is [nCont x nRP]

testTimes = rp>0.0005 # (in seconds) 
#only test for refractory period durations greater than 0.5 ms

maxConfidenceAt10Cont = max(confMatrix(cont==10, testTimes)); 

[ii,jj] = find(confMatrix(:,testTimes)>90); 
[minI, idx] = min(ii); 
minContWith90Confidence = cont(minI);  
if isempty(minContWith90Confidence); minContWith90Confidence = NaN; end

[~,minRP] = max(confMatrix(minI,testTimes)); 
timeOfLowestCont = rp(minRP+find(testTimes,1));
if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end

nSpikesBelow2 = sum(nACG(1:find(rp>0.002,1)));


return maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
    nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate
    
    
    
def computeMatrix(spikeTimes, params): 
    
    return confMatrix, cont, rp, nACG, firingRate
