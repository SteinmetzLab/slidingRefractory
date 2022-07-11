# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 11:34:59 2022

@author: Noam Roth

compute the metric for a single cluster (neuron) in a recording
"""


from phylib.stats import correlograms
from types import SimpleNamespace


def slidingRP(spikeTimes, params):

    '''     
    
        Parameters
    ----------
    spikeTimes : TYPE
        DESCRIPTION.
    params : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    
    
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
    
    # maxConfidenceAt10Cont = max(confMatrix(cont==10, testTimes)); 
    
    # [ii,jj] = find(confMatrix(:,testTimes)>90); 
    # [minI, idx] = min(ii); 
    # minContWith90Confidence = cont(minI);  
    # if isempty(minContWith90Confidence); minContWith90Confidence = NaN; end
    
    # [~,minRP] = max(confMatrix(minI,testTimes)); 
    # timeOfLowestCont = rp(minRP+find(testTimes,1));
    # if isempty(timeOfLowestCont); timeOfLowestCont = NaN; end
    
    # nSpikesBelow2 = sum(nACG(1:find(rp>0.002,1)));


    # return maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,...
    # nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate
    
    
    
def computeMatrix(spikeTimes, params): 
    '''
    

    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (ms)
    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize? 
        params.sampleRate : sample rate of the recording (Hz)
    DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    
    cont = np.arange(0.5, 35, 0.5)  # vector of contamination values to test
    rpBinSize = 1 / 30000  
    rpEdges = np.arange(0, 10/1000, rpBinSize) # in ms  
    rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
    
    #compute firing rate and spike count
    spikeCount = len(spikeTimes)
    #setup for acg
    clustersIds = [0] # call the cluster id 0 (not used, but required input for correlograms)
    spikeClusters = np.zeros(len(ts), dtype = 'int8') # each spike time gets cluster id 0 
    
    # compute an acg in 1s bins to compute the firing rate 
    nACG = correlograms(spikeTimes, spikeClusters, cluster_ids = clustersIds, bin_size = 1, sample_rate = params.sampleRate, window_size=2,symmetrize=False)[0][0] #compute acg
    firingRate = nACG[1] / spikeCount
        
    nACG = correlograms(spikeTimes, spikeClusters, cluster_ids = clustersIds, bin_size = params.binSizeCorr, sample_rate = params.sampleRate, window_size=2,symmetrize=False)[0][0] #compute acg

    confMatrix = np.empty((len(cont), len(rp)))
    confMatrix[:] = np.nan
    for rpIdx in range(len(rp)):
        
        # compute observed violations
        obsViol = sum(nACG[0:rpIdx+1]) #TODO this is off slightly (half-bin) from matlab...

        for cidx in range(len(cont)):

            confMatrix[cidx, rpIdx] = 100*computeViol(obsViol, firingRate, 
                          spikeCount, rp[rpIdx]+rpBinSize/2, cont[cidx]/100) #TODO FIX RP BIN
       
            

    return confMatrix, cont, rp, nACG, firingRate


def computeViol(obsViol, firingRate, spikeCount, refDur, contaminationProp):
    '''
    

    Parameters
    ----------
    obsViol : int
        the number of spikes observed within the refractory period duration.
    firingRate : float
        firing rate of the cluster (estimated from ACG) in spks/s
    spikeCount : int
        total spike count of cluster
    refDur : float
        refractory period duration in seconds
    contaminationProp : float
        the allowed contamination (proportion) i.e. 0.1 for 10% contamination

    Returns
    -------
    confidenceScore : float
        how confident we are that the cluster is less than contaminationProp 
        contaminated, given the observed violations and refractory period for 
        this cluster, under a Poisson assumption

    '''

    contaminationRate = firingRate * contaminationProp 
    expectedViol = contaminationRate* refDur * 2 * spikeCount 

    confidenceScore =  1 - stats.poisson.cdf(obsViol, expectedViol)
    
    return confidenceScore


#%% script testing

params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params = SimpleNamespace(**params) #convert to dot notation

mat = computeMatrix(spikeTimes, params)