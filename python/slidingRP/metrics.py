# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 11:34:59 2022

@author: Noam Roth

compute the metric for a single cluster (neuron) in a recording
"""


from phylib.stats import correlograms
from types import SimpleNamespace
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import time



def slidingRP_all(spikeTimes, spikeClusters, **params):
    '''
    
    Compute the metric for each cluster in a recording
    
    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (s)
    spikeClusters : numpy.ndarray
        array of spike cluster ids that corresponds to spikeTimes.
    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize? 
        params.sampleRate : sample rate of the recording (Hz)
        params.recordingDuration: length of the recording (s)
        params.2msNoSpikes : set to True to automatically accept neurons above FRthresh with no spikes below 2ms
        params.FRthresh : set a FR threshold to reject neurons below this number regardless of metric performance (spks/s)

    Returns
    -------
    rpMetrics: dict
        keys:
            maxConfidenceAt10Cont
            minContWith90Confidence 
            timeOfLowestCont
            nSpikesBelow2
            confMatrix (optional, if returnMatrix ==1)
        note: minContWith90Confidence, timeOfLowestCont will return np.nan 
        for neurons with too few spikes -- these neurons have "empty" 
        confidence and should be rejected. 
    cont: nd.array
        Vector of contamination values tested      
    rp: nd.array
        Vector of refractory period durations tested  

    '''

    if params and 'returnMatrix' in params:
        returnMatrix = params['returnMatrix'] 
    else:
        returnMatrix = False
    
    
    if params and 'verbose' in params:
        verbose = params['verbose']; 
    else:
        verbose = False

    if params:
        print('params already instantiated')
    else:
        params = {}
        #default params
        params['sampleRate'] = 30000
        params['binSizeCorr'] = 1/30000
        params['2msNoSpikesCondition'] = False
        params['recordingDuration']  = max(spikeTimes) - min(spikeTimes) #if no recording duration listed, compute a rough estimate
        # based on the spike times

    if 'recordingDuration' not in params:
        params['recordingDuration']  = max(spikeTimes) - min(spikeTimes) #if no recording duration listed, compute a rough estimate
        # based on the spike times


    cids = np.unique(spikeClusters)
    
    #initialize rpMetrics as dict
    rpMetrics = {}
    rpMetrics['cidx'] = []
    rpMetrics ['maxConfidenceAtContaminationThresh'] = []
    rpMetrics['minContAtConfidenceThresh'] = []
    rpMetrics['timeOfLowestCont'] = []
    rpMetrics['nSpikesBelow2'] = []
    rpMetrics['firingRate'] = []
    
    if verbose:
        print("Computing metrics for %d clusters \n" % len(cids))
     

    for cidx in range(len(cids)):
        st = spikeTimes[spikeClusters==cids[cidx]] 
        
        [maxConfidenceAtContaminationThresh, minContAtConfidenceThresh, timeOfLowestCont,
            nSpikesBelow2, confMatrix, cont, rp, nACG,
            firingRate, secondsElapsed] = slidingRP(st, params)
    
        rpMetrics['cidx'].append(cids[cidx]) 
        rpMetrics['maxConfidenceAtContaminationThresh'].append(maxConfidenceAtContaminationThresh)
        rpMetrics['minContAtConfidenceThresh'].append(minContAtConfidenceThresh)
        rpMetrics['timeOfLowestCont'].append(timeOfLowestCont)
        rpMetrics['nSpikesBelow2'].append(nSpikesBelow2)
        rpMetrics['firingRate'].append(firingRate)
       
        
        if returnMatrix:
            if 'confMatrix' not in rpMetrics:
                rpMetrics['confMatrix'] = []
            rpMetrics['confMatrix'].append(confMatrix)
            
        
        if 'value' not in rpMetrics:
            rpMetrics['value'] = []
        # Set value (pass or reject neuron):
        if params['2msNoSpikesCondition']: #In this case, reject neurons below FRthresh and accept if no spikes below 2
            if firingRate < params['FRthresh']:
                rpMetrics['value'].append(0)
            else:
                if rpMetrics['nSpikesBelow2'] == 0:
                    rpMetrics['value'].append(1)
                else:
                    if minContAtConfidenceThresh <= params['contaminationThresh']:
                        rpMetrics['value'].append(1)
                    else:
                        rpMetrics['value'].append(0)

        else: #regular behavior of the metric, disregarding whether neurons have spikes below 2ms
            if minContAtConfidenceThresh<= params['contaminationThresh']:
                rpMetrics['value'].append(1)
            else:
                rpMetrics['value'].append(0)
            
        if verbose:
            if minContAtConfidenceThresh<=params['contaminationThresh']:
                pfstring = 'PASS'
            else: 
                pfstring = 'FAIL'
            print('  %d: %s max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms, n below 2 ms = %d' % (cids[cidx], pfstring, maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont*1000, nSpikesBelow2))

    return rpMetrics


def slidingRP(spikeTimes, params = None):

    '''     
    Compute the metric for one cluster
    
    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (s) for one cluster

    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize? 
        params.sampleRate : sample rate of the recording (Hz)
        params.recordingDuration: length of the recording (s) (if not supplied, this will be estimated based on spikeTimes)

    Returns
    -------


    
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

    if params is None:
        params = {}
        params['sampleRate'] = 30000
        params['binSizeCorr'] = 1 / params['sampleRate']
        params['returnMatrix'] = True
        params['verbose'] = True
        params['cidx'] = [0]
        params['recordingDuration']  = max(spikeTimes) - min(spikeTimes) #if no recording duration listed, compute a rough estimate
        # based on the spike times of this neuron

    if 'recordingDuration' not in params:
        params['recordingDuration']  = max(spikeTimes) - min(spikeTimes) #if no recording duration listed, compute a rough estimate
        # based on the spike times of this neuron

    seconds_start = time.time()
    [confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params)
    # matrix is [nCont x nRP]
    
    testTimes = rp>0.0005 # (in seconds) 
    #only test for refractory period durations greater than 0.5 ms
    
    maxConfidenceAtContaminationThresh = max(confMatrix[cont==params['contaminationThresh'],testTimes]) #TODO check behavior if no max
    
    
    indsConf = np.row_stack(np.where(confMatrix[:,testTimes]>params['confidenceThresh']))
    ii = indsConf[0] #row inds
    jj = indsConf[1] #col inds
    

    try:
        minI = np.min(ii)
        idx = np.argmin(ii)
        minContAtConfidenceThresh = cont[minI]
        minRP = np.argmax(confMatrix[minI,testTimes])


    except:    
        minContAtConfidenceThresh = np.nan
    
        minRP = np.nan

    try:
        timeOfLowestCont = rp[minRP+np.where(testTimes)[0][0]+1]
    except: 
        timeOfLowestCont = np.nan
        
    
    nSpikesBelow2 = sum(nACG[0:np.where(rp>0.002)[0][0]+1])

    secondsElapsed = time.time()-seconds_start
    return maxConfidenceAtContaminationThresh, minContAtConfidenceThresh, timeOfLowestCont, nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate, secondsElapsed
    
    
    
def computeMatrix(spikeTimes, params): 
    """
    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (ms)
    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)
        params.sampleRate : sample rate of the recording (Hz)

    Returns
    -------
    None.
    """
    
    
    cont = np.arange(0.5, 35, 0.5)  # vector of contamination values to test

    #set bin size for ACG
    if params and 'sampleRate' in params:
        sampleRate = params['sampleRate']
    else:
        sampleRate = 30000
    if params and 'binSizeCorr' in params:
        binSizeCorr = params['binSizeCorr']
    elif params and 'sampleRate' in params:
        binSizeCorr = 1 / sampleRate
    else:
        binSizeCorr = 1 / 30000
    if params and 'recordingDuration' in params:
        recDur = params['recordingDuration']
    else:
        recDur = max(spikeTimes) - min(spikeTimes) #if no recording duration listed, compute a rough estimate
        # based on the spike times of this neuron

    rpEdges = np.arange(0, 10/1000, binSizeCorr)  # in s
    rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2 # vector of refractory period durations to test
    
    #compute firing rate and spike count
    n_spikes = len(spikeTimes)
    #setup for acg
    clustersIds = [0]  # call the cluster id 0 (not used, but required input for correlograms)
    spikeClustersACG = np.zeros(n_spikes, dtype=np.int8)  # each spike time gets cluster id 0

    # compute an acg in 1s bins to compute the firing rate  
    nACG = correlograms(spikeTimes, spikeClustersACG, cluster_ids = clustersIds, bin_size = 1, sample_rate = params['sampleRate'], window_size=2,symmetrize=False)[0][0] #compute acg
    firingRate = nACG[1] / n_spikes
        
    nACG = correlograms(spikeTimes, spikeClustersACG, cluster_ids=clustersIds, bin_size=params['binSizeCorr'], sample_rate=params['sampleRate'], window_size=2, symmetrize=False)[0][0]  # compute acg

    confMatrix = 100 * computeViol(np.cumsum(nACG[0:rp.size])[np.newaxis, :], firingRate, n_spikes, rp[np.newaxis, :] + binSizeCorr / 2, cont[:, np.newaxis] / 100, recDur)

    return confMatrix, cont, rp, nACG, firingRate


def computeViol(obsViol, firingRate, spikeCount, refDur, contaminationProp, recDur):
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
        the contamination (proportion) i.e. 0.1 for 10% contamination

    Returns
    -------
    confidenceScore : float
        how confident we are that the cluster is less than contaminationProp 
        contaminated, given the observed violations and refractory period for 
        this cluster, under a Poisson assumption

    '''

    # the rate (spk/s) of the contamination level we are testing (contaminationProp)
    contaminationRate = firingRate * contaminationProp

    # the number of violations (spikes) we expect to see under this contamination rate
    #expectedViol = contaminationRate * refDur * 2 * spikeCount
    #as computed *not* in the same way as originally taken from Hill (see above commented out)
    N_t = firingRate * recDur #total number of spikes
    N_c = contaminationRate * recDur  #total number of contaminating spikes you would expect under the inputted CR
    N_b = N_t - N_c # the "base" number of spikes under the inputted CR

    expectedViol = 2 * refDur * 1/recDur * N_c * (N_b + (N_c - 1)/2) #number of expected violations, as defined in Llobet et al.

    # the confidence that this neuron is contaminated at a level less than contaminationProp, given the number of true
    # observed violations and under the assumption of Poisson firing
    confidenceScore = 1 - stats.poisson.cdf(obsViol, expectedViol)

    return confidenceScore


def plotSlidingRP(spikeTimes, params = None):
    '''
    

    Parameters
    ----------
    spikeTimes : numpy.ndarray
        array of spike times (ms)
    params : dict
        params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize? 
        params.sampleRate : sample rate of the recording (Hz)

    Returns
    -------
    None.

    '''
    if params is None: 
        clusterlabel = False
    else:
        clusterlabel = params['clusterLabel']

    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
        nSpikesBelow2, confMatrix, cont, rp, nACG, 
        firingRate,xx]  = slidingRP(spikeTimes, params)
    
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize = (12,4))
    
    ax = axs[0]

    ax.bar(rp*1000, nACG[0:len(rp)], width = np.diff(rp)[0]*1000, color = 'k', edgecolor = (1, 0, 0, 0)) #TODO width??
    ax.set_xlim([0, 5]) 
    ax.set_xlabel('Time from spike (ms)')
    ax.set_ylabel('ACG count (spks)')
    if clusterlabel:
        t1 = ('Cluster #%d: FR=%.2f' %(params['cidx'][0], firingRate))
    ax.fill(np.array([0, 1, 1, 0])*0.5, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'k',alpha= 0.2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
    ax = axs[1]
    c = ax.imshow(confMatrix, extent = [rp[0]*1000, rp[-1]*1000, cont[0], cont[-1]], aspect = 'auto', vmin = 0, vmax = 100, origin = 'lower')
    ax.set_xlim((0,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    cbar = fig.colorbar(c, ax = ax, location = 'right')
    cbar.set_label('Confidence (%)')
    ax.invert_yaxis()
    ax.plot([rp[0]*1000, rp[-1]*1000], [10, 10], 'r', linewidth = 1)
    
    if ~np.isnan(timeOfLowestCont):
        ax.plot(timeOfLowestCont*1000*np.array([1, 1]), [cont[0], cont[-1]],'r', linewidth = 1)
    
        # compute the conf=90 contour
        #zeropad confMatrix 
        z = np.zeros((np.shape(confMatrix)[0]+1,np.shape(confMatrix)[1]))
        z[1:,:] = confMatrix
        
        ii = np.argmax(z>90, 0).astype(float)  
        ii[ii==0] = np.nan
        contContour = np.empty(np.shape(ii)); contContour[:] = np.nan
        contContour[~np.isnan(ii)] = cont[(ii[~np.isnan(ii)]-1).astype(int)]
        ax.plot(rp*1000, contContour, 'r', linewidth = 2)
    val = ax.get_ylim()[1]
    # ax.fill(np.array([0, 1, 1, 0])*0.5, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'k',alpha= 1) 
    ax.fill(np.array([0, 1, 1, 0])*0.5, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'k',alpha= 0.2)
    # ax.add_patch(patches.Rectangle((0,0), 0.5, val, fc = 'k') )
    ax.set_xlabel('Time from spike (ms)')
    ax.set_xlim([0, 5]) 
    ax.set_ylabel('Contamination (%)')
    # ax.set_ylim([max(cont),0]) #seems unnecessary?
    t2 = ('max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms'% (maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont*1000))
    
    
    if minContWith90Confidence >= 10:
        print('FAIL')
        axs[0].set_title(t1, color='r')                   
        axs[1].set_title(t2,color = 'r')
    elif nSpikesBelow2 == 0:
        print('NOSPIKES2')
        axs[0].set_title(t1, color='b')                   
        axs[1].set_title(t2,color = 'b')
    else:
        print('PASS')
        axs[0].set_title(t1, color='g')                   
        axs[1].set_title(t2,color = 'g')
        
        
    
    ax = axs[2]
    ax.plot(rp*1000, np.squeeze(confMatrix[cont==10,:]), 'k', linewidth = 2.0)
    ax.set_xlabel('Time from spike (ms)')
    ax.set_ylabel('Confidence of ≤10% contamination (%)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.plot([0, 5], [90, 90], 'r'); 
    ax.fill(np.array([0, 1, 1, 0])*0.5, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'k',alpha= 0.2)
    ax.set_xlim([0, 5]); 
    ax.set_ylim([0, 100]); 
    
    fig.tight_layout()
    
    if params['savefig']:
        plt.savefig(params['figpath'], dpi=300)

def fitSigmoidACG_All(spikeTimes, spikeClusters, brainRegions, spikeAmps, rp, params):

    cids = np.unique(spikeClusters)
    
    #initialize rpMetrics as dict
    rpFit = {}
    rpFit['cidx'] = []
    rpFit['rpEstimate']  = []
    rpFit['brainRegion'] = []
    rpFit['amp'] = []
    rpFit['fr'] = []

    for cidx in range(len(cids)):
        st = spikeTimes[spikeClusters==cids[cidx]] 
        brainRegion = brainRegions[cidx]
        fr = len(st)/(st[-1]-st[0])

        sa = spikeAmps[spikeClusters==cids[cidx]] 
        amp = np.nanmedian(sa)*1000000

        #estimate rp
        minFR = 1 #spks/s
        minAmp = 50 #uV
        
        try:
            estimatedRP, estimateIdx, xSigmoid, ySigmoid = fitSigmoidACG(st, rp, fr, amp, minFR, minAmp, params)

            rpFit['cidx'].append(cids[cidx]) 
            rpFit['rpEstimate'].append(estimatedRP)
            rpFit['brainRegion'].append(brainRegion)
            rpFit['fr'].append(fr)
            rpFit['amp'].append(amp)
            if verbose:

                print('Estimated RP for cluster %d is %.2f ms'%(cids[cidx],estimatedRP))
        except:
            continue

    return rpFit


def fitSigmoidACG(st, timeBins, fr, amp, minFR = 1, minAmp = 50, params = None):
    '''
    

    Parameters
    ----------
    acg : np.array
        heights of acg bins (probability of spike for each timebin)
    timeBins : np.array
        timeBins used when computing acg
    params : TYPE
        DESCRIPTION.

    Returns
    -------
    estimatedRP:
        
    x: 
        
    y: 

    '''
    if params and 'verbose' in params:
        verbose = params['verbose']; 
    else:
        verbose = False
        
    if params and 'sampleRate' in params:
        sampleRate = params['sampleRate']
    else:
        sampleRate = 30000
    if params and 'binSizeCorr' in params:
        binSizeCorr = params['binSizeCorr']
    else:
        binSizeCorr = 1/30000
    
    
    
    #setup for acg
    clustersIds = [0] # call the cluster id 0 (not used, but required input for correlograms)
    spikeClustersACG = np.zeros(len(st), dtype = 'int8') # each spike time gets cluster id 0 
    #compute acg
    acg = correlograms(st, spikeClustersACG, cluster_ids = clustersIds, bin_size = binSizeCorr, sample_rate = sampleRate, window_size=2,symmetrize=False)[0][0] #compute acg

    if len(acg) > len(timeBins):
        acg = acg[0:len(timeBins)]

    if 'numSpikesThresh' not in params.keys():
        numSpikesThresh = 20  # need at least this many total spikes to compute fit
    else:
        numSpikesThresh = params['numSpikesThresh']
    
    peakDistFromEndBin = 5 #todo: params!! 
        
    estimatedRP = np.nan #initialize as a nan for cases where it doesn't work
    if fr >= minFR and amp >= minAmp:    #(sum(acg)>numSpikesThresh):   
        #potential todo: insert a case here for if the acg is symmetric?
        # p0 = [np.mean(ydata),2,1,min(ydata)] #starting point for fit? 
        minSigmoid = np.mean(acg[timeBins<0.0005]) #first 0.5 ms of data #todo make this parameter
        peakIdx = np.argmax(acg)        
        peakVal = np.max(acg)
        
        #2ms around peak
        timeValuesMin = np.where(timeBins >= timeBins[peakIdx]-0.001)[0][0]
        timeValuesMax = np.where(timeBins <= timeBins[peakIdx]+0.001)[0][-1]
        maxSigmoid = np.mean(acg[timeValuesMin:(timeValuesMax+1)])
    
        #if the peak is well before the end of the acg, only fit data up to the peak + n bins
        if peakIdx < len(acg)-peakDistFromEndBin:
            acg = acg[0:peakIdx + peakDistFromEndBin]
            timeBins = timeBins[0:peakIdx + peakDistFromEndBin]
    
        #fit the sigmoid with max and min fixed
        try:
            popt, pcov = curve_fit(lambda x, x0, k: sigmoid(x, maxSigmoid, x0, k, minSigmoid ), timeBins, acg)      
            fitParams = [maxSigmoid, popt[0], popt[1], minSigmoid]
    
            xSigmoid = timeBins#np.linspace(0, timeBins[-1], timeBins[-]) #evenly spaced vector for plotting sigmoid
            ySigmoid = sigmoid(xSigmoid, *fitParams)
        
                #find RP            
            RPEstimateFromPercentageOfSlope = 0.10
            estimateIdx, _ = closest(ySigmoid, RPEstimateFromPercentageOfSlope*(maxSigmoid - minSigmoid) + minSigmoid) 
            estimatedRP = 1000* xSigmoid[estimateIdx] # in ms
    
        except:
            estimatedRP = np.nan
            estimateIdx = np.nan
            xSigmoid = np.nan
            ySigmoid = np.nan
    
                #print('fit error')

        

    else:
        #print('fit error')
        estimatedRP = np.nan
        estimateIdx = np.nan
        xSigmoid = np.nan
        ySigmoid = np.nan
        
    return estimatedRP, estimateIdx, xSigmoid, ySigmoid
    


def plotSigmoid(ax, acg, timeBins, ySigmoid, estimatedIdx, estimatedRP):
    if len(acg) > len(timeBins):
        acg = acg[0:len(timeBins)]
    ax.bar(timeBins*1000, acg, width = np.diff(timeBins*1000)[0],alpha = 0.5)
    ax.set_xlim(0,5)

    ax.plot(timeBins*1000, ySigmoid,'k')
    ax.plot(timeBins[estimatedIdx]*1000, ySigmoid[estimatedIdx],'rx')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Estimated RP:%.2f ms'%estimatedRP)
    ax.set_ylabel('Number of spikes')
    ax.set_xlabel('Time (s)')
    return acg
    
    
    
    
    
    #helper functions    
def find_nearest(array, value):
    array = np.asarray(array)
    subtracted = (array - value)
    valid_idx = np.where(subtracted >= 0)[0]
    if len(valid_idx)>0:

        out = valid_idx[subtracted[valid_idx].argmin()]

    else:
        out = np.nan
    
    return out


def sigmoid(x, L, x0, k, b): # L is max value; x0 is the midpoint x value; k is the steepness, b is the baseline shift
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return y


def closest(lst, K):
      
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return idx, lst[idx]