# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 11:34:59 2022

@author: Gaelle Chapuis ; Legacy code from Noam Roth commented below

compute the metric for a single cluster (neuron) in a recording
"""
import warnings
from scipy.optimize import OptimizeWarning
from phylib.stats import correlograms
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import scipy


def closest(lst, K):
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return idx, lst[idx]


def compute_timebins(acg, bin_size_secs):
    x = np.arange(len(acg))
    timeBins = x.dot(bin_size_secs)
    return timeBins


def sigmoid(x, L, x0, k, b):  # L is max value; x0 is the midpoint x value; k is the steepness, b is the baseline shift
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


def compute_rf(acg,
               min_sig=np.array([0.0, 0.0005]),  # 0-0.5 ms
               t_medfilter=0.00083,
               bin_size_secs=1 / 30_000,
               timeBins=None,
               RPEstimateFromPercentageOfSlope = 0.05,
               fr_percentage=10 / 100):
    '''
    Compute the refractory period (RP) from an auto-correlogram (ACG) array, by :
    - filtering the ACG (median filter)
    - finding the ACG first local peak above a certain firing rate threshold
    - fitting a sigmoid to the first portion of the filtered ACG
    - The RP is defined as the time at which a certain percentage of the sigmoid fit is reached.
    :param acg: the autocorrelogram, 1D numpy array containing spike number in bins of size bin_size_secs
    :param min_sig: time window in second to compute the minimum of the sigmoid fit
    :param t_medfilter: window for the median filter, in seconds
    :param bin_size_secs: size of the acg bins, in seconds
    :param timeBins: time of the bins, re-computed if None
    :param RPEstimateFromPercentageOfSlope: portion of the sigmoid fit at which to take the time as RP; value from 0-1
    :param fr_percentage: percentage of the max-min firing rate, used to find the first peak.
    :return:
    estimatedRP: the estimated RP in milliseconds
    estimateIdx: the index of the RP (which corresponds to the bin index of the ACG)
    xSigmoid, ySigmoid: the sigmoid fit (x-values: times, y-values: curve values)

    '''
    # This is a hack: we do not want to return the fit on an ACG where the fit is poor
    warnings.simplefilter("error", OptimizeWarning)

    if timeBins is None:  # Compute timebins
        timeBins = compute_timebins(acg, bin_size_secs)

    # Median filter
    med_filt = scipy.ndimage.median_filter(acg, size=int(np.round(t_medfilter / bin_size_secs)))
    # Find all the peaks on this filtered trace
    peaks_idx = scipy.signal.find_peaks(med_filt)[0]

    # Compute value (max) in short time window near 0
    minSigmoid = np.max(med_filt[(timeBins > min_sig[0]) & (timeBins < min_sig[1])])
    # Note: could be using either the max() or mean()
    # max forces the peak to be outside of this time window

    # To find the peak, make sure it is strictly higher than this value
    # and above a baseline percentage firing rate
    fr_baseline = fr_percentage * (med_filt.max() - med_filt.min())
    peaks_possible = np.where((med_filt[peaks_idx] > minSigmoid) &
                              (med_filt[peaks_idx] > fr_baseline))[0]

    # if no peak is found, abort and return NaN
    if len(peaks_possible) == 0:
        estimatedRP = np.nan
        estimateIdx = np.nan
        xSigmoid = np.array([])
        ySigmoid = np.array([])
        return estimatedRP, estimateIdx, xSigmoid, ySigmoid

    # Use first peak possible found (i.e. closest to 0 second on the ACG)
    peak_idx = peaks_idx[peaks_possible[0]]
    maxSigmoid = med_filt[peak_idx]

    # Truncate ACG and time bins according to max
    timeBins_fit = timeBins[0:peak_idx]
    acg_fit = med_filt[0:peak_idx]

    # fit the sigmoid with max and min fixed
    try:
        popt, pcov = curve_fit(lambda x, x0, k: sigmoid(x, maxSigmoid, x0, k, minSigmoid), timeBins_fit, acg_fit)
        fitParams = [maxSigmoid, popt[0], popt[1], minSigmoid]

        xSigmoid = timeBins
        ySigmoid = sigmoid(xSigmoid, *fitParams)

        # find RP
        estimateIdx, _ = closest(ySigmoid, RPEstimateFromPercentageOfSlope * (maxSigmoid - minSigmoid) + minSigmoid)

        # Compute the index of the first ACG bin with non-null firing rate
        first_index = np.where(acg != 0)[0][0]
        if estimateIdx < first_index:  # If the estimate RP is BEFORE the first ACG bin
            estimateIdx = first_index  # Replace

        estimatedRP = 1000 * xSigmoid[estimateIdx]  # in ms
        if np.max(ySigmoid) - np.min(ySigmoid) < 1:  # The fit is essentially flat
            raise OptimizeWarning
    except (OptimizeWarning, RuntimeWarning, RuntimeError):  # This is in the case the bins are too few
        # print('fit error')
        estimatedRP = np.nan
        estimateIdx = np.nan
        xSigmoid = np.array([])
        ySigmoid = np.array([])
    return estimatedRP, estimateIdx, xSigmoid, ySigmoid


def remove_lowrp_confmat(confMatrix, rp, rp_reject=0.0005):
    # We want to compute on the matrix only for RPs above a certain value
    # Remove those small RP values from the rp vector and conf matrix
    rp_idx_keep = rp > rp_reject
    rp = rp[rp_idx_keep]
    confMatrix = confMatrix[:, rp_idx_keep]
    return confMatrix, rp


def confidence_contamin(confMatrix, cont, rp, cont_thresh=10.0, rp_reject = 0.0005):
    '''
    For a level of contamination contamin_level given (default 10%), find the smallest confidence
    value for which the minimum value of the contamination curve is equal or lower to the
    contamin_level.
    In practice, this is equivalent to finding the maximum value of the confidence of
    the confidence matrix at the contamination level row.

    Uses the output of the function computeMatrix()

    :param confMatrix: the confidence matrix (contamination x RP, values: confidence, ranging from 0-1)
    :param cont: contamination vector at which the confidence is computed (ranges by default from 0-35)
    :param rp: refractory period vector at which the confidence is computed
    :param cont_level: level of contamination searched for, default is 10% (0.1)
    :return:
    '''
    # Find index in cont vector where there is cont_thresh or closest (higher) value
    idx_cont = np.where(cont >= cont_thresh)[0][0]
    cont_thresh = cont[idx_cont]  # Return actual level of contamination studied

    # We want to compute the curve of contamination only for RPs above a certain value
    # Remove those small RP values from the rp vector and conf matrix
    confMatrix, _ = remove_lowrp_confmat(confMatrix, rp, rp_reject=rp_reject)

    # At the contamination level studied, find the maximal value of confidence
    max_conf = np.max(confMatrix[idx_cont, :])  # Legacy name: 'maxConfidenceAt10Cont'
    return max_conf, idx_cont, cont_thresh


def pass_slidingRP_confmat(confMatrix, cont, rp, conf_thresh=90, cont_thresh=10, rp_reject=0.0005):
    '''
    Given a confidence matrix, a confidence threshold (default 90%) and a contamination threshold (default=10),
    assess whether the unit passes the sliding RP metric

    Uses the output of the function computeMatrix()

    :param confMatrix:
    :param cont:
    :param rp:
    :param conf_thresh:
    :param cont_thresh:
    :param rp_reject:
    :return:
    '''
    # We want to compute the curve of contamination only for RPs above a certain value
    # Remove those small RP values from the rp vector and conf matrix
    confMatrix, rp = remove_lowrp_confmat(confMatrix, rp, rp_reject=rp_reject)

    # Find matrix indices that are above or equal to the confidence threshold
    a = np.where(confMatrix >= conf_thresh)
    if len(a[0]) > 0:
        # Find minimum contamination value for this conf threshold
        min_idx = np.min(a[0])  # Min on rows axis = contamination axis
        min_cont = cont[min_idx]  # Legacy name: minContWith90Confidence
        # Find the smallest RP possible at the contamination level at the max confidence val
        minRP = np.argmax(confMatrix[min_idx, :])
        rp_min_val = rp[minRP]  # Legacy name: timeOfLowestCont
        # TODO this code could be deleted as rp_mini_val is not informative

        # Check if this unit passes the sliding RP metric
        pass_cont_thresh = min_cont <= cont_thresh
    else:
        pass_cont_thresh = False
        min_cont = np.nan
        rp_min_val = np.nan

    return pass_cont_thresh, min_cont, rp_min_val  # Legacy: value, minContWith90Confidence, timeOfLowestCont


def slidingRP(spikeTimes, conf_thresh=90, cont_thresh=10, rp_reject=0.0005,
              params=None):
    if params is None:
        params = {}
        params['sampleRate'] = 30000
        params['binSizeCorr'] = 1 / params['sampleRate']
        params['returnMatrix'] = True
        params['verbose'] = True
        params['cidx'] = [0]

    [confMatrix, cont, rp, nACG, firing_rate] = computeMatrix(spikeTimes, params)
    # Legacy: PASS, minContWith90Confidence, timeOfLowestCont
    pass_cont_thresh, min_cont, rp_min_val = \
        pass_slidingRP_confmat(confMatrix, cont, rp, conf_thresh, cont_thresh, rp_reject)

    # Legacy 'maxConfidenceAt10Cont'
    max_conf, _, _ = \
        confidence_contamin(confMatrix, cont, rp, cont_thresh, rp_reject)

    # Legacy
    n_spikes_below2 = sum(nACG[0:np.where(rp > 0.002)[0][0] + 1])

    # We apply this for IBL data, that has a duration of 1h on average hence the hardcoded FR>0.5:
    if (n_spikes_below2 == 0) and (firing_rate > 0.5) and (pass_cont_thresh is False):
        pass_forced = True
    else:
        pass_forced = False

    # Legacy:
    # return maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont, nSpikesBelow2,
    # confMatrix, cont, rp, nACG, firingRate, secondsElapsed
    return max_conf, min_cont, rp_min_val, \
        n_spikes_below2, firing_rate,\
        pass_cont_thresh, pass_forced


def slidingRP_all(spikeTimes, spikeClusters,
                  conf_thresh=90, cont_thresh=10, rp_reject=0.0005,
                  params=None):
    """
    :param spikeTimes:  array of spike times (s)
    :param spikeClusters:  array of spike cluster ids that corresponds to spikeTimes
    :param params:
    :return: dictionary of values
    """

    cids = np.unique(spikeClusters)

    # initialize rpMetrics as dict
    rpMetrics = {}
    rpMetrics['cidx'] = []
    rpMetrics['max_confidence'] = []
    rpMetrics['min_contamination'] = []
    rpMetrics['rp_min_val'] = []
    rpMetrics['n_spikes_below2'] = []
    rpMetrics['firing_rate'] = []
    rpMetrics['value'] = []
    rpMetrics['value_forced'] = []

    # Loop over clusters
    for cidx in range(len(cids)):
        st = spikeTimes[spikeClusters == cids[cidx]]

        [max_confidence, min_contamination, rp_min_val,
         n_spikes_below2, firing_rate,
         pass_cont_thresh, pass_forced] = slidingRP(st,
                                                    conf_thresh=conf_thresh, cont_thresh=cont_thresh,
                                                    rp_reject=rp_reject, params=params)

        rpMetrics['cidx'].append(cids[cidx])
        rpMetrics['max_confidence'].append(max_confidence)
        rpMetrics['min_contamination'].append(min_contamination)
        rpMetrics['rp_min_val'].append(rp_min_val)
        rpMetrics['n_spikes_below2'].append(n_spikes_below2)
        rpMetrics['firing_rate'].append(firing_rate)
        rpMetrics['value'].append(int(pass_cont_thresh))
        rpMetrics['value_forced'].append(int(pass_forced))

    return rpMetrics



## Code from OW


def computeMatrix(spikeTimes, params):
    """
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
    """

    cont = np.arange(0.5, 35, 0.5)  # vector of contamination values to test
    rpBinSize = 1 / 30000
    rpEdges = np.arange(0, 10 / 1000, rpBinSize)  # in s
    rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2  # vector of refractory period durations to test

    # compute firing rate and spike count
    n_spikes = spikeTimes.size
    # setup for acg
    clustersIds = [0]  # call the cluster id 0 (not used, but required input for correlograms)
    spikeClustersACG = np.zeros(n_spikes, dtype=np.int8)  # each spike time gets cluster id 0

    # compute an acg in 1s bins to compute the firing rate
    nACG = \
    correlograms(spikeTimes, spikeClustersACG, cluster_ids=clustersIds, bin_size=1, sample_rate=params['sampleRate'],
                 window_size=2, symmetrize=False)[0][0]  # compute acg
    firingRate = nACG[1] / n_spikes

    nACG = correlograms(spikeTimes, spikeClustersACG, cluster_ids=clustersIds, bin_size=params['binSizeCorr'],
                        sample_rate=params['sampleRate'], window_size=2, symmetrize=False)[0][0]  # compute acg

    # confMatrix = np.zeros((cont.size, rp.size)) * np.nan
    # for ir in np.arange(confMatrix.shape[1]):
    #     # compute observed violations
    #     obsViol = np.sum(nACG[0:ir + 1])  # TODO this is off slightly (half-bin) from matlab...
    #     for cidx in np.arange(confMatrix.shape[0]):
    #         confMatrix[cidx, ir] = 100 * computeViol(obsViol, firingRate, n_spikes, rp[ir] + rpBinSize / 2, cont[cidx] / 100)  # TODO FIX RP BIN

    confMatrix = 100 * computeViol(np.cumsum(nACG[0:rp.size])[np.newaxis, :], firingRate, n_spikes,
                                   rp[np.newaxis, :] + rpBinSize / 2, cont[:, np.newaxis] / 100)

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
    expectedViol = contaminationRate * refDur * 2 * spikeCount
    confidenceScore = 1 - stats.poisson.cdf(obsViol, expectedViol)

    return confidenceScore

# --------
##

#
# '''
# The below is legacy code, written prior to 22-Feb-2024 by N. Steinmetz and N. Roth
# '''
#
# def fit_sigmoid(acg, timeBins, min_sig=[0.0004, 0.0008], peakDistFromEndBin=5):
#     # remove first ms from the computation
#     # idx_time = (timeBins > min_sig[0])
#     # acg = acg[idx_time]
#     # timeBins = timeBins[idx_time]
#
#     # Force first ms to 0 in acg
#     # idx_time = (timeBins < min_sig[0])
#     # acg[idx_time] = 0
#
#     #Compute
#     minSigmoid = np.mean(acg[(timeBins > min_sig[0]) & (timeBins < min_sig[1])])  # first 0.4-0.8 ms of data
#     peakIdx = np.argmax(acg)
#     peakVal = np.max(acg)
#
#     # 2ms around peak
#     timeValuesMin = np.where(timeBins >= timeBins[peakIdx] - 0.001)[0][0]
#     timeValuesMax = np.where(timeBins <= timeBins[peakIdx] + 0.001)[0][-1]
#     maxSigmoid = np.mean(acg[timeValuesMin:(timeValuesMax + 1)])
#
#     # if the peak is well before the end of the acg, only fit data up to the peak + n bins
#     if peakIdx < len(acg) - peakDistFromEndBin:
#         acg = acg[0:peakIdx + peakDistFromEndBin]
#         timeBins = timeBins[0:peakIdx + peakDistFromEndBin]
#
#     # fit the sigmoid with max and min fixed
#     try:
#         popt, pcov = curve_fit(lambda x, x0, k: sigmoid(x, maxSigmoid, x0, k, minSigmoid), timeBins, acg)
#         fitParams = [maxSigmoid, popt[0], popt[1], minSigmoid]
#
#         xSigmoid = timeBins
#         ySigmoid = sigmoid(xSigmoid, *fitParams)
#
#         # find RP
#         RPEstimateFromPercentageOfSlope = 0.10
#         estimateIdx, _ = closest(ySigmoid, RPEstimateFromPercentageOfSlope * (maxSigmoid - minSigmoid) + minSigmoid)
#         estimatedRP = 1000 * xSigmoid[estimateIdx]  # in ms
#     except:
#         # print('fit error')
#         estimatedRP = np.nan
#         estimateIdx = np.nan
#         xSigmoid = np.nan
#         ySigmoid = np.nan
#     return estimatedRP, estimateIdx, xSigmoid, ySigmoid
#
#
# def slidingRP_all(spikeTimes, spikeClusters, params=None):
#     '''
#
#     Compute the metric for each cluster in a recording
#
#     Parameters
#     ----------
#     spikeTimes : numpy.ndarray
#         array of spike times (ms)
#     spikeClusters : numpy.ndarray
#         array of spike cluster ids that corresponds to spikeTimes.
#     params : dict
#         params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize?
#         params.sampleRate : sample rate of the recording (Hz)
#
#     Returns
#     -------
#     rpMetrics: dict
#         keys:
#             maxConfidenceAt10Cont
#             minContWith90Confidence
#             timeOfLowestCont
#             nSpikesBelow2
#             confMatrix (optional, if returnMatrix ==1)
#         note: minContWith90Confidence, timeOfLowestCont will return np.nan
#         for neurons with too few spikes -- these neurons have "empty"
#         confidence and should be rejected.
#     cont: nd.array
#         Vector of contamination values tested
#     rp: nd.array
#         Vector of refractory period durations tested
#
#     '''
#
#     if params and 'returnMatrix' in params:
#         returnMatrix = params['returnMatrix']
#     else:
#         returnMatrix = False
#
#     if params and 'verbose' in params:
#         verbose = params['verbose'];
#     else:
#         verbose = False
#
#     cids = np.unique(spikeClusters)
#
#     # initialize rpMetrics as dict
#     rpMetrics = {}
#     rpMetrics['cidx'] = []
#     rpMetrics['maxConfidenceAt10Cont'] = []
#     rpMetrics['minContWith90Confidence'] = []
#     rpMetrics['timeOfLowestCont'] = []
#     rpMetrics['nSpikesBelow2'] = []
#
#     if verbose:
#         print("Computing metrics for %d clusters \n" % len(cids))
#
#     for cidx in range(len(cids)):
#         st = spikeTimes[spikeClusters == cids[cidx]]
#
#         [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
#          nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate,
#          pass_cont_thresh, pass_forced] = slidingRP_2(st, params=params)
#
#         '''
#         [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
#          nSpikesBelow2, confMatrix, cont, rp, nACG,
#          firingRate, secondsElapsed] = slidingRP(st, params=params)
#         '''
#
#         rpMetrics['cidx'].append(cids[cidx])
#         rpMetrics['maxConfidenceAt10Cont'].append(maxConfidenceAt10Cont)
#         rpMetrics['minContWith90Confidence'].append(minContWith90Confidence)
#         rpMetrics['timeOfLowestCont'].append(timeOfLowestCont)
#         rpMetrics['nSpikesBelow2'].append(nSpikesBelow2)
#
#         if returnMatrix:
#             if 'confMatrix' not in rpMetrics:
#                 rpMetrics['confMatrix'] = []
#             rpMetrics['confMatrix'].append(confMatrix)
#
#         if 'value' not in rpMetrics:
#             rpMetrics['value'] = []
#         if minContWith90Confidence <= 10:
#             rpMetrics['value'].append(1)
#         else:
#             rpMetrics['value'].append(0)
#
#         if verbose:
#             if minContWith90Confidence <= 10:
#                 pfstring = 'PASS'
#             else:
#                 pfstring = 'FAIL'
#             print('  %d: %s max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms, n below 2 ms = %d' % (
#             cids[cidx], pfstring, maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont * 1000,
#             nSpikesBelow2))
#
#     return rpMetrics
#
#
#
# ## ------ LEGACY CODE BELOW BY NICK S. & NOAM R. -------
#
# def slidingRP(spikeTimes, params=None):
#     '''
#     Compute the metric for one cluster
#
#     Parameters
#     ----------
#     spikeTimes : numpy.ndarray
#         array of spike times (ms) for one cluster
#
#     params : dict
#         params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize?
#         params.sampleRate : sample rate of the recording (Hz)
#
#     Returns
#     -------
#
#
#
#     maxConfidenceAt10Cont:   Max confidence that you have <= 10% contamination
#     minContWith90Confidence: Minimum contamination for which you have >=90% confidence
#     timeOfLowestCont:        Time at which best score happens
#     nSpikesBelow2:           Number of observed spikes that occur before 2 ms
#     confMatrix:              Full confidence matrix of size nCont x nRP
#     cont:Vector of contamination values tested
#     rp: Vector of refractory period durations tested
#     nACG: the autocorrelogram of the neuron
#     firingRate: firing rate of the cluster, computed as the average acg value from 1-2 seconds
#     '''
#
#     if params is None:
#         params = {}
#         params['sampleRate'] = 30000
#         params['binSizeCorr'] = 1 / params['sampleRate']
#         params['returnMatrix'] = True
#         params['verbose'] = True
#         params['cidx'] = [0]
#
#     seconds_start = time.time()
#     [confMatrix, cont, rp, nACG, firingRate] = computeMatrix(spikeTimes, params)
#     # matrix is [nCont x nRP]
#
#     testTimes = rp > 0.0005  # (in seconds)
#     # only test for refractory period durations greater than 0.5 ms
#
#     maxConfidenceAt10Cont = max(confMatrix[cont == 10, testTimes])  # TODO check behavior if no max
#
#     indsConf90 = np.row_stack(np.where(confMatrix[:, testTimes] >= 90))
#     ii = indsConf90[0]  # row inds
#     jj = indsConf90[1]  # col inds
#
#     try:
#         minI = np.min(ii)
#         idx = np.argmin(ii)
#         minContWith90Confidence = cont[minI]
#         minRP = np.argmax(confMatrix[minI, testTimes])
#
#
#     except:
#         minContWith90Confidence = np.nan
#
#         minRP = np.nan
#
#     try:
#         timeOfLowestCont = rp[minRP + np.where(testTimes)[0][0] + 1]
#     except:
#         timeOfLowestCont = np.nan
#
#     nSpikesBelow2 = sum(nACG[0:np.where(rp > 0.002)[0][0] + 1])
#
#     secondsElapsed = time.time() - seconds_start
#     return maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont, nSpikesBelow2, confMatrix, cont, rp, nACG, firingRate, secondsElapsed
#
#
# def plotSlidingRP(spikeTimes, params=None):
#     '''
#
#
#     Parameters
#     ----------
#     spikeTimes : numpy.ndarray
#         array of spike times (ms)
#     params : dict
#         params.binSizeCorr : bin size for ACG, usually set to 1/sampleRate (s)    TODO: set this up somewhere as same as refDur binsize?
#         params.sampleRate : sample rate of the recording (Hz)
#
#     Returns
#     -------
#     None.
#
#     '''
#     if params is None:
#         clusterlabel = False
#
#     [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
#      nSpikesBelow2, confMatrix, cont, rp, nACG,
#      firingRate, xx] = slidingRP(spikeTimes, params)
#
#     fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
#
#     ax = axs[0]
#
#     ax.bar(rp * 1000, nACG[0:len(rp)], width=np.diff(rp)[0] * 1000, color='k', edgecolor=(1, 0, 0, 0))  # TODO width??
#     ax.set_xlim([0, 5])
#     ax.set_xlabel('Time from spike (ms)')
#     ax.set_ylabel('ACG count (spks)')
#     if clusterlabel:
#         t1 = ('Cluster #%d: FR=%.2f' % (params['cidx'][0], firingRate))
#     ax.fill(np.array([0, 1, 1, 0]) * 0.5, np.array([0, 0, 1, 1]) * ax.get_ylim()[1], 'k', alpha=0.2)
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#
#     ax = axs[1]
#     c = ax.imshow(confMatrix, extent=[rp[0] * 1000, rp[-1] * 1000, cont[0], cont[-1]], aspect='auto', vmin=0, vmax=100,
#                   origin='lower')
#     ax.set_xlim((0, 5))
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.spines['left'].set_visible(False)
#     ax.spines['bottom'].set_visible(False)
#     cbar = fig.colorbar(c, ax=ax, location='right')
#     cbar.set_label('Confidence (%)')
#     ax.invert_yaxis()
#     ax.plot([rp[0] * 1000, rp[-1] * 1000], [10, 10], 'r', linewidth=1)
#
#     if ~np.isnan(timeOfLowestCont):
#         ax.plot(timeOfLowestCont * 1000 * np.array([1, 1]), [cont[0], cont[-1]], 'r', linewidth=1)
#
#         # compute the conf=90 contour
#         # zeropad confMatrix
#         z = np.zeros((np.shape(confMatrix)[0] + 1, np.shape(confMatrix)[1]))
#         z[1:, :] = confMatrix
#
#         ii = np.argmax(z > 90, 0).astype(float)
#         ii[ii == 0] = np.nan
#         contContour = np.empty(np.shape(ii));
#         contContour[:] = np.nan
#         contContour[~np.isnan(ii)] = cont[(ii[~np.isnan(ii)] - 1).astype(int)]
#         ax.plot(rp * 1000, contContour, 'r', linewidth=2)
#     val = ax.get_ylim()[1]
#     # ax.fill(np.array([0, 1, 1, 0])*0.5, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'k',alpha= 1)
#     ax.fill(np.array([0, 1, 1, 0]) * 0.5, np.array([0, 0, 1, 1]) * ax.get_ylim()[1], 'k', alpha=0.2)
#     # ax.add_patch(patches.Rectangle((0,0), 0.5, val, fc = 'k') )
#     ax.set_xlabel('Time from spike (ms)')
#     ax.set_xlim([0, 5])
#     ax.set_ylabel('Contamination (%)')
#     # ax.set_ylim([max(cont),0]) #seems unnecessary?
#     t2 = ('max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms' % (
#     maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont * 1000))
#
#     if minContWith90Confidence >= 10:
#         axs[0].set_title(t1, color='r')
#         axs[1].set_title(t2, color='r')
#     elif nSpikesBelow2 == 0:
#         axs[0].set_title(t1, color='b')
#         axs[1].set_title(t2, color='b')
#     else:
#         axs[0].set_title(t1, color='g')
#         axs[1].set_title(t2, color='g')
#
#     ax = axs[2]
#     ax.plot(rp * 1000, np.squeeze(confMatrix[cont == 10, :]), 'k', linewidth=2.0)
#     ax.set_xlabel('Time from spike (ms)')
#     ax.set_ylabel('Confidence of \leq10% contamination (%)')
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#
#     ax.plot([0, 5], [90, 90], 'r');
#     ax.fill(np.array([0, 1, 1, 0]) * 0.5, np.array([0, 0, 1, 1]) * ax.get_ylim()[1], 'k', alpha=0.2)
#     ax.set_xlim([0, 5]);
#     ax.set_ylim([0, 100]);
#
#     fig.tight_layout()
#
#
# def fitSigmoidACG_All(spikeTimes, spikeClusters, brainRegions, spikeAmps, rp, params):
#     cids = np.unique(spikeClusters)
#
#     # initialize rpMetrics as dict
#     rpFit = {}
#     rpFit['cidx'] = []
#     rpFit['rpEstimate'] = []
#     rpFit['brainRegion'] = []
#     rpFit['amp'] = []
#     rpFit['fr'] = []
#
#     for cidx in range(len(cids)):
#         st = spikeTimes[spikeClusters == cids[cidx]]
#         brainRegion = brainRegions[cidx]
#         fr = len(st) / (st[-1] - st[0])
#
#         sa = spikeAmps[spikeClusters == cids[cidx]]
#         amp = np.nanmedian(sa) * 1000000
#
#         # estimate rp
#         minFR = 1  # spks/s
#         minAmp = 50  # uV
#
#         try:
#             estimatedRP, estimateIdx, xSigmoid, ySigmoid = fitSigmoidACG(st, rp, fr, amp, minFR, minAmp, params)
#
#             rpFit['cidx'].append(cids[cidx])
#             rpFit['rpEstimate'].append(estimatedRP)
#             rpFit['brainRegion'].append(brainRegion)
#             rpFit['fr'].append(fr)
#             rpFit['amp'].append(amp)
#             if verbose:
#                 print('Estimated RP for cluster %d is %.2f ms' % (cids[cidx], estimatedRP))
#         except:
#             continue
#
#     return rpFit
#
#
# def fitSigmoidACG(st, timeBins, fr, amp, minFR=1, minAmp=50, params=None, return_acg=False):
#     '''
#     Compute the autocorrelogram from spike train
#     Apply sigmoid fit if certain conditions are met and return it
#
#     Parameters
#     ----------
#     acg : np.array
#         heights of acg bins (probability of spike for each timebin)
#     timeBins : np.array
#         timeBins used when computing acg
#     params : TYPE
#         DESCRIPTION.
#
#     Returns
#     -------
#     estimatedRP:
#
#     x:
#
#     y:
#
#     '''
#     if params and 'verbose' in params.keys():
#         verbose = params['verbose']
#     else:
#         verbose = False
#     if params and 'sampleRate' in params.keys():
#         sampleRate = params['sampleRate']
#     else:
#         sampleRate = 30000
#     if params and 'binSizeCorr' in params.keys():
#         binSizeCorr = params['binSizeCorr']
#     else:
#         binSizeCorr = 1 / 30000
#     if params and 'numSpikesThresh' in params.keys():
#         numSpikesThresh = params['numSpikesThresh']
#     else:
#         numSpikesThresh = 20  # need at least this many total spikes to compute fit
#
#     # setup for acg
#     clustersIds = [0]  # call the cluster id 0 (not used, but required input for correlograms)
#     spikeClustersACG = np.zeros(len(st), dtype='int8')  # each spike time gets cluster id 0
#     # compute acg
#     acg = correlograms(st, spikeClustersACG, cluster_ids=clustersIds, bin_size=binSizeCorr, sample_rate=sampleRate,
#                        window_size=2, symmetrize=False)[0][0]  # compute acg
#
#     if len(acg) > len(timeBins):
#         acg = acg[0:len(timeBins)]
#
#     if fr >= minFR and amp >= minAmp and len(st) > numSpikesThresh:  # (sum(acg)>numSpikesThresh):
#         # potential todo: insert a case here for if the acg is symmetric?
#         estimatedRP, estimateIdx, xSigmoid, ySigmoid = fit_sigmoid(acg, timeBins)
#     else:
#         # Criteria not matched to run fit
#         estimatedRP = np.nan
#         estimateIdx = np.nan
#         xSigmoid = np.nan
#         ySigmoid = np.nan
#     if return_acg:
#         return estimatedRP, estimateIdx, xSigmoid, ySigmoid, acg
#     else:
#         return estimatedRP, estimateIdx, xSigmoid, ySigmoid
#
#


def plot_acg(ax, acg, timeBins, estimatedIdx=None):
    # Convert time in milliseconds for plotting
    timeBins = timeBins * 1000

    if len(acg) > len(timeBins):
        acg = acg[0:len(timeBins)]
        # TODO change this so the correct bins are taken for whatever values of timeBins
    ax.bar(timeBins, acg, width=np.diff(timeBins)[0], alpha=0.5)
    ax.set_xlim(0, timeBins[-1])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if estimatedIdx is not None:
        # Plot RP as black vertical line
        ax.plot([timeBins[estimatedIdx], timeBins[estimatedIdx]], [0, max(acg)], 'k-')

    ax.set_ylabel('Number of spikes')
    ax.set_xlabel('Time (ms)')

def plotSigmoid(ax, acg, timeBins, ySigmoid, estimatedIdx, estimatedRP):
    plot_acg(ax, acg, timeBins, estimatedIdx=estimatedIdx)
    # Plot on top of ACG the sigmoid
    ax.plot(timeBins[0:len(ySigmoid)] * 1000, ySigmoid, 'b')
    ax.plot(timeBins[estimatedIdx] * 1000, ySigmoid[estimatedIdx], 'rx')
    ax.set_title('Estimated RP:%.2f ms' % estimatedRP)

#     return acg
#
#     # helper functions
#
#
# def find_nearest(array, value):
#     array = np.asarray(array)
#     subtracted = (array - value)
#     valid_idx = np.where(subtracted >= 0)[0]
#     if len(valid_idx) > 0:
#
#         out = valid_idx[subtracted[valid_idx].argmin()]
#
#     else:
#         out = np.nan
#
#     return out
#
#

