# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:40:23 2022

@author: Noam Roth

code to run simulations for slidingRefractory metric
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
from slidingRP.metrics import slidingRP
from scipy import stats
import time
import colorcet as cc
from statsmodels.stats.proportion import proportion_confint as binofit




def genST(rate, duration, rp = 0, params=None):

    '''
    

    Parameters
    ----------
    rate : float
        desired firing rate of simulated neuron (spks/s).
    duration: float
        length of recording of simulated neuron (s)
    rp:
    params : dict, optional

    Returns
    -------
    st: np.array
        array of spike times in seconds. 

    '''
    # print('generating spike train...')

    rSim = rate / (1-(rp*rate)) # rate for simulation, taking account the refractory period of the neuron
    mu = 1 / rSim
    n = rate * duration
    isi = rp + np.random.exponential(mu, int(np.ceil(n * 2)))  # generate

    while sum(isi) < duration:
        isi = np.append(isi, np.random.exponential(mu))

    st = np.cumsum(isi) #convert from isi's to spike times
    if len(np.where(st < duration)[0]) > 0:
        st = st[0:np.where(st < duration)[0][(-1)]] #cut spike times after the recording duration
    else:
        st = []

    if params is None:
        print('No firing rate plot')
    elif params['checkFR']:
        print('plotting FR...')
        fig, ax = plt.subplots(1, 1)
        histWin = 0.5  # seconds
        binnedFR = firing_rate(st, hist_win=histWin, fr_win=10)
        ax.plot(np.arange(binnedFR.size) * histWin, binnedFR, 'k')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('FR (spks/s)')
        ax.set_title('FR for simulated neuron, mean %d spikes/second, duration %d seconds' % (rate, duration))

    return st

def genChangingST(rate,duration,params,delta):

    # delta is in fractional increase e.g. 0.5
    if delta ==0:
        nchunk=1
        rates = np.arange(rate)

    else:
        nchunk = 100  # parameter, how many chunks to split the data into
        rate_end = rate * (1 + delta)
        rate_diff = rate_end - rate
        rateDiffChunk = rate_diff / (nchunk - 1)
        rates = np.arange(rate, rate_end + rateDiffChunk, rateDiffChunk)

    chunkLength = duration / nchunk
    st = np.array([])

    for r in rates:
        stChunk = genST(r, chunkLength, params)
        if len(st) > 0:
            stChunkConcatenate = stChunk + st[-1]
        else:
            stChunkConcatenate = stChunk
        st = np.append(st, stChunkConcatenate)
    return st

def simulateContNeurons(params):

    #initialize matrices of percent neurons pass
    passPct = np.empty([len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPct[:] = np.nan

    passPct2MsNoSpikes = np.zeros(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPct2MsNoSpikes[:] = np.nan

    passPctHalfInactive = np.zeros(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctHalfInactive[:] = np.nan

    passPctHill15 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctHill15[:] = np.nan

    passPctHill2 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctHill2[:] = np.nan

    passPctHill3 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctHill3[:] = np.nan

    #check whether to simultaneously run comparison Llobet simulations
    if params and 'runLlobet' in params:
        runLlobet = params['runLlobet']
        print('running Llobet regular')
    else:
        runLlobet = False #default is false to decrease runtime if not testing simulations for analysis purposes

    if params and 'runLlobetPoiss' in params:
        runLlobetPoiss = params['runLlobetPoiss']
        print('running Llobet Poiss')
    else:
        runLlobetPoiss = False #default is false to decrease runtime if not testing simulations for analysis purposes

    #create empty passPct to output either way, if runLlobet False these will not get filled
    passPctLlobet15 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobet15[:] = np.nan

    passPctLlobet2 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobet2[:] = np.nan

    passPctLlobet3 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobet3[:] = np.nan

    #create empty passPct to output either way, if runLlobetPoiss False these will not get filled
    passPctLlobetPoiss15 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobetPoiss15[:] = np.nan

    passPctLlobetPoiss2 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobetPoiss2[:] = np.nan

    passPctLlobetPoiss3 = np.empty(
        [len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPctLlobetPoiss3[:] = np.nan



    # start time to time simulations
    start_time = time.time()
    for j, recDurScalar in enumerate(params['recDurs']):
        recDur = recDurScalar * 3600
        print('recording Duration %d' % recDur)
        params['recordingDuration'] = recDur
        for i, rp in enumerate(params['RPs']):
            print('refractory period duration %.3f' % rp)

            bidx = 0
            for baseRate in params['baseRates']:
                print('baseRate %.2f' % baseRate)
                cidx = 0
                for c in params['contRates']:
                    # contRate = baseRate * c   #contamination rate is the base firing rate times the fraction contamination

                    contRate = c * baseRate / (1-c)
                    #from Nick's code:
                    #contRate = contProp(c) * baseRate / (1 - contProp(c));

                    #initialize vectors to save passing neurons
                    passVec = np.empty(params['nSim'])
                    passVec[:] = np.nan

                    passVec2MsNoSpikes = np.empty(params['nSim'])
                    passVec2MsNoSpikes[:] = np.nan

                    passVecHalfInactive = np.empty(params['nSim'])
                    passVecHalfInactive[:] = np.nan

                    passVecHill15 = np.empty(params['nSim'])
                    passVecHill15[:] = np.nan

                    passVecHill2 = np.empty(params['nSim'])
                    passVecHill2[:] = np.nan

                    passVecHill3 = np.empty(params['nSim'])
                    passVecHill3[:] = np.nan

                    # set runLlobet to True to also run the standard Llobet metric for comparison
                    if runLlobet:
                        passVecLlobet15 = np.empty(params['nSim'])
                        passVecLlobet15[:] = np.nan

                        passVecLlobet2 = np.empty(params['nSim'])
                        passVecLlobet2[:] = np.nan

                        passVecLlobet3 = np.empty(params['nSim'])
                        passVecLlobet3[:] = np.nan

                    # set runLlobetPoiss to True to also run a Poisson version of the Llobet metric for comparison
                    if runLlobetPoiss:
                        passVecLlobetPoiss15 = np.empty(params['nSim'])
                        passVecLlobetPoiss15[:] = np.nan

                        passVecLlobetPoiss2 = np.empty(params['nSim'])
                        passVecLlobetPoiss2[:] = np.nan

                        passVecLlobetPoiss3 = np.empty(params['nSim'])
                        passVecLlobetPoiss3[:] = np.nan

                    for n in range(params['nSim']):
                        if n % 20 == 0:
                            print(':', end="")
                        if c == (params['contRates'][-1]) and n == (params['nSim'] - 1):
                            print(' ')

                        st = genST(baseRate, recDur, rp, params)  # generate a spike train with the current base rate

                        if c > 0:
                            contST = genST(contRate, recDur, 0, params) #add contaminating neuron
                        else:
                            contST = []

                        combST = np.sort(np.concatenate(
                            (st, contST)))  # put spike times from both spike trains together (and sort chronologically)


                        #run the slidingRP metric on combST - the contaminated neuron spike train.
                        [maxConfidenceAtContThresh, minContWithConfidenceThresh, timeOfLowestCont,
                         nSpikesBelow2, confMatrix, cont, rpVec, nACG,
                         firingRate, secondsElapsed] = slidingRP(combST, params)

                        #the neuron passes if the minimum contamination (at our confidence threshold) is smaller than
                        # # our contamination threshold (in percent)
                        passVec[n] = int(maxConfidenceAtContThresh>params['confidenceThresh'])
                        # if minContWithConfidenceThresh <= params['contaminationThresh']:
                        #     passVec[n] = 1
                        # else:
                        #     passVec[n] = 0
                        # testTimes = rpVec>0.0005
                        # passVec[n] = max(confMatrix())
                        # if testSliding
                        #     % normal
                        #     sliding
                        #     rp
                        #     test
                        #     result
                        #     simRes(n) = max(confMatrix(rpTestVals > 0.0005)) > conf;

                        # to evaluate the 2ms no spikes condition, check whether there are no spikes below 2ms.
                        passVec2MsNoSpikes[n] = passVec[n]
                        if nSpikesBelow2 == 0:
                            passVec2MsNoSpikes[n] = 1

                        # To compare the performance of our metric on neurons that stop firing halfway through the
                        # recording duration, take combST (our contaminated neuron's spike train) and stop it halfway
                        frHalfInactive = combST[0:np.where(combST > recDur / 2)[0][0]]
                             # remove spikes from second half of the recording

                        # Here we only need to save minCont to determine pass/fail
                        minContWith90ConfidenceHalfInactive = slidingRP(frHalfInactive, params)[1]

                        if minContWith90ConfidenceHalfInactive <= 10:
                            passVecHalfInactive[n] = 1
                        else:
                            passVecHalfInactive[n] = 0


                        # Hill comparison with 1.5ms, 2ms, or 3 ms
                        fpRate15 = HillMetric(firingRate, recDur, nACG, rpVec, refDur=0.0015, minISI=0)
                        fpRate2 = HillMetric(firingRate, recDur, nACG, rpVec, refDur = 0.002, minISI=0)
                        fpRate3 = HillMetric(firingRate, recDur, nACG, rpVec, refDur = 0.003, minISI=0)

                        # add these false positive rates to passVec with a threshold of 10% contamination
                        if fpRate15 <= params['contaminationThresh']/100:
                            passVecHill15[n] = 1
                        else:
                            passVecHill15[n] = 0


                        if fpRate2 <= params['contaminationThresh']/100:
                            passVecHill2[n] = 1
                        else:
                            passVecHill2[n] = 0

                        if fpRate3 <= params['contaminationThresh']/100:
                            passVecHill3[n] = 1
                        else:
                            passVecHill3[n] = 0



                        if runLlobet:
                            testedCont = params['contaminationThresh']/100
                            # Llobet comparison with 1.5ms, 2ms, or 3 ms
                            fpLlobet15, confLlobet15 = LlobetMetric(firingRate, recDur, nACG, rpVec, testedCont, refDur=0.0015, minISI=0)
                            fpLlobet2, confLlobet2 = LlobetMetric(firingRate, recDur, nACG, rpVec, testedCont, refDur=0.002, minISI=0)
                            fpLlobet3, confLlobet3 = LlobetMetric(firingRate, recDur, nACG, rpVec, testedCont,  refDur=0.003, minISI=0)

                            # add these false positive rates to passVec with a threshold of 10% contamination
                            if fpLlobet15 <= params['contaminationThresh']/100:
                                passVecLlobet15[n] = 1
                            else:
                                passVecLlobet15[n] = 0

                            if fpLlobet2 <= params['contaminationThresh']/100:
                                passVecLlobet2[n] = 1
                            else:
                                passVecLlobet2[n] = 0

                            if fpLlobet3 <= params['contaminationThresh']/100:
                                passVecLlobet3[n] = 1
                            else:
                                passVecLlobet3[n] = 0
                        # add these false positive rates to passVec with a threshold of 10% contamination
                        if confLlobet15 <= params['contaminationThresh'] / 100:
                            passVecLlobetPoiss15[n] = 1
                        else:
                            passVecLlobetPoiss15[n] = 0

                        if confLlobet2 <= params['contaminationThresh'] / 100:
                            passVecLlobetPoiss2[n] = 1
                        else:
                            passVecLlobetPoiss2[n] = 0

                        if confLlobet3 <= params['contaminationThresh'] / 100:
                            passVecLlobetPoiss3[n] = 1
                        else:
                            passVecLlobetPoiss3[n] = 0

                    passPct[j, i, bidx, cidx] = sum(passVec) / params['nSim'] * 100

                    passPct2MsNoSpikes[j, i, bidx, cidx] = sum(passVec2MsNoSpikes) / params['nSim'] * 100

                    passPctHalfInactive[j, i, bidx, cidx] = sum(passVecHalfInactive) / params['nSim'] * 100

                    passPctHill15[j, i, bidx, cidx] = sum(passVecHill15) / params['nSim'] * 100
                    passPctHill2[j, i, bidx, cidx] = sum(passVecHill2) / params['nSim'] * 100
                    passPctHill3[j, i, bidx, cidx] = sum(passVecHill3) / params['nSim'] * 100

                    if runLlobet:
                        passPctLlobet15[j, i, bidx, cidx] = sum(passVecLlobet15) / params['nSim'] * 100
                        passPctLlobet2[j, i, bidx, cidx] = sum(passVecLlobet2) / params['nSim'] * 100
                        passPctLlobet3[j, i, bidx, cidx] = sum(passVecLlobet3) / params['nSim'] * 100
                    if runLlobetPoiss:
                        passPctLlobetPoiss15[j, i, bidx, cidx] = sum(passVecLlobetPoiss15) / params['nSim'] * 100
                        passPctLlobetPoiss2[j, i, bidx, cidx] = sum(passVecLlobetPoiss2) / params['nSim'] * 100
                        passPctLlobetPoiss3[j, i, bidx, cidx] = sum(passVecLlobetPoiss3) / params['nSim'] * 100

                    cidx += 1

                bidx += 1

    current_time = time.time()
    elapsed_time = current_time - start_time
    print('Loop with %f iterations takes %f seconds' % (params['nSim'], elapsed_time))
    return passPct, passPct2MsNoSpikes, passPctHalfInactive, passPctHill15, passPctHill2, passPctHill3, passPctLlobet15, passPctLlobet2, passPctLlobet3, passPctLlobetPoiss15, passPctLlobetPoiss2, passPctLlobetPoiss3


def simulateChangingContNeurons(params):
    # initialize matrices of percent neurons pass
    delta = params['delta']
    passPct = np.empty([len(params['recDurs']), len(params['RPs']), len(params['baseRates']), len(params['contRates'])])
    passPct[:] = np.nan

    # start time to time simulations
    start_time = time.time()
    for j, recDurScalar in enumerate(params['recDurs']):
        recDur = recDurScalar * 3600
        print('recording Duration %d' % recDur)
        for i, rp in enumerate(params['RPs']):
            print('refractory period duration %.3f' % rp)


            bidx = 0
            for baseRate in params['baseRates']:
                print('baseRate %.2f' % baseRate)
                cidx = 0
                for c in params['contRates']:
                    contRate = baseRate * c

                    passVec = np.empty(params['nSim'])
                    passVec[:] = np.nan

                    for n in range(params['nSim']):
                        if n % 20 == 0:
                            print('-', end="")
                        if c == (params['contRates'][-1]) and n == (params['nSim'] - 1):
                            print(' ')
                        st = genChangingST(baseRate, recDur, params,delta)  # generate a spike train with the current base rate
                        isi = np.diff(np.insert(st, 0, 0))
                        isi = np.delete(isi, np.where(isi < rp)[
                            0])  # get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
                        st = np.cumsum(isi)

                        if c > 0:
                            contST = genChangingST(contRate, recDur, params,delta)
                        else:
                            contST = []
                        combST = np.sort(np.concatenate(
                            (st, contST)))  # put spike times from both spike trains together (and sort chronologically)

                        [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
                         nSpikesBelow2, confMatrix, cont, rpVec, nACG,
                         firingRate, secondsElapsed] = slidingRP(combST, params)

                        if minContWith90Confidence <= params['contaminationThresh']:
                            passVec[n] = 1
                        else:
                            passVec[n] = 0


                    passPct[j, i, bidx, cidx] = sum(passVec) / params['nSim'] * 100

                    cidx += 1

                bidx += 1

    current_time = time.time()
    elapsed_time = current_time - start_time
    print('Loop with %f iterations takes %f seconds' % (params['nSim'], elapsed_time))


    return passPct


def HillMetric(firingRate, recDur, nACG, rpVec, refDur, minISI=0):

    # number of violations between minISI and refDur
    nViol = np.sum(nACG[np.where(rpVec > minISI)[0][0] : np.where(rpVec > refDur)[0][0] + 1])

    # time for violations to occur
    N_t = firingRate * recDur #total number of spikes

    #original Hill metric fpRate
    # false positive rate (f_1^p in Hill paper)
    #this is the case of one contaminating neuron, see Llobet et al.
    fpRate = 1/2 *  (1- np.sqrt(1- ( (2*nViol*recDur) / (N_t**2 * (refDur - minISI)) ) ) )

    # For cases where this is greater than 1, set to 1
    if fpRate > 1:
        fpRate = 1

    return fpRate


def LlobetMetric(firingRate, recDur, nACG, rpVec, testedCont, refDur, minISI=0):
    #compute false positive rate as defined in Llobet et al.
    #This is essentially the same as our (revised) metric, but only tested at one time point instead of all possible bins.

    # number of violations between minISI and refDur
    nViol = np.sum(nACG[np.where(rpVec > minISI)[0][0]: np.where(rpVec > refDur)[0][0] + 1])
    contaminationRate = firingRate * testedCont

    # time for violations to occur
    N_t = firingRate * recDur  # total number of spikes
    N_c = contaminationRate * recDur  #total number of contaminating spikes you would expect under the tested contamination value
    N_b = N_t - N_c # the "base" number of spikes under the inputted (tested) contamination value


    violationTime = 2 * N_t * (refDur - minISI)
    # this violationTime was previously written equivalently as  2 * firingRate * recDur * (refDur - minISI)

    # rate of violations
    violationRate = nViol / violationTime

    # false positive rate (called f_1^p in Hill paper)
    fpRate = 1 - np.sqrt(1 - ((nViol * recDur)/ (N_t**2 * (refDur - minISI))))

    expectedViol = 2 * refDur * 1/recDur * N_c * (N_b + N_c/2) #number of expected violations, as defined in Llobet et al., this gets rid of the minus 1, not sure what this difference is about

    # the confidence that this neuron is contaminated at a level less than contaminationProp, given the number of true
    # observed violations and under the assumption of Poisson firing
    confidenceScore = 1 - stats.poisson.cdf(nViol, expectedViol)

    # For cases where this is greater than 1, set to 1
    if fpRate > 1:
        fpRate = 1

    return fpRate, confidenceScore

def plotSimulations(pc, params, savefile, rp_valFig1 = 0.0025,frPlot = [0.5,1,2,5,10], input_color=cc.linear_protanopic_deuteranopic_kbw_5_95_c34,
                    subplot1 = False, subplot2 = False, subplot3 = False, subplot4 = False,
                    sensSpecPlot=False, plotType='paper', zoomCont = False, addPCflag = 0, highCont=False):
    # plot type = 'full' or plot type = 'paper'
    # compute confidence intervals

    if plotType == 'paper':
        spinesSetting = False
    else:
        spinesSetting = True

    # exception for concatenated params
    if False:
        count = np.empty(np.shape(pc))
        count[:, :, [0, 1, 2, 4, 8, 10, 11], :] = pc[:, :, [0, 1, 2, 4, 8, 10, 11], :] / 100 * 500
        count[:, :, [3, 5, 6, 7, 9], :] = pc[:, :, [3, 5, 6, 7, 9], :] / 100 * 100

        CI_scaled1 = binofit(count[:, :, [0, 1, 2, 4, 8, 10, 11], :], 500)
        CI_scaled2 = binofit(count[:, :, [3, 5, 6, 7, 9], :], 100)

        CI_scaled = np.empty((2,) + np.shape(pc))
        CI_scaled[:, :, :, [0, 1, 2, 4, 8, 10, 11], :] = CI_scaled1
        CI_scaled[:, :, :, [3, 5, 6, 7, 9], :] = CI_scaled2


    else:
        count = pc / 100 * params['nSim']  # number of correct trials
        CI_scaled = binofit(count, params['nSim'])
    CI = [x * 100 for x in CI_scaled]
    colors = matplotlib.cm.Blues(np.linspace(0.2, 1, len(params['baseRates'])))  # never go below .2, hard to see

    if subplot1:
        if plotType == 'paper':
            print('in plot type paper figure 1')
            # for this figure, plotType = 'full' or plotType = 'paper' are the same.
            numrows = 1  # len(params['recDurs'])
            numcols = 1  # len(params['RPs'])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax = axs

            fig2, axs2 = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax2 = axs2
            titlexval = 0.5
            titleyvals = [1, 0.77, 0.52, 0.27]

            # plot just recDur == 2:
            recDurs = params['recDurs']
            recDursInd = np.where(recDurs == 2)[0]

            # plot just RP = 2:
            rps = params['RPs']
            rpInd = np.where(rps == rp_valFig1)[0]

            # plot just baseRates = [0.5, 1, 5, 10]
            # fr_plot = [0.05,0.25,0.45,0.75,1,10]
            fr_plot = [0.5, 1, 2, 5, 10]
            # fr_plot = [0.05,0.25,0.45,0.75,1,2, 5,10]
            # fr_plot = params['baseRates']
            frs = params['baseRates']
            frInd = [x for x in range(len(frs)) if frs[x] in frPlot]

            c = input_color  # cc.linear_protanopic_deuteranopic_kbw_5_95_c34
            c = c[::-1]  # if using linear_blue37 or protanopic, flip the order
            colors = [c[x] for x in np.round(np.linspace(0.2, 0.75, len(params['baseRates'])) * 255).astype(int)]
            # colors = [c[x] for x in np.round(np.linspace(0.25, 1, len(fr_plot))*255).astype(int)]

            pltcnt = 0
            for j, recDur in enumerate(recDurs[recDursInd]):
                print(j)
                print(recDur)
                for i, rp in enumerate(rps[rpInd]):
                    pltcnt += 1
                    # different base rates get different colors
                    for b, baseRate in enumerate(frInd):
                        fr_use = frs[frInd[b]]
                        print('here in loop')
                        if zoomCont:
                            lowerCI = CI[0][recDursInd, rpInd, frInd[b], 6:15][0]
                            upperCI = CI[1][recDursInd, rpInd, frInd[b], 6:15][0]
                            x = params['contRates'][6:15] * 100
                            y = pc[recDursInd, rpInd, frInd[b], 6:15][0]
                            # ax.xaxis.set_ticks([6,7, 8, 9, 10, 11, 12, 13,14])


                        else:
                            lowerCI = CI[0][recDursInd, rpInd, frInd[b], :][0]
                            upperCI = CI[1][recDursInd, rpInd, frInd[b], :][0]
                            x = params['contRates'] * 100
                            y = pc[recDursInd, rpInd, frInd[b], :][0]

                        if highCont:
                            ax.plot(x[:-1], y[:-1], '.-', color=colors[frInd[b]],
                                label=fr_use)  # NR changed 12/30 frInd[b] to [b] becaues changed how indexing colors
                            ax.fill_between(x[:-1], lowerCI[:-1], upperCI[:-1], color=colors[frInd[b]], alpha=.5)

                        else:
                            ax.plot(x, y, '.-', color=colors[frInd[b]],
                                    label=fr_use)  # NR changed 12/30 frInd[b] to [b] becaues changed how indexing colors
                            ax.fill_between(x, lowerCI, upperCI, color=colors[frInd[b]], alpha=.5)



                        if addPCflag:
                            # find 10% point
                            tenPercentPoint = np.where(x == 10)[0][0]
                            percentCorrect = np.concatenate(
                                (y[0:(tenPercentPoint + 1)], 100 - y[(tenPercentPoint + 1):]))
                            ax2.plot(x[:-1], percentCorrect[:-1], color=colors[frInd[b]])
                            # ax2.plot(25, percentCorrect[-1], '.', color= colors[b])

                            ax2.set_ylabel('Percent Correct')
                        if highCont:
                            ax.plot(x[:-1], y[:-1], '.-', color=colors[frInd[b]],
                                    label=fr_use)  # NR changed 12/30 frInd[b] to [b] becaues changed how indexing colors

                            ax.fill_between(x[:-1], lowerCI[:-1], upperCI[:-1], color=colors[frInd[b]], alpha=.5)
                        # ax.vlines(x=25, ymin=lowerCI[-1], ymax=upperCI[-1],
                        #           colors=colors[frInd[b]])
                        if (i == 0):
                            ax.set_ylabel('Percent pass')
                        else:
                            ax.yaxis.set_ticklabels([])
                        if (pltcnt > (numplots - numcols)):
                            ax.set_xlabel('Contamination (%)')
                        if (pltcnt == 1):
                            ax.set_title('True RP: %.1f ms' % (rp * 1000))
                        else:
                            ax.set_title('%.1f ms' % (rp * 1000))
                        ax.set_ylim(-10, 110)
                        ax.spines.right.set_visible(spinesSetting)
                        ax.spines.top.set_visible(spinesSetting)
                        # ax.xaxis.set_ticks([0, 10, 20])
                        ax.vlines(10, -10, 110, color='k')
                # fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
                if recDur != 1:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hours' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
                else:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hour' % recDur, ha="center",
                                va="top", fontsize=14, color="k")

            print('now here')
            handles, labels = ax.get_legend_handles_labels()
            # fig.tight_layout(pad=1, w_pad=1.1, h_pad=1.3)
            # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            plt.subplots_adjust(hspace=0.6, wspace=0.5)
            fig.tight_layout()
            fig.legend(handles, labels, title='Firing rate (spk/s)', loc='upper right', bbox_to_anchor=(1, .9))
            handles2, labels2 = ax2.get_legend_handles_labels()
            fig2.legend(handles, labels, title='Firing rate (spk/s)', loc='upper right', bbox_to_anchor=(1.1, 1))

            fig.savefig(savefile + '_Main.svg', dpi=500)
            fig.savefig(savefile + '_Main.png', dpi=500)
            print(savefile)
            if addPCflag:
                fig2.savefig(savefile + '_MainPC.svg', dpi=500)
                fig2.savefig(savefile + '_MainPC.png', dpi=500)
        if plotType == 'heatmap':
            print('in plot type heatmap figure 1')
            # for this figure, plotType = 'full' or plotType = 'paper' are the same.
            numrows = 1  # len(params['recDurs'])
            numcols = 1  # len(params['RPs'])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax = axs

            fig2, axs2 = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax2 = axs2
            titlexval = 0.5
            titleyvals = [1, 0.77, 0.52, 0.27]

            # plot just recDur == 2:
            recDurs = params['recDurs']
            recDursInd = np.where(recDurs == 2)[0]

            # plot just RP = 2:
            rps = params['RPs']
            rpInd = np.where(rps == 0.002)[0]

            # plot just baseRates = [0.5, 1, 5, 10]
            # fr_plot = [0.5, 1, 5, 10] #changed on 12/26 to look at all
            fr_plot = params['baseRates']
            frs = params['baseRates']
            frInd = [x for x in range(len(frs)) if frs[x] in fr_plot]

            pltcnt = 0
            for j, recDur in enumerate(recDurs[recDursInd]):
                print(j)
                print(recDur)
                for i, rp in enumerate(rps[rpInd]):
                    pltcnt += 1
                    # different base rates get different colors
                    # for b, baseRate in enumerate(frInd):
                    # fr_use = frs[frInd[b]]
                    # lowerCI = CI[0][recDursInd, rpInd, frInd[b], :][0]
                    # upperCI = CI[1][recDursInd, rpInd, frInd[b], :][0]
                    x = params['contRates'] * 100
                    y = pc[recDursInd, rpInd, :, :][0]

                    # ax.plot(x, y, '.-', color=colors[frInd[b]], label=fr_use)
                    ax.imshow(y, cmap="plasma", aspect="auto")
                    if addPCflag:
                        # find 10% point
                        tenPercentPoint = np.where(x == 10)[0][0]
                        percentCorrect = np.concatenate(
                            (y[0:(tenPercentPoint + 1)], 100 - y[(tenPercentPoint + 1):]))
                        ax2.imshow(percentCorrect, cmap="plasma", aspect="auto")

                        # ax2.plot(x, percentCorrect, color=colors[frInd[b]], label=fr_use)
                        # ax2.set_ylabel('Percent Correct')
                        # ax.fill_between(x, lowerCI, upperCI, color=colors[frInd[b]], alpha=.5)
                    if (i == 0):
                        ax.set_ylabel('Percent pass')
                    else:
                        ax.yaxis.set_ticklabels([])
                    if (pltcnt > (numplots - numcols)):
                        ax.set_xlabel('Contamination (%)')
                    if (pltcnt == 1):
                        ax.set_title('True RP: %.1f ms' % (rp * 1000))
                    else:
                        ax.set_title('%.1f ms' % (rp * 1000))
                    ax.set_ylim(-10, 110)
                    ax.spines.right.set_visible(spinesSetting)
                    ax.spines.top.set_visible(spinesSetting)
                    ax.xaxis.set_ticks([0, 10, 20])
                    ax.vlines(10, -10, 110, color='k')
                # fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
                if recDur != 1:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hours' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
                else:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hour' % recDur, ha="center",
                                va="top", fontsize=14, color="k")

            handles, labels = ax.get_legend_handles_labels()
            # fig.tight_layout(pad=1, w_pad=1.1, h_pad=1.3)
            # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            plt.subplots_adjust(hspace=0.6, wspace=0.5)
            fig.tight_layout()
            fig.legend(handles, labels, title='Firing rate (spk/s)', loc='upper right', bbox_to_anchor=(1, .9))
            handles2, labels2 = ax2.get_legend_handles_labels()
            fig2.legend(handles, labels, title='Firing rate (spk/s)', loc='upper right', bbox_to_anchor=(1.1, 1))

            fig.savefig(savefile + '_Main.svg', dpi=500)
            fig2.savefig(savefile + '_MainPC.svg', dpi=500)
            fig.savefig(savefile + '_Main.png', dpi=500)
            fig2.savefig(savefile + '_MainPC.png', dpi=500)

        if plotType == 'full':
            # for this figure, plotType = 'full' or plotType = 'paper' are the same.
            numrows = len(params['recDurs'])
            numcols = len(params['RPs'])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(20, 20))
            titlexval = 0.5
            titleyvals = [0.9, 0.7, 0.5, 0.3]

            pltcnt = 0
            for j, recDur in enumerate(params['recDurs']):
                for i, rp in enumerate(params['RPs']):
                    pltcnt += 1
                    if len(params['recDurs']) > 1 and len(params['RPs']) > 1:
                        ax = axs[j, i]
                    else:
                        ax = axs[(j + 1) * i]
                    # different base rates get different colors
                    for b, baseRate in enumerate(params['baseRates']):
                        lowerCI = CI[0][j, i, b, :]
                        upperCI = CI[1][j, i, b, :]
                        x = params['contRates'] * 100
                        y = pc[j, i, b, :]
                        ax.plot(x, y, '.-', color=colors[b], label=baseRate)

                        ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                        if (i == 0):
                            ax.set_ylabel('Percent pass')
                        if (pltcnt > (numplots - numcols)):
                            ax.set_xlabel('Percent contamination')
                        if (pltcnt == 1):
                            ax.set_title('True RP: %.1f ms' % (rp * 1000))
                        else:
                            ax.set_title('%.1f ms' % (rp * 1000))
                        ax.set_ylim(-10, 110)
                        ax.spines.right.set_visible(spinesSetting)
                        ax.spines.top.set_visible(spinesSetting)
                # fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
                if recDur != 1:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hours' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
                else:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hour' % recDur, ha="center",
                                va="top", fontsize=14, color="k")

            handles, labels = ax.get_legend_handles_labels()
            # fig.tight_layout(pad=1, w_pad=1.1, h_pad=1.3)
            # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            plt.subplots_adjust(hspace=0.3)

            fig.legend(handles, labels, loc='upper right')  # ,bbox_to_anchor=(1.1, 1))

            fig.savefig(savefile + '_Main.png', dpi=500)
            # fig.savefig(savefile + '_Main.svg', dpi=500)

    if subplot2:
        if plotType == 'full':
            numrows = len(params['contRates'][::2])
            numcols = len(params['RPs'])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(20, 20))
            pltcnt = 0
            for j, contRate in enumerate(params['contRates'][::2]):
                for i, rp in enumerate(params['RPs']):
                    pltcnt += 1
                    if len(params['contRates'][::2]) > 1 and len(params['RPs']) > 1:
                        ax = axs[j, i]

                    else:
                        ax = axs[(j + 1) * i]

                    # different base rates get different colors
                    for b, baseRate in enumerate(params['baseRates']):

                        lowerCI = CI[0][:, i, b, j * 2]
                        upperCI = CI[1][:, i, b, j * 2]
                        x = params['recDurs']
                        y = pc[:, i, b, j * 2]
                        ax.plot(x, y, '.-', color=colors[b], label=baseRate)

                        ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                        if (i == 0):
                            ax.set_ylabel('Percent pass')
                        if (pltcnt > (numplots - numcols)):
                            ax.set_xlabel('Recording duration (hours)')
                        if (pltcnt == 1):
                            ax.set_title('True RP: %.1f ms; contamination: %d' % (rp * 1000, contRate * 100) + '%')
                        else:
                            ax.set_title('%.1f ms; %d' % (rp * 1000, contRate * 100) + '%')
                        ax.set_ylim(-10, 110)
                        ax.spines.right.set_visible(spinesSetting)
                        ax.spines.top.set_visible(spinesSetting)
                # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
                fig.suptitle('Proportion contamination: %.2f' % contRate, x=.5, y=1.1)
            handles, labels = ax.get_legend_handles_labels()

            # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            fig.tight_layout()
            fig.legend(handles, labels, loc='upper right')
            fig.savefig(savefile + '_recDur_full.svg', dpi=500)

        if plotType == 'paper':
            numrows = 1
            numcols = 1
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax = axs  # for the case of just one subplot
            # plot just contRates 0.08 to 0.12:
            cr = params['contRates'];

            # plot just RP = 2:
            rps = params['RPs']
            rpInd = np.where(rps == 0.002)[0]

            # plot just fr = 5:
            frs = np.array(params['baseRates'])
            frInd = np.where(frs == 5)[0][0]

            # colors = matplotlib.cm.Set1(np.linspace(0, 1, 10))
            c = cc.b_linear_bgyw_15_100_c67#input_color  # cc.linear_protanopic_deuteranopic_kbw_5_95_c34
            c = c[::-1]  # if using linear_blue37 or protanopic, flip the order
            colors = [c[x] for x in np.round(np.linspace(0.2, 0.75, len(params['recDurs'])) * 255).astype(int)]



            pltcnt = 0
            linewidths = [1,1,1,1]#[0.5, 1, 2, 3]
            for j, recDur in enumerate(params['recDurs']):
                for i, rp in enumerate(rps[rpInd]):
                    pltcnt += 1

                    # different base rates get different colors
                    # fix baseRate at frs[frInd]: previously did for b, baseRate in enumerate(frs[frInd]):
                    b = 0
                    baseRate = frs[frInd]# Todo change x axis
                    lowerCI = CI[0][j, rpInd[0], frInd, :]
                    upperCI = CI[1][j, rpInd[0], frInd, :]
                    x = cr * 100
                    y = pc[j, rpInd[0], frInd, :]
                    ax.plot(x, y, '.-', color=colors[j], linewidth=linewidths[j], label=recDur)

                    ax.fill_between(x, lowerCI, upperCI, color=colors[j], alpha=.3)
                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('Contamination (%)')
                    # ax.set_title('True RP: %.1f ms; contamination: %d' % (rp * 1000, contRate * 100) + '%')
                    ax.set_ylim(-10, 110)
                    ax.spines.right.set_visible(spinesSetting)
                    ax.spines.top.set_visible(spinesSetting)
                # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
                # fig.suptitle('Proportion contamination: %.2f' % contRate*100, x=.5, y=1.1)
            handles, labels = ax.get_legend_handles_labels()

            # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            fig.tight_layout()
            fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1), title='Recording duration (hrs)')
            fig.savefig(savefile + '_recDur.svg', dpi=500)
            fig.savefig(savefile + '_recDur.png', dpi=500)


        if plotType == 'paper_full':
            numrows = len(params['baseRates'])
            numcols = len(params['RPs'])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(20, 20))

            # plot just contRates 0.08 to 0.12:
            cr = params['contRates'];

            # plot just RP = 2:
            rps = params['RPs']
            rpInd = np.arange(len(params['RPs']))  # np.where(rps == 0.002)[0]

            # plot just fr = 5:
            frs = params['baseRates']
            frInd = np.arange(len(params['baseRates']))  # np.where(frs == 5)[0]

            pltcnt = 0
            linewidths = [0.5, 1, 2, 3]
            for j, recDur in enumerate(params['recDurs']):
                for i, rp in enumerate(rps[rpInd]):
                    test_12 = pc[j, i, :, 12]
                    print(test_12)
                    pltcnt += 1

                    # different base rates get different colors
                    print(' ')
                    for b, baseRate in enumerate(frs[frInd]):
                        ax = axs[b, i]  # for the case of just one subplot
                        # Todo change x axis
                        lowerCI = CI[0][j, i, b, :]
                        upperCI = CI[1][j, i, b, :]
                        x = cr * 100
                        y = pc[j, i, b, :]
                        ax.plot(x, y, '.-', color=colors[frInd[-1]], linewidth=linewidths[j], label=recDur)

                        ax.fill_between(x, lowerCI, upperCI, color=colors[frInd[-1]], alpha=.3)
                        ax.set_ylabel('Percent pass')
                        ax.set_xlabel('Contamination (%)')
                        # ax.set_title('True RP: %.1f ms; contamination: %d' % (rp * 1000, contRate * 100) + '%')
                        ax.set_ylim(-10, 110)
                        ax.spines.right.set_visible(spinesSetting)
                        ax.spines.top.set_visible(spinesSetting)
                        ax.set_title('RP %.2f; FR %.2f' % (rp * 1000, baseRate))
                        ymin, ymax = ax.get_ylim()
                        ax.vlines(10, ymin, ymax, color='k')
                        ax.vlines(12, ymin, ymax, color='b', linestyles='--')
                # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
                # fig.suptitle('Proportion contamination: %.2f' % contRate*100, x=.5, y=1.1)
            handles, labels = ax.get_legend_handles_labels()

            # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            fig.tight_layout()
            fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1), title='Recording duration (hrs)')
            fig.savefig(savefile + '_recDur_paperfull.svg', dpi=500)

    if subplot3:
        if plotType == 'full':
            numrows = len(params['recDurs'])
            numcols = len(params['contRates'][::2])
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(20, 20))
            titlexval = 0.5
            titleyvals = [1, 0.75, 0.5, 0.25]
            pltcnt = 0
            for j, recDur in enumerate(params['recDurs']):
                for i, contRate in enumerate(params['contRates'][::2]):
                    pltcnt += 1
                    if len(params['recDurs']) > 1 and len(params['contRates']) > 1:
                        ax = axs[j, i]
                    else:
                        ax = axs[(j + 1) * i]
                    # different base rates get different colors
                    for b, baseRate in enumerate(params['baseRates']):
                        lowerCI = CI[0][j, :, b, i * 2]
                        upperCI = CI[1][j, :, b, i * 2]
                        x = params['RPs']
                        y = pc[j, :, b, i * 2]
                        ax.plot(x, y, '.-', color=colors[b], label=baseRate)
                        ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                        if (i == 0):
                            ax.set_ylabel('Percent pass')
                        if (pltcnt > (numplots - numcols)):
                            ax.set_xlabel('True RP')
                        if (pltcnt == 1):
                            ax.set_title('Contamination: %.2f' % (contRate) + '%')
                        else:
                            ax.set_title('%.2f' % (contRate) + '%')
                        ax.set_ylim(-10, 110)
                        ax.spines.right.set_visible(spinesSetting)
                        ax.spines.top.set_visible(spinesSetting)

                # fig.text(0.65, 0.9-(.17*j), 'Recording Duration: %d hours'%recDur)
                if recDur != 1:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hours' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
                else:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hour' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
            # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            handles, labels = ax.get_legend_handles_labels()
            fig.tight_layout()
            plt.subplots_adjust(hspace=0.3)

            fig.legend(handles, labels, loc='upper right')
            fig.savefig(savefile + '_RP.svg', dpi=500)
        if plotType == 'paper':
            numrows = 1
            numcols = 1
            numplots = numrows * numcols
            fig, axs = plt.subplots(numrows, numcols, figsize=(4, 5))
            ax = axs  # for the case of just one subplot
            # plot just contRates 0.08 to 0.12:
            cr = params['contRates'];

            # plot just RP = 2:
            # rps = params['RPs']
            # rpInd = np.where(rps == 0.002)[0]
            #plot just  recDur = 2:
            recDurs = params['recDurs']
            rdInd = np.where(recDurs == 2)[0]

            # plot just fr = 5:
            frs = np.array(params['baseRates'])
            frInd = np.where(frs == 2)[0][0]
            print('Firing rate is 2')

            # colors = matplotlib.cm.Set1(np.linspace(0, 1, 10))
            c = cc.linear_bmw_5_95_c89#input_color  # cc.linear_protanopic_deuteranopic_kbw_5_95_c34
            c = c[::-1]  # if using linear_blue37 or protanopic, flip the order
            colors = [c[x] for x in np.round(np.linspace(0.2, 0.75, len(params['RPs'])) * 255).astype(int)]


            pltcnt = 0
            linewidths = [1,1,1,1,1,1]#[0.5, 1, 2, 3]
            for j, recDur in enumerate(recDurs[rdInd]):
                for i, rp in enumerate(params['RPs']):

                    # different base rates get different colors
                    # fix baseRate at frs[frInd]: previously did for b, baseRate in enumerate(frs[frInd]):
                    b = 0
                    baseRate = frs[frInd]# Todo change x axis
                    if zoomCont:
                        lowerCI = CI[0][rdInd[0], i, frInd, 7:14]
                        upperCI = CI[1][rdInd[0], i, frInd, 7:14]
                        x = cr[7:14] * 100
                        y = pc[rdInd[0], i, frInd, 7:14]

                    else:
                        lowerCI = CI[0][rdInd[0], i, frInd, :]
                        upperCI = CI[1][rdInd[0], i, frInd, :]
                        x = cr * 100
                        y = pc[rdInd[0], i, frInd, :]



                    ax.plot(x, y, '.-', color=colors[i], linewidth=linewidths[i], label=rp*1000)

                    ax.fill_between(x, lowerCI, upperCI, color=colors[j], alpha=.3)
                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('Contamination (%)')
                    # ax.set_title('True RP: %.1f ms; contamination: %d' % (rp * 1000, contRate * 100) + '%')
                    ax.set_ylim(-10, 110)
                    ax.spines.right.set_visible(spinesSetting)
                    ax.spines.top.set_visible(spinesSetting)
                # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
                # fig.suptitle('Proportion contamination: %.2f' % contRate*100, x=.5, y=1.1)
            handles, labels = ax.get_legend_handles_labels()

            # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
            fig.tight_layout()
            fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1), title='Refractory Period (ms)')
            fig.savefig(savefile + '_RP.svg', dpi=500)
            fig.savefig(savefile + '_RP.png', dpi=500)

    if subplot4: #zoom
        fig, axs = plt.subplots(1, 1, figsize=(3, 5))
        ax = axs  # for the case of just one subplot
        # plot just contRates 0.08 to 0.12:
        cr = params['contRates'];
        crInds = np.where((cr >= 0.07) & (cr <= 0.13))[0]
        print(crInds)
        # plot just recDur = 1
        rd = params['recDurs']
        rdInd = np.where(rd == 2)[0]
        # plot just RP = 2:
        rp = params['RPs']
        rpInd = np.where(rp == 0.0025)

        for j, recDur in enumerate(rd[rdInd]):
            for i, contRate in enumerate(rp[rpInd]):

                # different base rates get different colors
                for b, baseRate in enumerate(params['baseRates']):
                    lowerCI = CI[0][rdInd, rpInd, b, crInds][0]
                    upperCI = CI[1][rdInd, rpInd, b, crInds][0]
                    x = cr[crInds] * 100
                    y = pc[rdInd, rpInd, b, crInds][0]
                    ax.plot(x, y, '.-', color=colors[b], label=baseRate)
                    ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('Contamination (%)')
                    ax.set_ylim(-10, 110)
                    ax.spines.right.set_visible(spinesSetting)
                    ax.spines.top.set_visible(spinesSetting)
                    ax.set_title('True RP: %.1f ms; recording duration: %d hour' % (rp[rpInd] * 1000, rd[rdInd]))
                    ax.xaxis.set_ticks([7, 8, 9, 10, 11, 12, 13])
        # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
        handles, labels = ax.get_legend_handles_labels()
        plt.subplots_adjust(hspace=0.3)

        fig.savefig(savefile + '_zoom.svg', dpi=500)
        fig.savefig(savefile + '_zoom.pdf', dpi=500)


def plotSimulationsOverlay(pcDict,params,savefile, rpPlot=2.5,frPlot = 5,recDurPlot = 2,legendLabels = None,colorflag = False):
    '''


    Parameters
    ----------
    pcDict - a dictionary where each value is a matrix of percent correct for the simulations across all parameters
    params - the parameters used for running the simulations
    savefile - path to save the figures
    rpPlot - RP to use for plot, default is 2.5 ms.
    frPlot - FR to use for plot, default is 5 spks/s.
    recDurPlot - recording duration to use for plot, default is 2 hours.
    legendLabels - optional, np.array containing labels (str) for each set of simulations. e.g. different confidence
        levels or metrics run. Title of the plot is the last element of this array.
        default values here are the keys of pcDict.
    colorflag - optional, only set to True if needed for an error in the number of colors (per runs of the simulation)

    Returns
    saves the plot in 'savefile'

    '''

    #set up figure and plotting parameters
    h=0;l=0;po=0;k=0
    spinesSetting = False
    fig, axs = plt.subplots(1, 1, figsize=(4, 5))
    ax = axs  # for the case of just one subplot
    ax.vlines(10,0,100,'k','dashed',alpha=0.5)
    ax.hlines([0,20,40,60,80,100],0,20,'k','solid',alpha = 0.2)

    #each key of pcDict is the name of the run of the simulation (e.g. 'slidingRP', 'Llobet', 'slidingRP 70% conf')
    pcDict_keys = list(pcDict.keys())
    #this checks the type of variable in the keys
    #   to assign diferent colors for string names (e.g. Llobet, Hill) vs. confidence values (e.g. '70')
    varTypes = [type(pc_key) for pc_key in pcDict_keys]
    numStrings = len([v for v in varTypes if v is str])
    nConfs = len([v for v in varTypes if v == int])


    for p,pc_key in enumerate(pcDict_keys): #loop through the different simulation matrices

        # pc matrix to use for plotting simulations:
        pc = pcDict[pc_key]

        #compute error bars
        count = []
        count = pc / 100 * params['nSim']
        CI_scaled = binofit(count, params['nSim'])
        CI = [x * 100 for x in CI_scaled]

        # contamination rates to plot:
        cr = params['contRates'];

        # plot just RP = rpPlot:
        rps = params['RPs']
        rpInd = np.where(rps == rpPlot/1000)[0] #rpPlot in ms, convert to s here

        #plot just  recDur = recDurPlot:
        recDurs = params['recDurs']
        rdInd = np.where(recDurs == recDurPlot)[0]

        # plot just fr = frPlot:
        frs = np.array(params['baseRates'])
        frInd = np.where(frs == frPlot)[0][0]
        print('Firing rate is', str(frPlot))



        if type(pc_key) is str and pc_key[0]=='H': # Hill
            colors = np.flip(matplotlib.cm.Reds(np.linspace(0.3, 1, 3)),0)  #red for Hill
            color = colors[h]
            h+=1
        elif type(pc_key) is str and len(pc_key)>7 and pc_key[7] == 'P':  # purple for Llobet Poiss
            colors = np.flip(matplotlib.cm.Purples(np.linspace(0.3, 1, 3)),0)
            color = colors[po]
            po += 1
        elif type(pc_key) is str and pc_key[0] == 'L':  # green for Llobet
            colors = np.flip(matplotlib.cm.Greens(np.linspace(0.3, 1, 3)),0)
            color = colors[l]
            l+=1

        else: #not hill (conf)
            colors = np.flip(matplotlib.cm.Blues(np.linspace(0.3, 1, 3)),0) #blue for regular version of metric
            k+=1
            if colorflag:
                color = colors[p+2]
            else:
                color = colors[p]

        pltcnt = 0
        linewidths = [1,1,1,1,1,1]
        for j, recDur in enumerate(recDurs[rdInd]):
            for i, rp in enumerate(rps[rpInd]):

                # fix baseRate at frs[frInd]
                b = 0
                baseRate = frs[frInd]# Todo change x axis

                #confidence interval for this recDur, RP, FR
                lowerCI = CI[0][rdInd[0], rpInd[0], frInd, :]
                upperCI = CI[1][rdInd[0], rpInd[0], frInd, :]

                x = cr * 100 #convert contamination rate to percent (from fractional)
                y = pc[rdInd[0], rpInd[0], frInd, :] #percent correct values for this recDur, RP, FR

                #plot simulations, set labels and spines
                ax.plot(x, y, '.-', color=color, linewidth=linewidths[i], label=rp*1000)
                ax.fill_between(x, lowerCI, upperCI, color=color, alpha=.3)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('Contamination (%)')
                ax.set_ylim(-10, 110)
                ax.spines.right.set_visible(spinesSetting)
                ax.spines.top.set_visible(spinesSetting)

    handles, xx = ax.get_legend_handles_labels()

    if legendLabels is None:
        labels = pcDict_keys
        title = 'Contamination Metric'
    else:
        labels = legendLabels[0:-1]
        title = legendLabels[-1]

    #figure legend/title
    fig.tight_layout()
    fig.legend(handles,labels, bbox_to_anchor = (1.07, 1.0),title=title)
    fig.suptitle('RP = {}    FR = {}'.format(rpPlot,frPlot))

    #save and display figure
    fig.savefig(savefile + '_RP.svg', dpi=500)
    fig.savefig(savefile + '_RP.png', dpi=500)
    fig.show()


def plotDriftOverlay(pcDict,paramsDict,savefile, rpPlot=2.5,frPlotInput = 5,driftDir = 'Inc'):
    # (pcSliding,pcHill15,pcHill2,pcHill3,params,savefile, rpPlot=2.5):
    spinesSetting = False
    fig, axs = plt.subplots(1, 1, figsize=(6, 8))
    ax = axs  # for the case of just one subplot
    ax.vlines(10,0,100,'k','dashed',alpha=0.5)
    ax.hlines([0,20,40,60,80,100],0,20,'k','solid',alpha = 0.2)

    if driftDir == 'Inc':
        pcDict_keys = list(pcDict.keys())
        pcDict_keys.pop(0) #get rid of dec
    elif driftDir =='Dec':
        pcDict_keys = list(pcDict.keys()) #get rid of inc
        pcDict_keys.pop(2) #get rid of dec


    for p,pc_key in enumerate(pcDict_keys):#,pcHill3]):
        if pc_key == 0:
            if driftDir == 'Inc':
                frPlot = frPlotInput*(2+0.5)/2 #mean firing rate across session
                frMean = frPlot
            elif driftDir == 'Dec':
                frPlot = frPlotInput*(2-0.5)/2 #mean firing rate across session
                frMean = frPlot
        else:
            frPlot = frPlotInput
            frMean = frPlot
        pc = pcDict[pc_key]
        params = paramsDict[pc_key]
        print(pc_key)
        print(pc.shape)
        count = []
        print(params['nSim'])
        count = pc / 100 * params['nSim']  # number of correct trials

        CI_scaled = binofit(count, params['nSim'])
        CI = [x * 100 for x in CI_scaled]

        # plot just contRates 0.08 to 0.12:
        cr = params['contRates'];

        # plot just RP = rpPlot:
        rps = params['RPs']
        rpInd = np.where(rps == rpPlot/1000)[0] #rpPlot in ms, convert to s here
        #plot just  recDur = 2:
        recDurs = params['recDurs']
        rdInd = np.where(recDurs == 2)[0]

        # plot just fr = 5:
        frs = np.array(params['baseRates'])
        print(frs)
        print(frPlot)
        frInd = np.where(frs == frPlot)[0][0]
        print('Firing rate is 2')

        # colors = matplotlib.cm.Set1(np.linspace(0, 1, 10))
        c = cc.linear_bmw_5_95_c89#input_color  # cc.linear_protanopic_deuteranopic_kbw_5_95_c34
        c = c[::-1]  # if using linear_blue37 or protanopic, flip the order
        # if p==0:
        #     color = [c[x] for x in np.round(np.linspace(0.2, 0.75, len(params['RPs'])) * 255).astype(int)][5]
        # else:
        if pc_key !=0 and driftDir == 'Dec': # Hill
            colors = matplotlib.cm.Reds(np.linspace(0.3, 1, 3))
            color = colors[p-1]
        elif pc_key !=0 and driftDir == 'Inc': # Hill
            colors = matplotlib.cm.Greens(np.linspace(0.3, 1, 3))
            color = colors[p-1]

        else: #not hill (conf)
            colors = matplotlib.cm.Blues(np.linspace(0.3, 1, 3))
            color = colors[p-1]


        pltcnt = 0
        linewidths = [1,1,1,1,1,1]#[0.5, 1, 2, 3]
        for j, recDur in enumerate(recDurs[rdInd]):
            for i, rp in enumerate(rps[rpInd]):

                # different base rates get different colors
                # fix baseRate at frs[frInd]: previously did for b, baseRate in enumerate(frs[frInd]):
                b = 0
                baseRate = frs[frInd]# Todo change x axis


                if pc.ndim == 3:
                    print(CI[0].shape)
                    lowerCI = CI[0][rpInd[0], frInd, :]
                    upperCI = CI[1][rpInd[0], frInd, :]

                    x = cr * 100
                    y = pc[rpInd[0], frInd, :]
                else:
                    lowerCI = CI[0][rdInd[0], rpInd[0], frInd, :]
                    upperCI = CI[1][rdInd[0], rpInd[0], frInd, :]
                    x = cr * 100
                    y = pc[rdInd[0], rpInd[0], frInd, :]

                ax.plot(x, y, '.-', color=color, linewidth=linewidths[i], label=rp*1000)

                ax.fill_between(x, lowerCI, upperCI, color=color, alpha=.3)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('Contamination (%)')
                # ax.set_title('True RP: %.1f ms; contamination: %d' % (rp * 1000, contRate * 100) + '%')
                ax.set_ylim(-10, 110)
                ax.spines.right.set_visible(spinesSetting)
                ax.spines.top.set_visible(spinesSetting)
            # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
            # fig.suptitle('Proportion contamination: %.2f' % contRate*100, x=.5, y=1.1)
    handles, xx = ax.get_legend_handles_labels()

    labels = ['Flexible RP metric', 'Hill metric; threshold = 1.5 ms','Hill metric; threshold = 2 ms','Hill metric; threshold = 3 ms']
    # labels = ['50%', '75%','90%']#,'99%']
    # labels = ['10%', '10%','10%','10%']
    labels = pcDict_keys
    if driftDir == 'Inc':
        labels = ['Stable FR', driftDir + ' FR 50%']
    elif driftDir == 'Dec':
        labels = [ driftDir + ' FR 50%', 'Stable FR']

    #tsting just one
    # labels=labels[0]
    # handles=handles[0]
    # title = 'Confidence threshold'
    # title = 'Refractory Period (ms)'
    title = 'Drift'

    # fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)


    fig.tight_layout()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1),title=title )
    fig.legend(handles,labels, bbox_to_anchor = (1.05, 1.0), loc = 'upper left',title=title)
    fig.suptitle('Mean FR = {}'.format(frMean))
    fig.savefig(savefile + '_RP.svg', dpi=500)
    fig.savefig(savefile + '_RP.png', dpi=500)




def plotSensitivitySpecificity(pcOrig, pc2Ms, params, savefile, plusMinusThresh):
    """
    function to plot sensitivity and specificity functions (data exploration, not in paper)

        crRange = range(len(params['contRates']))
        crInd = np.where(params['contRates'] == 0.1])[0] #the cutoff point where contamination is exactly 10%
       #version using all contRates
        TP = pc[:, :, :, np.arange(0,crInd)]
        FP = pc[:,:,:,  np.arange(crInd+1, len(params['contRates']))]
        TN = 100 - pc[:,:,:,  np.arange(crInd+1, len(params['contRates']))]
        FN = 100- pc[:, :, :, np.arange(0,crInd)]

    #version using 5/15
    """
    # plusMinusThresh = 4
    minusThresh = 10 - plusMinusThresh
    plusThresh = 10 + plusMinusThresh

    recDurUse = 1
    RPuse = 2
    useAllParams = 1
    if useAllParams:
        figsize = (25, 20)
        recDurInds = np.arange(len(params['recDurs']))
        rpInds = np.arange(len(params['RPs']))
    else:
        figsize = (5, 5)
        recDurInds = np.where(params['recDurs'] == recDurUse)[0]
        rpInds = np.where(params['RPs'] == RPuse / 1000)[0]

    fig, axs = plt.subplots(len(recDurInds), len(rpInds), figsize=figsize)

    for p, pc in enumerate([pcOrig, pc2Ms]):
        if p == 0:
            linestyle = '-'
        else:
            linestyle = '--'
        for j, recDur in enumerate(params['recDurs'][recDurInds]):
            for i, rp in enumerate(params['RPs'][rpInds]):
                print(i)
                print(rp)
                # fig,axs = plt.subplots(1,1,figsize = (5,5))
                if np.size(axs) == 1:
                    ax = axs  # for the case of just one subplot
                else:
                    ax = axs[j, i]

                if False:
                    recDurChosen = 1
                    rdInd = np.where(params['recDurs'] == recDurChosen)[0][0]

                    rpChosen = 2
                    rpInd = np.where(params['RPs'] == rpChosen / 1000)[0][0]
                recDurChosen = recDur
                rdInd = j

                rpChosen = rp * 1000
                rpInd = i

                TP = pc[rdInd, rpInd, :, minusThresh]
                FP = pc[rdInd, rpInd, :, plusThresh]
                TN = 100 - pc[rdInd, rpInd, :, plusThresh]
                FN = 100 - pc[rdInd, rpInd, :, minusThresh]

                sensitivity = TP / (TP + FN)  # ability to correctly detect neurons that should pass
                specificity = TN / (TN + FP)  # ability to correctly reject neurons that should fail

                PerCorr = ((pc[rdInd, rpInd, :, minusThresh] + (100 - pc[rdInd, rpInd, :, plusThresh])) / 2) / 100

                ax.plot(params['baseRates'], sensitivity, linestyle=linestyle, marker='.', color='green',
                        label='sensitivity')
                ax.plot(params['baseRates'], specificity, linestyle=linestyle, marker='.', color='purple',
                        label='specificity')
                ax.plot(params['baseRates'], PerCorr, linestyle=linestyle, marker='.', color='red',
                        label='percent corrects')
                ax.set_xscale('log')
                ax.set_xlabel('Base firing rate (spk/s)')
                if j == 4:
                    ax.set_xlabel('Base firing rate (spk/s)')
                if i % 7 == 0:
                    ax.set_ylabel('Sensitivity or specificity (proportion)')
                ax.set_title('recDur %d hr, RP %d ms  |  %d / %d ' % (recDurChosen, rpChosen, minusThresh, plusThresh))
        # fig.legend(bbox_to_anchor=(.9,.8))
    fig.show()
    print('plotting sensitivity specificity')
    fig.savefig(savefile + 'SpecSens%dRD_%dRP_%d.svg' % (recDurChosen, rpChosen, plusMinusThresh), dpi=500)
    # fig.savefig(r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs\sensitivity_specificity_%dRD_%dRP_%d.png'%(recDurChosen,rpChosen,plusMinusThresh), dpi=500)