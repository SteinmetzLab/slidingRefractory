import numpy as np
import matplotlib.pyplot as plt
from slidingRP.metrics import slidingRP, plotSlidingRP
from slidingRP.simulations import *
from scipy import stats

#%%

#first, simulate a spike train for a fake example neuron
firingRate = 10 #spks/s
recDur = 3600 #seconds, i.e. a 1 hour recording
rp = 2.5 /1000 #seconds, i.e. a 2.5 ms refractory period

st = genST(firingRate, recDur, rp, params=None)

#now add contamination
contaminationProp = .05 #fractional contamination, i.e. 10% contamination
contRate = firingRate * contaminationProp #rate of the contaminating source
contST = genST(contRate, recDur, 0, params=None)  # add contaminating source, has 0 rp because may come from multiple neurons
combST = np.sort(np.concatenate((st, contST)))  # put spike times from both spike trains together (and sort chronologically)


#plot ACG
params = { 'contaminationThresh': 10,
           'clusterLabel': False,
           'sampleRate': 30000,
           'binSizeCorr': 1/30000,
           'confidenceThresh': 90,
           'savefig': False
           }
plotSlidingRP(combST,params)

fig,axs = plt.subplots(1,1)

tauR = 1.5 #ms

plotFig2(combST,tauR,params)
#%%
def plotFig2(spikeTimes, tauR,params):
    refDur = tauR / 1000
    n_spikes = len(spikeTimes)
    clustersIds = [0]

    spikeClustersACG = np.zeros(n_spikes, dtype=np.int8)  # each spike time gets cluster id 0

    # compute an acg in 1s bins to compute the firing rate
    nACG = correlograms(spikeTimes, spikeClustersACG, cluster_ids=clustersIds, bin_size=1, sample_rate=params['sampleRate'],
                        window_size=2, symmetrize=False)[0][0]  # compute acg
    firingRate = nACG[1] / n_spikes
    contaminationRate = firingRate * contaminationProp

    nACG = correlograms(spikeTimes, spikeClustersACG, cluster_ids=clustersIds, bin_size=params['binSizeCorr'],
                        sample_rate=params['sampleRate'], window_size=2, symmetrize=False)[0][0]  # compute acg
    #
    # confMatrix = 100 * computeViol(np.cumsum(nACG[0:rp.size])[np.newaxis, :], firingRate, n_spikes,
    #                                rp[np.newaxis, :] + binSizeCorr / 2, cont[:, np.newaxis] / 100, recDur)

    # the number of violations (spikes) we expect to see under this contamination rate
    # expectedViol = contaminationRate * refDur * 2 * spikeCount
    # as computed *not* in the same way as originally taken from Hill (see above commented out)
    N_t = firingRate * recDur  # total number of spikes
    N_c = contaminationRate * recDur  # total number of contaminating spikes you would expect under the inputted CR
    N_b = N_t - N_c  # the "base" number of spikes under the inputted CR
    refDur = tauR/1000

    rpEdges = np.arange(0, 10 / 1000, params['binSizeCorr'])  # in s
    rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2  # vector of refractory period durations to test

    #find rp closest to tauR
    rpInd = np.argmin(abs(rp-tauR/1000))
    if rpInd == len(rp)-1 : #edge case: last rp tested
        rpInd = rpInd-1
    print(rp[rpInd])


    # expectedViol = 2 * refDur * 1/recDur * N_c * (N_b + (N_c - 1)/2) #number of expected violations, as defined in Llobet et al.
    expectedViol = 2 * refDur * 1 / recDur * N_c * (
                N_b + N_c / 2)  # number of expected violations, as defined in Llobet et al.

    obsViol = np.cumsum(nACG[0:rpInd])[-1]




    fig, axs = plt.subplots(2,1)

    ax = axs[0]

    ax.bar(rp*1000, nACG[0:len(rp)], width = np.diff(rp)[0]*1000, color = 'k', edgecolor = (1, 0, 0, 0))
    ax.set_xlim([0, 5])
    ax.set_xlabel('Time from spike (ms)')
    ax.set_ylabel('ACG count (spks)')
    ylim = ax.get_ylim()
    ax.axvline(x=rp[rpInd]*1000,ymin=0,ymax=1,color='b')
    t1 = ('Cluster FR = %.2f' %firingRate)
    ax.set_title('Observed Violations for tauR = {0}: {1}'.format(tauR,obsViol) +
                 '\n Expected Violations for this neuron: {0:.{1}f} '.format(expectedViol,2))


    ax.fill(np.array([0, 1, 1, 0])*tauR, np.array([0, 0, 1, 1])*ax.get_ylim()[1], 'b',alpha= 0.2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.show()









    # got rid of the -1, not sure whether that's right but wanted to be consistent with llobet

# the confidence that this neuron is contaminated at a level less than contaminationProp, given the number of true
# observed violations and under the assumption of Poisson firing
# confidenceScore = 1 - stats.poisson.cdf(obsViol, expectedViol)

nACG = plotFig2(combST,1.5,params)
