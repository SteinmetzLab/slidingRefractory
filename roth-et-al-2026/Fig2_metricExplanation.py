#imports
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pickle
from slidingRP.metrics import slidingRP, plotSlidingRP, computeACG, computeViol
from slidingRP.simulations import genST



def runSaveFig2(figsavepath,resultsBasePath,savedSimNeuronFlag=False):
    #run and save the plots for the main figure 2

    # Simulate neuron and run the metric on it, plotting a schematic of the steps of the metric

    #possible to save combST, if savedSimNeuronFlag = True
    if savedSimNeuronFlag: #load a saved simulated neuron
        #path to load the simulated neuron
        resultsPath = resultsBasePath + '\\savedExampleSimulatedNeuron.pickle'

        file = open(resultsPath, 'rb')
        combST = pickle.load(file)
        file.close()
    else:
        #simulate a new contaminated neuron:

        #first, simulate a spike train for an uncontaminated example neuron
        firingRate = 10 #spks/s
        recDur = 3600 #seconds, i.e. a 1 hour recording
        rp = 2.5 /1000 #seconds, i.e. a 2.5 ms refractory period

        st = genST(firingRate, recDur, rp, params=None)

        #now add contamination
        contaminationProp = .1 #fractional contamination, i.e. 10% contamination
        contRate = firingRate * contaminationProp #rate of the contaminating source
        contST = genST(contRate, recDur, 0, params=None)  # add contaminating source, has 0 rp because may come from multiple neurons
        combST = np.sort(np.concatenate((st, contST)))  # put spike times from both spike trains together (and sort chronologically)


        # for this example neuron, save it to file to be able to reload
        exampleNeuronFilename = resultsBasePath + '\\savedExampleSimulatedNeuron2.pickle'
        with open(exampleNeuronFilename, 'wb') as handle:
            pickle.dump(combST, handle, protocol=pickle.HIGHEST_PROTOCOL)



    #parameters needed for running metric and plotting ACG
    params = { 'contaminationThresh': 10,
               'clusterLabel': False,
               'sampleRate': 30000,
               'binSizeCorr': 1/30000,
               'confidenceThresh': 90,
               'savefig': False,
               'recDur':3600
               }

    # create a figure with 2 rows and nColumns columns (set to 3 examples for paper)
    nColumns = 3
    fig, axs = plt.subplots(2,nColumns)
    verbose = False
    nACG = plotFig2(combST,1,axs,params,columnValue=0,color = 'steelblue',color2='navy',verbose=verbose)
    nACG = plotFig2(combST,1.5,axs,params,columnValue=1,color='forestgreen',color2 = 'darkgreen',verbose=verbose)
    nACG = plotFig2(combST,2.2,axs,params,columnValue=2,color='orchid',color2='darkorchid',verbose=verbose)
    fig.show()

    print(figsavepath)
    fig.savefig(os.path.join(figsavepath, 'shiftSchematic.svg'), dpi=300)
    fig.savefig(os.path.join(figsavepath, 'shiftSchematic.pdf'), dpi=300)
    fig.savefig(os.path.join(figsavepath, 'shiftSchematic.png'), dpi=300)


    #add subplots b,c, run on this same simulated neuron

    #mark X values corresponding to the example tauR's above
    XsToPlot = [1,1.5,2.2] #tauR's to plot
    Xcolors = ['steelblue','forestgreen','orchid'] #corresponding colors for each tauR value
    plotXs = [XsToPlot,Xcolors] #re-format for input to plotSlidingRP below

    #generate a figure object and format inputs to plotSlidingRP to save the figure
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize = (10,4))
    params['savefig']= True
    params['figpath'] = os.path.join(figsavepath, 'traceAndMatrix')

    plotSlidingRP(combST,params,plotXs = plotXs,inputAxes=(fig,axs),plotExtraContours=True)



def plotFig2(spikeTimes, tauR,axs,params,columnValue=0,color = 'b',color2 = 'darkblue',verbose=False):
    # inputs:

    #spikeTimes is an np.array of spike times (in ms)
    #tauR is a float, time for which to plot the ACG/distribution (ms)
    #axs is a created axes object, which axis within n columns will be indexed by columnValue if specified.
    # params should include keys 'binSizeCorr' for the ACG and 'sampleRate'. Otherwise this is set by default to every bin at 30kHz
    #color is for the main label text and lines in the plots
    #color2 is for the shaded areas and text to label shading
    #verbose indicates whether to have extra text and labels, or make this figure formally for the paper

    refDur = tauR / 1000 #convert to seconds
    spikeTimes = np.asarray(spikeTimes)
    n_spikes = len(spikeTimes)

    recDur = params['recDur']
    rpBinSize = params['binSizeCorr']

    # ACG via the same histdiff-equivalent the metric uses (no phylib needed);
    # firing rate from the explicit recording duration.
    rpEdges = np.arange(0, 10 / 1000, rpBinSize)  # in s
    rp = rpEdges + rpBinSize / 2  # refractory period durations to test (bin centres)
    nACG = computeACG(spikeTimes, rpBinSize, rp.size)
    firingRate = n_spikes / recDur

    # values for plot labels: expected violations under the contamination
    # threshold (Llobet), and the observed violations in the ACG up to tauR.
    contaminationProp = params['contaminationThresh'] / 100
    N_t = n_spikes
    N_c = contaminationProp * N_t           # contaminating spikes under the threshold
    N_b = N_t - N_c                         # base-neuron spikes

    #find rp closest to tauR
    rpInd = np.argmin(abs(rp - tauR / 1000))
    if rpInd == len(rp) - 1:  # edge case: last rp tested
        rpInd = rpInd - 1

    expectedViol = 2 * refDur / recDur * N_c * (N_b + (N_c - 1) / 2)  # Llobet et al.
    obsViol = int(np.sum(nACG[rp < refDur]))


    #plot this in two plots (forming one column of the axes):

    # top plot of the column is a bar plot of the ACG
    ax = axs[0,columnValue]
    ax.bar(rp*1000, nACG[0:len(rp)], width = np.diff(rp)[0]*1000, color = 'k', edgecolor = (1, 0, 0, 0))
    ax.set_xlim([0, 5])
    ax.set_xlabel('Time from spike (ms)')
    ax.set_ylabel('ACG count (spks)')
    ax.axvline(x=rp[rpInd]*1000,ymin=0,ymax=1,color=color2)
    t1 = ('Cluster FR = %.2f' %firingRate) #label with firing Rate

    if verbose:
        ax.set_title('Observed Violations for tauR = {0}: {1}'.format(tauR,obsViol) +
                 '\n Expected Violations for this neuron if it had 10% contamination: {0:.{1}f} '.format(expectedViol,2))
    else:
        ax.text(0.25,32,'{}'.format(obsViol),color = color)
        ax.text(tauR,35,str(tauR),color=color2)
        ax.text(3.5, 30, ' {0:.{1}f} '.format(expectedViol,2), wrap=True, horizontalalignment='center', fontsize=12,color=color2)

    ax.fill(np.array([0, 1, 1, 0])*tauR, np.array([0, 0, 1, 1])*ax.get_ylim()[1], color=color,alpha= 0.2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # bottom plot of the column is the probability distribution for this tauR
    ax = axs[1,columnValue]
    xvalues = np.arange(0, 3*obsViol) #plot for value that sets all 3 columns the same

    computedConf = 1 - stats.poisson.cdf(obsViol, expectedViol)

    ax.plot(xvalues,stats.poisson.pmf(xvalues,expectedViol),color=color2)
    ax.fill_between(xvalues[0:obsViol],np.zeros(len(xvalues[0:obsViol])),stats.poisson.pmf(xvalues[0:obsViol],expectedViol),color=color,alpha=0.2)
    ax.set_xlim(50,250)
    ax.set_ylim(0,.05)
    if verbose:
        ax.set_title('Probability of observing x spikes before tauR ms: {0:.{1}}'.format(computedConf,2))
    else:
        ax.set_title('{0:.{1}}'.format(computedConf,2))
    ax.set_xlabel('Number of spikes')
    ax.set_ylabel('Probability')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()



