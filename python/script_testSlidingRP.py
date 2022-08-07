# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:58:20 2022

@author: Noam Roth

Script for running and testing slidingRP functions
"""

#imports

import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde

import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\paper-reproducible-ephys')
sys.path.append(r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\metrics\slidingRP')
sys.path.append(r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python')
#from reproducible_ephys_functions import query
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
from one.api import ONE
from slidingRP import *
import pickle
one = ONE()

#%% load data

params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
# params = SimpleNamespace(**params) #convert to dot notation
params['returnMatrix'] = True
params['verbose'] = True

#%%
# mat = computeMatrix(spikeTimes, params)
#load sample data
datapath = r'E:\Hopkins_CortexLab'
nickname = 'Hopkins'
spikeTimes = np.load(datapath + '\\spike_times.npy').flatten() 
if nickname == 'Hopkins':
    #convert from samples to seconds
    spikeTimes = spikeTimes/params['sampleRate']

spikeClusters = np.load(datapath + '\\spike_clusters.npy')

#%%
#run slidingRP for one cluster
params['cidx'] = [0]
st = spikeTimes[spikeClusters == params['cidx'][0]]
print(len(st))
[maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
nSpikesBelow2, confMatrix, cont, rp, nACG,
firingRate,xx] = slidingRP(st, params)




 #%%
#run plotting code for one cluster
#for RE presentation: ran below plus plotSlidingRP with a temporary "bug"

params['cidx'] = [0]
st = spikeTimes[spikeClusters == params['cidx'][0]]
plotSlidingRP(st, params)
from brainbox.singlecell import firing_rate, acorr

hist_win = 0.05
fr_win = 1
if(1):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,4))
    st = st[:-500]
    fr = firing_rate(st, hist_win = hist_win, fr_win=fr_win)
    x = np.arange(fr.size) * hist_win
    ax.plot(x, fr, color='k', linewidth=0.4)
    ax.set_title('Firing Rate')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Rate (spks/s)')



#%%
#run slidingRP for the loaded recording

[rpMetrics, cont, rp, sc, frrd] = slidingRP_all(spikeTimes, spikeClusters, params = params)
#%%


fig,ax = plt.subplots(1,1)
ax.scatter(sc[0:10],frrd[0:10])
ax.set_xlabel('spikecount')
ax.set_ylabel('firingrate * recDur')
xlims = ax.get_xlim()
ax.plot( [0,xlims[1]],[0,xlims[1]] )

#%%
#look into neurons where spike count is greater than firing rate * recDur 
# these should be neurons that now pass (spike count) and previously failed  (fr * recDur)
candidateNeurs = np.where((sc-frrd)>0)[0]

oldVal = np.array([rpMetricsUpdated['oldMetricValue'][x] for x in candidateNeurs])
newVal = np.array([rpMetricsUpdated['value'][x] for x in candidateNeurs])

print(newVal-oldVal) # there is only one such neuron!

#%%
# try it the other way:
oldVal = np.array([rpMetricsUpdated['oldMetricValue'][x] for x in range(len(rpMetricsUpdated['oldMetricValue']))])
newVal = np.array([rpMetricsUpdated['value'][x] for x in range(len(rpMetricsUpdated['oldMetricValue']))])

idxFailPass = np.where(newVal-oldVal > 0 )[0]
fig,ax = plt.subplots(1,1)
ax.scatter(sc[idxFailPass], frrd[idxFailPass])
ax.set_xlabel('spikecount')
ax.set_ylabel('firingrate * recDur')
xlims = ax.get_xlim()
ax.plot( [0,xlims[1]],[0,xlims[1]] )


#%% script for testing sigmoid fits
estimatedRP, estimateIdx, xSigmoid, ySigmoid = fitSigmoidACG(nACG, rp, params)

#%%
fig,axs = plt.subplots(1,1,figsize = (6,5))
ax = axs
acg = plotSigmoid(ax, nACG, xSigmoid, ySigmoid, estimateIdx, estimatedRP)


#%%
rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters, rp, params)