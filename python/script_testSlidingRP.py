# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:58:20 2022

@author: Noam Roth

Script for running and testing slidingRP functions
"""


#%% load data

params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
# params = SimpleNamespace(**params) #convert to dot notation
params['returnMatrix'] = True
params['verbose'] = True

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

[maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
nSpikesBelow2, confMatrix, cont, rp, nACG,
firingRate,xx] = slidingRP(st, params)


 #%%
#run plotting code for one cluster

params['cidx'] = [0]
st = spikeTimes[spikeClusters == params['cidx'][0]]
plotSlidingRP(st, params)


#%%
#run slidingRP for the loaded recording

slidingRP_all(spikeTimes, spikeClusters, params = params)




