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