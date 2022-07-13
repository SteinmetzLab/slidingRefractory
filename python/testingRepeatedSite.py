# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 17:29:46 2022

@author: Noam Roth
run slidingRP on repeated site and compare to previous version of metrics
"""


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


from phylib.stats import correlograms
from RP_plotting import *
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 #sets fonts as editable

#pull down RS data
from reproducible_ephys_functions import filter_recordings, labs, BRAIN_REGIONS, query, get_insertions
from reproducible_ephys_functions import combine_regions, BRAIN_REGIONS, get_insertions, save_data_path, save_dataset_info

from brainbox.io.one import SpikeSortingLoader
from ibllib.atlas import AllenAtlas
ba = AllenAtlas()

insertions = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSmetrics\\'
#set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = True


all_recordings = []
for iIns, ins in enumerate(insertions):

    try:
        print(f'processing {iIns + 1}/{len(insertions)}')
        eid = ins['session']['id']
        probe = ins['probe_name']
        pid = ins['probe_insertion']
    
        data = {}
    
        # Load in spikesorting
        sl = SpikeSortingLoader(eid=eid, pname=probe, one=one, atlas=ba)
        spikes, clusters, channels = sl.load_spike_sorting()
        clusters = sl.merge_clusters(spikes, clusters, channels)

    
  
    
        clusters['rep_site_acronym'] = combine_regions(clusters['acronym'])
        # Find clusters that are in the repeated site brain regions
        cluster_idx = np.sort(np.where(np.isin(clusters['rep_site_acronym'], BRAIN_REGIONS))[0])
        
        #save the 'old' (computed in pipeline) slidingRP metric value
        oldMetricValue = [clusters['slidingRP_viol'][x] for x in cluster_idx]
        
        #now rerun code in steinmetzlabrepo
        print('Loading spike times and clusters...')
        cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
        spikeTimes = spikes.times[cluInds]
        spikeClusters = spikes.clusters[cluInds]
        
        print('Running RP metric...')
        rpMetrics, cont, rp = slidingRP_all(spikeTimes, spikeClusters, params)
        
        rpMetrics['oldMetricValue'] = oldMetricValue
        
        
        with open(savefile + ins['session']['subject'] + '.pickle', 'wb') as handle:
            pickle.dump(rpMetrics, handle)
        
    except:
        print('error')

#%%

# compare old metric value to new for one session 
print('Comparing old and  new RP values...')
passpass = 0
failpass = 0
failfail = 0
passfail = 0
nClusters = len(rpMetrics['oldMetricValue'])
for c in range(nClusters):
    if rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==1:
        passpass +=1
    elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==1:
        failpass +=1
    elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==0:
        failfail +=1
    elif rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==0:
        passfail +=1
        
        
plt.bar([1,2,3,4], np.array([passpass,failfail,passfail, failpass])/nClusters)
plt.title('Session #1: 639 clusters: 86.5% remained and 13.5% changed')
plt.ylabel('Proportion of clusters')
plt.xticks([1, 2,3,4], ['Pass/Pass', 'Fail/Fail', 'Pass/Fail','Fail/Pass'],
       rotation=20)


#%% 
# compare old metric value to new one for all sessions 
nSess = len(insertions)
plotEach = False
passpassAll = 0
failfailAll = 0
failpassAll = 0
passfailAll = 0
nClustersAll = 0
for s in range(nSess):
    subject = insertions[s]['session']['subject']
    #load saved rpMetrics
    file = open(savefile + subject + '.pickle','rb')
    rpMetrics = pickle.load(file)
    file.close()
    
    print('Comparing old and  new RP values for session %d / %d'%(s, nSess))
    passpass = 0
    failpass = 0
    failfail = 0
    passfail = 0
    nClusters = len(rpMetrics['oldMetricValue'])
    nClustersAll += nClusters
    for c in range(nClusters):
        if rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==1:
            passpass +=1
            passpassAll +=1
        elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==1:
            failpass +=1
            failpassAll +=1
        elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==0:
            failfail +=1
            failfailAll +=1
        elif rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==0:
            passfail +=1
            passfailAll +=1
            
      
    if plotEach:
        fig,ax = plt.subplots(1,1)
        ax.bar([1,2,3,4], np.array([passpass,failfail,passfail, failpass])/nClusters)
        pRemain = (passpass + failfail) / nClusters *100 
        pChange = (passfail + failpass) / nClusters *100
        ax.set_title('Session %s; %s;  %d clusters; %.2f%% remained and %.2f%% changed'%(s, subject, nClusters, pRemain, pChange ))
        ax.set_ylabel('Proportion of clusters')
        ax.set_xticks([1, 2,3,4], ['Pass/Pass', 'Fail/Fail', 'Pass/Fail','Fail/Pass'],
               rotation=20)
        
    fig,ax = plt.subplots(1,1)
    ax.bar([1,2,3,4], np.array([passpassAll, failfailAll, passfailAll, failpassAll])/nClustersAll)
    pRemainAll = (passpassAll + failfailAll) / nClustersAll *100 
    pChangeAll = (passfailAll + failpassAll) / nClustersAll *100
    ax.set_title('All sessions:  %d clusters; %.2f%% remained and %.2f%% changed'%(nClustersAll, pRemainAll, pChangeAll ))
    ax.set_ylabel('Proportion of clusters')
    ax.set_xticks([1, 2,3,4], ['Pass/Pass', 'Fail/Fail', 'Pass/Fail','Fail/Pass'],
           rotation=20)
        
    
#%%

#Look into why neurons go from pass to fail and why some neurons go from fail to pass
nSess = len(insertions)

for s in range(nSess)[0:1]:
    subject = insertions[s]['session']['subject']
    #load saved rpMetrics
    file = open(savefile + subject + '.pickle','rb')
    rpMetrics = pickle.load(file)
    file.close()
    
    print(f'processing {s + 1}/{len(insertions)}')
    eid = insertions[s]['session']['id']
    probe = insertions[s]['probe_name']
    pid = insertions[s]['probe_insertion']


    # Load in spikesorting
    sl = SpikeSortingLoader(eid=eid, pname=probe, one=one, atlas=ba)
    spikes, clusters, channels = sl.load_spike_sorting()
    clusters = sl.merge_clusters(spikes, clusters, channels)


  

    clusters['rep_site_acronym'] = combine_regions(clusters['acronym'])
    # Find clusters that are in the repeated site brain regions
    cluster_idx = np.sort(np.where(np.isin(clusters['rep_site_acronym'], BRAIN_REGIONS))[0])
    
    #save the 'old' (computed in pipeline) slidingRP metric value
    oldMetricValue = [clusters['slidingRP_viol'][x] for x in cluster_idx]
    
    #now rerun code in steinmetzlabrepo
    print('Loading spike times and clusters...')
    cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
    spikeTimes = spikes.times[cluInds]
    spikeClusters = spikes.clusters[cluInds]
    
    #%%
    
from brainbox.singlecell import firing_rate, acorr
from phylib.stats import correlograms

subject = insertions[s]['session']['subject']
#load saved rpMetrics
file = open(savefile + subject + '.pickle','rb')
rpMetrics = pickle.load(file)
file.close()    
 
timeBins = np.arange(1/sampleRate, 1, 1/sampleRate)
nClusters = len(rpMetrics['oldMetricValue'])
for cluster_id in range(nClusters)[0:10]:   
    c = rpMetrics['cidx'][cluster_id]
    st = spikeTimes[spikeClusters==c] 
    spikeCount = len(st)
    sampleRate = 30000
    
    [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
    nSpikesBelow2, confMatrix, cont, rp, nACG,
    firingRate,xx] = slidingRP(st, params)
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize = (12,4))
    
    ax = axs[0]
    histWin = 0.5 #seconds
    binnedFR = firing_rate(st,hist_win = histWin, fr_win =  10)
    ax.plot(np.arange(0,recDur,histWin),binnedFR, 'k')
    recDur = st[-1]-st[0]
    firingRateSC = spikeCount / recDur
    rpOld = bool(rpMetrics['oldMetricValue'][cluster_id])
    rpNew = bool(rpMetrics['value'][cluster_id])
    
    if rpOld:
        rpOld = 'PASS'
    else:
        rpOld = 'FAIL'
    if rpNew:
        rpNew= 'PASS'
    else:
        rpNew = 'FAIL' 
        
    ax.set_title('FR from ACG = %.2f   FR from SC = %.2f'%(firingRate, firingRateSC))
    ax.set_ylabel('FR (spks/s)')

    ax = axs[1]
    ax.plot(timeBins,nACG[1:-1],'k')
    ax.set_title('old:%s, new:%s'%(rpOld, rpNew))
    ax.set_ylabel('ACG count (spks)')
    ax.set_xlabel('Time (s)')



#%%
#TODO: also check what about passing/failing in general!! i.e. including other 2 single unit quality metrics!

# compare old metric value to new one for all sessions 
nSess = len(insertions)
plotEach = False
passpassAll = 0
failfailAll = 0
failpassAll = 0
passfailAll = 0
nClustersAll = 0
for s in range(nSess):
    subject = insertions[s]['session']['subject']
    #load saved rpMetrics
    file = open(savefile + subject + '.pickle','rb')
    rpMetricsUpdated = pickle.load(file)
    file.close()
    
    print('Comparing old and  new RP values for session %d / %d'%(s, nSess))
    passpass = 0
    failpass = 0
    failfail = 0
    passfail = 0
    nClusters = len(rpMetrics['oldMetricValue'])
    nClustersAll += nClusters
    for c in range(nClusters):
        if rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==1:
            passpass +=1
            passpassAll +=1
        elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==1:
            failpass +=1
            failpassAll +=1
        elif rpMetrics['oldMetricValue'][c] == 0 and rpMetrics['value'][c] ==0:
            failfail +=1
            failfailAll +=1
        elif rpMetrics['oldMetricValue'][c] == 1 and rpMetrics['value'][c] ==0:
            passfail +=1
            passfailAll +=1
            
      
    if plotEach:
        fig,ax = plt.subplots(1,1)
        ax.bar([1,2,3,4], np.array([passpass,failfail,passfail, failpass])/nClusters)
        pRemain = (passpass + failfail) / nClusters *100 
        pChange = (passfail + failpass) / nClusters *100
        ax.set_title('Session %s; %s;  %d clusters; %.2f%% remained and %.2f%% changed'%(s, subject, nClusters, pRemain, pChange ))
        ax.set_ylabel('Proportion of clusters')
        ax.set_xticks([1, 2,3,4], ['Pass/Pass', 'Fail/Fail', 'Pass/Fail','Fail/Pass'],
               rotation=20)
        
    fig,ax = plt.subplots(1,1)
    ax.bar([1,2,3,4], np.array([passpassAll, failfailAll, passfailAll, failpassAll])/nClustersAll)
    pRemainAll = (passpassAll + failfailAll) / nClustersAll *100 
    pChangeAll = (passfailAll + failpassAll) / nClustersAll *100
    ax.set_title('All sessions:  %d clusters; %.2f%% remained and %.2f%% changed'%(nClustersAll, pRemainAll, pChangeAll ))
    ax.set_ylabel('Proportion of clusters')
    ax.set_xticks([1, 2,3,4], ['Pass/Pass', 'Fail/Fail', 'Pass/Fail','Fail/Pass'],
           rotation=20)
        



