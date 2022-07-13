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
nClusters = len(oldMetricValue)
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
print('  %d: %s max conf = %.2f%%, min cont = %.1f%%, time = %.2f ms, n below 2 ms = %d' % (cids[cidx], pfstring, maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont*1000, nSpikesBelow2))

