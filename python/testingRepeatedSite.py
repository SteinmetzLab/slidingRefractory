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

#set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = True


all_recordings = []
for iIns, ins in enumerate(insertions[0:1]):

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
        

        with open('filename.pickle', 'wb') as handle:
            pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    except:
        print('error')

#%%
print('Comparing old and  new RP values...')
#now compare old metric value to new
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
        data['cluster_ids'] = clusters['cluster_id'][cluster_idx]
        data['RPmetric'] = clusters['cluster']
        
        #Todo: better way to do this?
        x = []
        for cluster in cluster_idx:
            x.append(list(spikes['times'][np.where(spikes['clusters'] == cluster)[0]]))
        clusters['ts'] = x
        # # Find spikes that are from the clusterIDs
        # spike_idx = np.isin(spikes['clusters'], data['cluster_ids'])
        # if np.sum(spike_idx) == 0:
        #     continue
        all_recordings.append(clusters)
            
    except Exception as err:
        print(f'{pid} errored: {err}')