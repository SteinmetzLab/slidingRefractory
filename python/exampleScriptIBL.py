# -*- coding: utf-8 -*-
"""
Created on Wed Aug 3 16:29:20 2022

@author: Noam Roth

example script to compute sliding refractory metric for one IBL session 
"""

import numpy as np
import matplotlib.pyplot as plt
from reproducible_ephys_functions import BRAIN_REGIONS, get_insertions, combine_regions
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from one.api import ONE
from slidingRP import *
import pickle
one = ONE()

from brainbox.io.one import SpikeSortingLoader
from ibllib.atlas import AllenAtlas
ba = AllenAtlas()

insertionsAll = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')
saveLocally = False
# savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSmetrics\\'

#set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = True
#%%
insertions = insertionsAll[0:1] #for this example, just run for a sample session 
all_recordings = []
for iIns, ins in enumerate(insertions):

    print(f'processing {iIns + 1}/{len(insertions)}')
    eid = ins['session']['id']
    probe = ins['probe_name']
    pid = ins['probe_insertion']

    # Load in spikesorting
    sl = SpikeSortingLoader(eid=eid, pname=probe, one=one, atlas=ba)
    spikes, clusters, channels = sl.load_spike_sorting()
    clusters = sl.merge_clusters(spikes, clusters, channels)

    #Find and load clusters from this insertion (in RS regions)
    clusters['rep_site_acronym'] = combine_regions(clusters['acronym'])
    # Clusters that are in the repeated site brain regions
    cluster_idx = np.sort(np.where(np.isin(clusters['rep_site_acronym'], 
                                           BRAIN_REGIONS))[0])

    print('Loading spike times and clusters...')
    cluInds = np.array(list(x for x in range(len(spikes.clusters))
                            if spikes.clusters[x] in cluster_idx)) 
    spikeTimes = spikes.times[cluInds]
    spikeClusters = spikes.clusters[cluInds]
    
    print('Running RP metric...')
    rpMetrics, cont, rp = slidingRP_all(spikeTimes, spikeClusters, params)
    
    
    if saveLocally:
        with open(savefile + ins['session']['subject'] + '.pickle', 'wb') as handle:
            pickle.dump(rpMetrics, handle)
        