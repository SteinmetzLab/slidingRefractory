# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 22:10:45 2022

@author: Noam Roth

fit estimate histograms

"""


import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde


#%% IBL 

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
# from RP_plotting import *
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 #sets fonts as editable

#pull down RS data
from reproducible_ephys_functions import filter_recordings, labs, BRAIN_REGIONS, query, get_insertions
from reproducible_ephys_functions import combine_regions, BRAIN_REGIONS, get_insertions, save_data_path, save_dataset_info



from brainbox.io.one import SpikeSortingLoader
from ibllib.atlas import AllenAtlas
ba = AllenAtlas()

insertions = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSFits\\'
#set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = False


cont = np.arange(0.5, 35, 0.5)  # vector of contamination values to test
rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 

all_recordings = []
for iIns, ins in enumerate(insertions):

    try:
        print(f'processing {iIns + 1}/{len(insertions)}')
        eid = ins['session']['id']
        probe = ins['probe_name']
        pid = ins['probe_insertion']
    
    
        # Load in spikesorting
        sl = SpikeSortingLoader(eid=eid, pname=probe, one=one, atlas=ba)
        spikes, clusters, channels = sl.load_spike_sorting()
        clusters = sl.merge_clusters(spikes, clusters, channels)
    
    
      
    
        clusters['rep_site_acronym'] = combine_regions(clusters['acronym'])
        # Find clusters that are in the repeated site brain regions
        cluster_idx = np.where(np.isin(clusters['rep_site_acronym'], BRAIN_REGIONS))[0]
        brainRegions  = [clusters['rep_site_acronym'][i] for i in cluster_idx]

        #save the 'old' (computed in pipeline) slidingRP metric value
        # oldMetricValue = [clusters['slidingRP_viol'][x] for x in cluster_idx]
        
        #now rerun code in steinmetzlabrepo
        print('Loading spike times and clusters...')
        cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
        spikeTimes = spikes.times[cluInds]
        spikeClusters = spikes.clusters[cluInds]
        spikeAmps = spikes.amps[cluInds]
        
        print('Fitting RP Estimate...')
        
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegions, spikeAmps, rp, params)
        
        with open(savefile + ins['session']['subject'] + '.pickle', 'wb') as handle:
            pickle.dump(rpEstimates, handle)
    
    except:
        print('error')
        
        
# Nick data

from one.api import ONE
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits\\'

one = ONE(cache_dir='E:\Steinmetz2019\Steinmetz_et_al_2019_9974357\9974357') # The location of the unarchived data
sessions = one.search()

BRAIN_REGIONS = ['ACA', 'AUD','ILA' , 'MOp', 'MOs',  'OLF', 'ORB', 'ORBm',
                     'PIR', 'PL', 'RSP', 'SSp','SSs',  'VISa', 'VISam', 'VISl',
                     'VISp', 'VISpm', 'VISrl',
                 'CA', 'CA1', 'CA2', 'CA3','DG', 'POST', 'SUB'
                 'TH', 'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG','PO', 'POL', 
                     'PT','RT','SPF','VAL', 'VPL', 'VPM', 
                 'SNr','APN', 'IC','MB','MRN', 'NB','PAG','RN','SCig', 'SCm', 
                     'SCs', 'SCsg']
rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}

for e,eid in enumerate(sessions):
    
        print(f'processing {e + 1}/{len(sessions)}')
        
        print('Loading spike times and clusters...')
        # Load data
        try:
            spikes, clusters, channels =pickle.load(open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "rb"))
        except FileNotFoundError:
            channels = one.load_object(eid,'channels')
            spikes = one.load_object(eid,'spikes')
            clusters = one.load_object(eid, 'clusters')
            pickle.dump((spikes, clusters, channels),(open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "wb")))    
            

        #find all descendant acronyms of brain area of interest
        cluster_areas =[channels.brainLocation.allen_ontology[int(clusters.peakChannel[i])-1] for i in range(len(clusters.peakChannel))]
        cluster_idx = np.where(np.isin(cluster_areas, BRAIN_REGIONS))[0]
        brainRegions  = [cluster_areas[i] for i in cluster_idx]

        cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
        spikeTimes = spikes.times[cluInds]
        spikeClusters = spikes.clusters[cluInds]
        spikeAmps = spikes.amps[cluInds]
        
        print('Fitting RP Estimate...')
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegions, spikeAmps, rp, params)

        with open(savefile + eid + '.pickle', 'wb') as handle:
            pickle.dump(rpEstimates, handle)

#%% load allen data

#code to save all allen data in a dataframe is here: Desktop\ecephys_data_access-Copy1

#TODO: edit this filename and rest of this code to run on recently saved (all!) data
#make sure you know that spike times are in useconds i think? check that

#try to rerun for LP only and see what happens
#TODO plotting: make the histograms smooth like in Allen log firing rate plot 



import pickle
import pandas as pd


savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedAllenFits\\'

rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}

minFR = 1
minAmp = 50

#find session names (from the saved files)
from os import listdir
from os.path import isfile, join
mypath = 'E:\AllenBrainObservatory\saved_units'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
sessions = [i[11:20] for i in onlyfiles]

for j, session_id in enumerate(sessions):
    file_name = 'E:\AllenBrainObservatory\saved_units\Allen_units%s.pkl'%session_id
    print('loading data for session %d out of %d'%(j+1, len(sessions)))

    #load the dataframe for this session
    with open(file_name, 'rb') as f:
        units = pickle.load(f)


    #initialize rpEstimates dict
    rpEstimates = {}
    rpEstimates['cidx'] = []
    rpEstimates['rpEstimate']  = []
    rpEstimates['brainRegion'] = []
    rpEstimates['parentRegion'] = [] 
    rpEstimates['amp'] = []
    rpEstimates['fr'] = []
    
    for i, (unit_id, row) in enumerate(units.iterrows()):
        st = row.spike_times
        timeBins = rp
        fr = row.firing_rate
        amp = row.waveform_amplitude
        brainRegion = row.subregion
        parentRegion = row.region

        # print('Estimating RP Estimate...')
        # compute estimated RP for all neurons
        estimatedRP, estimateIdx, xSigmoid, ySigmoid = fitSigmoidACG(st, rp, fr, amp, minFR, minAmp, params)

        rpEstimates['cidx'].append(unit_id) 
        rpEstimates['rpEstimate'].append(estimatedRP)
        rpEstimates['brainRegion'].append(brainRegion)
        rpEstimates['parentRegion'].append(parentRegion)
        rpEstimates['fr'].append(fr)
        rpEstimates['amp'].append(amp)
        
        

    with open(savefile + session_id + '.pickle', 'wb') as handle:
        pickle.dump(rpEstimates, handle)


#%%
# file_name = savefile + session_id + '.pickle'
# with open(file_name, 'rb') as f:
#     rpEstimates = pickle.load(f)

#look into neuron 2 of session 2 -- huge amplitude but very short rp?? 

