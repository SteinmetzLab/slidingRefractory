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
#%%
#IBL 
insertions = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSFits\\'
#set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = False


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
        
        print('Fiting RP Estimate...')
        
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegions, rp, params)

        with open(savefile + ins['session']['subject'] + '.pickle', 'wb') as handle:
            pickle.dump(rpEstimates, handle)
        
    except:
        print('error')
        
        
#%% Nick data
from one.api import ONE
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits\\'

one = ONE(cache_dir='E:\Steinmetz2019\Steinmetz_et_al_2019_9974357\9974357') # The location of the unarchived data
sessions = one.search()#dataset='trials') # search for all sessions that have a `trials` object
# session = sessions[0] # take the first session
# trials = one.load_object(session, 'trials')




#first, find areas:
cluster_areas = []
for e,eid in enumerate(sessions[0:1]):

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
    cluster_areas.append([channels.brainLocation.allen_ontology[int(clusters.peakChannel[i])-1] for i in range(len(clusters.peakChannel))])

#%%
#TODO find a way to do this programatically....
clusterAreasFlat = [item for sublist in cluster_areas for item in sublist]
# uniqueAreas = np.unique(clusterAreasFlat)
rerunAncestors = False
#takes 8 min to run!!
if rerunAncestors:
    ancestors = [br.ancestors(br.acronym2id(acr)).acronym for acr in clusterAreasFlat]
    with open('SteinmetzClusterAreas.pickle', 'wb') as handle:
        pickle.dump(ancestors, handle)
else:
     file = open('SteinmetzClusterAreas.pickle','rb')
     ancestors = pickle.load(file)
     file.close()   
#%%
# ancestors = [br.ancestors(br.acronym2id(acr)).acronym for acr in uniqueAreas]

clusterRelevantAreas =np.empty(len(ancestors), dtype=object)
for a in range(len(ancestors)):
    clusterRelevantAreas[a] = [ancestors[a][i] for i in range(len(ancestors[a])) if ancestors[a][i] in bigAreas]

clusterRelevantAreasFlat = [item for sublist in clusterRelevantAreas for item in sublist]


with open('SteinmetzClusterAreas.pickle', 'wb') as handle:
    pickle.dump(ancestors, handle)
        



bigAreas = ['Isocortex', 'HIP','MB','TH']

relevantAreas =np.empty(len(ancestors), dtype=object)
for j in range(len(ancestors)):
    xx =  [ancestors[j][i] for i in range(len(ancestors[j])) if ancestors[j][i] in bigAreas]
    if len(xx)>0:
        relevantAreas[j] = xx[0]

[x for i in range(len(clusterAreasFlat)) if clusterAreasFlat[i] in bigAreas
res = [([idx for idx, val in enumerate(relevantAreas) if val == sub] if sub in relevantAreas else [None])
      for sub in clusterAreasFlat]
  
res = [([idx for idx, val in enumerate(relevantAreas) if val == sub] if sub in relevantAreas else [None])
      for sub in clusterAreasFlat]

    #%%
from one.api import ONE
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits\\'

one = ONE(cache_dir='E:\Steinmetz2019\Steinmetz_et_al_2019_9974357\9974357') # The location of the unarchived data
sessions = one.search()#dataset='trials') # search for all sessions that have a `trials` object
# session = sessions[0] # take the first session
# trials = one.load_object(session, 'trials')

#find sessions that have at least one channel in Cortex, Thalamus
# BRAIN_REGIONS = ['FRP','MO','SS','AUD','VIS','ACA','PL','ILA','ORB','RSP','TEa',
#                'CA1','CA2','CA3','DG',
#                'MB',
#                'DORsm','DORpm'] #change from Isocortext to 'CTX'??

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

for e,eid in enumerate(sessions[1:]):
    
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
        
        
        print('Fiting RP Estimate...')
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegions, rp, params)

        with open(savefile + eid + '.pickle', 'wb') as handle:
            pickle.dump(rpEstimates, handle)

#%% load allen data
# from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache

# data_directory = 'E:\AllenBrainObservatory' # must be updated to a valid directory in your filesystem

# manifest_path = os.path.join(data_directory, "manifest.json")

# cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)
# units_all = cache.get_units(amplitude_cutoff_maximum = np.inf,
#                         presence_ratio_minimum = -np.inf,
#                         isi_violations_maximum = np.inf)
# units = cache.get_units()

# num_units = len(units)
import pickle
#now save for use in iblenv
filename = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\metrics\slidingRP\ts_region_dict_all.pkl'

file = open(filename,'rb')

ts_dict = pickle.load(file)
file.close()


savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedAllenFits\\'

rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}

for b,brainRegion in enumerate(ts_dict.keys()):
    
        print(f'processing {b + 1}/{len(ts_dict.keys())}')
        
        print('Loading spike times and clusters...')
        # Load data
        # try:
        #     spikes, clusters, channels =pickle.load(open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "rb"))
        # except FileNotFoundError:
        #     channels = one.load_object(eid,'channels')
        #     spikes = one.load_object(eid,'spikes')
        #     clusters = one.load_object(eid, 'clusters')
        #     pickle.dump((spikes, clusters, channels),(open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "wb")))    
        clusters = []    
        spikes = []
        for i in range(len(ts_dict[brainRegion])):
           numspikes = len(ts_dict[brainRegion][i])
           clusters.append(np.repeat(i,numspikes))
           spikes.append(ts_dict[brainRegion][i])
        
        clustersFlat =  np.array([item for sublist in clusters for item in sublist])
        spikesFlat = np.array([item for sublist in spikes for item in sublist])

        #find all descendant acronyms of brain area of interest
        # cluster_areas =[channels.brainLocation.allen_ontology[int(clusters.peakChannel[i])-1] for i in range(len(clusters.peakChannel))]
        # cluster_idx = np.where(np.isin(cluster_areas, BRAIN_REGIONS))[0]
        brainRegions  = [brainRegion for i in np.unique(clustersFlat)]

        # cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
        # spikeTimes = spikes.times[cluInds]
        # spikeClusters = spikes.clusters[cluInds]
        
        
        print('Fiting RP Estimate...')
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikesFlat, clustersFlat,brainRegions, rp, params)

        with open(savefile + brainRegion + '.pickle', 'wb') as handle:
            pickle.dump(rpEstimates, handle)


