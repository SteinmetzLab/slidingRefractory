# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:31:29 2022

@author: Noam Roth
scratch code to figure out whether amplitudes/frs matter for low fits
"""


from one.api import ONE
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits_new\\'

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
        cluster_areas =[channels.brainLocation.allen_ontology[int(clusters.peakChannel[i])-1] for i in range(len(clusters.peakChannel))]
        cluster_idx = np.where(np.isin(cluster_areas, BRAIN_REGIONS))[0]
        brainRegions  = [cluster_areas[i] for i in cluster_idx]

        cluInds = np.array(list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in cluster_idx)) #this takes forever??
        spikeTimes = spikes.times[cluInds]
        spikeClusters = spikes.clusters[cluInds]
        
        
        print('Fitting RP Estimate...')
        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegions, rp, params)

        #add amp and fr info
        rpEstimates = collectAmpFR(spikes,cluInds,rpEstimates)


        # with open(savefile + eid + '.pickle', 'wb') as handle:
        #     pickle.dump(rpEstimates, handle)
#%%
fig,axs = plt.subplots(2,1)
axs[0].scatter(rpEstimates['rpEstimate'],rpEstimates['fr'])
axs[1].scatter(rpEstimates['rpEstimate'],rpEstimates['med_amp'])

#%%
def collectAmpFR(spikes,cluInds,rpEstimate):
    spikeTimes = spikes.times[cluInds]
    spikeClusters = spikes.clusters[cluInds]
    spikeAmps = spikes.amps[cluInds]
    
    
    cids = np.unique(spikeClusters)
    
    #add to dict
    rpEstimate['med_amp'] = []
    rpEstimate['fr']  = []

    for cidx in range(len(cids)):
        if cids[cidx] in rpEstimate['cidx']:
            st = spikeTimes[spikeClusters==cids[cidx]] 
            fr = len(st)/(st[-1]-st[0])
            amp = np.median(spikeAmps[spikeClusters==cids[cidx]] )
            rpEstimate['med_amp'].append(amp)
            rpEstimate['fr'].append(fr)

    
    return rpEstimate
    