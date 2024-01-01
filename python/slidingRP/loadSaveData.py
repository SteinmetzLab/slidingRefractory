# A script to load and save datasets to file for sliding refractory paper

import os
import pickle
import numpy as np
from phylib.stats import correlograms
from slidingRP.metrics import fitSigmoidACG,fitSigmoidACG_All,slidingRP, slidingRP_all, LlobetMetric

import pandas as pd

#establish paths:
dataBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP'


#load Horowitz data and reorganize
dataPathHorowitz = dataBasePath + '\\Horwitz_data'

spikeClustersPath = dataPathHorowitz + '\\spike_clusters.npy'
spikeTimesPath = dataPathHorowitz + '\\spike_times.npy' #in samples

#load datasets
spikeClusters = np.load(spikeClustersPath)
spikeTimes = np.load(spikeTimesPath)
#convert spikeTimes from samples to seconds:
samplingRate = 30000 #30kHz sampling rate
binSizeCorr = 1/samplingRate
spikeTimes = np.array([spikeTimes[i][0]/samplingRate for i in range(len(spikeTimes))]) #seconds
recordingDuration = spikeTimes[-1]-spikeTimes[0] #seconds
#set up rp vector:
rpEdges = np.arange(0, 10 / 1000, binSizeCorr)  # in s
rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2  #refractory period durations to test


#for each neuron:
# compute firing rate and spike count
clusterIds = np.unique(spikeClusters)

nSpikes = np.empty(len(clusterIds))
for i in range(len(clusterIds)):
    cl = clusterIds[i]
    idx = np.where(spikeClusters==cl)[0]
    nSpikes[i] = len(spikeTimes[idx])

# setup for acg

# compute an acg in 1s bins to compute the firing rate
window_size = 2 #seconds
nACGFR = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=1,
                    sample_rate=samplingRate, window_size=2, symmetrize=False)  # compute acg

firingRate = [nACGFR[i,i,1]/nSpikes[i] for i in range(len(clusterIds))] #divide each mean spike count by the total number of spikes for that neuron

window_size = 0.02 #seconds (20 ms)
nACG = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=binSizeCorr,
                    sample_rate=samplingRate, window_size=window_size, symmetrize=False)  # compute acg

acgs = np.empty((len(clusterIds),len(nACG[0][0])))
for i in range(len(clusterIds)):
    acgs[i,:]=nACG[i][i]

#set up for dataframe
dataframeColumns =['source', 'sessionID','brainRegion','recordingDuration', 'totalSpikeCount', 'ACG','estimatedRP', 'maxConfAtContThresh', 'Llobet2','Llobet3']

source = np.array(['Horwitz']*len(clusterIds))
sessionID = np.array([1]*len(clusterIds))
brainRegion = np.array(['VISp']*len(clusterIds))
recordingDuration = np.array([recordingDuration]*len(clusterIds))
totalSpikeCount = nSpikes


#run RP estimate and slidingRP:
#this dataset has no spike amps, so run without this by setting them all to 100:
spikeAmps = 100*np.ones(len(spikeTimes))
params={}
params['binSizeCorr']=binSizeCorr
params['runLlobet'] = True
params['contaminationThresh']=10

rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegion, spikeAmps, rp, params)
estimatedRP = rpEstimates['rpEstimate']

rpMetrics = slidingRP_all(spikeTimes, spikeClusters, **params)
maxConfAtContThresh = rpMetrics['maxConfidenceAtContaminationThresh']
Llobet2 = rpMetrics['Llobet2']
Llobet3 = rpMetrics['Llobet3']






dataList = list(zip(source,sessionID,brainRegion,recordingDuration,totalSpikeCount,acgs, estimatedRP,maxConfAtContThresh,Llobet2,Llobet3))
df_Horwitz = pd.DataFrame(dataList,columns = dataframeColumns)



#now do the same for IBL, Steinmetz,Allen:
#%%Steinmetz data

from one.api import ONE
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits\\'

one = ONE(cache_dir='E:\Steinmetz2019\Steinmetz_et_al_2019_9974357\9974357')  # The location of the unarchived data
sessions = one.search()

BRAIN_REGIONS = ['ACA', 'AUD', 'ILA', 'MOp', 'MOs', 'OLF', 'ORB', 'ORBm',
                 'PIR', 'PL', 'RSP', 'SSp', 'SSs', 'VISa', 'VISam', 'VISl',
                 'VISp', 'VISpm', 'VISrl',
                 'CA', 'CA1', 'CA2', 'CA3', 'DG', 'POST', 'SUB',
                  'TH', 'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'POL',
                 'PT', 'RT', 'SPF', 'VAL', 'VPL', 'VPM',
                 'SNr', 'APN', 'IC', 'MB', 'MRN', 'NB', 'PAG', 'RN', 'SCig', 'SCm',
                 'SCs', 'SCsg']
rpBinSize = 1 / 30000
rpEdges = np.arange(0, 10 / 1000, rpBinSize)  # in s
rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2  # vector of refractory period durations to test
params = {}

for e, eid in enumerate(sessions):

    print(f'processing {e + 1}/{len(sessions)}')

    print('Loading spike times and clusters...')
    # Load data
    try:
        spikes, clusters, channels = pickle.load(
            open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "rb"))
    except FileNotFoundError:
        channels = one.load_object(eid, 'channels')
        spikes = one.load_object(eid, 'spikes')
        clusters = one.load_object(eid, 'clusters')
        pickle.dump((spikes, clusters, channels),
                    (open("E:\Steinmetz2019\loaded_data\data_{eid}.p".format(eid=eid), "wb")))

        # find all descendant acronyms of brain area of interest
    cluster_areas = [channels.brainLocation.allen_ontology[int(clusters.peakChannel[i]) - 1] for i in
                     range(len(clusters.peakChannel))]
    clusterIds = np.where(np.isin(cluster_areas, BRAIN_REGIONS))[0] #this step discards neurons not in BRAIN_REGIONS
    brainRegions = [cluster_areas[i] for i in clusterIds]

    cluInds = np.array(
        list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in clusterIds))  # this takes forever??
    spikeTimes = spikes.times[cluInds]
    spikeClusters = spikes.clusters[cluInds]
    spikeAmps = spikes.amps[cluInds]
    recordingDuration = spikeTimes[-1] - spikeTimes[0]  # seconds

    #compute total number of spikes
    nSpikes = np.empty(len(clusterIds))
    for i in range(len(clusterIds)):
        cl = clusterIds[i]
        idx = np.where(spikeClusters == cl)[0]
        nSpikes[i] = len(spikeTimes[idx])

    #compute ACG:
    # compute an acg in 1s bins to compute the firing rate
    window_size = 2  # seconds
    nACGFR = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=1,
                          sample_rate=samplingRate, window_size=2, symmetrize=False)  # compute acg

    firingRate = [nACGFR[i, i, 1] / nSpikes[i] for i in
                  range(len(clusterIds))]  # divide each mean spike count by the total number of spikes for that neuron

    window_size = 0.02  # seconds (20 ms)
    nACG = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=binSizeCorr,
                        sample_rate=samplingRate, window_size=window_size, symmetrize=False)  # compute acg

    acgs = np.empty((len(clusterIds), len(nACG[0][0])))
    for i in range(len(clusterIds)):
        acgs[i, :] = nACG[i][i]

    # compute estimated RP for all neurons
    rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters, brainRegions, spikeAmps, rp, params)


    source = np.array(['Steinmetz']*len(clusterIds))
    sessionID = np.array([eid]*len(clusterIds))
    recordingDuration = np.array([recordingDuration]*len(clusterIds))
    totalSpikeCount = nSpikes

    #run RP estimate and slidingRP:
    params={}
    params['binSizeCorr']=binSizeCorr
    params['runLlobet'] = True
    params['recordingDuration'] = recordingDuration[0]
    params['contaminationThresh'] = 10

    # rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegion, spikeAmps, rp, params)
    estimatedRP = rpEstimates['rpEstimate']

    rpMetrics = slidingRP_all(spikeTimes, spikeClusters, **params)
    maxConfAtContThresh = rpMetrics['maxConfidenceAtContaminationThresh']
    Llobet2 = rpMetrics['Llobet2']
    Llobet3 = rpMetrics['Llobet3']


    dataList = list(zip(source,sessionID,brainRegion,recordingDuration,totalSpikeCount,acgs, estimatedRP,maxConfAtContThresh,Llobet2,Llobet3))
    df = pd.DataFrame(dataList,columns = dataframeColumns)

    if e==0:
        df_Steinmetz = df
    else:
        df_Steinmetz= pd.concat([df_Steinmetz,df])


#%%

rpBinSize = 1 / 30000
rpEdges = np.arange(0, 10 / 1000, rpBinSize)  # in s
rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2  # vector of refractory period durations to test
params = {}

minFR = 1
minAmp = 50

# find session names (from the saved files)
from os import listdir
from os.path import isfile, join

mypath = 'E:\AllenBrainObservatory\saved_units'
sessionFiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
sessionNames = [i[11:20] for i in sessionFiles] #get the session name from the filename string

for j, session_id in enumerate(sessionNames):
    file_name = 'E:\AllenBrainObservatory\saved_units\Allen_units%s.pkl' % session_id
    print('loading data for session %d out of %d' % (j + 1, len(sessionNames)))

    # load the dataframe for this session
    # with open(file_name, 'rb') as f:
    #     units = pickle.load(f)
    units = pd.read_pickle(file_name)

    # initialize rpEstimates dict
    rpEstimates = {}
    rpEstimates['cidx'] = []
    rpEstimates['rpEstimate'] = []
    rpEstimates['brainRegion'] = []
    rpEstimates['parentRegion'] = []
    rpEstimates['amp'] = []
    rpEstimates['fr'] = []

    recordingDurationVec  =[]
    nSpikesVec = []
    maxConfVec = []
    Llobet2Vec = []
    Llobet3Vec = []

    for i, (unit_id, row) in enumerate(units.iterrows()):
        st = row.spike_times
        recordingDuration = st[-1] - st[0]  # seconds
        nSpikes = len(st)
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

        # run RP estimate and slidingRP:
        params = {}
        params['binSizeCorr'] = binSizeCorr
        params['runLlobet'] = True
        params['recordingDuration'] = recordingDuration
        params['contaminationThresh'] = 10
        testedCont = params['contaminationThresh']/100

        [maxConfAtContThresh, minContAtConfidenceThresh, timeOfLowestCont,
            nSpikesBelow2, confMatrix, cont, rpVec, nACG,
            firingRate, secondsElapsed] = slidingRP(st, params)

        fp, confLlobet2 = LlobetMetric(firingRate, params['recordingDuration'], nACG, rpVec, testedCont, 0.002,
                                       minISI=0)
        fp, confLlobet3 = LlobetMetric(firingRate, params['recordingDuration'], nACG, rpVec, testedCont, 0.003,
                                       minISI=0)
        confLlobet2 = confLlobet2 * 100
        confLlobet3 = confLlobet3 * 100

        nACG = nACG[0:len(rpVec)+1]
        if i ==0:
            acgs = np.empty((len(units),len(nACG)))

        acgs[i, :] = nACG

        recordingDurationVec.append(recordingDuration)
        nSpikesVec.append(nSpikes)
        maxConfVec.append(maxConfAtContThresh)
        Llobet2Vec.append(confLlobet2)
        Llobet3Vec.append(confLlobet3)

    source = np.array(['Allen'] * len(units))
    sessionID = np.array([session_id] * len(units))
    estimatedRP = rpEstimates['rpEstimate']
    brainRegion = rpEstimates['brainRegion']
    recordingDuration = recordingDurationVec
    totalSpikeCount = nSpikesVec
    maxConfAtContThresh = maxConfVec
    Llobet3 = Llobet3Vec
    Llobet2 = Llobet2Vec



    dataList = list(
        zip(source, sessionID, brainRegion, recordingDuration, totalSpikeCount, acgs, estimatedRP, maxConfAtContThresh,
            Llobet2, Llobet3))
    df = pd.DataFrame(dataList, columns=dataframeColumns)

    if j==0:
        df_Allen = df
    else:
        df_Allen = pd.concat([df_Allen, df])

#%%
#IBL Repeated site:
from reproducible_ephys_functions import combine_regions, BRAIN_REGIONS, get_insertions, save_data_path, \
    save_dataset_info

from brainbox.io.one import SpikeSortingLoader
from ibllib.atlas import AllenAtlas

ba = AllenAtlas()

insertions = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSFits\\'
# set up params for slidingRP_all
params = {}
params['sampleRate'] = []
params['sampleRate'] = 30000
params['binSizeCorr'] = 1 / params['sampleRate']
params['returnMatrix'] = True
params['verbose'] = False

cont = np.arange(0.5, 35, 0.5)  # vector of contamination values to test
rpBinSize = 1 / 30000
rpEdges = np.arange(0, 10 / 1000, rpBinSize)  # in s
rp = rpEdges + np.mean(np.diff(rpEdges)[0]) / 2

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
        clusterIds = np.where(np.isin(clusters['rep_site_acronym'], BRAIN_REGIONS))[0]
        brainRegions = [clusters['rep_site_acronym'][i] for i in clusterIds]


        print('Loading spike times and clusters...')
        cluInds = np.array(
            list(x for x in range(len(spikes.clusters)) if spikes.clusters[x] in clusterIds))  # this takes forever??
        spikeTimes = spikes.times[cluInds]
        spikeClusters = spikes.clusters[cluInds]
        spikeAmps = spikes.amps[cluInds]
        brainRegions = [cluster_areas[i] for i in clusterIds]
        recordingDuration = spikeTimes[-1] - spikeTimes[0]  # seconds

        # compute total number of spikes
        nSpikes = np.empty(len(clusterIds))
        for i in range(len(clusterIds)):
            cl = clusterIds[i]
            idx = np.where(spikeClusters == cl)[0]
            nSpikes[i] = len(spikeTimes[idx])

        # compute ACG:
        # compute an acg in 1s bins to compute the firing rate
        window_size = 2  # seconds
        nACGFR = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=1,
                              sample_rate=samplingRate, window_size=2, symmetrize=False)  # compute acg

        firingRate = [nACGFR[i, i, 1] / nSpikes[i] for i in
                      range(
                          len(clusterIds))]  # divide each mean spike count by the total number of spikes for that neuron

        window_size = 0.02  # seconds (20 ms)
        nACG = correlograms(spikeTimes, spikeClusters, cluster_ids=clusterIds, bin_size=binSizeCorr,
                            sample_rate=samplingRate, window_size=window_size, symmetrize=False)  # compute acg

        acgs = np.empty((len(clusterIds), len(nACG[0][0])))
        for i in range(len(clusterIds)):
            acgs[i, :] = nACG[i][i]

        # compute estimated RP for all neurons
        rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters, brainRegions, spikeAmps, rp, params)
        source = np.array(['IBL']*len(clusterIds))
        sessionID = np.array([eid]*len(clusterIds))
        recordingDuration = np.array([recordingDuration]*len(clusterIds))
        totalSpikeCount = nSpikes

        #run RP estimate and slidingRP:
        params={}
        params['binSizeCorr']=binSizeCorr
        params['runLlobet'] = True
        params['recordingDuration'] = recordingDuration[0]
        params['contaminationThresh'] = 10

        # rpEstimates = fitSigmoidACG_All(spikeTimes, spikeClusters,brainRegion, spikeAmps, rp, params)
        estimatedRP = rpEstimates['rpEstimate']

        rpMetrics = slidingRP_all(spikeTimes, spikeClusters, **params)
        maxConfAtContThresh = rpMetrics['maxConfidenceAtContaminationThresh']
        Llobet2 = rpMetrics['Llobet2']
        Llobet3 = rpMetrics['Llobet3']


        dataList = list(zip(source,sessionID,brainRegion,recordingDuration,totalSpikeCount,acgs, estimatedRP,maxConfAtContThresh,Llobet2,Llobet3))
        df = pd.DataFrame(dataList,columns = dataframeColumns)

        if e==0:
            df_IBL= df
        else:
            df_IBL= pd.concat([df_IBL,df])



#%%
#concatenate all to a dataframe, and save:

listDFs = [df_Horwitz, df_Steinmetz, df_IBL, df_Allen]

df_AllData = pd.concat(listDFs)

#save to tsv:
saveFileName = dataBasePath + '\\AllRPNeurons.tsv'

df_AllData.to_csv(saveFileName, sep = "\t")