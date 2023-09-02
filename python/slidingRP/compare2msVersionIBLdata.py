'''
Codce to look into IBL metrics for:
regular slidingRP
2ms slidingRP for neurons with FR:
>0.33 Hz,
<0.5 Hz,
0.75 Hz.

'''

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from one.api import ONE
from brainbox.io.one import SessionLoader
from brainbox.io.one import SpikeSortingLoader
from brainwidemap import bwm_query, load_good_units, load_trials_and_mask, filter_regions, filter_sessions, \
    download_aggregate_tables
from slidingRP.metrics import slidingRP_all
import pickle
#%%
# one = ONE()
print('instantiate one')
one = ONE()

#download all IBL frozen BWM data
print('download frozen df')
bwm_df = bwm_query(one = one)

"""
Map probes to regions and filter regions
"""
# First, download the latest clusters table
clusters_table = download_aggregate_tables(one, overwrite = False) #type='clusters', local_path=local_path)


# Map probes to regions (here: Beryl mapping) and filter according to QC, number of units per region and number of
# probes per region. These are the default values explicitly specified only for demonstration purposes.
region_df = filter_regions(bwm_df['pid'], clusters_table=clusters_table, mapping='Beryl', min_qc=1,
                           min_units_region=10, min_probes_region=2)

"""
Filter sessions on min trial number
"""
# First, download the latest trials table. Note that currently, the bwm_include column is currently based on the
# whether a trial is included based on the default parameters of the load_trials_and_mask function (see below).
# If you need other criteria you can add a column to the trials table, or ask for help
trials_table = download_aggregate_tables(one, type='trials')
eids = filter_sessions(bwm_df['eid'], trials_table=trials_table, min_trials=200)


"""
Remove probes and sessions based on those filters
"""
bwm_df = bwm_df[bwm_df['pid'].isin(region_df['pid'].unique())]
bwm_df = bwm_df[bwm_df['eid'].isin(eids)]

"""
Load spike sorting, good units only
"""
# Load cluster information and spike trains for all good units.
# The first time you run this, it will take a long time because it needs to download all data


overall_struct = {}
overall_struct['cidx'] = []
overall_struct['maxConfidenceAt10Cont'] = []
overall_struct['minContWith90Confidence'] = []
overall_struct['timeOfLowestCont'] = []
overall_struct['nSpikesBelow2'] = []
overall_struct['firingRate'] = []
overall_struct['pid'] = []
overall_struct['eid'] = []
overall_struct['subject'] = []
overall_struct['previous_pass']=[]
overall_struct['pt5_pass'] = []
overall_struct['pt33_pass'] = []
overall_struct['pt75_pass'] = []
overall_struct['dynamic_pass'] = []
overall_struct['dynamic_threshold'] = []
overall_struct['prev_RP']=[]
overall_struct['label']=[]

#%%define a closest function for finding the nearest recording duration and firing rate
def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]
#%%
# run the loop
for i, pid in enumerate(bwm_df['pid']):

    #for each trajectory, load spikes
    try:
        pid = bwm_df['pid'][i]
        eid = bwm_df['eid'][i]
        pname = bwm_df['probe_name'][i]
        subject = bwm_df['subject'][i]

        spike_loader = SpikeSortingLoader(pid=pid, one=one, eid=eid, pname=pname)
        spikes, clusters, channels = spike_loader.load_spike_sorting()
        clusters_labeled = SpikeSortingLoader.merge_clusters(
            spikes, clusters, channels).to_df()
    except:
        print('spike loading error')
        continue
    print('Running code for pid number ' + str(i))
    print(subject)

    return_struct = slidingRP_all(spikes.times, spikes.clusters)



    '''
    goal here is to add:
    for each recording, pull out the recording length
    then look up on the heat map, what is the corresponding value for
     1) 90% reject on 50% contamijnation,
     2)  50% reject on 20% contamination
    I guess pick the shortest of these two (but also save them, so you can plot this later)
    and then for this value, put it in a new struct, tailored, where you run the inclusion criteria based on this value

    '''
    # first spike to last spike (across all recorded neurons), in hours
    recordingDuration = (spikes.times[-1]-spikes.times[0]) / 3600

    #load dict of firing rates from contour plots
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\contFRContourDict.pickle'
    file = open(savefile, 'rb')
    contContourDict = pickle.load(file)
    file.close()

    #also load corresponding params for the corresponding recording durations:
    # load data
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_122.pickle'
    file = open(savefile, 'rb')
    results = pickle.load(file)
    file.close()
    params = results[2]

    recDurClosest = closest(params['recDurs'], recordingDuration) #closest rec dur for which I've computed contours
    recDurInd = np.where(params['recDurs'] == recDurClosest)


    #The rule is: 20% contamination we want at least 50% reject, for 50% contamination we want at least 90% reject
    lowContFR = contContourDict[(0.2, 50)][recDurInd]
    highContFR = contContourDict[(0.5, 90)][recDurInd]

    frThreshRecDur = min(lowContFR, highContFR)


    for n in range(len(return_struct['cidx'])):

        #log if this neuron passes generally
        try:
            cid = return_struct['cidx'][n]
            label = clusters.label[np.where(clusters.cluster_id == cid)[0][0]]
            prev_RP =clusters.slidingRP_viol[np.where(clusters.cluster_id == cid)[0][0]]
        except:
            label = np.nan
            prev_RP = np.nan

        if return_struct['minContWith90Confidence'][n] <= 10:
            passfail = 1
        else:
            passfail = 0

        frThresh = frThreshRecDur
        dynamic_pass = passfail #by default, set the value to the same as the normal metric
        if return_struct['nSpikesBelow2'][n]==0 and return_struct['firingRate'][n]>frThresh:
            dynamic_pass = 1 #but, if there are no spikes and fr is higher than the threshold, set to pass
        if return_struct['firingRate'][n] < frThresh:
            dynamic_pass = 0 #if the firing rate is less than the threshold, set to fail


        frThresh = 0.5
        pt5_pass = passfail #by default, set the value to the same as the normal metric
        if return_struct['nSpikesBelow2'][n]==0 and return_struct['firingRate'][n]>frThresh:
            pt5_pass = 1 #but, if there are no spikes and fr is higher than the threshold, set to pass
        if return_struct['firingRate'][n] < frThresh:
            pt5_pass = 0 #if the firing rate is less than the threshold, set to fail

        frThresh = 0.33
        pt33_pass = passfail #by default, set the value to the same as the normal metric
        if return_struct['nSpikesBelow2'][n]==0 and return_struct['firingRate'][n]>frThresh:
            pt33_pass = 1 #but, if there are no spikes and fr is higher than the threshold, set to pass
        if return_struct['firingRate'][n] < frThresh:
            pt33_pass = 0 #if the firing rate is less than the threshold, set to fail

        frThresh = 0.75
        pt75_pass = passfail #by default, set the value to the same as the normal metric
        if return_struct['nSpikesBelow2'][n]==0 and return_struct['firingRate'][n]>frThresh:
            pt75_pass = 1 #but, if there are no spikes and fr is higher than the threshold, set to pass
        if return_struct['firingRate'][n] < frThresh:
            pt75_pass = 0 #if the firing rate is less than the threshold, set to fail






        overall_struct['cidx'].append(return_struct['cidx'][n])

        overall_struct['prev_RP'].append(prev_RP)
        overall_struct['label'].append(label)

        overall_struct['maxConfidenceAt10Cont'].append(return_struct['maxConfidenceAt10Cont'][n])
        overall_struct['minContWith90Confidence'].append(return_struct['minContWith90Confidence'][n])
        overall_struct['timeOfLowestCont'].append(return_struct['timeOfLowestCont'][n])
        overall_struct['nSpikesBelow2'].append(return_struct['nSpikesBelow2'][n])
        overall_struct['firingRate'].append(return_struct['firingRate'][n])
        overall_struct['previous_pass'].append(passfail)
        overall_struct['pt5_pass'].append(pt5_pass)
        overall_struct['pt33_pass'].append(pt33_pass)
        overall_struct['pt75_pass'].append(pt75_pass)
        overall_struct['dynamic_pass'].append(dynamic_pass)
        overall_struct['dynamic_threshold'].append(frThreshRecDur)

        overall_struct['pid'].append(pid)
        overall_struct['eid'].append(eid)
        overall_struct['subject'].append(subject)

    #compute fr and 2ms spikes pass version from metric code (For all neurons do 2ms version!!)
    #put in a struc with session and subject and neuron number
    #struct of neurons for 0.33 thresh, struct of neurons for 0.5 thresh, struct of neurons of 0.75 hz
    #now go from pass to fail for neurons < threshold, and put those in each of the structs



savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\saved_2msIBL_v2.pkl'
with open(savefile, 'wb') as f:
    pickle.dump(overall_struct, f)
#%%
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\saved_2msIBL_v2.pkl'

with open(savefile, 'rb') as f:
    overall_struct = pickle.load(f)
#%%
#here, only look at those which prev_RP=0, recDurRP=1, and label = 0.66

prevRP = np.array(overall_struct['prev_RP'])
dynamic_pass = np.array(overall_struct['dynamic_pass'])
label = np.array(overall_struct['label'])
previous = len(np.where(label==1)[0]) #neurons that pass all 3 metrics currently!
added_with_dynamic2ms = np.where((prevRP==0) & (dynamic_pass == 1) & (label == 2/3))[0]
dynamic = len(added_with_dynamic2ms) + previous #the number of units (added w dynamic plus pervious pass)

#compute percent
total_neurons = len(overall_struct['pt5_pass'])
previous_pct = previous/total_neurons *100
dynamic_pct = dynamic/total_neurons *100

labels = ['Original metric', 'Proposed adjustment']

fig, ax1 = plt.subplots(figsize = (5.5,5))
# ax2 = ax1.twinx()

ax1.bar(labels, [previous,dynamic], color = '0.6')
ax1.set_ylabel('Number of passing neurons')
# ax2.bar(labels, [previous_pct,dynamic_pct], color = '0.6')
# ax2.set_ylabel('Percent of passing neurons')
plt.tight_layout()
plt.show()

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\2msComparison'
fig.savefig( savefile+ '_recDur.png', dpi=500)


#%%
previous = sum(overall_struct['previous_pass'])
pt33 = sum(overall_struct['pt33_pass'])
pt5 = sum(overall_struct['pt5_pass'])
pt75 = sum(overall_struct['pt75_pass'])
dynamic = sum(overall_struct['dynamic_pass'])
labels = ['Original metric', '0.33 spks/s', '0.5 spks/s', '0.75 spks/s', 'Based on Rec Dur']
title = 'Number of passing neurons with different FR thresholds'

fig,ax = plt.subplots()
ax.bar(labels, [previous, pt33,pt5,pt75,dynamic], color = '0.6')
ax.set_ylim([80000, 130000])
ax.set_ylabel('Number of neurons')
ax.set_title(title)
plt.show()

#%%
title = 'Percent of passing neurons with different FR thresholds'
total_neurons = len(overall_struct['pt5_pass'])
fig,ax = plt.subplots()
previous_pct = previous/total_neurons *100
pt33_pct = pt33/total_neurons * 100
pt5_pct = pt5/total_neurons * 100
pt75_pct = pt75/total_neurons * 100
dynamic_pct = dynamic/total_neurons *100
ax.bar(labels, [previous_pct, pt33_pct,pt5_pct,pt75_pct,dynamic_pct], color = '0.6')
ax.set_ylabel('Percent of total IBL neurons')
ax.set_title(title)
ax.set_ylim([30,50])
plt.show()

#%%
fig, ax1 = plt.subplots(figsize = (5.5,5))
ax2 = ax1.twinx()

previous = sum(overall_struct['previous_pass'])
# pt33 = sum(overall_struct['pt33_pass'])
# pt5 = sum(overall_struct['pt5_pass'])
# pt75 = sum(overall_struct['pt75_pass'])
dynamic = sum(overall_struct['dynamic_pass'])

total_neurons = len(overall_struct['pt5_pass'])

previous_pct = previous/total_neurons *100
# pt33_pct = pt33/total_neurons * 100
# pt5_pct = pt5/total_neurons * 100
# pt75_pct = pt75/total_neurons * 100
dynamic_pct = dynamic/total_neurons *100


labels = ['Original metric', 'Proposed adjustment']

ax1.bar(labels, [previous,dynamic], color = '0.6')
ax1.set_ylabel('Number of passing neurons')
ax2.bar(labels, [previous_pct,dynamic_pct], color = '0.6')
ax2.set_ylabel('Percent of passing neurons')
plt.tight_layout()
plt.show()



#%%

#plot minFR as a function of recDur
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\contFRContourDict.pickle'
file = open(savefile, 'rb')
contContourDict = pickle.load(file)
file.close()

#need to also load in the corresponding params file for the correct recording durations
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_122.pickle'
file = open(savefile, 'rb')
results = pickle.load(file)
file.close()

params = results[-1] #the last in the results
pc = results[0] #normal metric results (not 2ms cond)

frThreshRecDur = []
recDurVec = params['recDurs']
for recordingDuration in recDurVec:
    recDurClosest = closest(params['recDurs'], recordingDuration)  # closest rec dur for which I've computed contours
    recDurInd = np.where(params['recDurs'] == recDurClosest)
    # The rule is: 20% contamination we want at least 50% reject, for 50% contamination we want at least 90% reject
    lowContFR = contContourDict[(0.2, 50)][recDurInd]
    highContFR = contContourDict[(0.5, 90)][recDurInd]
    print(lowContFR)
    frThreshRecDur.append(max(lowContFR, highContFR)[0])
#%%
fig,ax = plt.subplots()
ax.plot(recDurVec,frThreshRecDur,'.-',color = 'k')
ax.set_xlabel('Recording duration (hrs)')
ax.set_ylabel('Minimum passing FR (spks/s)')

ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
# ax.set_yticks()
plt.show()

#todo: make spines go away, make dots instead of line, make yticks fewer

#%% rerun to compute lowest FR for recDur

#run simulations:
sampleRate = 30000
params = {
    'recDurs': np.array([0.5 , 0.75, 1.  , 1.25, 1.5 , 1.75, 2.  , 2.25, 2.5 , 2.75, 3.  ]),  #recording durations (hours)
    'RPs': np.array([0.002]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates': np.array([0]),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
    'nSim': 200,
    'contaminationThresh': 10,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'confidenceThresh': 90,
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True,
    'savePCfile': True

}



#%% run and save

print('in simulations, minFR')
[pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
    params['nSim']) + 'iterTestingMinFR.pickle'

results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, params]
if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)

#%% load and plot
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
    params['nSim']) + 'iterTestingMinFR.pickle'
file = open(savefile,'rb')
results = pickle.load(file)
file.close()
pc = results[0]

x = params['baseRates']
y=[]
for r, recDur in enumerate(params['recDurs']):
    firstNonZero = np.nonzero(pc[r,0,:,0])
    if len(firstNonZero[0]):
        val = x[firstNonZero[0][0]]
    else:
        val = np.nan
    y.append(val)



rd = params['recDurs']


fig, ax = plt.subplots(1, 1, figsize=(6, 4))

ax.set_xlabel('Recording Duration (hours)')
ax.set_ylabel('First FR with nonzero PC ')
ax.plot(rd,y,'k.-')
ax.set_title('Minimum passing FR (uncontaminated)')

spinesSetting = False
ax.spines.right.set_visible(spinesSetting)
ax.spines.top.set_visible(spinesSetting)
fig.show()


