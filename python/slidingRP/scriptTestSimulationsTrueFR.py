import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')


import pickle
import numpy as np
from phylib.stats import correlograms

from scipy import stats
from slidingRP.simulations import *
import datetime

date_now = datetime.datetime.now().strftime('_%m_%d')
#%%
#run simulations:
sampleRate = 30000
params = {
    'recDurs': np.array([2]),  #recording durations (hours)
    'RPs': np.array([0.002]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.25, 0.5,0.75, 1,1.5, 2,3, 4,5,6, 7.5,10,15],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates': np.arange(0.00,0.21, 0.02),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
    'nSim': 500,
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
drift_values = [-0.5,0, 0.5]
drift_strings = ['Dec','Same','Inc']
for drift,string in zip(drift_values,drift_strings):
    if drift ==0:
        params['baseRates']=[0.75, 1, 1.25, 1.5, 2, 2.5, 3.75, 5, 6.25, 7.5,10, 12.5]
    else:
        params['baseRates'] = [1,2,5,10]
    print(params['baseRates'])
    params['delta'] = drift
    print('in simulations, {0} drift'.format(drift))
    [pc] = simulateChangingContNeurons(params)

    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
        params['nSim']) + 'iter' + date_now + 'delta' + string + '.pickle'


    results = [pc, params]
    if params['savePCfile']:
        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)



#%% rerun just 0 because there was some issue?
[pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
    params['nSim']) + 'iter' + date_now + 'delta' + string + '.pickle'

results = [pc, params]
if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)

#%% plot
pcDict = {}
paramsDict = {}
date_now = '_08_02'

for drift, string in zip(drift_values,drift_strings):
    print('loading sim results {0} drift'.format(drift))
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
        params['nSim']) + 'iter' + date_now + 'delta' + string + '.pickle'

    file = open(savefile,'rb')
    results = pickle.load(file)
    file.close()


    pcDict[drift] = results[0]
    paramsDict[drift] = results[1]

#now also add hill
# HillRPs = [15,2,3]
# for HillRP in HillRPs:
#     print('loading sim results {0} conf'.format(HillRP))
#     savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(500) + 'iter' + date_now + str(conf) + '.pickle'
#
#     file = open(savefile,'rb')
#     results = pickle.load(file)
#     file.close()
#
#     pcDict[conf] = results[0]
#%%
date_now = datetime.datetime.now().strftime('_%m_%d')
# date_now = '_07_25'
rp = params['RPs'][0]*1000
for driftDir in ['Inc','Dec']:
    for frPlot in [1,2,5,10]:
        figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsDrift' + date_now + 'Drift' + driftDir + '_500iterFR' + str(frPlot)
        plotDriftOverlay(pcDict,paramsDict,figsavefile, rpPlot=rp, frPlotInput = frPlot,driftDir = driftDir)



#%%


#%%

#now plot pc as a function of true rp for hill and ours

def plotHillvsSRP(pcDict,params,savefile):
    # (pcSliding,pcHill15,pcHill2,pcHill3,params,savefile, rpPlot=2.5):
    spinesSetting = False
    fig, axs = plt.subplots(1, 1, figsize=(6, 8))
    ax = axs  # for the case of just one subplot
    ax.vlines(10,0,100,'k','dashed',alpha=0.5)
    ax.hlines([0,20,40,60,80,100],0,20,'k','solid',alpha = 0.2)

    pcDict_keys = list(pcDict.keys())

# [type(list(pcDict.keys())[i]) for i in range(len(pcDict.keys()))]



#%%


spinesSetting = False
fig, axs = plt.subplots(1, 2, figsize=(16, 6))

# ax.vlines(10,0,100,'k','dashed',alpha=0.5)
# ax.hlines([0,20,40,60,80,100],0,20,'k','solid',alpha = 0.2)
pcDictKeys = list(pcDict.keys())
pcDictKeysTypes = [type(pcDictKeys[i]) for i in range(len(pcDictKeys))]
numInts = sum([pcDictKeysTypes[i]==int for i in range(len(pcDictKeys))])
pcDictKeys.pop(pcDictKeys.index('Hill 1.5ms')) #remove 1.5 ms from plot

crVec = [8,12] #contamination rates to plot
for c,crPlot in enumerate(crVec):
    ax = axs[c]
    for p, pc_key in enumerate(pcDictKeys):
        pc = pcDict[pc_key]
        print(pc_key)
        count = []
        count = pc / 100 * params['nSim']  # number of correct trials

        CI_scaled = binofit(count, params['nSim'])
        CI = [x * 100 for x in CI_scaled]

        # plot just contRates 0.08 to 0.12:
        cr = params['contRates']
        crInd = np.where(cr == crPlot / 100)[0]  # rpPlot in fraction, convert to percent here

        # plot just RP = rpPlot:
        # rps = params['RPs']
        # rpInd = np.where(rps == rpPlot / 1000)[0]  # rpPlot in ms, convert to s here
        # plot just  recDur = 2:
        recDurPlot = 2
        recDurs = params['recDurs']
        rdInd = np.where(recDurs == recDurPlot)[0]

        # plot just fr = 2:
        frPlot = 2
        frs = np.array(params['baseRates'])
        frInd = np.where(frs == frPlot)[0][0]

        # colors = matplotlib.cm.Set1(np.linspace(0, 1, 10))
        c = cc.linear_bmw_5_95_c89  # input_color  # cc.linear_protanopic_deuteranopic_kbw_5_95_c34
        c = c[::-1]  # if using linear_blue37 or protanopic, flip the order
        # if p==0:
        #     color = [c[x] for x in np.round(np.linspace(0.2, 0.75, len(params['RPs'])) * 255).astype(int)][5]
        # else:
        if type(pc_key) is str:  # Hill
            colors = matplotlib.cm.Reds(np.linspace(0.3, 1, 2))
            color = colors[p-numInts]
        else:  # not hill (conf)
            colors = matplotlib.cm.Blues(np.linspace(0.3, 1, numInts))
            color = colors[p]

        pltcnt = 0
        linewidths = [1, 1, 1, 1, 1, 1]  # [0.5, 1, 2, 3]



        lowerCI = CI[0][rdInd[0], :, frInd, crInd[0]]
        upperCI = CI[1][rdInd[0], :, frInd, crInd[0]]
        x = rps * 1000
        y = pc[rdInd[0], :, frInd,crInd[0]]

        ax.plot(x, y, '.-', color=color, label=rp * 1000)

        ax.fill_between(x, lowerCI, upperCI, color=color, alpha=.3)
        ax.set_ylabel('Percent pass')
        ax.set_xlabel('True RP')
        ax.set_title('Contamination = {0}%'.format(crPlot))
        ax.set_ylim(0,100)
        spinesSetting = False
        ax.spines.right.set_visible(spinesSetting)
        ax.spines.top.set_visible(spinesSetting)


    handles, xx = ax.get_legend_handles_labels()
    labels = pcDictKeys
    labels = ['Flexible RP metric; 70% confidence threshold', 'Flexible RP metric; 75% confidence threshold', 'Flexible RP metric; 80% confidence threshold', 'Hill metric; threshold = 2ms','Hill metric; threshold = 3ms']
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1))

figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\ConfidenceHillOverlay' + date_now + 'Conf_all'
fig.savefig(figsavefile + '_RP.svg', dpi=500)
fig.savefig(figsavefile + '_RP.png', dpi=500)
print('hi')
#now match the confidence

