import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')


import pickle
import numpy as np
from phylib.stats import correlograms

from scipy import stats
from slidingRP.simulations import *
import datetime
from statsmodels.stats.proportion import proportion_confint as binofit


date_now = datetime.datetime.now().strftime('_%m_%d')
#%%
#run simulations:
sampleRate = 30000
params = {
    'recDurs':np.array([0.5, 1 , 2 , 3 ]),  #recording durations (hours) np.array([0.5, 1 , 2 , 3 ])
    'RPs': np.array([0.0015,0.002,0.006])#np.array([0.001,0.0015, 0.002,0.0025, 0.003,0.004,0.005,0.006]),# , np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.5,1,2,5,10],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
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
    'runLlobet': True

}



#%% run and save
confidence_values = [70,80,90]
for conf in confidence_values:
    params['confidenceThresh'] = conf
    print('in simulations, {0} conf'.format(conf))
    [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)

    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
        params['nSim']) + 'iter' + date_now + str(conf) + '.pickle'

    results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, params]
    if params['savePCfile']:
        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)


#%% test plot just one
date_now = datetime.datetime.now().strftime('_%m_%d')
date_now = '_09_04'
nIter = 1000
conf=80
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(nIter) + 'iter' + date_now + str(conf) + '.pickle'
         #%%
file = open(savefile,'rb')
results = pickle.load(file)
params = results[-1]
file.close()


for rp in np.arange(1,7):
    figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlayConf' + str(rp) + str(conf) + date_now

    plotHillOverlay(results[0],results[0],results[0],results[0],params,figsavefile, rpPlot=rp)



#%%plot just one
pcDict = {}
date_now = '_09_04'
nIter = 1000

frPlot = 2
rpPlot = 2

conf=80
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_08_2290.pickle'
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(nIter) + 'iter' + date_now + str(conf) + '.pickle'


file = open(savefile,'rb')
results = pickle.load(file)
file.close()

pcDict[0] = results[0]
params = results[-1]

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_07_1990.pickle'
conf = 70
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(nIter) + 'iter' + date_now + str(conf) + '.pickle'


file = open(savefile,'rb')
results = pickle.load(file)
file.close()

pcDict[1] = results[0]

pcDict['Hill 1.5ms'] = results[3]
pcDict['Hill 2ms'] = results[4]
pcDict['Hill 3ms'] = results[5]

figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlayConf7080_newcalcCompare'

plotHillOverlay(pcDict, params, figsavefile, rpPlot=rpPlot, frPlot = frPlot, legendLabels=['80','70','Hill 1.5','Hill 2', 'Hill3','Confidence'])

#%% plot
pcDict = {}
confidence_values = [50,60,70,75,80,85,90,99]
dates = ['_07_19','_07_25', '_07_25', '_07_19','_07_25', '_07_25', '_07_19', '_07_19']
#
# confidence_values = [70,75,80]
# dates = [ '_07_25', '_07_19','_07_25']

confidence_values = [90]

#plot just 90
# confidence_values = [90]
# dates = ['_07_19']

# date_now = datetime.datetime.now().strftime('_%m_%d')

for conf, date_now in zip(confidence_values,dates):
    print('loading sim results {0} conf'.format(conf))
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(500) + 'iter' + date_now + str(conf) + '.pickle'

    file = open(savefile,'rb')
    results = pickle.load(file)
    file.close()

    pcDict[conf] = results[0]


simulationsPC500iter_07_1990

pcDict['Hill 1.5ms'] = results[3]
pcDict['Hill 2ms'] = results[4]
pcDict['Hill 3ms'] = results[5]

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
for rp in np.arange(1,7):
    figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlay' + date_now + 'Conf' + str(rp) + '500iter'

    # plotHillOverlay(pcDict[confidence_values[0]],pcDict[confidence_values[1]],pcDict[confidence_values[2]],pcDict[confidence_values[3]],params,figsavefile, rpPlot=rp)
    # plotHillOverlay(pcDict[confidence_values[0]],pcDict[confidence_values[1]],pcDict[confidence_values[2]],params,figsavefile, rpPlot=rp)

    plotHillOverlay(pcDict,params,figsavefile, rpPlot=rp)



#%%

#determine AUC

contAvg = False #use the average across all contaminations on each side of contThresh (uncontaminated, or contaminated)
    #otherwise, use a fixed point
if contAvg ==False :
    threshDist = 0.02 #2%

rpPlot = 2
pcDict_keys = list(pcDict.keys())

pcDict_keys = pcDict_keys[0:-2]


threshDists = np.arange(0.02,0.1,0.02)
from matplotlib.pyplot import cm
color = iter(cm.rainbow(np.linspace(0, 1, len(threshDists))))

fig, axs = plt.subplots(1, 1, figsize=(6, 6))
ax = axs  # for the case of just one subplot
for threshDist in np.arange(0.02,0.1,0.02):
    c = next(color)
    print(c)
    TPR = []
    FPR = []
    for p,pc_key in enumerate(pcDict_keys):#,pcHill3]):
        pc = pcDict[pc_key]
        count = []
        count = pc / 100 * params['nSim']  # number of correct trials

        CI_scaled = binofit(count, params['nSim'])
        CI = [x * 100 for x in CI_scaled]

        # plot just contRates 0.08 to 0.12:
        cr = params['contRates'];

        # plot just RP = rpPlot:
        rps = params['RPs']
        rpInd = np.where(rps == rpPlot/1000)[0] #rpPlot in ms, convert to s here
        #plot just  recDur = 2:
        recDurs = params['recDurs']
        rdInd = np.where(recDurs == 2)[0]

        # plot just fr = 5:
        frs = np.array(params['baseRates'])
        frInd = np.where(frs == 2)[0][0]
        print('Firing rate is 2')


        lowerCI = CI[0][rdInd[0], rpInd[0], frInd, :]
        upperCI = CI[1][rdInd[0], rpInd[0], frInd, :]

        x = cr * 100
        y = pc[rdInd[0], rpInd[0], frInd, :]

        contThresh = params['contaminationThresh'] / 100  # convert to percent as in contamination rate
        if contAvg:
            uncontPerf = np.nanmean(y[cr < contThresh])
            contPerf = np.nanmean(y[cr > contThresh])
        else:
            uncontPerf = y[cr == np.round((contThresh - threshDist),2)]
            contPerf = y[cr == np.round((contThresh + threshDist),2)]
        TPR.append(uncontPerf)
        FPR.append(contPerf)

    ax.plot(FPR[0:-1],TPR[0:-1],'.', color = c,label = threshDist)
    ax.plot(FPR[-1],TPR[-1],'x', color = c)

    print(c)

ax.set_ylabel('True Positive Rate')
ax.set_xlabel('False Positive Rate')
ax.set_title('ROC')
ax.set_ylim(0,100)
ax.set_xlim(0,100)


handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(.7, .4), title = 'distance from contThresh')


ax.plot(0, 0, '.', color='k', label='slidingRP')
ax.plot(0,0, 'x', color='k', label='Hill')
handles, labels = ax.get_legend_handles_labels()
handles = handles[-2:]
labels = labels[-2:]
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(.9, .4), title = 'metric')


fig.show()


#%%
#compute AUC




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
fig, axs = plt.subplots(1, 2, figsize=(12, 4))

# ax.vlines(10,0,100,'k','dashed',alpha=0.5)
# ax.hlines([0,20,40,60,80,100],0,20,'k','solid',alpha = 0.2)
pcDictKeys = list(pcDict.keys())

pcDictKeys.pop(pcDictKeys.index('Hill 1.5ms')) #remove 1.5 ms from plot
pcDictKeys.pop(pcDictKeys.index(50)) #remove 1.5 ms from plot
pcDictKeys.pop(pcDictKeys.index(60)) #remove 1.5 ms from plot
pcDictKeys.pop(pcDictKeys.index(85)) #remove 1.5 ms from plot
pcDictKeys.pop(pcDictKeys.index(90)) #remove 1.5 ms from plot
pcDictKeys.pop(pcDictKeys.index(99)) #remove 1.5 ms from plot
pcDictKeysTypes = [type(pcDictKeys[i]) for i in range(len(pcDictKeys))]
numInts = sum([pcDictKeysTypes[i]==int for i in range(len(pcDictKeys))])

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

