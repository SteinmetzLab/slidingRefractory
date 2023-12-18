#goal here is to test different versions of Hill/Llobet/sliding to check
#what the performance is for 10% or greater contamination
#hypothesis is that at 10% or greater contamination, we see 10% percent accept, but this
#seems to not be happening

import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')


import pickle
import numpy as np
from phylib.stats import correlograms

from scipy import stats

import datetime
from statsmodels.stats.proportion import proportion_confint as binofit

from slidingRP.simulations import *
date_now = datetime.datetime.now().strftime('_%m_%d')
#%%
#run simulations:
sampleRate = 30000
params = {
    'recDurs':np.array([2 ]),  #recording durations (hours) np.array([0.5, 1 , 2 , 3 ])
    'RPs': np.array([0.0015, 0.002,0.0025]),# , np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.5,1,2,5,10],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates': np.arange(0.00,0.21, 0.02),#np.arange(0.00,0.21, 0.02),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
    'nSim': 500,
    'contaminationThresh': 10,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,
    'confidenceThresh': 90,
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True,
    'savePCfile': True,
    'runLlobet': True,
    'runLlobetPoiss': True
}



#%% run and save
confidence_values = [90]
for conf in confidence_values:
    params['confidenceThresh'] = conf
    print('in simulations, {0} conf'.format(conf))
    [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
     pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3] = simulateContNeurons(params)

    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
        params['nSim']) + 'iter' + date_now + str(conf) + '.pickle'

    results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
               pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3, params]
    if params['savePCfile']:
        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)


#%%plot just one
pcDict = {}
date_now = '_10_25'
# nIter = 500
nIter = params['nSim']

frPlot = 2
rpPlot = 2

#load
conf = 90
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(nIter) + 'iter' + date_now + str(conf) + '.pickle'

file = open(savefile,'rb')
results = pickle.load(file)
file.close()
pcDict[0] = results[0]

# pcDict['Hill 1.5ms'] = results[3]
# pcDict['Hill 2ms'] = results[4]
# pcDict['Hill 3ms'] = results[5]

pcDict['Llobet 1.5ms'] = results[6]
pcDict['Llobet 2ms'] = results[7]
pcDict['Llobet 3ms'] = results[8]

pcDict['Llobet Poiss 1.5ms'] = results[6]
pcDict['Llobet Pioss 2ms'] = results[7]
pcDict['Llobet Poiss 3ms'] = results[8]

figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlay_newCalcCompare' + date_now

plotHillOverlay(pcDict, params, figsavefile, rpPlot=rpPlot, frPlot = frPlot, legendLabels=['90',
                                                                                           # 'Hill 1.5', 'Hill 2', 'Hill 3',
                                                                                           'Llobet 1.5','Llobet 2','Llobet 3',
                                                                                           'LlobetP 1.5','LlobetP 2','LlobetP 3','Confidence/Metric'])
