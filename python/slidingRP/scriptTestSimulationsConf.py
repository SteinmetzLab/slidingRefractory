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
    'recDurs': np.array([0.5, 1 , 2 , 3 ]),  #recording durations (hours)
    'RPs': np.array([0.001,0.002,0.003,0.004,0.005,0.006]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
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

}



#%% run and save
confidence_values = [50,75,90,99]
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

conf=75
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(10) + 'iter' + date_now + str(conf) + '.pickle'

file = open(savefile,'rb')
results = pickle.load(file)
file.close()


for rp in np.arange(1,7):
    figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlayConf' + str(rp) + str(conf)

    plotHillOverlay(results[0],results[0],results[0],results[0],params,figsavefile, rpPlot=rp)






#%% plot
pcDict = {}
confidence_values = [50,75,90]
date_now = '_07_19'
# date_now = datetime.datetime.now().strftime('_%m_%d')

for conf in confidence_values:
    print('loading sim results {0} conf'.format(conf))
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(500) + 'iter' + date_now + str(conf) + '.pickle'

    file = open(savefile,'rb')
    results = pickle.load(file)
    file.close()

    pcDict[conf] = results[0]
#%%
date_now = datetime.datetime.now().strftime('_%m_%d')
for rp in np.arange(1,7):
    figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPCHillOverlay' + date_now + 'Conf' + str(rp) + '500iter'

    # plotHillOverlay(pcDict[confidence_values[0]],pcDict[confidence_values[1]],pcDict[confidence_values[2]],pcDict[confidence_values[3]],params,figsavefile, rpPlot=rp)
    plotHillOverlay(pcDict[confidence_values[0]],pcDict[confidence_values[1]],pcDict[confidence_values[2]],params,figsavefile, rpPlot=rp)



#%%
#to figure out next: why are results with 500 iter not


#%%

#now plot
