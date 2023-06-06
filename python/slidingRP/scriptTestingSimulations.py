# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:29:37 2022

@author: Noam Roth

script for running and plotting various simulations
"""

#%%
import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')


import pickle
import numpy as np
from phylib.stats import correlograms
from scipy import stats
from slidingRP.simulations import *
#%%
#%%
#script testing (to be moved later for import purposes)

# from simulations import genST  # this is not working, for some reason it's importing a different simulations gensT
baseRate = 10 #spks/s
recDur = 3600 #seconds
params = {}
params['checkFR']  = 1
st = genST(meanRate, recDur, xx = None, params = params)

isi = np.diff(np.insert(st,0,0)) #add a zero to the beginning to includ the "first" isi?
isi = np.delete(isi,np.where(isi<rp)[0])
st = np.cumsum(isi)
if c>0:
    contST = genST(contRate,recDur)
else:
    contST=[]
combST = np.sort(np.concatenate((st, contST)))



#%% run (and save) simulations for only one recDur
sampleRate = 30000
params = {
    'recDurs': np.array([0.5, 1 , 2 , 3 ]),  #recording durations (hours)
    'RPs': np.array([0.002]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.5,1,2,5,10],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates': np.arange(0.00,0.21, 0.02),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
    'nSim': 500,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True,
    'savePCfile': True

}

#%% run and save condition with changing firing rates
for delta in [0.5,-0.5,0.1,-0.1]:
    params['delta'] = delta
    print('in simulations:, delta = ', delta)
    pc= simulateChangingContNeurons(params)
    import datetime
    date_now  = datetime.datetime.now().strftime('_%m_%d')
    version = '2' #adjust if running more than once in the same day
    if delta < 0:
        deltaprint = 'neg' + str(int(delta*10))
    else:
        deltaprint = str(int(delta*10))
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPC' + str(params['nSim']) + 'iter' + date_now + version + 'delta' + deltaprint +  '.pickle'
    results = [pc, params]
    if params['savePCfile']:
        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)
#%%
plotfile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPC' + str(params['nSim']) + 'iter' + date_now + version + 'delta' + deltaprint
plotSimulations(pc,params,plotfile)
#%%
#rerun for 5 and 6 ms
sampleRate = 30000
params = {
    'recDurs': np.array([0.5, 1 , 2 , 3 ]),  #recording durations (hours)
    'RPs': np.array([0.001,0.002,0.003,0.004,0.005,0.006]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': [0.5,1,2,5,10],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates': np.arange(0.00,0.21, 0.01),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
    'nSim': 10,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'confidence': 90,
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True,
    'savePCfile': True

}


#%%
print('in simulations')
[pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)
import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')
version = '1' #adjust if running more than once in the same day
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + version +  '.pickle'


results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15,pcHill2, pcHill3, params]

if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)



#also run for different confidences:
params['confidence'] = 80
print('in simulations. 80 conf')
[pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)
import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')
version = '80' #adjust if running more than once in the same day
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + version +  '.pickle'


results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15,pcHill2, pcHill3, params]

if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)




params['confidence'] = 99
print('in simulations, 99 conf')
[pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3] = simulateContNeurons(params)
import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')
version = '99' #adjust if running more than once in the same day
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + version +  '.pickle'


results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15,pcHill2, pcHill3, params]

if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)


#%%

import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')

#from simulationsFunctions import *
import pickle
import numpy as np
from phylib.stats import correlograms
from scipy import stats
from slidingRP.simulations import *

#load data
# savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_11_03_21.pickle'
# savefile= r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_12_26.pickle'
# savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'

file = open(savefile,'rb')
results = pickle.load(file)

file.close()

pc = results[0]
pc2MsNoSpikes = results[1]
pcHalfInactive = results[2]
pcHill2 = results[3]
pcHill3 = results[4]
params = results[5]


#%%
import datetime
#prep for figure saving
date_now  = datetime.datetime.now().strftime('_%m_%d')
#make RP plots for regular and Hill
nSim = 500
figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now
figsavefile2 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes'
figsavefile3 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill2'
figsavefile4 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill3'

plotSimulations(pc, params,figsavefile1)
plotSimulations(pc2MsNoSpikes, params, figsavefile2)
plotSimulations(pcHill2,params,figsavefile3)
plotSimulations(pcHill3,params,figsavefile4)
#NEXT, try this with zoomed in on 7-13 contamination ish (or 8-12)

#%%
if False:
    #put together two sets of results:
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_11_03_21.pickle'
    file = open(savefile, 'rb')
    results = pickle.load(file)
    file.close()

    pc1 = results[0]
    pc2MsNoSpikes1 = results[1]
    params1 = results[2]

    savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC100iter_12_11.pickle'
    file = open(savefile, 'rb')
    results = pickle.load(file)
    file.close()

    pc2 = results[0]
    pc2MsNoSpikes2 = results[1]
    params2 = results[2]

    pcc1 = np.concatenate((pc1[:,:,0:3,:], pc2[:,:,0:1,:]),2)
    pcc2 = np.concatenate((pcc1, pc1[:,:,3:4,:]),2)
    pcc3 = np.concatenate((pcc2,pc2[:,:,1:4,:]),2)
    pcc4 = np.concatenate((pcc3,pc1[:,:,4:5,:]),2)
    pcc5 = np.concatenate((pcc4,pc2[:,:,4:7,:]),2)
    pcc6 = np.concatenate((pcc5,pc1[:,:,5:7,:]),2)
    pc = pcc6


    pc1 = pc2MsNoSpikes1
    pc2 = pc2MsNoSpikes2


    pcc1 = np.concatenate((pc1[:,:,0:3,:], pc2[:,:,0:1,:]),2)
    pcc2 = np.concatenate((pcc1, pc1[:,:,3:4,:]),2)
    pcc3 = np.concatenate((pcc2,pc2[:,:,1:4,:]),2)
    pcc4 = np.concatenate((pcc3,pc1[:,:,4:5,:]),2)
    pcc5 = np.concatenate((pcc4,pc2[:,:,4:7,:]),2)
    pcc6 = np.concatenate((pcc5,pc1[:,:,5:7,:]),2)
    pc2MsNoSpikes = pcc6

    params = params1
    params['baseRates'] = [ 0.1 ,  0.25,  0.5 ,0.75,  1.  , 2, 3, 4, 5.  , 7.5, 10.  , 20.  ]  # [0.75, 2, 3, 4, 7.5]

if False:
    #put together two sets of results:
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_12_26.pickle'
    file = open(savefile, 'rb')
    results = pickle.load(file)
    file.close()

    pc2 = results[0]
    pc2MsNoSpikes2 = results[1]
    params2 = results[2]

    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_12_262.pickle'
    file = open(savefile, 'rb')
    results = pickle.load(file)
    file.close()

    pc1 = results[0]
    pc2MsNoSpikes1 = results[1]
    params1 = results[2]

    pcc1 = np.concatenate((pc1, pc2[:,:,1:,:]),2)
    pc = pcc1

    pc1 = pc2MsNoSpikes1
    pc2 = pc2MsNoSpikes2
    pcc1 = np.concatenate((pc1, pc2[:,:,1:,:]),2)
    pc2MsNoSpikes = pcc1


    params = params1
    params['baseRates'] = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1, 2,  5,  10, 20]  # [0.75, 2, 3, 4, 7.5]


RP = 0.0015
figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_15msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now
figsavefile2 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_15msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes'
figsavefile3 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_15msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill2'
figsavefile4 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_15msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill3'
figsavefile5 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_15msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'FRHalf'

import colorcet as cc
plotSimulations(pc, params,figsavefile1,rp_valFig1=RP)
plotSimulations(pc2MsNoSpikes, params, figsavefile2,rp_valFig1=RP)
plotSimulations(pcHill2,params,figsavefile3, rp_valFig1=RP,input_color=cc.linear_tritanopic_krw_5_95_c46)
plotSimulations(pcHill3,params,figsavefile4,rp_valFig1 = RP,input_color=cc.linear_tritanopic_krw_5_95_c46)
plotSimulations(pcHalfInactive, params, figsavefile5,rp_valFig1 = RP)

RP = 0.004
figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_4msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now
figsavefile2 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_4msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes'
figsavefile3 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_4msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill2'
figsavefile4 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_4msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill3'
figsavefile5 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\results_122_4msRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'FRHalf'

plotSimulations(pc, params,figsavefile1,rp_valFig1=RP)
plotSimulations(pc2MsNoSpikes, params, figsavefile2,rp_valFig1=RP)
plotSimulations(pcHill2,params,figsavefile3, rp_valFig1=RP,input_color=cc.linear_tritanopic_krw_5_95_c46)
plotSimulations(pcHill3,params,figsavefile4,rp_valFig1 = RP,input_color=cc.linear_tritanopic_krw_5_95_c46)
plotSimulations(pcHalfInactive, params, figsavefile5,rp_valFig1 = RP)


#%%
#replotRecDur
figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now

plotSimulations(pc, params,figsavefile1,rp_valFig1=RP)

#%%
#replot RP
RP = 0.015
figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now

plotSimulations(pc, params,figsavefile1,rp_valFig1=RP)

#%%



figsavefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(params['nSim']) + 'iter' + date_now
# plotSensitivitySpecificity(pc,pc2MsNoSpikes,params, figsavefile,5)


#%%


import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')

import pickle
import numpy as np
from phylib.stats import correlograms
from scipy import stats
from slidingRP.simulations import *

#put all RP runs together and plot RP plot (magenta plot)
# put together two sets of results:
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'
file = open(savefile, 'rb')
results = pickle.load(file)
file.close()

pc1 = results[0]
pc2MsNoSpikes1 = results[1]
pcHalfInactive1 = results[2]
pcHill21 = results[3]
pcHill31 = results[4]
params1 = results[5]

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_241.pickle' #this one has 1, 5, 6
file = open(savefile, 'rb')
results = pickle.load(file)
file.close()

pc2 = results[0]
pc2MsNoSpikes2 = results[1]
pcHalfInactive2 = results[2]
pcHill22 = results[3]
pcHill32 = results[4]
params2 = results[5]

pcc1 = np.concatenate((pc1, pc2[:, 1:, :, :]), 1)
pc = pcc1

pc1 = pc2MsNoSpikes1
pc2 = pc2MsNoSpikes2
pcc1 = np.concatenate((pc1, pc2[:, 1:, :, :]), 1)
pc2MsNoSpikes = pcc1

pc1 = pcHill21
pc2 = pcHill22
pcc1 = np.concatenate((pc1, pc2[:, 1:, :, :]), 1)
pcHill2 = pcc1

pc1 = pcHill31
pc2 = pcHill32
pcc1 = np.concatenate((pc1, pc2[:, 1:, :, :]), 1)
pcHill3 = pcc1

params = params1
params['RPs'] = np.array( [0.0015, 0.002 , 0.003 , 0.004 , 0.005,0.006])
#%%

import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')

figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now
figsavefile2 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes'
figsavefile3 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill2'
figsavefile4 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill3'

plotSimulations(pc, params,figsavefile1)
plotSimulations(pc2MsNoSpikes, params, figsavefile2)
plotSimulations(pcHill2,params,figsavefile3)
plotSimulations(pcHill3,params,figsavefile4)


figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'zoom'
figsavefile2 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes' + 'zoom'
figsavefile3 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill2'+ 'zoom'
figsavefile4 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\RPmagenta\simulationsPC' + str(params['nSim']) + 'iter' + date_now + 'Hill3'+ 'zoom'

plotSimulations(pc, params,figsavefile1,zoomCont=True)
plotSimulations(pc2MsNoSpikes, params, figsavefile2,zoomCont=True)
plotSimulations(pcHill2,params,figsavefile3,zoomCont=True)
plotSimulations(pcHill3,params,figsavefile4,zoomCont=True)



#%%
#plot main plot for increase decrease
nSim = 500
import datetime
#prep for figure saving
date_now  = datetime.datetime.now().strftime('_%m_%d')

version = '2'  # adjust if running more than once in the same day

for delta in [0.5,-0.5,0.1,-0.1]:

    if delta < 0:
        deltaprint = 'neg' + str(int(delta*10))
    else:
        deltaprint = str(int(delta*10))
    print(deltaprint)
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPC500iter_01_242delta' + deltaprint + '.pickle'
    file = open(savefile,'rb')
    results = pickle.load(file)
    file.close()
    pc = results[0]
    params = results[1]
    #make RP plots for regular and Hill
    figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPCdelta' + deltaprint + str(params['nSim']) + 'iter' + date_now
    plotSimulations(pc, params,figsavefile1, Fig1=True, Fig2=False, Fig3=False, Fig4=False,    input_color = cc.linear_green_5_95_c69)






#%%
#plot just regular first one and zoom
import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')

import pickle
import numpy as np
from phylib.stats import correlograms
from scipy import stats
from slidingRP.simulations import *

#put all RP runs together and plot RP plot (magenta plot)
# put together two sets of results:
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'
file = open(savefile, 'rb')
results = pickle.load(file)
file.close()

pc = results[0]

params = results[5]
import datetime
date_now  = datetime.datetime.now().strftime('_%m_%d')

figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs\simulationsPCAll' + str(params['nSim']) + 'iter' + date_now
plotSimulations(pc, params,figsavefile1,frPlot = [0.5,1,5,10],zoomCont=False, Fig1=True, Fig2=False, Fig3=False, Fig4=False)

figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs\simulationsPC10' + str(params['nSim']) + 'iter' + date_now
plotSimulations(pc, params,figsavefile1,frPlot = [10],zoomCont=False, Fig1=True, Fig2=False, Fig3=False, Fig4=False)

#now zoom

figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs\simulationsPCAll' + str(params['nSim']) + 'iter' + date_now +'zoom'
plotSimulations(pc, params,figsavefile1,frPlot = [0.5,1,5,10],zoomCont=True, Fig1=True, Fig2=False, Fig3=False, Fig4=False)

figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs\simulationsPC10' + str(params['nSim']) + 'iter' + date_now + 'zoom'
plotSimulations(pc, params,figsavefile1,frPlot = [10],zoomCont=True, Fig1=True, Fig2=False, Fig3=False, Fig4=False)




#%%
#run to look at




#%%
# savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pickle'
#
# sampleRate = 30000
# params = {
#     'recDurs': np.array([0.25, 0.5, 1, 2, 5]),  #recording durations (hours)
#     'RPs': np.array([0.001,0.0015, 0.002, 0.003, 0.004, 0.005]), #true RP (s)
#     'baseRates': np.array([ 0.5, 1, 2, 5, 10, 20 ]), #F1, 2, 5, 10 , 20 R (spk/s)
#     'contRates':  np.arange(0.00,0.225, 0.025), #contamination levels (proportion) #.025
#     'nSim': 20,
#     'threshold': 0.1,
#     'binSize': 1 / sampleRate,
#     'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
#     'checkFR': False,
#     'binSizeCorr': 1 / sampleRate,
#     'returnMatrix': True,
#     'verbose': True ,
#     'savePCfile': True
#
# }
#
#
# file = open(savefile,'rb')
# pc = pickle.load(file)
# file.close()
#
# figsavefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pdf'
# plotSimulations(pc, params, figsavefile)
#
#
#
# #%%
# #test a closed form expression for probability of passing
# P = 0.1 #true contamination
# R = 1
# r = 0.002 #true frefractory period
# D = 2*3600; #recording duration
#
#
# B = S*R*P*r
# S = R*D
#
# lam = P* R * D * S * S
# C = 1-stats.poisscdf(B, lam)
#
# #%%
#
# #define params
# contRate = 0.1
# rp = 0.002
# baseRate = 4
#
# recDur = 2*3600
# nSim = 1000
#
#
# rpBinSize = 1 / 30000
# rpEdges = np.arange(0, 10/1000, rpBinSize) # in s
# rpVec = rpEdges + np.mean(np.diff(rpEdges)[0])/2 #
#
# obsViol = np.empty(nSim)
# obsViol[:] = np.nan
# expectedViol = np.empty(nSim)
# expectedViol[:] = np.nan
# for i in range(nSim):
#     #generate contaminated spiketrain
#     st = genST(baseRate,recDur,params) #generate a spike train with the current base rate
#     isi = np.diff(np.insert(st,0,0))
#     isi = np.delete(isi,np.where(isi<rp)[0]) #get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
#     st = np.cumsum(isi)
#     if contRate>0: #contamination greater than 0
#         contST = genST(contRate*baseRate,recDur,params)
#     else:
#         contST=[]
#     combST = np.sort(np.concatenate((st, contST))) # put spike times from both spike trains together (and sort chronologically)
#
#     st = combST
#     spikeCount = len(st)
#
#
#     #compute firing rate
#     clustersIds = [0] # call the cluster id 0 (not used, but required input for correlograms)
#     spikeClustersACG = np.zeros(len(st), dtype = 'int8') # each spike time gets cluster id 0
#     nACGFR = correlograms(st, spikeClustersACG, cluster_ids = clustersIds, bin_size = 1, sample_rate = 30000, window_size=2,symmetrize=False)[0][0] #compute acg
#     firingRate = nACGFR[1] / spikeCount
#
#     #compute ACG
#     nACG = correlograms(st, spikeClustersACG, cluster_ids = clustersIds, bin_size = params['binSizeCorr'], sample_rate = params['sampleRate'], window_size=2,symmetrize=False)[0][0] #compute acg
#
#     # compute observed violations
#     rpIdx = np.where(rpVec > rp)[0][0]
#     obsViol[i] = sum(nACG[0:rpIdx])
#     # print('obsViol: ', obsViol)
#
#
#     expectedViol[i] = firingRate * contRate * rp * 2 * spikeCount
#     # print('expectedViol: ', expectedViol)
#
# #%%
#
# #%% test the 'expected count'
#
# baseRates = np.logspace(-0.6,1.5,80)
# recDur = 7200;
# rp = 0.002;
# contProp = 0.1;
# # params = struct(); params.cont = 10;
# obsViol = np.empty(len(baseRates))
# obsViol[:] = np.nan
# expViol = np.empty(len(baseRates))
# expViol[:] = np.nan
# for bidx in range(len(baseRates)):
#     baseRate = baseRates[bidx];
#     contRate = baseRate*contProp;
#
#     st = genST(baseRate,recDur,params);
#     isi = np.diff(np.insert(st,0,0))
#     isi = np.delete(isi,np.where(isi<rp)[0]) #get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
#     st = np.cumsum(isi)
#
#
#
#     contST = genST(contRate,recDur,params)
#
#     combST = np.sort(np.concatenate((st, contST))) # put spike times from both spike trains together (and sort chronologically)
#
#     [ confMatrix, cont, rpTestVals, nACG, firingRate]= computeMatrix(combST, params)
#
#     nACG = nACG[0:len(rpTestVals)]
#     contaminationRate = firingRate*contProp;
#     expectedViol = contaminationRate*rp*2*len(combST);
#
#     obsViol[bidx] = sum(nACG[rpTestVals<=0.002]);
#     expViol[bidx] = expectedViol;
#
# #%%
# fig,ax= plt.subplots(1,1)
# ax.plot(baseRates, obsViol, 'r.-',label = 'observed violations');
# ax.plot(baseRates, expViol, 'k.-',label = 'expected violations');
# ax.set_xlabel('Base firing rate (spks/s)')
# ax.set_ylabel('Number of violations');
# ax.spines.right.set_visible(False)
# ax.spines.top.set_visible(False)
# fig.legend()
# #%%
#
#
#
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
# from scipy.stats import poisson
#
# bins = np.arange(25) - 0.5
# entries, bin_edges, patches = plt.hist(obsViol, bins=bins, density=True, label='Data')
#
# # calculate bin centres
# bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#
#
# def fit_function(k, lamb):
#     '''poisson function, parameter lamb is the fit parameter'''
#     return poisson.pmf(k, lamb)
#
#
# # fit with curve_fit
# parameters, cov_matrix = curve_fit(fit_function, bin_middles, entries)
# x_plot = np.arange(0, 25)
#
# plt.plot(
#     x_plot,
#     fit_function(x_plot, *parameters),
#     marker='o', linestyle='',
#     label='Fit result',
# )
# plt.legend()
# plt.title('Observed Viols (mean: %.2f) over %d simulations. Mean of ExpViol: %.2f'%(parameters[0],nSim, np.mean(expectedViol)))
#
# plt.show()
#
#
# #%%
# plt.plot([0.5,1,2,3,4],[1.66,3.18,6.21,9.88,13.57],'r.-',label = 'observed violations')
# plt.plot([0.5,1,2,3,4],[1.03, 3.47,12.61,27.37,47.66], 'k.-', label = 'expected violations')
# plt.xlabel('Firing Rate')
# plt.ylabel('Number of spikes (averaged across simulations)')
# plt.title('Number of violations: contRate = 0.1; rp = 0.002; recDur = 2hr')
# plt.legend()
#
#
# #%%
#
# fig,axs = plt.subplots(1,3)
# axs[0].hist(obsViol)
# axs[0].set_title('obsViol')
# axs[1].hist(expectedViol)
# axs[1].set_title('expectedViol')
#
# conf = 1 - stats.poisson.cdf(obsViol, expectedViol)
# axs[2].hist(conf)
# axs[2].set_title('Confidence')
#
#
#
#
# percentPass = sum(conf>confThresh)/len(conf)
#
#
# #now test Nick's code
# PLimit = 0.1
# PTrue = contRate
# R = firingRate
# r = rp
# S = spikeCount
# confLimit = 0.1
#
# S = firingRate*recDur#R*D; % spike count
# violLimit = PLimit*R*r*2*S; # mean count of a poisson process at the limit of acceptable contamination
# violTrue = PTrue*R*r*2*S; # mean count of poisson process at the true contamination --> this one is not right
#
#
# # % what is the value of the limit contamination's CDF that would exceed the
# # % acceptable confidence? I.e. when does the limit CDF = 10%? this is the
# # % max acceptable.
# k = stats.poisson.ppf(1-confLimit, violLimit);
#
# # % now what is the value of the true CDF just below that value? that's the
# # % proportion of simulations that should pass
# propPass = stats.poisson.cdf(k-1, violTrue); % k can be negative, in which case it returns zero