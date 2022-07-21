# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:29:37 2022

@author: Noam Roth

script for running and plotting various simulations
"""

import simulations
import pickle

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



#%%
sampleRate = 30000
params = {
    'recDurs': np.array([1, 2, 5]),#np.array([0.25, 0.5, 1, 2, 5]),  #recording durations (hours)
    'RPs': np.array([0.001, 0.002, 0.003]), #np.array([0.001,0.0015, 0.002, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': np.array([ 0.25, 0.5, 1, 5, 10]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates':  np.arange(0.00,0.225, 0.05), #contamination levels (proportion) #.025
    'nSim': 10,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True ,
    'savePCfile': True
    
}
#%%
[pc, pc2MsNoSpikes] = simulateContNeurons(params)
import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pickle'


results = [pc, pc2MsNoSpikes, params]
if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(pc, handle)
     
        
     
        #%%
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pickle'
#load data
file = open(savefile,'rb')
results = pickle.load(file)
file.close()

#prep for figure saving
figsavefile1 = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pdf'
figsavefile2 = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '2MsNoSpikes.pdf'


plotSimulations(pc, params, figsavefile1)
plotSimulations(pc2MsNoSpikes, params, figsavefile2)
#%%
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pickle'

sampleRate = 30000
params = {
    'recDurs': np.array([0.25, 0.5, 1, 2, 5]),  #recording durations (hours)
    'RPs': np.array([0.001,0.0015, 0.002, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': np.array([ 0.5, 1, 2, 5, 10, 20 ]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates':  np.arange(0.00,0.225, 0.025), #contamination levels (proportion) #.025
    'nSim': 20,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True ,
    'savePCfile': True
    
}


file = open(savefile,'rb')
pc = pickle.load(file)
file.close()

figsavefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pdf'
plotSimulations(pc, params, figsavefile)




