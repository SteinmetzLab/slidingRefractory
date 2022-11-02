# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:29:37 2022

@author: Noam Roth

script for running and plotting various simulations
"""
import sys
sys.path.append(r'C:\Users\Steinmetz Lab User\int-brain-lab\phylib')

#from simulationsFunctions import *
import pickle
import numpy as np
from phylib.stats import correlograms
from types import SimpleNamespace
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import time
from slidingRP.simulations import *

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
    'recDurs': np.array([1]),#%0.25, 0.5, 1, 2,3,4, 5]),  #recording durations (hours)
    'RPs': np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': np.array([0.1, 0.25, 0.5, 1, 5, 10,20]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates':  np.arange(0.00,0.225, 0.01),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.225, 0.01), #contamination levels (proportion) #.025
    'nSim': 50,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True ,
    'savePCfile': True
    
}

[pc, pc2MsNoSpikes] = simulateContNeurons(params)
import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pickle'


results = [pc, pc2MsNoSpikes, params]
if params['savePCfile']:
    with open(savefile, 'wb') as handle:
        pickle.dump(results, handle)
     
        #%% load saved results and plot
# savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pickle'
from slidingRP import simulations
from simulations import plotSimulations
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter_10_14.pickle'
#load data

import datetime

date_now  = datetime.datetime.now().strftime('_%m_%d')

file = open(savefile,'rb')
results = pickle.load(file)
file.close()

pc = results[0]
pc2MsNoSpikes = results[1]
params = results[2]


#prep for figure saving
figsavefile1 = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.svg'
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



#%%
#test a closed form expression for probability of passing
P = 0.1 #true contamination
R = 1
r = 0.002 #true frefractory period
D = 2*3600; #recording duration


B = S*R*P*r
S = R*D

lam = P* R * D * S * S 
C = 1-stats.poisscdf(B, lam)

#%%

#define params
contRate = 0.1
rp = 0.002
baseRate = 4

recDur = 2*3600
nSim = 1000


rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rpVec = rpEdges + np.mean(np.diff(rpEdges)[0])/2 #

obsViol = np.empty(nSim)
obsViol[:] = np.nan
expectedViol = np.empty(nSim)
expectedViol[:] = np.nan
for i in range(nSim):
    #generate contaminated spiketrain 
    st = genST(baseRate,recDur,params) #generate a spike train with the current base rate
    isi = np.diff(np.insert(st,0,0)) 
    isi = np.delete(isi,np.where(isi<rp)[0]) #get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
    st = np.cumsum(isi)
    if contRate>0: #contamination greater than 0 
        contST = genST(contRate*baseRate,recDur,params)
    else:
        contST=[]
    combST = np.sort(np.concatenate((st, contST))) # put spike times from both spike trains together (and sort chronologically)
    
    st = combST
    spikeCount = len(st)
    
    
    #compute firing rate 
    clustersIds = [0] # call the cluster id 0 (not used, but required input for correlograms)
    spikeClustersACG = np.zeros(len(st), dtype = 'int8') # each spike time gets cluster id 0 
    nACGFR = correlograms(st, spikeClustersACG, cluster_ids = clustersIds, bin_size = 1, sample_rate = 30000, window_size=2,symmetrize=False)[0][0] #compute acg
    firingRate = nACGFR[1] / spikeCount
      
    #compute ACG
    nACG = correlograms(st, spikeClustersACG, cluster_ids = clustersIds, bin_size = params['binSizeCorr'], sample_rate = params['sampleRate'], window_size=2,symmetrize=False)[0][0] #compute acg
  
    # compute observed violations  
    rpIdx = np.where(rpVec > rp)[0][0]
    obsViol[i] = sum(nACG[0:rpIdx])
    # print('obsViol: ', obsViol)
    
    
    expectedViol[i] = firingRate * contRate * rp * 2 * spikeCount 
    # print('expectedViol: ', expectedViol)
    
#%%

#%% test the 'expected count'

baseRates = np.logspace(-0.6,1.5,80)
recDur = 7200; 
rp = 0.002;
contProp = 0.1;
# params = struct(); params.cont = 10;
obsViol = np.empty(len(baseRates))
obsViol[:] = np.nan
expViol = np.empty(len(baseRates))
expViol[:] = np.nan
for bidx in range(len(baseRates)):
    baseRate = baseRates[bidx];
    contRate = baseRate*contProp;
    
    st = genST(baseRate,recDur,params);  
    isi = np.diff(np.insert(st,0,0)) 
    isi = np.delete(isi,np.where(isi<rp)[0]) #get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
    st = np.cumsum(isi)
    
    
    
    contST = genST(contRate,recDur,params)
    
    combST = np.sort(np.concatenate((st, contST))) # put spike times from both spike trains together (and sort chronologically)
    
    [ confMatrix, cont, rpTestVals, nACG, firingRate]= computeMatrix(combST, params)
    
    nACG = nACG[0:len(rpTestVals)]
    contaminationRate = firingRate*contProp; 
    expectedViol = contaminationRate*rp*2*len(combST);
    
    obsViol[bidx] = sum(nACG[rpTestVals<=0.002]); 
    expViol[bidx] = expectedViol;

#%%
fig,ax= plt.subplots(1,1)
ax.plot(baseRates, obsViol, 'r.-',label = 'observed violations'); 
ax.plot(baseRates, expViol, 'k.-',label = 'expected violations');
ax.set_xlabel('Base firing rate (spks/s)')
ax.set_ylabel('Number of violations'); 
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
fig.legend()
#%%



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson

bins = np.arange(25) - 0.5
entries, bin_edges, patches = plt.hist(obsViol, bins=bins, density=True, label='Data')

# calculate bin centres
bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])


def fit_function(k, lamb):
    '''poisson function, parameter lamb is the fit parameter'''
    return poisson.pmf(k, lamb)


# fit with curve_fit
parameters, cov_matrix = curve_fit(fit_function, bin_middles, entries)
x_plot = np.arange(0, 25)

plt.plot(
    x_plot,
    fit_function(x_plot, *parameters),
    marker='o', linestyle='',
    label='Fit result',
)
plt.legend()
plt.title('Observed Viols (mean: %.2f) over %d simulations. Mean of ExpViol: %.2f'%(parameters[0],nSim, np.mean(expectedViol)))

plt.show()


#%%
plt.plot([0.5,1,2,3,4],[1.66,3.18,6.21,9.88,13.57],'r.-',label = 'observed violations')
plt.plot([0.5,1,2,3,4],[1.03, 3.47,12.61,27.37,47.66], 'k.-', label = 'expected violations')
plt.xlabel('Firing Rate')
plt.ylabel('Number of spikes (averaged across simulations)')
plt.title('Number of violations: contRate = 0.1; rp = 0.002; recDur = 2hr')
plt.legend()


#%%

fig,axs = plt.subplots(1,3)
axs[0].hist(obsViol)
axs[0].set_title('obsViol')
axs[1].hist(expectedViol)
axs[1].set_title('expectedViol')

conf = 1 - stats.poisson.cdf(obsViol, expectedViol)
axs[2].hist(conf)
axs[2].set_title('Confidence')




percentPass = sum(conf>confThresh)/len(conf)


#now test Nick's code
PLimit = 0.1
PTrue = contRate
R = firingRate
r = rp
S = spikeCount
confLimit = 0.1

S = firingRate*recDur#R*D; % spike count
violLimit = PLimit*R*r*2*S; # mean count of a poisson process at the limit of acceptable contamination
violTrue = PTrue*R*r*2*S; # mean count of poisson process at the true contamination --> this one is not right


# % what is the value of the limit contamination's CDF that would exceed the
# % acceptable confidence? I.e. when does the limit CDF = 10%? this is the
# % max acceptable.
k = stats.poisson.ppf(1-confLimit, violLimit);

# % now what is the value of the true CDF just below that value? that's the
# % proportion of simulations that should pass
propPass = stats.poisson.cdf(k-1, violTrue); % k can be negative, in which case it returns zero