# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 13:31:20 2022

@author: Noam Roth
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import time
from simulations import genST
import scipy
import multiprocessing   
from slidingRP import *
 
# def func(n): 
#     for i in range(1000): 
#         for j in range(1000): 
#             s=j*i 
#     print(n) 
 
# if __name__ == '__main__': 
#     pool = multiprocessing.Pool(processes=4) 
#     pool.map(func, range(10)) 
#     pool.close() 
#     pool.join()    
#     print('done') 
    
    
def simulateContNeurons(n):

    sampleRate = 30000
    params = {
    'recDurs': np.array([0.25, 0.5]),  #recording durations (hours)
    'RPs': np.array([0.001, 0.005]), #true RP (s)
    'baseRates': np.array([ 0.5, 20 ]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates':  np.arange(0.00,0.225, 0.1), #contamination levels (proportion) #.025
    'nSim': 5,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True ,
    'savePCfile': True
    
    }
    #initialize matrices of percent neurons pass
    passPct = np.empty([len(params['recDurs']), len(params['RPs']),len(params['baseRates']), len(params['contRates'])])
    passPct[:] = np.nan
    passPct2MsNoSpikes = np.zeros([len(params['recDurs']), len(params['RPs']),len(params['baseRates']), len(params['contRates'])])
    passPct2MsNoSpikes[:] = np.nan
    
    # passPctOld = np.zeros([len(recDurScalarVec), len(rpvec),len(baseRates), len(contPct)])
    start_time = time.time()
    
    
    #set up plotting
    
    for j, recDurScalar in enumerate(params['recDurs']):
        recDur = recDurScalar*3600
        print('recording Duration %d'%recDur)   
        for i, rp in enumerate(params['RPs']):
            print('refractory period duration %.3f'%rp)
            thresh = params['threshold']
    
            bidx=0
            for baseRate in params['baseRates']:
                print('baseRate %.2f'%baseRate)
                cidx=0
                for c in params['contRates']:
                    contRate = baseRate*c
                    print('contRate %.2f'%c)

                    # rpidx=np.where(b<rp)[0][-1]+1 #???????
                    # mfunc =np.vectorize(genST)
                    #### stopped above here. 
    
    #in the old version, I used max acceptable: computed a max acceptable for this set of parameters
    # and then (in the loop below nSim) simulated a bunch of neurons and compared each to this same max acceptable
    
                    # m = mfunc(baseRate,b[1:-1],recDur,baseRate/10,thresh)
        
                    
                    # simRes = np.zeros([nSim,len(b[1:-1])])
                    # simResOld = np.zeros([nSim])
                    passVec = np.empty(params['nSim'])
                    passVec[:] = np.nan
                    passVec2MsNoSpikes = np.empty(params['nSim'])
                    passVec2MsNoSpikes[:] = np.nan
                    for n in range(params['nSim']):
                        st = genST(baseRate,recDur,params) #generate a spike train with the current base rate
                        isi = np.diff(np.insert(st,0,0)) 
                        # print(rp)
                        isi = np.delete(isi,np.where(isi<rp)[0]) #get rid of already contaminating spikes (does this make sense??? why are there already contaminating spikes)
                        st = np.cumsum(isi)
                        if c>0:
                            contST = genST(contRate,recDur,params)
                        else:
                            contST=[]
                        combST = np.sort(np.concatenate((st, contST))) # put spike times from both spike trains together (and sort chronologically)
                        
                        [maxConfidenceAt10Cont, minContWith90Confidence, timeOfLowestCont,
                             nSpikesBelow2, confMatrix, cont, rpVec, nACG,
                             firingRate, secondsElapsed] = slidingRP(combST, params)
                        
                        if minContWith90Confidence <=10:
                            passVec[n] = 1
                        else:
                            passVec[n] = 0
                        
                        passVec2MsNoSpikes[n] = passVec[n]
                        if nSpikesBelow2 ==0:
                            passVec2MsNoSpikes[n] = 1
                        # c0 = correlograms(combST,np.zeros(len(combST),dtype='int8'),cluster_ids=[0],bin_size=binSize,sample_rate=20000,window_size=.05,symmetrize=False)
                        # simRes[n,:] = np.cumsum(c0[0,0,1:(len(b)-1)])
                        # len(simRes)
                        # simResOld[n] = contamination_alt(np.asarray(combST))
                        # #pass_vec_old.append(int(ce<0.1))
                    passPct[j, i, bidx,cidx]=sum(passVec)/params['nSim']*100
                    passPct2MsNoSpikes[j, i, bidx,cidx]=sum(passVec2MsNoSpikes)/params['nSim']*100

                    # passPctOld[j, i, bidx,cidx]= sum([ce<0.1 for ce in simResOld])/nSim*100
                    cidx+=1
    
                bidx+=1
    
    current_time = time. time()
    elapsed_time = current_time - start_time
    print('Loop with %f iterations takes %f seconds' % (params['nSim'], elapsed_time))    
    print(n)
    return passPct, passPct2MsNoSpikes


if __name__ == '__main__': 
    pool = multiprocessing.Pool(processes=4) 
    pool.map(simulateContNeurons, range(3)) 
    pool.close() 
    pool.join()    
    print('done') 
    