# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:40:23 2022

@author: Noam Roth

code to run simulations for slidingRefractory metric
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import slidingRP
import time
# import slidingRP
# from slidingRP import slidingRP



def genST(rate, duration, params = None):
    '''
    

    Parameters
    ----------
    rate : float
        firing rate of simulated neuron (spks/s).
    duration: float
        length of recording of simulated neuron (s)
    params : dict, optional
        Todo: here include any parameters like bursty, drifty? The default is None.

    Returns
    -------
    st: np.array
        array of spike times in seconds. 

    '''
    # print('generating spike train...')
    mu = 1 / rate
    n = rate * duration
    isi = np.random.exponential(mu, int(np.ceil(n * 2))) #generate 
    while sum(isi) < duration:
        np.append(isi, np.random.exponential(mu))

    st = np.cumsum(isi)
    if len(np.where(st<duration)[0])>0:
        st = st[0:np.where(st < duration)[0][(-1)]]
    else:
        st =[]

    if params['checkFR']:
        print('plotting FR...')
        fig,ax = plt.subplots(1,1)
        histWin = 0.5 #seconds
        binnedFR = firing_rate(st,hist_win = histWin, fr_win =  10)
        ax.plot(np.arange(binnedFR.size) * histWin, binnedFR, 'k')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('FR (spks/s)')
        ax.set_title('FR for simulated neuron, mean %d spikes/second, duration %d seconds'%(rate, duration) )

    return st


def simulateContNeurons(params):

    
    #initialize matrices of percent neurons pass
    passPct = np.empty([len(params['recDurs']), len(params['RPs']),len(params['baseRates']), len(params['contRates'])])
    passPct[:] = np.nan
    passPct2MsNoSpikes = np.zeros([len(params['recDurs']), len(params['RPs']),len(params['baseRates']), len(params['contRates'])])
    passPct2MsNoSpikes[:] = np.nan
    
    # passPctOld = np.zeros([len(recDurScalarVec), len(rpvec),len(baseRates), len(contPct)])
    start_time = time. time()
    
    
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
                             firingRate, secondsElapsed] = slidingRP.slidingRP(combST, params)
                        
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
    return passPct, passPct2MsNoSpikes
         
    
def plotSimulations(pc,params, savefile):
 
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(params['baseRates'])))
    fig,axs = plt.subplots(len(params['recDurs']),len(params['RPs']), figsize = (12,3*len(params['recDurs'])))

    for j, recDur in enumerate(params['recDurs']):
        for i, rp in enumerate(params['RPs']):
            
            if len(params['recDurs']) > 1 and len(params['RPs'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['contRates'], pc[j, i, b,:], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('Proportion contamination')
                ax.set_title('True RP %d ms'%(rp*1000))
        fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
   
    handles, labels = ax.get_legend_handles_labels()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')


    fig.savefig(savefile, dpi = 500)
    
    
    fig,axs = plt.subplots(len(params['contRates'][::2]),len(params['RPs']), figsize = (12*2,3*len(params['recDurs'])))
    for j, contRate in enumerate(params['contRates'][::2]):
        print(j)
        for i, rp in enumerate(params['RPs']):
            
            if len(params['recDurs']) > 1 and len(params['RPs'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['recDurs'], pc[:, i, b,j*2], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('recording Duration')
                ax.set_title('True RP %d ms'%(rp*1000))
        fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
    fig.subplots_adjust(left=0.5, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')
    fig.show()
    
    
    
    fig,axs = plt.subplots(len(params['recDurs']), len(params['contRates'][::2]), figsize = (12*2,3*len(params['recDurs'])))
    for j, recDur in enumerate(params['recDurs']):
        print(j)
        for i, contRate in enumerate(params['contRates'][::2]):
            
            if len(params['recDurs']) > 1 and len(params['contRates'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['RPs'], pc[j,:, b,i*2], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('True RP')
                ax.set_title('contamination %.2f '%contRate)
        fig.text(0.65, 0.9-(.17*j), 'Recording Duration: %d hours'%recDur)
    fig.subplots_adjust(left=0.5, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')
    fig.show()