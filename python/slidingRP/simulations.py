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
from slidingRP.metrics import slidingRP 
import time

from statsmodels.stats.proportion import proportion_confint as binofit
    #(I had to pip install statsmodels; TODO check for a better package? )



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

                    passVec = np.empty(params['nSim'])
                    passVec[:] = np.nan
                    passVec2MsNoSpikes = np.empty(params['nSim'])
                    passVec2MsNoSpikes[:] = np.nan
                    for n in range(params['nSim']):
                        st = genST(baseRate,recDur,params) #generate a spike train with the current base rate
                        isi = np.diff(np.insert(st,0,0)) 
                        print(n,end="")
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

                    passPct[j, i, bidx,cidx]=sum(passVec)/params['nSim']*100
                    passPct2MsNoSpikes[j, i, bidx,cidx]=sum(passVec2MsNoSpikes)/params['nSim']*100

                    cidx+=1
    
                bidx+=1
    
    current_time = time. time()
    elapsed_time = current_time - start_time
    print('Loop with %f iterations takes %f seconds' % (params['nSim'], elapsed_time))       
    return passPct, passPct2MsNoSpikes
         
    
def plotSimulations(pc,params, savefile, Fig1 = False, Fig2 = False, Fig3 = False, Fig4 = True):
 
    #compute confidence intervals
    count = pc / 100 * params['nSim'] #number of correct trials 
    CI_scaled = binofit(count,params['nSim'])
    CI = [x*100 for x in CI_scaled]
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(params['baseRates'])))

    if Fig1:
        fig,axs = plt.subplots(len(params['recDurs']),len(params['RPs']), figsize = (12,3*len(params['recDurs'])))

        for j, recDur in enumerate(params['recDurs']):
            for i, rp in enumerate(params['RPs']):

                if len(params['recDurs']) > 1 and len(params['RPs'])>1:
                    ax = axs[j,i]
                else:
                    ax = axs[(j+1)*i]
                #different base rates get different colors
                for b, baseRate in enumerate(params['baseRates']):
                    lowerCI = CI[0][j, i, b,:]
                    upperCI = CI[1][j, i, b,:]
                    x = params['contRates']
                    y =  pc[j, i, b,:]
                    ax.plot(x, y, '.-',color = colors[b], label = baseRate)

                    ax.fill_between(x, lowerCI,upperCI, color=colors[b], alpha=.3)


                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('Prop. cont.')
                    ax.set_title('True RP %d ms'%(rp*1000))
            # fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
            fig.suptitle('Recording duration: %d hour'%recDur, x=.5, y=1.1)

        handles, labels = ax.get_legend_handles_labels()
        fig.tight_layout(pad=1, w_pad=1.1, h_pad=1.0)
        # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

        fig.legend(handles, labels, loc='upper right',bbox_to_anchor=(1.1, 1))


        fig.savefig(savefile + '_Main.svg', dpi = 500)
    
    if Fig2:
        fig,axs = plt.subplots(len(params['contRates'][::2]),len(params['RPs']), figsize = (12*2,3*len(params['recDurs'])))
        for j, contRate in enumerate(params['contRates'][::2]):
            for i, rp in enumerate(params['RPs']):

                if len(params['contRates'][::2]) > 1 and len(params['RPs'])>1:
                    ax = axs[j,i]
                else:
                    ax = axs[(j+1)*i]

                #different base rates get different colors
                for b, baseRate in enumerate(params['baseRates']):

                    lowerCI = CI[0][:, i, b,j*2]
                    upperCI = CI[1][:, i, b,j*2]
                    x = params['recDurs']
                    y =  pc[:, i, b,j*2]
                    ax.plot(x, y, '.-',color = colors[b], label = baseRate)

                    ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)

                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('recording Duration')
                    ax.set_title('True RP %d ms'%(rp*1000))
            # fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
            fig.suptitle('Proportion contamination: %.2f'%contRate, x=.5, y=1.1)
        handles, labels = ax.get_legend_handles_labels()

        fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

        fig.legend(handles, labels, loc='upper right')
        fig.savefig(savefile + '_recDur.svg', dpi = 500)

    
    if Fig3:
        fig,axs = plt.subplots(len(params['recDurs']), len(params['contRates'][::2]), figsize = (12*2,3*len(params['recDurs'])))
        for j, recDur in enumerate(params['recDurs']):
            for i, contRate in enumerate(params['contRates'][::2]):

                if len(params['recDurs']) > 1 and len(params['contRates'])>1:
                    ax = axs[j,i]
                else:
                    ax = axs[(j+1)*i]
                #different base rates get different colors
                for b, baseRate in enumerate(params['baseRates']):
                    lowerCI = CI[0][j,:, b,i*2]
                    upperCI = CI[1][j,:, b,i*2]
                    x = params['RPs']
                    y =  pc[j,:, b,i*2]
                    ax.plot(x, y, '.-',color = colors[b], label = baseRate)
                    ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('True RP')
                    ax.set_title('contamination %.2f '%contRate)
            # fig.text(0.65, 0.9-(.17*j), 'Recording Duration: %d hours'%recDur)
            fig.suptitle('Recording Duration: %d hour'%recDur, x=.5, y=1.1)
        fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
        handles, labels = ax.get_legend_handles_labels()

        fig.legend(handles, labels, loc='upper right')
        print('hi')
        fig.savefig(savefile + '_RP.svg', dpi = 500)

    if Fig4:
        fig, axs = plt.subplots(1, 1, figsize=(12 * 2, 3))
        ax = axs #for the case of just one subplot
        #plot just contRates 0.08 to 0.12:
        cr = params['contRates'];
        crInds = np.where((cr >= 0.08) & (cr <= 0.12))[0]
        #plot just recDur = 1
        rd = params['recDurs']
        rdInd = np.where(rd==1)[0]

        for j, recDur in enumerate(rd[rdInd]):
            for i, contRate in enumerate(cr[crInd]):

                # different base rates get different colors
                for b, baseRate in enumerate(params['baseRates']):
                    lowerCI = CI[0][j, :, b, i * 2]
                    upperCI = CI[1][j, :, b, i * 2]
                    x = params['RPs']
                    y = pc[j, :, b, i * 2]
                    ax.plot(x, y, '.-', color=colors[b], label=baseRate)
                    ax.fill_between(x, lowerCI, upperCI, color=colors[b], alpha=.3)
                    ax.set_ylabel('Percent pass')
                    ax.set_xlabel('True RP')
                    ax.set_title('contamination %.2f ' % contRate)
            # fig.text(0.65, 0.9-(.17*j), 'Recording Duration: %d hours'%recDur)
        fig.subplots_adjust(left=0.7, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)
        handles, labels = ax.get_legend_handles_labels()

        fig.legend(handles, labels, loc='upper right')
        print('hi1')
        fig.savefig(savefile + '_individual.svg', dpi=500)