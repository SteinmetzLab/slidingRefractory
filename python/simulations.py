# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:40:23 2022

@author: Noam Roth

code to run simulations for slidingRefractory metric
"""
import numpy as np


def genST(rate, duration, xx = None,params = None):
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
    print('generating spike train...')
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
passPct = np.zeros([len(recDurScalarVec), len(rpvec),len(baseRates), len(contPct)])
passPctOld = np.zeros([len(recDurScalarVec), len(rpvec),len(baseRates), len(contPct)])
start_time = time. time()


#set up plotting
colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(params['baseRates'])))

for j, recDurScalar in enumerate(params['recDurs']):
    recDur = recDurScalar*3600
    print('recording Duration %d'%recDur)   
    for i, rp in enumerate(params['RPs']):
        print('refractory period duration %d'%rp)
        thresh = params['threshold']

        bidx=0
        for baseRate in baseRates:
            cidx=0
            for c in params['contRates']:
                
                contRate = baseRate*c
                rpidx=np.where(b<rp)[0][-1]+1
                mfunc =np.vectorize(genST)
                #### stopped above here. 

                m = mfunc(baseRate,b[1:-1],recDur,baseRate/10,thresh)
    
                
                simRes = np.zeros([nSim,len(b[1:-1])])
                simResOld = np.zeros([nSim])
                for n in range(nSim):
                    st = genST(baseRate,recDur)
                    isi = np.diff(np.insert(st,0,0)) 
                    isi = np.delete(isi,np.where(isi<rp)[0]) 
                    st = np.cumsum(isi)
                    if c>0:
                        contST = genST(contRate,recDur)
                    else:
                        contST=[]
                    combST = np.sort(np.concatenate((st, contST)))
                    c0 = correlograms(combST,np.zeros(len(combST),dtype='int8'),cluster_ids=[0],bin_size=binSize,sample_rate=20000,window_size=.05,symmetrize=False)
                    simRes[n,:] = np.cumsum(c0[0,0,1:(len(b)-1)])
                    len(simRes)
                    simResOld[n] = contamination_alt(np.asarray(combST))
                    #pass_vec_old.append(int(ce<0.1))
                passPct[j, i, bidx,cidx]=sum(np.any(np.less_equal(simRes[:,0:],m),axis=1))/nSim*100
                passPctOld[j, i, bidx,cidx]= sum([ce<0.1 for ce in simResOld])/nSim*100
                cidx+=1

            bidx+=1

current_time = time. time()
elapsed_time = current_time - start_time
print('Loop with %f iterations takes %f seconds' % (nSim, elapsed_time))       
         

with open('percentPass%0.0f.pickle'%round(nSim),'wb') as f:
    pickle.dump([passPct, passPctOld],f)





    
def plotSimulations(params):
    
    fig,axs = plt.subplots(1,params['nPlots'])
    

#%%
#script testing (to be moved later for import purposes)

# from simulations import genST  # this is not working, for some reason it's importing a different simulations gensT
meanRate = 10 #spks/s
recDur = 3600 #seconds
params = {}
params['checkFR']  = 1
st = genST(meanRate, recDur, xx = None, params = params)

#%%
sampleRate = 30000
params = {
    'recDurs': np.array([0.25, 0.5, 1, 2, 5]),  #recording durations (hours)
    'RPs': np.array([0.001,0.0015, 0.002, 0.003, 0.004, 0.005]),
    'baseRates': np.array([ 0.5, 1, 2, 5, 10, 20]),
    'contRates':  np.arange(0.00,0.2,0.02),
    'nSim': 10,
    'threshold': 0.1
    'binSize': 1 / sampleRate
}


