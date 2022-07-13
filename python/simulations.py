# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:40:23 2022

@author: Noam Roth

code to run simulations for slidingRefractory metric
"""
import numpy as np


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
    return st

    if params['checkFR']:
        fig,ax = plt.subplots(1,1)
        histWin = 0.5 #seconds
        binnedFR = firing_rate(st,hist_win = histWin, fr_win =  10)
        ax.plot(np.arange(binnedFR.size) * histWin, binnedFR, 'k')
        ax.set_xlabel('Time (s)'))
        ax.set_ylabel('FR (spks/s)')

    return st