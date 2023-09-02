import numpy as np
import matplotlib.pyplot as plt
from brainbox.singlecell import firing_rate
from slidingRP.simulations import *


def plotFR(ax, st,color,histWin=0.5,fr_win =10):
    print('plotting FR...')
    histWin = 0.5  # seconds
    binnedFR = firing_rate(st, hist_win=histWin, fr_win=10)
    if ax == None:
        fig, ax = plt.subplots(1, 1)

    ax.plot(np.arange(binnedFR.size) * histWin, binnedFR, color=color)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('FR (spks/s)')
    # ax.set_title('FR for simulated neuron, mean %d spikes/second, duration %d seconds' % (rate, duration))
    spinesSetting = False
    ax.spines.right.set_visible(spinesSetting)
    ax.spines.top.set_visible(spinesSetting)

rate = 10; duration = 3600; params = None
delta = 0.5 #fractional increase
duration = 3600 #1 hour

params = {}
params['checkFR']=False
st = genST(rate+(0.25*rate),duration = duration, params=params)
driftST = genChangingST(rate,duration=duration,params=params,delta=delta)
fig, ax = plt.subplots(1, 1)
plotFR(ax,st,'b')
plotFR(ax,driftST,'g')

plt.show()


#%%
st = genST(rate,duration = duration, params=params)
dieST = genST(rate,duration=duration/2,params=params)
fig, ax = plt.subplots(1, 1)
plotFR(ax,st,'blue')
plotFR(ax,dieST,'blueviolet')

plt.show()

