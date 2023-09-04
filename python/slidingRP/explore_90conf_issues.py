from slidingRP.simulations import genST
from slidingRP.metrics import computeMatrix
import numpy as np
import pickle
import matplotlib.pyplot as plt

recDur = 2 * 3600
rp = 0.002
contLevels = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

baseRates = np.arange(1,35,1)
b = [] #fit error

fig, axs = plt.subplots(2, 1, figsize=(7, 10))
ax1 = axs[0]
ax2 = axs[1]
for contLevel in contLevels :
    obs = []
    exp = []
    for baseRate in baseRates:

        contRate = baseRate * contLevel

        #load some params
        savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_122.pickle'
        file = open(savefile, 'rb')
        results = pickle.load(file)
        file.close()
        params = results[-1] #the last element saved in the results

        st = genST(baseRate, recDur, rp,params)
        contST = genST(contRate, recDur, 0,params)
        combST = np.sort(np.concatenate((st, contST)))

        [confMatrix, cont, rpTestVals, nACG, firingRate] = computeMatrix(combST, params)

        expectedViol = contRate * rp * 2 * (baseRate + contRate) * recDur
        #         expectedViol = contRate * (1-contLevel/3) * rp * 2 * (baseRate + contRate) * recDur
        obs.append(sum(nACG[np.where(rpTestVals <=rp)[0]]))
        exp.append(expectedViol)


    ratioOE = np.divide(obs,exp)
    ax1.plot(baseRates,ratioOE,'k.')
    ax1.hlines(1,0,baseRates[-1])
    ax1.hlines(0.95,0,baseRates[-1],'k',linestyle='dashed')
    ax1.set_ylabel('ratio of obs FR to exp')
    ax1.set_xlabel('baseRate')
    # ax1.set_title('recording Duration = {} hours   rp = {} s   contLevel = {}'.format((str(int(recDur/3600))), str(rp),str(contLevel)))
    ax1.set_title('recording Duration = {} hours   rp = {} s ; all contLevels'.format((str(int(recDur/3600))), str(rp)))

    ax2.plot(baseRates, obs, 'k.', label = 'observed')
    ax2.plot(baseRates,exp,'r.', label = 'expected')
    ax2.set_ylabel('Num violations')
    ax2.set_xlabel('baseRate')


handles, labels = ax2.get_legend_handles_labels()
fig.legend(handles[0:2], labels[0:2], title='Estimated error', loc='upper right', bbox_to_anchor=(.4, .4))
fig.show()





    #find intercept
    b.append(1-np.mean(ratioOE))

#%%
fig, ax = plt.subplots()
byEyeError = [1- i for i in [0.95,0.91,0.88,0.85,0.83,0.81,0.79,0.77,0.75,0.74]]

ax.plot(contLevels,byEyeError,'k.-', label = 'by eye')
ax.plot(contLevels,b,'r.-', label = 'mean value')
ax.set_xlabel('contamination level')
ax.set_ylabel('estimated error')
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title='Estimated error', loc='upper left', bbox_to_anchor=(.2, .8))
fig.show()


#%%
fig, ax = plt.subplots()
byEyeError = [1- i for i in [0.95,0.91,0.88,0.85,0.83,0.81,0.79,0.77,0.75,0.74]]
ax.plot(contLevels,byEyeError,'k.-')
ax.set_xlabel('contamination level')
ax.set_ylabel('by eye estimated error')
fig.show()


#%%
fig, ax = plt.subplots()
ax.plot(baseRates,obs,'r.-',label = 'observed')
ax.plot(baseRates,exp,'k.-',label = 'expected')
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title='Firing rate (spk/s)', loc='upper left', bbox_to_anchor=(.2, .8))
fig.show()


#%% spike collisions?
len(np.where(np.diff(combST) < (params['binSize']))


#%%

baseRates = np.arange(0.5,35,1)
genFR =[]
for baseRate in baseRates:
    st = genST(baseRate, recDur, rp, params)
    FR = len(st) / recDur
    genFR.append(FR)

genFR/baseRates
plt.plot(genFR/baseRates,'k.')
plt.hlines(1,0,baseRates[-1])
plt.ylabel('ratio of genST FR to baseRate input')
plt.xlabel('baseRate')
plt.show()

#%%
baseRates = np.arange(0.5,35,1)
genFR =[]
for baseRate in baseRates:
    st = genST(baseRate, recDur, rp, params)
    FR = len(st) / recDur
    genFR.append(FR)

plt.plot(genFR/baseRates,'k.')
plt.hlines(1,0,baseRates[-1])
plt.ylabel('ratio of genST FR to baseRate input')
plt.xlabel('baseRate')
plt.show()
