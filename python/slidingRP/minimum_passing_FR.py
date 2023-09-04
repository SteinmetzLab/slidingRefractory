#plot minimium firing rate curve by just running computeviol for different firing rates and recDur



from slidingRP.simulations import *
from slidingRP.metrics import computeViol

contaminationProp = 0.1
confThresh = 0.9
obsViol = 0 #looking at uncontaminated neurons only

refDur = 0.002

firingRates = np.arange(0.1,5,0.01)
recDursHours = np.arange(0.5,5.5,0.25)
recDurs = [i*3600 for i in recDursHours]

confMat = np.empty([len(firingRates), len(recDurs)])
for f, firingRate in enumerate(firingRates):
    for r, recDur in enumerate(recDurs):
        spikeCount = firingRate * recDur

        conf = computeViol(obsViol, firingRate, spikeCount, refDur, contaminationProp)

        confMat[f,r] = conf


print(confMat)
confMat>.9
minFR = []
for i in range(len(recDurs)):
    res = next(x for x, val in enumerate(confMat[:,i]) if val > confThresh)
    print(firingRates[res])
    minFR.append(firingRates[res])

fig, ax = plt.subplots(1, 1, figsize=(6, 4))

ax.set_xlabel('Recording Duration (hours)')
ax.set_ylabel('Lowest passing firing rate (spk/s)')
ax.plot(recDursHours,minFR,'k.-')
ax.set_title('Minimum passing FR (uncontaminated)')

spinesSetting = False
ax.spines.right.set_visible(spinesSetting)
ax.spines.top.set_visible(spinesSetting)
fig.show()



#%%

def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


#plot minFR as a function of recDur
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\contFRContourDict.pickle'
file = open(savefile, 'rb')
contContourDict = pickle.load(file)
file.close()

#need to also load in the corresponding params file for the correct recording durations
savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_122.pickle'
file = open(savefile, 'rb')
results = pickle.load(file)
file.close()

params = results[-1] #the last in the results
pc = results[0] #normal metric results (not 2ms cond)

frThreshRecDur = []
recDurVec = params['recDurs']
for recordingDuration in recDurVec:
    recDurClosest = closest(params['recDurs'], recordingDuration)  # closest rec dur for which I've computed contours
    recDurInd = np.where(params['recDurs'] == recDurClosest)
    # The rule is: 20% contamination we want at least 50% reject, for 50% contamination we want at least 90% reject
    lowContFR = contContourDict[(0.2, 50)][recDurInd]
    highContFR = contContourDict[(0.5, 90)][recDurInd]
    print(lowContFR)
    frThreshRecDur.append(max(lowContFR, highContFR)[0])


#%%
fig,ax = plt.subplots()
ax.plot(recDurVec,frThreshRecDur,'.-',color = 'k')
ax.set_xlabel('Recording duration (hrs)')
ax.set_ylabel('Minimum passing FR (spks/s)')

ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
# ax.set_yticks()
fig.show()
