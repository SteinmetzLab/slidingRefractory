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