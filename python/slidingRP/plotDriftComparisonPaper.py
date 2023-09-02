import pickle
import matplotlib.pyplot as plt
import numpy as np
from phylib.stats import correlograms
from scipy import stats
from slidingRP.simulations import *


def plotDriftComp(pcReg,pcDrift, savefile, recDurPlot = 2, rpPlot = 2,frPlot=2,highCont=False):

    #preferences:
    titleflag = 0
    spinesSetting = False

    numrows = 1
    numcols = 1
    numplots = numrows * numcols
    fig, axs = plt.subplots(numrows, numcols, figsize=(4, 5))
    ax = axs

    if titleflag:
        titlexval = 0.5
        titleyvals = [1, 0.77, 0.52, 0.27]

    # plot just recDur == recDurPlot:
    recDurs = params['recDurs']
    recDursInd = np.where(recDurs == recDurPlot)[0]

    # plot just RP = rpPlot:
    rps = params['RPs']
    rpInd = np.where(rps == rpPlot/1000)[0]

    #plot this set of firing rates (baseRates)
    frs = np.array(params['baseRates'])
    print(frs)
    frInd = np.where(frs == frPlot)[0]# [x for x in range(len(frs)) if frs[x] in frPlot]
    print(frInd)

    frIndDrift = np.where(frs+(0.5*frs*delta) == frPlot)[0]
    print('drift Ind', frIndDrift)



    cReg = cc.linear_protanopic_deuteranopic_kbw_5_95_c34
    cReg = cReg[::-1]  # if using linear_blue37 or protanopic, flip the order
    colorsReg = [cReg[x] for x in np.round(np.linspace(0.2, 0.75, len(frs)) * 255).astype(int)]

    cDrift = cc.linear_green_5_95_c69
    cDrift = cDrift[::-1]  # if using linear_blue37 or protanopic, flip the order
    colorsDrift = [cDrift[x] for x in np.round(np.linspace(0.2, 0.75, len(frs)) * 255).astype(int)]

    #error bars
    count = pcReg / 100 * params['nSim']  # number of correct trials
    CI_scaled = binofit(count, params['nSim'])
    CIReg = [x * 100 for x in CI_scaled]

    count = pcHill/ 100 * params['nSim']  # number of correct trials
    CI_scaled = binofit(count, params['nSim'])
    CIHill = [x * 100 for x in CI_scaled]

    print(recDurs[recDursInd])
    print(rps[rpInd])
    print(frs[frInd])
    for j, recDur in enumerate(recDurs[recDursInd]):
        for i, rp in enumerate(rps[rpInd]):
            # different base rates get different colors
            print('here1')
            for b, baseRate in enumerate(frs[frInd]):
                print('here2')
                fr_use = frs[frInd[b]]

                #regular plot
                lowerCI = CIReg[0][recDursInd, rpInd, frInd, :][0]
                upperCI = CIReg[1][recDursInd, rpInd, frInd, :][0]
                x = params['contRates'] * 100
                y = pcReg[recDursInd, rpInd, frInd, :][0]
                print(y)

                ax.plot(x, y, '.-', color=colorsReg[frInd[b]], label='Our metric')
                ax.fill_between(x, lowerCI, upperCI, color=colorsReg[frInd[b]], alpha=.5)

                #Hill plot
                lowerCI = CIHill[0][recDursInd, rpInd, frInd, :][0]
                upperCI = CIHill[1][recDursInd, rpInd, frInd, :][0]
                x = params['contRates'] * 100
                y = pcHill[recDursInd, rpInd, frInd, :][0]
                ax.plot(x, y, '.-', color=colorsHill[frInd[b]], label='Hill metric')
                ax.fill_between(x, lowerCI, upperCI, color=colorsHill[frInd[b]], alpha=.5)


                #plot points at 50% cont
                if highCont:
                    ax.plot(25,y[-1])
                    ax.vlines(x=25, ymin=lowerCI[-1], ymax=upperCI[-1],
                              colors=colors[b])



                ax.set_ylabel('Percent pass')
                ax.set_xlabel('Contamination (%)')

                ax.set_ylim(-10, 110)
                ax.spines.right.set_visible(spinesSetting)
                ax.spines.top.set_visible(spinesSetting)
                ax.xaxis.set_ticks([0, 10, 20])

                #vertical line at 10 cont
                ax.vlines(10, -10, 110, color='k')
                #light gray horizontal lines
                xmin, xmax = ax.get_xlim()
                # ax.hlines([0,20,40,60,80,100], xmin, xmax,color = '0.8',alpha = 0.5)
                ax.legend(frameon=False)

        if titleflag:
            if recDur != 1:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hours' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
                else:
                    plt.figtext(0.5, titleyvals[j], 'Recording duration: %.1f hour' % recDur, ha="center",
                                va="top", fontsize=14, color="k")
            ax.set_title('True RP: %.1f ms' % (rpPlot))

    handles, labels = ax.get_legend_handles_labels()

    plt.subplots_adjust(hspace=0.6, wspace=0.5)
    fig.tight_layout()
    # fig.legend(handles, labels, title='Firing rate (spk/s)', loc='upper right', bbox_to_anchor=(1, .9))
    # fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, .9))

    fig.savefig(savefile + 'FR' + str(frPlot) + '_Main.svg', dpi=500)
    fig.savefig(savefile + 'FR' + str(frPlot) + '_Main.png', dpi=500)


#%%load data

#first load original data
loadfile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'
file = open(loadfile,'rb')
results = pickle.load(file)
file.close()

pc = results[0]
pc2MsNoSpikes = results[1]
pcHill2 = results[3]
pcHill3 = results[4]
params = results[5]






import datetime
#prep for figure saving
date_now  = datetime.datetime.now().strftime('_%m_%d')

version = '2'  # adjust if running more than once in the same day

for delta in [0.5,-0.5,0.1,-0.1]:

    if delta < 0:
        deltaprint = 'neg' + str(int(delta*10))
    else:
        deltaprint = str(int(delta*10))
    print(deltaprint)
    savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPC500iter_01_242delta' + deltaprint + '.pickle'
    file = open(savefile,'rb')
    results = pickle.load(file)
    file.close()
    pcDrift = results[0]
    paramsDrift = results[1]
    #make RP plots for regular and Hill
    figsavefile1 = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\resultsIncreaseDecrease\simulationsPCdelta' + deltaprint + str(params['nSim']) + 'iter' + date_now
    for fr in params['baseRates']:
        plotDriftComp(pc, pcDrift, savefile, frPlot=fr)







loadfile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC500iter_01_221.pickle'
file = open(loadfile,'rb')
results = pickle.load(file)
file.close(



pc = results[0]
pc2MsNoSpikes = results[1]
pcHill2 = results[3]
pcHill3 = results[4]
params = results[5]


#%%

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\HillComp'
for fr in params['baseRates']:
    plotHillComp(pc,pcHill2,savefile,frPlot = fr)
#%%

savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\HillComp2MsNoSpikes'
for fr in params['baseRates']:
    plotHillComp(pc2MsNoSpikes,pcHill2,savefile,frPlot = fr)


