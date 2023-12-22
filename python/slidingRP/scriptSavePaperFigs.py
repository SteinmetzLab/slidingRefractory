# A script to run and save figures to file (pdf) for sliding refractory paper

#establish paths:
figSaveBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs_Nick'
resultsBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationResults'

import os #import to create folders for each figure


#%%
#Set rerun flag:
rerunAllSims = False

if rerunAllSims = True:
    rerunFig3 = True
    rerunFig4 = True
    rerunFig5 = True
    rerunFig6 = True
    rerunFig7 = True



#%% Figure 1: #TODO

#%% Figure 2: Metric algorithm
from slidingRP.Fig2_ACGMetric import * #import code to make this figure
figsavepath = figSaveBasePath + '\\figure2' #destination folder for all subfigures
os.makedirs(figsavepath,exist_ok=True) #if this folder doesn't yet exist, create it
savedSimNeuronFlag=True #set to true to load a previously saved simulated contaminated neuron
runSaveFig2(figsavepath,resultsBasePath,savedSimNeuronFlag=savedSimNeuronFlag) #run code and save figure to figsavepath



#%% Figure 3: Simulations
from slidingRP.Fig3_simulations import *
figsavepath = figSaveBasePath + '\\figure3'
os.makedirs(figsavepath,exist_ok=True)


rerunFig3 = False
runSaveFig3(figsavepath,resultsBasePath, rerunFig3 = rerunFig3)



#%% Figure 4:Compare with Llobet
from slidingRP.Fig4_HillCompConfidence import *
figsavepath = figSaveBasePath + '\\figure4'
os.makedirs(figsavepath,exist_ok=True)

rerunFig4 = False
print('Running Figure 4')
runSaveFig4(figsavepath,resultsBasePath, rerunFig4 = rerunFig4)
print('Figure 4 saved')



#%% Figure 5: Confidence and ROC
from slidingRP.Fig5_ConfROC import *
figsavepath = figSaveBasePath + '\\figure5'
os.makedirs(figsavepath,exist_ok=True)

rerunFig5 = False
# runSaveFig5_ab(figsavepath)
runSaveFig5_ROC(figsavepath,resultsBasePath, rerunFig5 = rerunFig5)



#%% Figure 6:
# from slidingRP.Fig6_drift import *
# figsavepath = figSaveBasePath + '\\figure6'
# os.makedirs(figsavepath,exist_ok=True)

rerunFig6 = False
# runSaveFig6(figsavepath, resultsBasePath, rerunFig6 = rerunFig6)

# Figure 7:

rerunFig7 = False

