# A script to run and save figures to file (pdf) for sliding refractory paper

#establish paths:
figSaveBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs_Nick'
resultsBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationResults'

import os #import to create folders for each figure

#%% Figure 1: #TODO

#%% Figure 2: Metric algorithm
from slidingRP.Fig2_ACGMetric import * #import code to make this figure
figsavepath = figSaveBasePath + '\\figure2' #destination folder for all subfigures
os.makedirs(figsavepath,exist_ok=True) #if this folder doesn't yet exist, create it
runSaveFig2(figsavepath) #run code and save figure to figsavepath

#%% Figure 3: Simulations
from slidingRP.Fig3_simulations import *
figsavepath = figSaveBasePath + '\\figure3'
os.makedirs(figsavepath,exist_ok=True)
runSaveFig3(figsavepath,resultsBasePath)

#%% Figure 4:Compare with Llobet
from slidingRP.Fig4_HillCompConfidence import *
figsavepath = figSaveBasePath + '\\figure4'
os.makedirs(figsavepath,exist_ok=True)
print('Running Figure 4')
runSaveFig4(figsavepath,resultsBasePath)
print('Figure 4 saved')
#%% Figure 5: Confidence and ROC
from slidingRP.Fig5_ConfROC import *
figsavepath = figSaveBasePath + '\\figure5'
os.makedirs(figsavepath,exist_ok=True)
# runSaveFig5_ab(figsavepath)
runSaveFig5_ROC(figsavepath,resultsBasePath)
#%% Figure 6:
# from slidingRP.Fig6_drift import *
# figsavepath = figSaveBasePath + '\\figure6'
# os.makedirs(figsavepath,exist_ok=True)
# runSaveFig6(figsavepath)

# Figure 7:



