# A script to run and save figures to file (pdf) for sliding refractory paper

#establish paths:
basepath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\paper_figs_Nick'

#imports:
import os



# Figure 1: #TODO
#a left:  a schematic
#a right: example neurons #TODO
# b:

# c:
#%% Figure 2: Metric algorithm
from slidingRP.plotFig2_ACG import * #import code to make this figure
figsavepath = basepath + '\\figure2' #destination folder for all subfigures
os.makedirs(figsavepath,exist_ok=True) #if this folder doesn't yet exist, create it
runSaveFig2(figsavepath) #run code and save figure to figsavepath

#%% Figure 3: Simulations
from slidingRP.Fig3_simulations import *
figsavepath = basepath + '\\figure3'
os.makedirs(figsavepath,exist_ok=True)
runSaveFig3(figsavepath)

#%% Figure 4:Compare with Llobet
from slidingRP.Fig4_HillCompConfidence import *
figsavepath = basepath + '\\figure4'
os.makedirs(figsavepath,exist_ok=True)
print('Running Figure 4')
runSaveFig4(figsavepath)
print('Figure 4 saved')
#%% Figure 5: Confidence and ROC
from slidingRP.Fig5_ConfROC import *
figsavepath = basepath + '\\figure5'
os.makedirs(figsavepath,exist_ok=True)
# runSaveFig5_ab(figsavepath)
runSaveFig5_ROC(figsavepath)
#%% Figure 6:
# from slidingRP.Fig


# Figure 7:

