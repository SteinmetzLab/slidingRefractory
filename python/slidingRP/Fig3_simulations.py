import pickle

from slidingRP.simulations import *

def runSaveFig3(figSavePath, resultsBasePath):

    #path to load results

    #load simulation results
    nIter = 500 #number of simulation iterations run
    date = '_06_29' #date simulations were run
    version = '1'
    resultsPath = resultsBasePath + '\\simulationsPC' + str(nIter) + 'iter' + date + version + '.pickle'

    file = open(resultsPath, 'rb')
    results = pickle.load(file)
    file.close()
    pc = results[0]
    params = results[-1]


    #make RP plots for regular and Hill
    nSim = 500

    figsavefile = figSavePath + '\\fig3'

    plotSimulations(pc, params,figsavefile, subplot1=True)
    plotSimulations(pc, params,figsavefile, subplot2=True)
    plotSimulations(pc, params,figsavefile, subplot3=True)
    plotSimulations(pc, params,figsavefile, subplot4=True)
