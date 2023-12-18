import pickle

from slidingRP.simulations import *

def runSaveFig3(figsavepath):

    #path to load results
    resultsBasePath = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC'

    #load simulation results
    nIter = 500 #number of simulation iterations run
    date = '_06_29' #date simulations were run
    version = '1'
    resultsPath = resultsBasePath + str(nIter) + 'iter' + date + version + '.pickle'

    file = open(resultsPath, 'rb')
    results = pickle.load(file)
    file.close()
    pc = results[0]
    params = results[-1]


    #make RP plots for regular and Hill
    nSim = 500

    figsavefile = figsavepath + '\\fig3'

    plotSimulations(pc, params,figsavefile, subplot1=True)
    plotSimulations(pc, params,figsavefile, subplot2=True)
    plotSimulations(pc, params,figsavefile, subplot3=True)
    plotSimulations(pc, params,figsavefile, subplot4=True)
