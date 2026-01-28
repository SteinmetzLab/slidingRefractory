import pickle
from slidingRP.simulations import *
import datetime


def runSaveFig6(figSavePath, resultsBasePath, rerunFig6 = False):


    #parameters for different drift settings:
    drift_values = [-0.5, 0, 0.5]
    drift_strings = ['Dec', 'Same', 'Inc']

    if rerunFig6:
        #set up parameters to rerun simulations
        sampleRate = 30000
        params = {
            'recDurs': np.array([2]),  #recording durations (hours)
            'RPs': np.array([0.002]),#np.array([0.0015,0.002,0.003,0.004]),#np.array([0.001,0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]), #true RP (s)
            # 'baseRates': [0.25, 0.5,0.75, 1,1.5, 2,3, 4,5,6, 7.5,10,15],#np.arange(0.05, 1, 0.05) ,#   [0.05, np.arange(0.05, 1.4, 0.1)[:],2,4,5,10,20] #np.array([0.75,2,3,4,7.5]), #F1, 2, 5, 10 , 20 R (spk/s)
            'contRates': np.arange(0.00,0.21, 0.02),#np.array([.2, .5]),#%np.array([0.09,0.095,0.1,0.105,0.11]),#np.arange(0.00,0.21, 0.01), #contamination levels (proportion) #.025
            'nSim': 500,
            'contaminationThresh': 10,
            'binSize': 1 / sampleRate,
            'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
            'confidenceThresh': 90,
            'checkFR': False,
            'binSizeCorr': 1 / sampleRate,
            'returnMatrix': True,
            'verbose': True,
            'savePCfile': True
        }

        # get date for filename
        date_now = datetime.datetime.now().strftime('_%m_%d')

        #run simulation for each of the drift parameters:
        for drift, string in zip(drift_values, drift_strings):
            if drift == 0: #simulate many baserates to compare to overall rates of neurons with drift
                params['baseRates'] = [0.75, 1, 1.25, 1.5, 2, 2.5, 3.75, 5, 6.25, 7.5, 10, 12.5]
            else:
                params['baseRates'] = [1, 2, 5, 10]
            params['delta'] = drift

            print('in simulations, {0} drift'.format(drift))
            [pc] = simulateChangingContNeurons(params)

            savefile = resultsBasePath + '\\simulationsPC' + str(
                params['nSim']) + 'iter' + date_now + 'delta' + string + '.pickle'

            results = [pc, params]
            if params['savePCfile']:
                with open(savefile, 'wb') as handle:
                    pickle.dump(results, handle)


    else: # load previously saved simulation results
        # initialize dictionary to load simulation results for multiple drift levels
        pcDict = {}
        # initialize dictionary to load parameters for multiple runs of the simulation
        paramsDict = {}

        #parameters for loading file
        date_now = '_08_02'
        nIter = 500

        for drift, string in zip(drift_values, drift_strings):
            print('loading sim results {0} drift'.format(drift))
            savefile = resultsBasePath + '\\simulationsPC' + str(nIter)+ 'iter' + date_now + 'delta' + string + '.pickle'

            file = open(savefile, 'rb')
            results = pickle.load(file)
            file.close()

            pcDict[string] = results[0]
            paramsDict[string] = results[1]



    #now plot drift overlay for increasing and decreasing firing rates (subplots a-b)
    rpPlot = 2
    frPlot=5
    for driftDir in ['Inc', 'Dec']:
        params = paramsDict[driftDir]
        for frPlot in [1, 2, 5, 10]:
            figsavefile = resultsBasePath + '\\simulationsDrift' + driftDir + str(params['nSim']) + 'iterFR' + str(frPlot)
            plotDriftOverlay(pcDict, paramsDict, figsavefile, rpPlot=rpPlot, frPlotInput=frPlot, driftDir=driftDir)


