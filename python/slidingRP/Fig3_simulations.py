import pickle
import datetime

from slidingRP.simulations import *

def runSaveFig3(figSavePath, resultsBasePath,rerunFig3 = False):


    if rerunFig3:
        #set up parameters to rerun simulations
        sampleRate = 30000
        params = {
            'recDurs': np.array([0.5, 1, 2, 3]),  # recording durations (hours)
            'RPs': np.array([0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006]), #true RP (s)
            'baseRates': [0.5, 1, 2, 5, 10],#   (spk/s)
            'contRates': np.arange(0.00, 0.21, 0.01), #contamination levels (proportion)
            'nSim': 1000,
            'contaminationThresh': 10,
            'binSize': 1 / sampleRate,
            'sampleRate': sampleRate,
            'confidenceThresh': 90,
            'checkFR': False,
            'binSizeCorr': 1 / sampleRate,
            'returnMatrix': True,
            'verbose': True,
            'runLlobet': True,
            'runLlobetPoiss': True
        }


        # params = {
        #     'recDurs': np.array([ 2]),  # recording durations (hours)
        #     'RPs': np.array([0.0025]), #true RP (s)
        #     'baseRates': [0.5, 1, 2, 5, 10],#   (spk/s)
        #     'contRates': np.arange(0.00, 0.21, 0.01), #contamination levels (proportion)
        #     'nSim': 1000,
        #     'contaminationThresh': 10,
        #     'binSize': 1 / sampleRate,
        #     'sampleRate': sampleRate,
        #     'confidenceThresh': 90,
        #     'checkFR': False,
        #     'binSizeCorr': 1 / sampleRate,
        #     'returnMatrix': True,
        #     'verbose': True,
        #     'runLlobet': True,
        #     'runLlobetPoiss': True
        # }

        #run simulations with parameters above (params)
        [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
         pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3] = simulateContNeurons(params)

        #save all results to variable
        results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
                   pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3, params]

        #save results to pickle file
        date_now = datetime.datetime.now().strftime('_%m_%d')
        savefile = resultsBasePath + '\\simulationsPC' + str(params['nSim']) + 'iter' + date_now +'.pickle'

        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)

    else:
        #load previously saved simulation results
        nIter = 1000 #number of simulation iterations run
        date_now = '_12_29' #date simulations were run


        resultsPath = resultsBasePath + '\\simulationsPC' + str(nIter) + 'iter' + date_now + '.pickle'


        file = open(resultsPath, 'rb')
        results = pickle.load(file)
        file.close()

    #whether these were just loaded or rerun, load percent correct matrix and parameters for plotting
    pc = results[0]
    params = results[-1]


    #make RP plots for regular and Hill
    figsavefile = figSavePath + '\\fig3'

    plotSimulations(pc, params,figsavefile, subplot1=True)
    plotSimulations(pc, params,figsavefile, subplot2=True)
    plotSimulations(pc, params,figsavefile, subplot3=True)
    plotSimulations(pc, params,figsavefile, subplot4=True)
