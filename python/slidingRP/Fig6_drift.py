import pickle
from slidingRP.simulations import *
import datetime


def runSaveFig6(figSavePath, resultsBasePath, rerunFig6 = False):




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


        drift_values = [-0.5,0, 0.5]
        drift_strings = ['Dec','Same','Inc']


        date_now = datetime.datetime.now().strftime('_%m_%d')
        # run simulations with parameters above (params)
        [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
         pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3] = simulateContNeurons(params)

        # save all results to variable
        results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
                   pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3, params]

        # save results to pickle file
        date_now = datetime.datetime.now().strftime('_%m_%d')
        savefile = resultsBasePath + '\\simulationsPC' + str(params['nSim']) + 'iter' + date_now + '.pickle'

        with open(savefile, 'wb') as handle:
            pickle.dump(results, handle)

    else:
    # load previously saved simulation results