#%%
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

def runSaveFig5_ROC(figSavePath,resultsBasePath, rerunFig5=False):


    if rerunFig5:
        # set up parameters to rerun simulations
        sampleRate = 30000
        params = {
            'recDurs': np.array([0.5, 1, 2, 3]),  # recording durations (hours)
            'RPs': np.array([0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006]),  # true RP (s)
            'baseRates': [0.5, 1, 2, 5, 10],  # (spk/s)
            'contRates': np.arange(0.00, 0.21, 0.02),  # contamination levels (proportion)
            'nSim': 1000,
            'contaminationThresh': 10,
            'binSize': 1 / sampleRate,
            'sampleRate': sampleRate,
            'confidenceThresh': 90, #default confidence value (%)
            'checkFR': False,
            'binSizeCorr': 1 / sampleRate,
            'returnMatrix': True,
            'verbose': True,
            'runLlobet': True,
            'runLlobetPoiss': True
        }

        #set date for saving results to file:
        date_now = datetime.datetime.now().strftime('_%m_%d')
        #set different confidence values to run:
        confidence_values = [60, 70, 80, 90, 99] #confidence (%)
        #initialize dictionary to save simulation results for multiple confidence levels
        pcDict = {}
        for conf in confidence_values:
            params['confidenceThresh'] = conf
            print('in simulations, {0} conf'.format(conf))

            # run simulations with parameters above (params) and confidence level (conf)
            [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
             pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3] = simulateContNeurons(params)

            # save all results to variable
            results = [pc, pc2MsNoSpikes, pcHalfInactive, pcHill15, pcHill2, pcHill3, pcLlobet15, pcLlobet2, pcLlobet3,
                       pcLlobetPoiss15, pcLlobetPoiss2, pcLlobetPoiss3, params]

            #save results to file indicating date and confidence level
            savefile = r'C:\Users\noamroth\int-brain-lab\slidingRefractory\python\slidingRP\simulationsPC' + str(
                params['nSim']) + 'iter' + date_now + str(conf) + '.pickle'


            with open(savefile, 'wb') as handle:
                pickle.dump(results, handle)

            pcDict[str(conf)] = results[0]

        # We do not currently plot Hill results, so this loading is commented out (but can add these comparisons if wanted.)
        # pcDict['Hill 1.5ms'] = results[3]
        # pcDict['Hill 2ms'] = results[4]
        # pcDict['Hill 3ms'] = results[5]

        pcDict['Llobet 1.5ms'] = results[6]
        pcDict['Llobet 2ms'] = results[7]
        pcDict['Llobet 3ms'] = results[8]

    else:
        #load datasets (each set is a different confidence level and concatenate them into a dictionary
        pcDictAllConf = {} #initialize dictionary (percent correct matrix for each confidence value)
        confidence_values = [50, 60, 70, 75, 80, 85, 90, 99] #percent
        dates = ['_07_19', '_07_25', '_07_25', '_07_19','_07_25', '_07_25', '_07_19', '_07_19'] #datasets saved on different days
        nIter = 500

        for conf, date in zip(confidence_values, dates):
            print('loading sim results {0} conf'.format(conf))
            resultsPath = resultsBasePath + '\\simulationsPC' + str(nIter) + 'iter' + date + str(conf) + '.pickle'

            file = open(resultsPath, 'rb')
            results = pickle.load(file)
            file.close()

            pcDictAllConf[conf] = results[0] #save into pcDict: key = confidence, value = pc matrix

        #also load parameter file to have all parameters used to run these simulations (they are the same across confidences)
        params = results[-1]

        #load corresponding Hill results for comparison
        pcDictAllConf['Hill 1.5ms'] = results[3]
        pcDictAllConf['Hill 2ms'] = results[4]
        pcDictAllConf['Hill 3ms'] = results[5]

    #set data parameters for plotting:
    frPlot = 2 # spks/s, base firing rate of the simulated neurons, parameter picked for plot
    rpPlot = 2 #ms, refractory period of the simulated neurons, parameter picked for plot
    recDurPlot = 2 #hours, recording duration of the simulated neurons, parameter picked for plot

    #set contamination values for which to compute ROC
    contAvg = False #use the average across all contaminations on each side of contThresh (uncontaminated, or contaminated)
    fixedCont = False #use only one distance from the contamination threshold
    if fixedCont == True : #then use one fixed point, make it an array to be compatible with code below
        threshDists = np.array([0.02]) #2% contamination (added or subtracted from contThresh) is default threshold to test
    else:
        threshDists = np.arange(0.02, 0.1, 0.02) #different distances from the contamination threshold to test and plot

    #get list of confidence values to plot
    pcDict_keys = list(pcDictAllConf.keys())
    pcDict_keys = pcDict_keys[0:-2] #remove last 2 values, this ensures we only plot one of the Hill comparison results

    color = iter(cm.viridis(np.linspace(0, 1, len(threshDists)))) #pick color for each trace on ROC plot

    fig, ax = plt.subplots(1, 1, figsize=(4, 4)) #only one subplot for ROC plot
    for threshDist in threshDists:
        c = next(color)

        TPR = [] #initialize true positive rate vector
        FPR = [] #initialize true positive rate vector
        for p, pc_key in enumerate(pcDict_keys): #loop through different levels of confidence
            pc = pcDictAllConf[pc_key]

            # all possible contamination rates
            cr = params['contRates'];

            # plot just RP = rpPlot:
            rps = params['RPs']
            rpInd = np.where(rps == rpPlot/1000)[0] # rpPlot in ms, convert to s here

            #plot just  recDur = recDurPlot:
            recDurs = params['recDurs']
            rdInd = np.where(recDurs == recDurPlot)[0]

            # plot just base firing rate = frPlot
            frs = np.array(params['baseRates'])
            frInd = np.where(frs == frPlot)[0][0]

            x = cr * 100 # convert contamination rate to percent contamination
            y = pc[rdInd[0], rpInd[0], frInd, :] #percent correct values for the indicated recDur,RP,andFR

            contThresh = params['contaminationThresh'] / 100  # convert to percent as in contamination rate

            if contAvg: #in this condition, average across all contaminations to one side of contThresh
                uncontPerf = np.nanmean(y[cr < contThresh])
                contPerf = np.nanmean(y[cr > contThresh])
            else: #otherwise, use threshDist as defined above
                uncontPerf = y[cr == np.round((contThresh - threshDist),2)]
                contPerf = y[cr == np.round((contThresh + threshDist),2)]

            TPR.append(uncontPerf) # add performance of metric on uncontaminated neurons to true positive rate
            FPR.append(contPerf) # add performance of metric on contaminated neurons to false positive rate


        ax.plot(FPR[0:-1],TPR[0:-1],'.', color = c,label = threshDist) #adds a point to ROC curve for this distance from the threshold
        ax.plot(FPR[-1],TPR[-1],'x', color = c) #adds the point from the Hill comparison to ROC curve for this distance from the threshold


    ax.set_ylabel('True Positive Rate')
    ax.set_xlabel('False Positive Rate')
    ax.set_title('ROC')
    ax.set_ylim(0,102) #set slightly outside of 100 to see points on the border
    ax.set_xlim(0,102)

    #add legend for different distances from contThresh
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(.7, .4), title = 'Distance from threshold',frameon=False)

    #now add a legend for slidingRP vs. Hill metric (by plotting points off of the figure boundaries and labeling them)
    ax.plot(-1, -1, '.', color='k', label='slidingRP')
    ax.plot(-1,-1, 'x', color='k', label='Hill')
    handles, labels = ax.get_legend_handles_labels()
    handles = handles[-2:]
    labels = labels[-2:]
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(.9, .4), title = 'Metric',frameon=False)

    spinesSetting = True #keep outer bounding box for ROC plot (unlike most other plots)
    ax.spines.right.set_visible(spinesSetting)
    ax.spines.top.set_visible(spinesSetting)

    #show and save figure
    fig.show()
    fig.savefig(figSavePath + '\ROC.svg', dpi=500)
    fig.savefig(figSavePath + '\ROC.pdf', dpi=500)
