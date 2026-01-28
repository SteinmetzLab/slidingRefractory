import pickle
from slidingRP.simulations import *
import datetime

def runSaveFig4(figsavepath, resultsBasePath,rerunFig4 = False):




    if rerunFig4:
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


    else: #load simulation results from file

        pcDict = {}

        #parameters for loading file
        date_load = '_12_29' #date of simulations
        nIter = 1000 #number of iterations run for this run of simulations

        # set parameters for plotting
        frPlot = 2 #spks/s
        rpPlot = 2 #ms

        #load results separately for each confidence level
        for c, conf in enumerate([70, 80, 90]):

            #path to load results based on the confidence level:
            resultsPath = resultsBasePath + '\\simulationsPC' + str(nIter) + 'iter' + date_load + str(conf) + '.pickle'

            file = open(resultsPath, 'rb')
            results = pickle.load(file)
            file.close()
            pcDict[str(conf)] = results[0]

        params = results[-1] #also load parameters that these simulations were run under

        #We do not currently plot Hill results, so this loading is commented out (but can add these comparisons if wanted.)
        # pcDict['Hill 1.5ms'] = results[3]
        # pcDict['Hill 2ms'] = results[4]
        # pcDict['Hill 3ms'] = results[5]

        pcDict['Llobet 1.5ms'] = results[6]
        pcDict['Llobet 2ms'] = results[7]
        pcDict['Llobet 3ms'] = results[8]


    #subplots get plotted out of order here, we first plot subplot 3 which includes all confidence values, and then we plot
    # subplots 1 and 2 which only consider 90% confidence level.


    #Plot subplot 3: compare performance with different confidence parameters with a fixed RP:

    figsavefile = figsavepath + '\SubPlot3_confCompRP{}FR{}'.format(rpPlot, frPlot)

    #run function that plots standard simulations plot (pc vs. cont) for different simulations overlaid
    plotSimulationsOverlay(pcDict, params, figsavefile, rpPlot=rpPlot, frPlot=frPlot, legendLabels=['70', '80', '90',
                                                                                               'Llobet 1.5','Llobet 2','Llobet 3',
                                                                                              'Confidence'])

    #subplots 4/5:
    #%%
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    crVec = [8,12] #contamination rates to plot, each contamination rate plots as a separate subplot
    figsavefile = figsavepath + '\SubPlots4_5_confCompFixedContFR{}'.format(frPlot)
    plotConfRP(pcDict,crVec,fig, axs,params,figsavefile)



    #%% Now plot and save subplots 1 and 2, just showing standard rp metric simulations plot with rp 1 and 6
    rpPlot = 1 # first plot for RP=1 ms
    frPlot=2 # base firing rate for plotting, spks/s

    pcDict.pop('70') #remove 70% confidence from dictionary
    pcDict.pop('80') #remove 80% confidence from dictionary
    figsavefile = figsavepath + '\SubPlot1_RP{}FR{}'.format(rpPlot,frPlot)

    plotSimulationsOverlay(pcDict, params, figsavefile, rpPlot=rpPlot, frPlot = frPlot, legendLabels=['sliding Refractory',
                                                                                               'Llobet 1.5','Llobet 2','Llobet 3',
                                                                                              'Metric'], colorflag=True)

    rpPlot = 6
    frPlot=2

    figsavefile = figsavepath + '\SubPlot2_RP{}FR{}'.format(rpPlot,frPlot)

    plotSimulationsOverlay(pcDict, params, figsavefile, rpPlot=rpPlot, frPlot = frPlot, legendLabels=['sliding Refractory',
                                                                                               'Llobet 1.5','Llobet 2','Llobet 3',
                                                                                              'Metric'], colorflag=True)



#%%function for subplots 4 and 5 (compare conf as a function of RP, for a specific contamination value)

def plotConfRP(pcDict,crVec, fig, axs,params, figsavefile):
    pcDictKeys = list(pcDict.keys()) #load keys as different confidences / metrics for comparison
    recDurPlot = 2 # recording duration of the simulated neurons to plot, hours
    frPlot = 2 # base firing rate of the simulated neurons to plot, spks/s

    for c, crPlot in enumerate(crVec): #loop through all contamination rates
        ax = axs[c] #plot values for each contamination rate in a separate subplot

        for p, pc_key in enumerate(pcDictKeys): #loop through different confidences / metrics for comparison

            pc = pcDict[pc_key]
            print(pc_key)

            #compute error bars
            count = []
            count = pc / 100 * params['nSim']  # number of iterations of the simulation
            CI_scaled = binofit(count, params['nSim']) #compute confidence intervals
            CI = [x * 100 for x in CI_scaled] # convert from fractional to percent correct

            # find index to plot the contamination rate in crPlot
            cr = params['contRates']
            crInd = np.where(cr == crPlot / 100)[0]  # crPlot is fractional, convert to percent here


            # plot just  recDur = recDurPlot:
            recDurs = params['recDurs']
            rdInd = np.where(recDurs == recDurPlot)[0]

            # plot just fr = frPlot:
            frs = np.array(params['baseRates'])
            frInd = np.where(frs == frPlot)[0][0]

            # set up colors for plotting
            if type(pc_key) is str and pc_key[0] == 'H':  # if plotting, plot Hill results in red
                colors = matplotlib.cm.Reds(np.linspace(0.3, 1, 6))
                print(p)
                color = colors[p - 3]
            elif type(pc_key) is str and pc_key[0] == 'L':  # plot Llobet in green
                colors = matplotlib.cm.Greens(np.linspace(0.3, 1, 7))
                color = colors[p - 6]
            else:  # not hill (conf)
                colors = matplotlib.cm.Blues(np.linspace(0.3, 1, 3)) #plot slidingRP in blue
                color = colors[p]


            lowerCI = CI[0][rdInd[0], :, frInd, crInd[0]] #define lower and upper CI's for the picked plotting params
            upperCI = CI[1][rdInd[0], :, frInd, crInd[0]]

            rps = params['RPs'] #all possible RPs
            x = rps * 1000 # convert from ms to s for the plot
            y = pc[rdInd[0], :, frInd, crInd[0]] #percent correct values for the picked plotting params

            #plot percent correct as a function of RP
            ax.plot(x, y, '.-', color=color)
            ax.fill_between(x, lowerCI, upperCI, color=color, alpha=.3) #shaded error bars
            ax.set_ylabel('Percent pass')
            ax.set_xlabel('True RP')
            ax.set_title('Contamination = {0}%'.format(crPlot))
            ax.set_ylim(0,100)

            #remove top and right axes (aesthetic)
            spinesSetting = False
            ax.spines.right.set_visible(spinesSetting)
            ax.spines.top.set_visible(spinesSetting)


        handles, xx = ax.get_legend_handles_labels() #set up legend
        labels = ['Flexible RP metric; 70% confidence threshold', 'Flexible RP metric; 75% confidence threshold', 'Flexible RP metric; 80% confidence threshold', 'Hill metric; threshold = 2ms','Hill metric; threshold = 3ms']

    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1)) #only plot legend one time (outside of subplots loop)

    fig.show()
    fig.savefig(figsavefile + '_RP.svg', dpi=500)
    fig.savefig(figsavefile + '_RP.pdf', dpi=500)
