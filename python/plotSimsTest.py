# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 09:20:53 2022

@author: Noam Roth
plot simulations
"""


def plotSimulations(pc,params, savefile):
 
    colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(params['baseRates'])))
    fig,axs = plt.subplots(len(params['recDurs'][2:3]),len(params['RPs']), figsize = (12,3*len(params['recDurs'][2:3])))

    for j, recDur in enumerate(params['recDurs'][2:3]):
        for i, rp in enumerate(params['RPs']):
            
            if len(params['recDurs'][2:3]) > 1 and len(params['RPs'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['contRates'], pc[j, i, b,:], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('Proportion contamination')
                ax.set_title('True RP %d ms'%(rp*1000))
        fig.text(0.425, 0.9-(.17*j), 'Recording duration: %d hours'%recDur)
   
    handles, labels = ax.get_legend_handles_labels()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')


    fig.savefig(savefile, dpi = 500)
    
    
    fig,axs = plt.subplots(len(params['contRates'][::2]),len(params['RPs']), figsize = (12*2,3*len(params['recDurs'])))
    for j, contRate in enumerate(params['contRates'][::2]):
        print(j)
        for i, rp in enumerate(params['RPs']):
            
            if len(params['recDurs']) > 1 and len(params['RPs'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['recDurs'], pc[:, i, b,j*2], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('recording Duration')
                ax.set_title('True RP %d ms'%(rp*1000))
        fig.text(0.65, 0.9-(.17*j), 'Proportion contamination: %.2f'%contRate)
    fig.subplots_adjust(left=0.5, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')
    fig.show()
    
    
    
    fig,axs = plt.subplots(len(params['recDurs']), len(params['contRates'][::2]), figsize = (12*2,3*len(params['recDurs'])))
    for j, recDur in enumerate(params['recDurs']):
        print(j)
        for i, contRate in enumerate(params['contRates'][::2]):
            
            if len(params['recDurs']) > 1 and len(params['contRates'])>1:
                ax = axs[j,i]
            else:
                ax = axs[(j+1)*i]
            #different base rates get different colors
            for b, baseRate in enumerate(params['baseRates']):
                ax.plot(params['RPs'], pc[j,:, b,i*2], '.-',color = colors[b], label = baseRate)
                ax.set_ylabel('Percent pass')
                ax.set_xlabel('True RP')
                ax.set_title('contamination %.2f '%contRate)
        fig.text(0.65, 0.9-(.17*j), 'Recording Duration: %d hours'%recDur)
    fig.subplots_adjust(left=0.5, bottom=None, right=None, top=None, wspace=0.5, hspace=1.2)

    fig.legend(handles, labels, loc='upper right')
    fig.show()


savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pickle'

sampleRate = 30000
params = {
    'recDurs': np.array([0.25, 0.5, 1, 2, 5]),  #recording durations (hours)
    'RPs': np.array([0.001,0.0015, 0.002, 0.003, 0.004, 0.005]), #true RP (s)
    'baseRates': np.array([ 0.5, 1, 2, 5, 10, 20 ]), #F1, 2, 5, 10 , 20 R (spk/s)
    'contRates':  np.arange(0.00,0.225, 0.025), #contamination levels (proportion) #.025
    'nSim': 20,
    'threshold': 0.1,
    'binSize': 1 / sampleRate,
    'sampleRate': 30000,  #TODO figure out a way to refer to this in binsize?
    'checkFR': False,
    'binSizeCorr': 1 / sampleRate,
    'returnMatrix': True,
    'verbose': True ,
    'savePCfile': True
    
}


file = open(savefile,'rb')
pc = pickle.load(file)
file.close()

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC8iter_08_25.pickle'

file = open(savefile,'rb')
pc2 = pickle.load(file)
file.close()

figsavefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\simulationsPC20iter.pdf'
plotSimulations(pc, params, figsavefile)
