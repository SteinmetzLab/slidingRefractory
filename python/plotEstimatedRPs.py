# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 12:02:56 2022

@author: Noam Roth

Code to plot histograms of estimated RPs, as computed in computeEstimatedRPs.py

Runs for 3 datasets: IBL repeated site; Steinmetz 2019; Allen 

"""
#%% 

#set filter parameters for firing rate and amplitude across all 3 datasets
minFR = 20; minAmp = 100




#%%
#plot IBL data
one = ONE()
insertions = get_insertions(level=2, one=one, freeze='biorxiv_2022_05')

savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedRSFits\\'

nSess = len(insertions)
plotEach = False
ca1All = []
dgAll = []
poAll = []
lpAll = []
visaAll = []



for s in range(nSess):
    subject = insertions[s]['session']['subject']
    #load saved rpMetrics
    try:
        file = open(savefile + subject + '.pickle','rb')
    except:
        continue
    rpEstimates = pickle.load(file)
    file.close()

    ca1 = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'CA1' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    po  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'PO' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    dg  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'DG' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    visa  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'PPC' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    lp= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'LP' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    
    ca1[ca1<.05] =np.nan
    dg[dg<.05] = np.nan
    po[po<.05] =np.nan
    lp[lp<.05] = np.nan
    visa[visa<.05] =np.nan
    
    if plotEach:
        fig,axs = plt.subplots(2,2,figsize = (12,10))
        ax = axs[0,0]
        ax.hist(np.concatenate((ca1, dg)), 100)
        ax.set_title('hippocampus median: %.2f'%np.nanmedian(np.concatenate((ca1, dg))))
        ax = axs[0,1]
        ax.hist(np.concatenate((po, lp)), 100)
        ax.set_title('thalamus median: %.2f'%np.nanmedian(np.concatenate((po, lp))))
        
        ax = axs[1,0]
        ax.hist(visa, 100)
        ax.set_title('cortex median: %.2f'%np.nanmedian(visa))
    
    if len(ca1All) == 0:
        ca1All = ca1
        dgAll = dg
        lpAll = lp
        poAll = po
        visaAll = visa
    else:    
        ca1All = np.concatenate((ca1All,ca1))
        dgAll = np.concatenate((dgAll,dg))
        lpAll = np.concatenate((lpAll,lp))
        poAll = np.concatenate((poAll,po))
        visaAll = np.concatenate((visaAll,visa))
        
        

cortexAllRS = visaAll
thalamusAllRS = np.concatenate((poAll, lpAll))
hippocampusAllRS = np.concatenate((ca1All, dgAll))

#%%    plot all for rs 


fig,axs = plt.subplots(2,2,figsize = (8,6))
ax = axs[0,0]
ax.hist(cortexAllRS, 100)
ax.set_title('cortex (median rp: %.2f)'%np.nanmedian(cortexAllRS))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[0,1]
ax.hist(thalamusAllRS, 100)
ax.set_title('thalamus (median rp: %.2f)'%np.nanmedian(thalamusAllRS))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[1,0]
ax.hist(hippocampusAllRS, 100)
ax.set_title('hippocampus (median rp: %.2f)'%np.nanmedian(hippocampusAllRS))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

axs[1, 1].set_axis_off()
plt.tight_layout()
plt.suptitle('IBL (repeated site)',y=1.1)

fig.show()



#%% plot steinmetz
savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedSteinmetzFits\\'

one = ONE(cache_dir='E:\Steinmetz2019\Steinmetz_et_al_2019_9974357\9974357') # The location of the unarchived data
sessions = one.search()#dataset='trials') # search for all sessions that have a `trials` object


BRAIN_REGIONS = ['ACA', 'AUD','ILA' , 'MOp', 'MOs',  'OLF', 'ORB', 'ORBm',
                     'PIR', 'PL', 'RSP', 'SSp','SSs',  'VISa', 'VISam', 'VISl',
                     'VISp', 'VISpm', 'VISrl',
                 'CA', 'CA1', 'CA2', 'CA3','DG', 'POST', 'SUB'
                 'TH', 'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG','PO', 'POL', 
                     'PT','RT','SPF','VAL', 'VPL', 'VPM', 
                 'SNr','APN', 'IC','MB','MRN', 'NB','PAG','RN','SCig', 'SCm', 
                     'SCs', 'SCsg']
rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}
nSess = len(sessions)
plotEach = False
cortexAll = []
hippocampusAll = []
thalamusAll = []
midbrainAll = []

for e,eid in enumerate(sessions):

    #load saved rpMetrics
    try:
        file = open(savefile + eid + '.pickle','rb')
    except:
        continue
    rpEstimates = pickle.load(file)
    file.close()



    #find neurons in each of the target brain regions that pass the firing rate and amplitude criteria
    hippocampus = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['CA', 'CA1', 'CA2', 'CA3','DG', 'POST', 'SUB'] and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    cortex  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['ACA', 'AUD','ILA' , 'MOp', 'MOs',  'OLF', 'ORB', 'ORBm',
                     'PIR', 'PL', 'RSP', 'SSp','SSs',  'VISa', 'VISam', 'VISl',
                     'VISp', 'VISpm', 'VISrl'] and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    thalamus= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['TH', 'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG','PO', 'POL', 
                     'PT','RT','SPF','VAL', 'VPL', 'VPM' ] and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    midbrain = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['SNr','APN', 'IC','MB','MRN', 'NB','PAG','RN','SCig', 'SCm',  'SCs', 'SCsg'] and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    
    
    #get rid of very low values #todo -- change this once filtered by neurons code is done!!
    hippocampus[hippocampus<.05] =np.nan
    thalamus[thalamus<.05] = np.nan
    cortex[cortex<.05] =np.nan
    midbrain[midbrain<.05] =np.nan
    
    #if this flag is True, plot the histogram for each session individually (sanity check)
    if plotEach:
        fig,axs = plt.subplots(2,2,figsize = (12,10))
        ax = axs[0,0]
        ax.hist(np.concatenate((ca1, dg)), 100)
        ax.set_title('hippocampus median: %.2f'%np.nanmedian(np.concatenate((ca1, dg))))
        ax = axs[0,1]
        ax.hist(np.concatenate((po, lp)), 100)
        ax.set_title('thalamus median: %.2f'%np.nanmedian(np.concatenate((po, lp))))
        
        ax = axs[1,0]
        ax.hist(visa, 100)
        ax.set_title('cortex median: %.2f'%np.nanmedian(visa))
    
    
    #concatenate each session's values for each region to a vector of all values for that region
    if len(hippocampusAll) == 0:
        hippocampusAll = hippocampus
    else:    
        hippocampusAll = np.concatenate((hippocampusAll, hippocampus))
        
    if len(cortexAll) == 0:
        cortexAll = cortex
    else:
        cortexAll = np.concatenate((cortexAll, cortex))

    if len(midbrainAll) == 0:
        midbrainAll = midbrain
    else:
        midbrainAll = np.concatenate((midbrainAll, midbrain))
        
    if len(thalamusAll) == 0:
        thalamusAll = thalamus
    else:
        thalamusAll = np.concatenate((thalamusAll, thalamus))

        
cortexAllSteinmetz = cortexAll
thalamusAllSteinmetz = thalamusAll
hippocampusAllSteinmetz = hippocampusAll
midbrainAllSteinmetz = midbrainAll
#%%   
#plot estimated RPs for all sessions from the Steinmetz 2019 dataset     
fig,axs = plt.subplots(2,2,figsize = (8,6))
ax = axs[0,0]
ax.hist(cortexAllSteinmetz, 100)
ax.set_title('cortex (median rp: %.2f)'%np.nanmedian(cortexAllSteinmetz))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[1,0]
ax.hist(hippocampusAllSteinmetz, 100)
ax.set_title('hippocampus (median rp: %.2f)'%np.nanmedian(hippocampusAllSteinmetz))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[0,1]
ax.hist(thalamusAllSteinmetz, 100)
ax.set_title('thalamus (median rp: %.2f)'%np.nanmedian(thalamusAllSteinmetz))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)


ax = axs[1,1]
ax.hist(midbrainAllSteinmetz, 100)
ax.set_title('midbrain (median rp: %.2f)'%np.nanmedian(midbrainAllSteinmetz))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

plt.tight_layout()
plt.suptitle('Steinmetz 2019',y=1.1)

fig.show()

 #%% polot allen       
rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}
nSess = len(sessions)
plotEach = False
cortexAll = []
hippocampusAll = []
thalamusAll = []
midbrainAll = []

#find session names (from the saved files)
from os import listdir
from os.path import isfile, join
mypath = 'E:\AllenBrainObservatory\saved_units'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
sessions = [i[11:20] for i in onlyfiles]

for j, session_id in enumerate(sessions):
    file_name = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedAllenFits\%s.pickle'%session_id
    print('loading data for session %d out of %d'%(j+1, len(sessions)))

    #load the dataframe for this session
    with open(file_name, 'rb') as f:
        rpEstimates = pickle.load(f)


    hippocampus = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'hippocampus' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    cortex  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'cortex' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    thalamus= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'thalamus' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    midbrain = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'midbrain' and  rpEstimates['amp'][i]>minAmp and rpEstimates['fr'][i]>minFR])
    
    hippocampus[hippocampus<.05] =np.nan
    thalamus[thalamus<.05] = np.nan
    cortex[cortex<.05] =np.nan
    midbrain[midbrain<.05] =np.nan
    
    if plotEach:
        fig,axs = plt.subplots(2,2,figsize = (12,10))
        ax = axs[0,0]
        ax.hist(np.concatenate((ca1, dg)), 100)
        ax.set_title('hippocampus median: %.2f'%np.nanmedian(np.concatenate((ca1, dg))))
        ax = axs[0,1]
        ax.hist(np.concatenate((po, lp)), 100)
        ax.set_title('thalamus median: %.2f'%np.nanmedian(np.concatenate((po, lp))))
        
        ax = axs[1,0]
        ax.hist(visa, 100)
        ax.set_title('cortex median: %.2f'%np.nanmedian(visa))
    
    if len(hippocampusAll) == 0:
        hippocampusAll = hippocampus
    else:    
        hippocampusAll = np.concatenate((hippocampusAll, hippocampus))
        
    if len(cortexAll) == 0:
        cortexAll = cortex
    else:
        cortexAll = np.concatenate((cortexAll, cortex))

    if len(midbrainAll) == 0:
        midbrainAll = midbrain
    else:
        midbrainAll = np.concatenate((midbrainAll, midbrain))
        
    if len(thalamusAll) == 0:
        thalamusAll = thalamus
    else:
        thalamusAll = np.concatenate((thalamusAll, thalamus))

        
cortexAllAllen = cortexAll
thalamusAllAllen = thalamusAll
hippocampusAllAllen = hippocampusAll
midbrainAllAllen = midbrainAll



#%%        
fig,axs = plt.subplots(2,2,figsize = (8,6))
ax = axs[0,0]
ax.hist(cortexAllAllen, 100)
ax.set_title('cortex (median rp: %.2f)'%np.nanmedian(cortexAllAllen))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[0,1]
ax.hist(thalamusAllAllen, 100)
ax.set_title('thalamus (median rp: %.2f)'%np.nanmedian(thalamusAllAllen))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[1,0]
ax.hist(hippocampusAllAllen, 100)
ax.set_title('hippocampus (median rp: %.2f)'%np.nanmedian(hippocampusAllAllen))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)


ax = axs[1,1]
ax.hist(midbrainAllAllen, 100)
ax.set_title('midbrain (median rp: %.2f)'%np.nanmedian(midbrainAllAllen))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

plt.tight_layout()
plt.suptitle('Allen',y=1.1)

fig.show()


#%%

#concatenate all datasets

cortexAllDatasets = np.concatenate((cortexAllRS, cortexAllSteinmetz, cortexAllAllen))
midbrainAllDatasets = np.concatenate((midbrainAllSteinmetz, midbrainAllAllen))
thalamusAllDatasets = np.concatenate((thalamusAllRS, thalamusAllSteinmetz, thalamusAllAllen))
hippocampusAllDatasets = np.concatenate((hippocampusAllRS, hippocampusAllSteinmetz,hippocampusAllAllen))





fig,axs = plt.subplots(2,2,figsize = (8,6))
ax = axs[0,0]
ax.hist(cortexAllDatasets, 100)
ax.set_title('cortex (median rp: %.2f)'%np.nanmedian(cortexAllDatasets))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[0,1]
ax.hist(thalamusAllDatasets, 100)
ax.set_title('thalamus (median rp: %.2f)'%np.nanmedian(thalamusAllDatasets))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax = axs[1,0]
ax.hist(hippocampusAllDatasets, 100)
ax.set_title('hippocampus (median rp: %.2f)'%np.nanmedian(hippocampusAllDatasets))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)


ax = axs[1,1]
ax.hist(midbrainAllDatasets, 100)
ax.set_title('midbrain (median rp: %.2f)'%np.nanmedian(midbrainAllDatasets))
ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

plt.tight_layout()
plt.suptitle('All Datasets',y=1.1)
fig.show()



#%%

#non-smoothed histograms, all data


fig,axs = plt.subplots(1,1,figsize = (5,3))
ax = axs
ax.hist(cortexAllDatasets, 100, histtype = 'step', color = 'blue', label = 'Cortex')
ax.hist(thalamusAllDatasets, 100, histtype = 'step', color = 'green', label = 'Thalamus')
ax.hist(hippocampusAllDatasets, 100,  histtype = 'step', color = 'purple', label = 'Hippocampus')

ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Number of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)


plt.legend()
plt.tight_layout()
plt.suptitle('All Datasets',y=1.1)
fig.show()






#%%
from scipy.ndimage.filters import gaussian_filter1d

nBins = 20
fig,axs = plt.subplots(1,1,figsize = (5,3))
ax = axs

def plot_hists(ax, data, color, linestyle, label= None):
    cc = data[~np.isnan(data)]
    xx = sum(~np.isnan(data))
    h, b = np.histogram(cc, bins=nBins, weights = np.ones(xx) / xx, density=False)
    x = b[:-1]
    y = gaussian_filter1d(h, 1)
    ax.plot(x, y, color=color,linestyle = linestyle, label = label)

plot_hists(ax, cortexAllRS, 'blue', '-', 'cortex')
plot_hists(ax, thalamusAllRS, 'green', '-', 'thalamus')
plot_hists(ax, hippocampusAllRS, 'purple', '-', 'hippocampus')

plot_hists(ax, cortexAllSteinmetz, 'blue', '--')
plot_hists(ax, thalamusAllSteinmetz, 'green', '--')
plot_hists(ax, hippocampusAllSteinmetz, 'purple', '--')

plot_hists(ax, cortexAllAllen, 'blue', ':')
plot_hists(ax, thalamusAllAllen, 'green', ':')
plot_hists(ax, hippocampusAllAllen, 'purple', ':')

ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Proportion of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

def add_median_arrows(ax, data, y_value, lenArrow, lenHead, wiArrow, color, linestyle):
    ind = np.nanmedian(data)
    ax.annotate('', xy=(ind, y_value), xytext=(ind,y_value + lenArrow + lenHead),
        arrowprops={'arrowstyle': '->','color':color, 'ls': linestyle})
    
#add arrows:
lenArrow = .05
lenHead = 0
wiArrow = 0
ind = np.nanmedian(cortexAllRS)
n = 0.2#len(cortexAll[(cortexAll > ind*0.95) & (cortexAll <  ind*1.05)])/2 #y value of overall histogram

add_median_arrows(ax, cortexAllRS, n, lenArrow, lenHead, wiArrow, 'blue', 'solid')
add_median_arrows(ax, cortexAllSteinmetz, n, lenArrow, lenHead, wiArrow, 'blue', 'dashed')
add_median_arrows(ax, cortexAllAllen, n, lenArrow, lenHead, wiArrow, 'blue', 'dotted')

add_median_arrows(ax, thalamusAllRS, n,  lenArrow, lenHead, wiArrow, 'green', 'solid')
add_median_arrows(ax, thalamusAllSteinmetz, n, lenArrow, lenHead, wiArrow, 'green', 'dashed')
add_median_arrows(ax, thalamusAllAllen, n, lenArrow, lenHead, wiArrow, 'green', 'dotted')

add_median_arrows(ax, hippocampusAllRS, n, lenArrow, lenHead, wiArrow, 'purple', 'solid')
add_median_arrows(ax, hippocampusAllSteinmetz, n, lenArrow, lenHead, wiArrow, 'purple', 'dashed')
add_median_arrows(ax, hippocampusAllAllen, n, lenArrow, lenHead, wiArrow, 'purple', 'dotted')




ax.plot(np.NaN, np.NaN, '-', color='black', label='IBL')
ax.plot(np.NaN, np.NaN, '--', color='black', label='Steinmetz 2019')
ax.plot(np.NaN, np.NaN, ':', color='black', label='Allen')

plt.legend(frameon=False)
plt.tight_layout()
fig.show()

plt.savefig(r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\estimatedRPs%dAmp_%dFR.pdf'%(minAmp,minFR), dpi=300, format='pdf')


#%% TODO tomorrow:
    
    # scatter plots: rpEstimate vs fr and rpEstimate vs amplitude
    # histograms with different subselections: change amp and fr and run this code for a few different seleections of the parameters
    # if time: look into low thalamus example neurons!