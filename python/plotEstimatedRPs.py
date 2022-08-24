# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 12:02:56 2022

@author: Noam Roth
"""
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

    ca1 = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'CA1'])
    po  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'PO'])
    dg  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'DG'])
    visa  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'PPC'])
    lp= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'LP'])
    
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

    hippocampus = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['CA', 'CA1', 'CA2', 'CA3','DG', 'POST', 'SUB']])
    cortex  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['ACA', 'AUD','ILA' , 'MOp', 'MOs',  'OLF', 'ORB', 'ORBm',
                     'PIR', 'PL', 'RSP', 'SSp','SSs',  'VISa', 'VISam', 'VISl',
                     'VISp', 'VISpm', 'VISrl']])
    thalamus= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['TH', 'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG','PO', 'POL', 
                     'PT','RT','SPF','VAL', 'VPL', 'VPM' ]])
    midbrain = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] in ['SNr','APN', 'IC','MB','MRN', 'NB','PAG','RN','SCig', 'SCm',  'SCs', 'SCsg']])
    
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

        
cortexAllSteinmetz = cortexAll
thalamusAllSteinmetz = thalamusAll
hippocampusAllSteinmetz = hippocampusAll
midbrainAllSteinmetz = midbrainAll
#%%        
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


    hippocampus = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'hippocampus'])
    cortex  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'cortex'])
    thalamus= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'thalamus'])
    midbrain = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['parentRegion'][i] == 'midbrain'])
    
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
nBins = 20
fig,axs = plt.subplots(1,1,figsize = (5,3))
ax = axs
cc = cortexAllRS[~np.isnan(cortexAllRS)]
xx = sum(~np.isnan(cortexAllRS))
ax.hist(cc, nBins, weights = np.ones(xx) / xx, histtype = 'step', color = 'blue', label = 'Cortex', linestyle = '-')

xx = sum(~np.isnan(thalamusAllRS))
tt = thalamusAllRS[~np.isnan(thalamusAllRS)]
ax.hist(tt, nBins,weights = np.ones(xx) / xx, histtype = 'step', color = 'green', label = 'Thalamus', linestyle = '-')

xx = sum(~np.isnan(hippocampusAllRS))
hh= hippocampusAllRS[~np.isnan(hippocampusAllRS)]

ax.hist(hh, nBins, weights = np.ones(xx) /xx, histtype = 'step', color = 'purple', label = 'Hippocampus', linestyle = '-')
ax.set_ylim(0,.3)
# plt.hist(data, weights=np.ones(len(data)) / len(data))



cc = cortexAllSteinmetz[~np.isnan(cortexAllSteinmetz)]
xx = sum(~np.isnan(cortexAllSteinmetz))
ax.hist(cc, nBins,  weights = np.ones(xx) / xx,histtype = 'step', color = 'blue', linestyle = '--')
xx = sum(~np.isnan(thalamusAllSteinmetz))
tt = thalamusAllSteinmetz[~np.isnan(thalamusAllSteinmetz)]
ax.hist(tt, nBins, weights = np.ones(xx) / xx, histtype = 'step', color = 'green', linestyle = '--')
xx = sum(~np.isnan(hippocampusAllSteinmetz))
hh= hippocampusAllSteinmetz[~np.isnan(hippocampusAllSteinmetz)]
ax.hist(hh, nBins, weights = np.ones(xx) / xx,  histtype = 'step', color = 'purple', linestyle = '--')

cc = cortexAllAllen[~np.isnan(cortexAllAllen)]
xx = sum(~np.isnan(cortexAllAllen))
ax.hist(cc, nBins,  weights = np.ones(xx) / xx,histtype = 'step', color = 'blue', linestyle = ':')
xx = sum(~np.isnan(thalamusAllAllen))
tt = thalamusAllAllen[~np.isnan(thalamusAllAllen)]
ax.hist(tt, nBins, weights = np.ones(xx) / xx, histtype = 'step', color = 'green', linestyle = ':')
xx = sum(~np.isnan(hippocampusAllAllen))
hh= hippocampusAllAllen[~np.isnan(hippocampusAllAllen)]
ax.hist(hh, nBins, weights = np.ones(xx) / xx,  histtype = 'step', color = 'purple', linestyle = ':')

ax.set_xlabel('Estimated RP (ms)')
ax.set_ylabel('Proportion of neurons')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)


#add arrows:
lenArrow = .2
lenHead = .01
wiArrow = .1
ind = np.nanmedian(cortexAllRS)
n = .25#len(cortexAll[(cortexAll > ind*0.95) & (cortexAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = '-')
ind = np.nanmedian(cortexAllSteinmetz)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = '--')
ind = np.nanmedian(cortexAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = ':')

ind = np.nanmedian(thalamusAllRS)
# n = len(thalamusAll[(thalamusAll > ind*0.95) & (thalamusAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = '-')
ind = np.nanmedian(thalamusAllSteinmetz) - .4
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = 'dashed')
ind = np.nanmedian(thalamusAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = ':')

ind = np.nanmedian(hippocampusAllRS)
# n = len(hippocampusAll[(hippocampusAll > ind*0.95) & (hippocampusAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = '-')
ind = np.nanmedian(hippocampusAllSteinmetz)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = '--')
ind = np.nanmedian(hippocampusAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = ':')







ax.plot(np.NaN, np.NaN, '-', color='black', label='IBL')
ax.plot(np.NaN, np.NaN, '--', color='black', label='Steinmetz 2019')
ax.plot(np.NaN, np.NaN, ':', color='black', label='Allen')

plt.legend(frameon=False)
plt.tight_layout()
# plt.suptitle('All Datasets',y=1.1)
fig.show()

plt.savefig(r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\estimatedRPs.pdf', dpi=300, format='pdf')
# cortexAllRS, cortexAllSteinmetz, cortexAllAllen

#%%

import matplotlib.pyplot as plt

from scipy.ndimage.filters import gaussian_filter1d
plt.rcParams.update({'font.size': 14})

def plot_metric(data, bins, x_axis_label, color, max_value=-1):
    
    h, b = np.histogram(data, bins=bins, density=True)

    x = b[:-1]
    y = gaussian_filter1d(h, 1)

    plt.plot(x, y, color=color)
    plt.xlabel(x_axis_label)
    plt.gca().get_yaxis().set_visible(False)
    [plt.gca().spines[loc].set_visible(False) for loc in ['right', 'top', 'left']]
    if max_value < np.max(y) * 1.1:
        max_value = np.max(y) * 1.1
    # plt.ylim([0, max_value])
    
    return max_value

bins = np.linspace(-3,2,100)
max_value = -np.inf

for idx, region in enumerate(region_dict.keys()):
    
    data = np.log10(units[units.ecephys_structure_acronym.isin(region_dict[region])]['firing_rate'])
    print('%s has %d neurons'%(region,len(data)))
    max_value = plot_metric(data, bins, 'log$_{10}$ firing rate (Hz)', color_dict[region], max_value)
    
_ = plt.legend(region_dict.keys())



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


#add arrows:
lenArrow = .03
lenHead = .01
wiArrow = .1
ind = np.nanmedian(cortexAllRS)
n = .25#len(cortexAll[(cortexAll > ind*0.95) & (cortexAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = '-')
ind = np.nanmedian(cortexAllSteinmetz)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = '--')
ind = np.nanmedian(cortexAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='blue', ec='none',linestyle = ':')

ind = np.nanmedian(thalamusAllRS)
# n = len(thalamusAll[(thalamusAll > ind*0.95) & (thalamusAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = '-')
ind = np.nanmedian(thalamusAllSteinmetz) - .4
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = 'dashed')
ind = np.nanmedian(thalamusAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='green', ec='none',linestyle = ':')

ind = np.nanmedian(hippocampusAllRS)
# n = len(hippocampusAll[(hippocampusAll > ind*0.95) & (hippocampusAll <  ind*1.05)])/2 #y value of overall histogram
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = '-')
ind = np.nanmedian(hippocampusAllSteinmetz)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = '--')
ind = np.nanmedian(hippocampusAllAllen)
ax.arrow(ind, n+lenArrow+lenHead, 0, -lenArrow, head_width=wiArrow*3, head_length=lenHead, width=wiArrow, fc='purple', ec='none',linestyle = ':')







ax.plot(np.NaN, np.NaN, '-', color='black', label='IBL')
ax.plot(np.NaN, np.NaN, '--', color='black', label='Steinmetz 2019')
ax.plot(np.NaN, np.NaN, ':', color='black', label='Allen')

plt.legend(frameon=False)
plt.tight_layout()
# plt.suptitle('All Datasets',y=1.1)
fig.show()

plt.savefig(r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\estimatedRPs.pdf', dpi=300, format='pdf')
# cortexAllRS, cortexAllSteinmetz, cortexAllAllen

#%%

xx = sum(~np.isnan(thalamusAllRS))
tt = thalamusAllRS[~np.isnan(thalamusAllRS)]
ax.hist(tt, nBins,weights = np.ones(xx) / xx, histtype = 'step', color = 'green', label = 'Thalamus', linestyle = '-')

xx = sum(~np.isnan(hippocampusAllRS))
hh= hippocampusAllRS[~np.isnan(hippocampusAllRS)]

ax.hist(hh, nBins, weights = np.ones(xx) /xx, histtype = 'step', color = 'purple', label = 'Hippocampus', linestyle = '-')
ax.set_ylim(0,.3)
# plt.hist(data, weights=np.ones(len(data)) / len(data))

