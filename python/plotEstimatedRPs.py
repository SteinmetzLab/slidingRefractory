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
import pickle
#now save for use in iblenv
filename = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\metrics\slidingRP\ts_region_dict_all.pkl'

file = open(filename,'rb')

ts_dict = pickle.load(file)
file.close()


savefile = r'C:\Users\Steinmetz Lab User\Documents\GitHub\analysis\slidingRefractory\python\savedAllenFits\\'

rpBinSize = 1 / 30000  
rpEdges = np.arange(0, 10/1000, rpBinSize) # in s  
rp = rpEdges + np.mean(np.diff(rpEdges)[0])/2 # vector of refractory period durations to test 
params = {}

for b, brainRegion in enumerate(ts_dict.keys()):
    #load saved rpMetrics
    try:
        file = open(savefile + brainRegion + '.pickle','rb')
    except:
        continue
    rpEstimates = pickle.load(file)
    file.close()

    if brainRegion == 'hippocampus':
        hippocampus  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'hippocampus'])
    elif brainRegion == 'cortex':
        cortex  = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'cortex'])
    elif brainRegion == 'thalamus':
        thalamus= np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'thalamus'])
    elif brainRegion == 'midbrain':
        midbrain = np.array([rpEstimates['rpEstimate'][i] for i in range(len(rpEstimates['brainRegion'])) if rpEstimates['brainRegion'][i] == 'midbrain'])#,'SNr','APN', 'IC','MB','MRN', 'NB','PAG','RN','SCig', 'SCm',  'SCs', 'SCsg']])
    
hippocampus[hippocampus<.05] =np.nan
cortex[cortex<.05] = np.nan
thalamus[thalamus<.05] =np.nan
midbrain[midbrain<.05] =np.nan


cortexAllAllen = cortex
thalamusAllAllen = thalamus
hippocampusAllAllen = hippocampus
midbrainAllAllen = midbrain

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

