'''
Figure 1A
Plot ACG for 3 example neurons
And tilted slice with their location as dot

To get the data, run:
slidingRefractory/python/slidingRP/elts/transform_acgs_into_rps.py
'''
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from iblatlas import atlas

data_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Data")

# Atlas init
ba = atlas.AllenAtlas(25)

# Load files (ACG and DF)
file = data_path.joinpath('ibl').joinpath(f'acgs.npy')
acgs_ibl = np.load(file, allow_pickle=True)

file = data_path.joinpath('ibl').joinpath(f'clusters_df.pqt')
df_clusters_good = pd.read_parquet(file)
##
# Find an insertion for which these conditions are met:
# - 1 neuron with high firing rate, in these 3 regions: HPF, TH and Isocortex

fr_min = 2  # Hz
region_set = ['HPF', 'TH', 'Isocortex']

# Keep only units with firing rate > fr_min
indx_good = np.where(df_clusters_good['firing_rate'] > fr_min)[0]
df_clusters_good = df_clusters_good.iloc[indx_good]
acgs_ibl = acgs_ibl[indx_good, :]

# Keep only units within region set
indx_good = np.where(df_clusters_good['Cosmos_acronym'].isin(region_set))[0]
df_clusters_good = df_clusters_good.iloc[indx_good]
acgs_ibl = acgs_ibl[indx_good, :]

##
# Group by PID and find insertion with of unit count in regions valid
# Group by PID and region
a = df_clusters_good.groupby(['pid', 'Cosmos_id']).count()
a = a.reset_index()
# Group by PID again, count the number of region and threshold
b = a.groupby(['pid']).count()
indx_pid = np.where(b['Cosmos_id'] == len(region_set))[0]
pids = b.iloc[indx_pid].reset_index()['pid']
##
# Keep only units within pids set
indx_good = np.where(df_clusters_good['pid'].isin(pids.to_list()))[0]
df_clusters_good = df_clusters_good.iloc[indx_good]
acgs_ibl = acgs_ibl[indx_good, :]

##
# Take random PID and plot tilted slice
pid = pids[5]

# Keep only units within 1 pid
indx_good = np.where(df_clusters_good['pid'] == pid)[0]
df_clusters_good = df_clusters_good.iloc[indx_good]
acgs_ibl = acgs_ibl[indx_good, :]

# Get channels xyz:
xyz = df_clusters_good[['x', 'y', 'z']].to_numpy()
##
# Plot
# Want to plot without the traj as line
plt.figure()
ax = plt.gca()
ax.axis('equal')
cmap = plt.get_cmap('bone')

tslice, width, height, depth = ba.tilted_slice(xyz, axis=1, volume='image')
width = width * 1e6
height = height * 1e6
depth = depth * 1e6

ax.imshow(tslice, extent=np.r_[width, height], cmap=cmap)
'''
# Use this for debugging
ba.plot_tilted_slice(xyz, axis=1)
'''
##
# Chose 3 random neurons

for region in region_set:
    indx_good = np.where(df_clusters_good['Cosmos_acronym'] == region)[0]
    df_neuron = df_clusters_good.iloc[indx_good[0]]  # select first neuron
    acgs_neuron = acgs_ibl[indx_good, :]

    # Plot neurons
    ax.plot(df_neuron['x'] * 1e6, df_neuron['y'] * 1e6, 'o')
