'''
Figure 1B
Plot ACG with sigmoiod fit for 1 example neuron
'''
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from slidingRP.metrics import plotSigmoid, compute_rf, compute_timebins
from slidingRP.Figures.Fig_1_format import palette

data_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Data")
fig_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Figure")

# Load files (ACG and DF)
file = data_path.joinpath('ibl').joinpath(f'acgs.npy')
acgs_ibl = np.load(file, allow_pickle=True)

file = data_path.joinpath('ibl').joinpath(f'clusters_df.pqt')
df_clusters_good = pd.read_parquet(file)

##
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
bin_size_secs = 1 / 30_000

# Pick one neuron at random
i_neuron = 0

acg = acgs_ibl[0, :]
df_neuron = df_clusters_good.iloc[i_neuron]

# Compute RF
estimatedRP, estimatedIdx, xSigmoid, ySigmoid = \
    compute_rf(acg, bin_size_secs=bin_size_secs)

# Plot ACG of one neuron
fig_sig = plt.figure()
ax = plt.gca()
timeBins = compute_timebins(acg, bin_size_secs)
plotSigmoid(ax, acg, timeBins, ySigmoid, estimatedIdx, estimatedRP)
ax.set_title('Estimated RP:%.2f ms' % estimatedRP,
             color=palette[df_neuron['Cosmos_acronym']])
