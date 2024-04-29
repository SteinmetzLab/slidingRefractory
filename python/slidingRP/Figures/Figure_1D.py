'''
Plot Horowitz and IBL (VIS) RP distributiions

To get data, use:
slidingRefractory/python/slidingRP/elts/VIS_transform_acgs_into_rps.py
'''

import numpy as np
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from slidingRP.Figures.Fig_1_format import palette, linestyles

fig_name = "Fig_1D"
fig_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Figure")
data_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Data")

datasets = ['horowitz', 'ibl']
regions = ['VIS']

fr_threshold = 2  # Set threshold to pass units in analysis

##
# Load data and create dataframe
df_all = pd.DataFrame()

for dataset in datasets:
    for region in regions:
        file = data_path.joinpath(dataset).joinpath(f'estimatedRP_{region}.npy')
        estimatedRP_array = np.load(file, allow_pickle=True)

        file = data_path.joinpath(dataset).joinpath(f'mean_fr_{region}.npy')
        mean_fr_array = np.load(file, allow_pickle=True)

        data = {'dataset': [dataset] * len(mean_fr_array), 'region': [region] * len(mean_fr_array),
                'mean_fr': np.squeeze(mean_fr_array), 'estimatedRP': np.squeeze(estimatedRP_array)}
        df = pd.DataFrame.from_dict(data)

        df_all = pd.concat([df_all, df])
##
# Drop estimate RP that are NaN
index_nonnan = np.where(~np.isnan(df_all['estimatedRP']))[0]
df_all = df_all.iloc[index_nonnan]
# Drop mean_fr below set firing rate
index_fr = np.where(df_all['mean_fr'] >= fr_threshold)
df_all = df_all.iloc[index_fr]
df_all = df_all.reset_index()

##
linestyles = ["--", "-"]
palette = {"VIS": "#118391"}

# Plot distributions and median

figstyle = 'kde'

if figstyle == 'hist':
    arrow_len = 2
    head_length = 0.4
elif figstyle == 'kde':
    arrow_len = 0.065
    head_length = 0.01

fig, ax = plt.subplots()
for dataset, linestyle in zip(datasets, linestyles):
    df_dataset = df_all.loc[df_all['dataset'] == dataset]
    if figstyle == 'hist':
        sns.histplot(df_dataset, x="estimatedRP", hue='region',
                    ax=ax, linestyle=linestyle, legend=False, palette=palette,
                    element="step", fill=False, kde=False, stat='percent')
    elif figstyle == 'kde':
        sns.kdeplot(df_dataset, x="estimatedRP", hue='region',
                    ax=ax, linestyle=linestyle, legend=False, palette=palette)
# Set x min axis limit to 0
axlim = ax.get_xlim()
# ax.set_xlim(0, axlim[1])
ax.set_xlim(0, 7.3)
# Count N units per categories
count = df_all.groupby(['dataset', 'region'])['estimatedRP'].count()
medianRP = df_all.groupby(['dataset', 'region'])['estimatedRP'].median()
# Create labels
labels = list()
for dataset in datasets:
    for region in regions:
        labels.append(f'{dataset} - {region} - {count[dataset][region]} - '
                      f' RP:%.2f ms' % medianRP[dataset][region])

plt.legend(title='Datasets',
           loc='upper right', labels=labels)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))  # move legend outside
# Set arrow limits based on axis limits
aylim = ax.get_ylim()
arrow_min = aylim[1]+10/100*(aylim[1]-aylim[0])  # Add 10% of range

# Need to plot median as arrows after otherwise taken into legends
for dataset, linestyle in zip(datasets, linestyles):
    for region in regions:
        plt.arrow(x=medianRP[dataset][region], dx=0, y=arrow_min, dy=-arrow_len, linestyle=linestyle,
                  head_width=0.1, head_length=head_length, color=palette[region],
                  linewidth=1.5)
        # Trick to get the arrow head filled (replot on top with small arrow tail)
        plt.arrow(x=medianRP[dataset][region], dx=0, y=arrow_min-arrow_len, dy=-0.0001, linestyle="-",
                  head_width=0.1, head_length=head_length, color=palette[region],
                  linewidth=1.5)

fig.tight_layout()
fig.set_size_inches([9.61, 4.32])
# Save
plt.savefig(fig_path.joinpath(f'{fig_name}_{fr_threshold}frmin.pdf'))
