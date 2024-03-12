from slidingRP.data_access.allen import AllenDataAccess
from slidingRP.data_access.steinmetz import SteinmetzDataAccess
from pathlib import Path
from slidingRP.metrics import compute_rf
import matplotlib.pyplot as plt
import numpy as np
from one.api import ONE
import ephys_atlas.data
from iblatlas.atlas import BrainRegions

pic_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Picture")
data_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Data")
is_plot = False

datasets = ['steinmetz', 'allen', 'ibl']
regions = ['Isocortex', 'TH', 'HPF']

# --- Datasets for Steinmetz and Allen were pre-downloaded from G-drive
# https://github.com/orgs/int-brain-lab/projects/11/views/6?sliceBy%5Bvalue%5D=GaelleChapuis&filterQuery=-status%3ANew%2C%22On+Hold%22%2C%22Info+needed%22+sprint%3A%40current+sprint%3A%22Sprint+6%22&pane=issue&itemId=56008576

# --- Load IBL dataset once (filter for good units and region later)
one = ONE(base_url="https://alyx.internationalbrainlab.org", mode='local')
LABEL = '2024_W04'
LOCAL_DATA_PATH = Path('/Users/gaelle/Documents/Work/EphysAtlas/')
df_raw_features, df_clusters, df_channels, df_probes = ephys_atlas.data.download_tables(
    label=LABEL, local_path=LOCAL_DATA_PATH, one=one, extended=True)
corr_rf = ephys_atlas.data.read_correlogram(LOCAL_DATA_PATH.joinpath(
    LABEL, 'clusters_correlograms_refractory_period.bin'), df_clusters.shape[0])

# Remap acronyms to Cosmos
br = BrainRegions()
df_clusters['cosmos_id'] = br.remap(df_clusters['atlas_id'], source_map='Allen', target_map='Cosmos')
df_clusters['cosmos_acronym'] = br.id2acronym(df_clusters['cosmos_id'])

# Keep only good units
df_clusters = df_clusters.drop(columns=['cluster_id']).reset_index()
indx_good = df_clusters.index[df_clusters['label'] == 1]
df_clusters_good = df_clusters.iloc[indx_good]
acgs_ibl = corr_rf[indx_good, :]

##
# --- For all datasets, we used a bin size of 1 / 30000 seconds to make the ACG
bin_size_secs = 1 / 30_000
for dataset in datasets:
    for region in regions:
        print(f'{dataset} - {region} : in process')

        pqt_path = data_path.joinpath(dataset).joinpath("clusters.pqt")
        data_dir = pqt_path.parent

        if dataset == "allen":
            ADA = AllenDataAccess(pqt_path, data_dir)
            # Load acgs for a given brain region
            acgs, df = ADA.load_acgs_region(region)
        elif dataset == "steinmetz":
            SDA = SteinmetzDataAccess(pqt_path, data_dir)
            acgs, df = SDA.load_acgs_region(region)
        elif dataset == 'ibl':
            # Filter for region
            indx_reg = df_clusters_good['Cosmos_acronym'] == region
            acgs = acgs_ibl[indx_reg, :]

        # Compute RF from ACG for single unit
        estimatedRP_array = np.empty((acgs.shape[0], 1))
        for i_unit in range(0, acgs.shape[0]):
            acg = acgs[:, i_unit]
            estimatedRP, estimateIdx, xSigmoid, ySigmoid = \
                compute_rf(acg, bin_size_secs=bin_size_secs)

            # Save into array
            estimatedRP_array[i_unit] = estimatedRP

            # Plot auto-corr
            if is_plot and not np.isnan(estimatedRP):
                x = np.arange(len(acg))
                fig, ax = plt.subplots()
                ax.bar(x, list(acg))
                xstep = 35
                ax.set_xticks(x[0:-1:xstep])
                # Plot in millisecond
                ax.set_xticklabels(np.around(x[0:-1:xstep].dot(bin_size_secs) * 1e3, decimals=2))
                ax.set_title(f'{dataset} - {region} - unit {i_unit}')
                ax.plot(np.array([estimatedRP, estimatedRP]) / (1000 * bin_size_secs), [acg.min(), acg.max()], 'k')
                plt.savefig(pic_path.joinpath(dataset).joinpath(f'{dataset}_{region}_{i_unit}.png'))
                plt.close(fig)

        # Save file
        file = data_path.joinpath(dataset).joinpath(f'{region}.npy')
        with open(file, 'wb') as f:
            np.save(file, estimatedRP_array, allow_pickle=True)
