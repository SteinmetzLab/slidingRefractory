from pathlib import Path
from slidingRP.metrics import compute_rf
import numpy as np
from iblatlas.atlas import BrainRegions
import pandas as pd

pic_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Picture")
data_path = Path("/Users/gaelle/Desktop/Reports/RefractoryPeriod/Data")

datasets = ['horowitz', 'ibl']
regions = ['VIS']

mapping = 'Beryl'
br = BrainRegions()

# --- Datasets for Horowitz were pre-downloaded from G-drive
# https://drive.google.com/drive/folders/1txYCRlk-2YV4IijV09Z0HtKE_gfeRfWR

##
# --- For all datasets, we used a bin size of 1 / 30000 seconds to make the ACG
bin_size_secs = 1 / 30_000
for dataset in datasets:
    for region in regions:
        print(f'{dataset} - {region} : in process')

        pqt_path = data_path.joinpath(dataset).joinpath("clusters.pqt")
        data_dir = pqt_path.parent

        if dataset == "horowitz":
            file = data_path.joinpath(dataset).joinpath(f'acgs.npy')
            acgs = np.load(file, allow_pickle=True)
        elif dataset == 'ibl':
            # Load files (ACG and DF)
            file = data_path.joinpath(dataset).joinpath(f'acgs.npy')
            acgs_ibl = np.load(file, allow_pickle=True)

            file = data_path.joinpath('ibl').joinpath(f'clusters_df.pqt')
            df_clusters_good = pd.read_parquet(file)

            # Remap (make sure atlas ID are int prior to remapping)
            df_clusters_good[mapping + '_id'] = br.remap(df_clusters_good['atlas_id'].values.astype(int),
                                                         source_map='Allen', target_map=mapping)
            df_clusters_good[mapping + '_acronym'] = br.id2acronym(df_clusters_good[mapping + '_id'])

            # Filter for region
            indx_reg = np.where(df_clusters_good[mapping + '_acronym'].str.startswith(region))[0]
            df_clusters_good = df_clusters_good.iloc[indx_reg]
            acgs = acgs_ibl[indx_reg, :]
            print(f'Regions used: {df_clusters_good.Beryl_acronym.unique()}')
            # ['VISa', 'VISpm', 'VISam', 'VISp', 'VISrl', 'VISpl', 'VISC',
            #        'VISpor', 'VISli', 'VISl', 'VISal']

        # Compute RF from ACG for single unit
        estimatedRP_array = np.empty((acgs.shape[0], 1))
        mean_fr_array = np.empty((acgs.shape[0], 1))
        for i_unit in range(0, acgs.shape[0]):
            acg = acgs[i_unit, :]
            mean_fr = sum(acg)/len(acg)
            estimatedRP, estimateIdx, xSigmoid, ySigmoid = \
                compute_rf(acg, bin_size_secs=bin_size_secs)

            # Save into array
            estimatedRP_array[i_unit] = estimatedRP
            mean_fr_array[i_unit] = mean_fr

        # Save file
        file = data_path.joinpath(dataset).joinpath(f'estimatedRP_{region}.npy')
        with open(file, 'wb') as f:
            np.save(file, estimatedRP_array, allow_pickle=True)

        file = data_path.joinpath(dataset).joinpath(f'mean_fr_{region}.npy')
        with open(file, 'wb') as f:
            np.save(file, mean_fr_array, allow_pickle=True)
