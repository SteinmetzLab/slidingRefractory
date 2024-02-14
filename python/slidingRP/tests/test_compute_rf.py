from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from one.api import ONE
import ephys_atlas.data
from ephys_atlas.plots import plt_unit_acg
from slidingRP.metrics import compute_rf, plotSigmoid
import scipy
from scipy import signal

one = ONE(base_url="https://alyx.internationalbrainlab.org", mode='local')

LABEL = '2024_W04'

LOCAL_DATA_PATH = Path('/Users/gaelle/Documents/Work/EphysAtlas/')
savepath = LOCAL_DATA_PATH.joinpath('Autocorrelograms').joinpath(LABEL)
if not savepath.parent.exists():
    savepath.parent.mkdir()
if not savepath.exists():
    savepath.mkdir()

df_raw_features, df_clusters, df_channels, df_probes = ephys_atlas.data.download_tables(
    label=LABEL, local_path=LOCAL_DATA_PATH, one=one, extended=True)

# corr_ts = ephys_atlas.data.read_correlogram(LOCAL_DATA_PATH.joinpath(
#     LABEL, 'clusters_correlograms_time_scale.bin'), df_clusters.shape[0])
corr_rf = ephys_atlas.data.read_correlogram(LOCAL_DATA_PATH.joinpath(
    LABEL, 'clusters_correlograms_refractory_period.bin'), df_clusters.shape[0])
bin_size_secs = 1 / 30_000

# Compute timebins
x = np.arange(corr_rf.shape[1])
timeBins = x.dot(bin_size_secs)
##

# Plot for 5 random units
idxplt = [2415, 80121, 35243, 328904, 369420, 70160]
estimated_RP_expected = [1.7, 2.8, 1.27, 0.83, 4.2, 1.93]
estimatedRP_computed = list()
# Plot to debug
is_plot = True

for i_cell in idxplt:
    acg = corr_rf[i_cell, :]
    # Compute
    estimatedRP, estimateIdx, xSigmoid, ySigmoid = \
        compute_rf(acg,
                   bin_size_secs=bin_size_secs,
                   timeBins=timeBins,
                   fr_percentage=10 / 100)

    estimatedRP_computed.append(estimatedRP)

    if is_plot:
        # Plot
        fig, axs = plt.subplots(1, 2)
        plt_unit_acg(i_cell, corr_rf, df_clusters, bin_size_secs, ax=axs[0], fig=fig)

        # TODO the function does not return the filtered traced, re-created here
        # med_filt = scipy.ndimage.median_filter(acg, size=25)
        # peaks = scipy.signal.find_peaks(med_filt)
        #
        # axs[0].plot(med_filt, 'g')
        # axs[0].plot(peaks[0], med_filt[peaks[0]], 'kx')
        axs[0].plot(np.array([estimatedRP, estimatedRP]) / (1000 * bin_size_secs),
                    [acg.min(), acg.max()], 'k')
        plotSigmoid(axs[1], acg, timeBins, ySigmoid, estimateIdx, estimatedRP)
        fig.set_size_inches([8.61, 3.94])
        plt.savefig(savepath.joinpath(f'{i_cell}_acg_TEST.png'))

    # Test
np.testing.assert_almost_equal(estimatedRP_computed, estimated_RP_expected, decimal=1)

