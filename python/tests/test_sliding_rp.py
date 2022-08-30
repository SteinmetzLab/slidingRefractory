import numpy as np
import slidingRP
from pathlib import Path

TEST_DATA_PATH = Path(slidingRP.__file__).parent.parent.joinpath("test-data")
EXPECTED = {
    167: (0.7968085162939342, np.nan, np.nan, 0),  # FAIL
    274: (100.0, 0.5, 0.0012166666666666667, 2),   # PASS
    275: (53.88397133550098, 15.0, 0.00055, 99),  # FAIL
}

def generate_test_data():
    #  167: FAIL max conf = 0.80%, min cont = nan%, time = nan ms, n below 2 ms = 0
    #  274: PASS max conf = 100.00%, min cont = 0.5%, time = 1.22 ms, n below 2 ms = 2
    #  275: FAIL max conf = 53.88%, min cont = 15.0%, time = 0.55 ms, n below 2 ms = 99
    from brainbox.io.one import SpikeSortingLoader
    from one.api import ONE
    pid = 'ce397420-3cd2-4a55-8fd1-5e28321981f4'
    one = ONE()
    spikes, clusters, channels = SpikeSortingLoader(pid, one=one)
    sel_c = np.array([167, 274, 275])
    ispi = np.isin(spikes.clusters, sel_c)
    np.save(TEST_DATA_PATH.joinpath("spikes.times.npy"), spikes.times[ispi])
    np.save(TEST_DATA_PATH.joinpath("spikes.clusters.npy"), spikes.clusters[ispi])


def test_single_cluster():
    spikes_times = np.load(TEST_DATA_PATH.joinpath('spikes.times.npy'))
    spikes_clusters = np.load(TEST_DATA_PATH.joinpath('spikes.clusters.npy'))
    params = {'sampleRate': 30000, 'binSizeCorr': 1 / 30000}
    for clu in np.unique(spikes_clusters):
        sel = spikes_clusters == clu
        out = slidingRP.slidingRP(spikes_times[sel], params)
        assert EXPECTED[clu] == out[:4]


def test_multi_clusters():
    spikes_times = np.load(TEST_DATA_PATH.joinpath('spikes.times.npy'))
    spikes_clusters = np.load(TEST_DATA_PATH.joinpath('spikes.clusters.npy'))
    params = {'sampleRate': 30000, 'binSizeCorr': 1 / 30000}
    table = slidingRP.slidingRP_all(spikes_times, spikes_clusters, **params)
    for i, clu in enumerate(table['cidx']):
        assert EXPECTED[clu] == (table['maxConfidenceAt10Cont'][i],
                                 table['minContWith90Confidence'][i],
                                 table['timeOfLowestCont'][i],
                                 table['nSpikesBelow2'][i])


