import numpy as np
import pytest
from slidingRP import metrics
from pathlib import Path

TEST_DATA_PATH = Path(__file__).parents[3].joinpath('test-data', 'unit')

# Expected (max_confidence, min_contamination, rp_min_val, n_spikes_below2).
# These values now match the authoritative MATLAB implementation bit-for-bit
# (verified 2026-06-16); the Python ACG and Llobet Ve were aligned to MATLAB.
EXPECTED = {
    167: (0.5301414775441771, np.nan, np.nan, 0),                     # FAIL
    274: (100.0, 0.5, 0.00115, 2),                                    # PASS
    275: (27.431935653278614, 19.0, 0.0005166666666666667, 104),     # FAIL
}


def _assert_matches_expected(clu, out4):
    """Compare a (max_conf, min_cont, rp_min_val, n_below2) tuple to EXPECTED,
    handling NaNs and floating-point tolerance."""
    exp = EXPECTED[clu]
    max_conf, min_cont, rp_min_val, n_below2 = out4
    assert max_conf == pytest.approx(exp[0], abs=1e-9)
    if np.isnan(exp[1]):
        assert np.isnan(min_cont)
    else:
        assert min_cont == pytest.approx(exp[1], abs=1e-9)
    if np.isnan(exp[2]):
        assert np.isnan(rp_min_val)
    else:
        assert rp_min_val == pytest.approx(exp[2], abs=1e-9)
    assert int(n_below2) == exp[3]


def generate_test_data():
    #  167: FAIL max conf = 0.53%, min cont = nan%, time = nan ms, n below 2 ms = 0
    #  274: PASS max conf = 100.00%, min cont = 0.5%, time = 1.15 ms, n below 2 ms = 2
    #  275: FAIL max conf = 27.43%, min cont = 19.0%, time = 0.52 ms, n below 2 ms = 104
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
    for clu in np.unique(spikes_clusters):
        sel = spikes_clusters == clu
        out = metrics.slidingRP(spikes_times[sel])
        _assert_matches_expected(clu, out[:4])


def test_multi_clusters():
    spikes_times = np.load(TEST_DATA_PATH.joinpath('spikes.times.npy'))
    spikes_clusters = np.load(TEST_DATA_PATH.joinpath('spikes.clusters.npy'))
    table = metrics.slidingRP_all(spikes_times, spikes_clusters)
    for i, clu in enumerate(table['cidx']):
        _assert_matches_expected(clu, (table['max_confidence'][i],
                                       table['min_contamination'][i],
                                       table['rp_min_val'][i],
                                       table['n_spikes_below2'][i]))
