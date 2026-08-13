"""Grouping spikes by cluster must not change any metric it feeds."""
import numpy as np

from slidingRP import metrics


def _fake_recording(n_clusters=12, seed=0):
    """Poisson-ish spikes for several clusters, interleaved in time as a sorter emits them."""
    rng = np.random.default_rng(seed)
    times, clusters = [], []
    for c in range(n_clusters):
        n = int(rng.integers(200, 2000))
        t = np.cumsum(rng.exponential(1.0 / rng.uniform(2, 20), n))
        times.append(t)
        clusters.append(np.full(n, c * 3))          # non-contiguous cluster ids
    times = np.concatenate(times)
    clusters = np.concatenate(clusters)
    order = np.argsort(times)                        # recordings arrive time-ordered
    return times[order], clusters[order]


def test_group_matches_boolean_masking():
    """The grouping returns exactly what the mask-per-cluster comprehension returned."""
    times, clusters = _fake_recording()
    cids, sts = metrics._group_by_cluster(times, clusters)

    expected_cids = np.unique(clusters)
    assert np.array_equal(cids, expected_cids)
    for cid, st in zip(cids, sts):
        assert np.array_equal(st, np.sort(times[clusters == cid]))


def test_group_handles_unsorted_and_single_cluster():
    times, clusters = _fake_recording(n_clusters=1)
    rng = np.random.default_rng(1)
    shuffle = rng.permutation(len(times))            # spikes not in time order
    cids, sts = metrics._group_by_cluster(times[shuffle], clusters[shuffle])
    assert len(cids) == 1
    assert np.array_equal(sts[0], np.sort(times))
    assert np.all(np.diff(sts[0]) >= 0)


def test_slidingRP_all_unchanged_by_grouping():
    """Every per-cluster metric matches computing it one cluster at a time."""
    times, clusters = _fake_recording()
    out = metrics.slidingRP_all(times, clusters, conf_thresh=90, cont_thresh=10)

    for i, cid in enumerate(out['cidx']):
        st = np.sort(times[clusters == cid])
        (max_conf, min_cont, rp_min_val, n_below2,
         firing_rate, pass_cont, pass_forced) = metrics.slidingRP(st)
        assert out['max_confidence'][i] == max_conf or (
            np.isnan(out['max_confidence'][i]) and np.isnan(max_conf))
        assert out['min_contamination'][i] == min_cont or (
            np.isnan(out['min_contamination'][i]) and np.isnan(min_cont))
        assert out['n_spikes_below2'][i] == n_below2
        assert out['firing_rate'][i] == firing_rate
        assert out['value'][i] == int(pass_cont)
        assert out['value_forced'][i] == int(pass_forced)
