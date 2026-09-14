"""Shared helpers for the J Neurophysiol resubmission analyses.

Builds and consumes the per-unit ACG table described in
``jNeurophysResubmission/01_fig1_rp_durations/plan.md`` (stage A): one row per
sorted unit with identifiers, region, spike count, recording duration and the
0-10 ms autocorrelogram at 1/30000 s resolution (300 bins).

Everything downstream (RP estimate, Sliding RP, Hill-Llobet at fixed RP,
tau_pass0, power flags) is computable from that row, so the raw spike files can
be deleted after the pass.

Nothing here changes the behaviour of the published package: the metric
functions below are re-expressions of ``slidingRP.metrics.slidingRP``'s fast
path that accept a precomputed ACG instead of spike times, and they are checked
against the package on real spike trains by ``test_acg_table.py``.
"""
from __future__ import annotations

import numpy as np

SAMPLE_RATE = 30000.0
BIN_SIZE = 1.0 / SAMPLE_RATE       # ACG bin width (s) = one sample at 30 kHz
N_BINS = 300                       # 0 to 10 ms
RP_REJECT = 0.0005                 # tau_min (s)

# tau_r vectors, matching slidingRP.metrics.slidingRP exactly:
#   rpEdges = arange(0, 10 ms, binSize);  rp = rpEdges + binSize/2  (bin centres)
#   refDur  = rp + binSize/2                                        (right edge)
RP_EDGES = np.arange(0, 10 / 1000, BIN_SIZE)
RP_CENTERS = RP_EDGES + BIN_SIZE / 2
REF_DUR = RP_CENTERS + BIN_SIZE / 2
TEST_TIMES = RP_CENTERS > RP_REJECT


def computeACG_samples(spike_samples: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """Exact integer-lag autocorrelogram from spike times in samples.

    Bin k counts ordered pairs (j < i) with ``s_i - s_j == k`` samples, for
    1 <= k <= n_bins - 1 (bin 0 holds exact coincidences, which are excluded,
    matching histdiff / computeACG).

    This is the float-free version of ``slidingRP.metrics.computeACG``: for
    sample-derived data the two agree only up to the seconds<->samples rounding
    (see the repo's code-review note on ACG binning), and this one is exact.
    """
    s = np.sort(np.asarray(spike_samples, dtype=np.int64))
    counts = np.zeros(n_bins, dtype=np.int64)
    n = s.size
    shift = 1
    while shift < n:
        d = s[shift:] - s[:-shift]
        m = (d > 0) & (d < n_bins)
        if not np.any(m):
            break
        counts += np.bincount(d[m], minlength=n_bins)[:n_bins]
        shift += 1
    counts[0] = 0
    return counts


def computeViol(obs_viol, spike_count, ref_dur, contamination_prop, rec_dur):
    """Poisson confidence that contamination is below ``contamination_prop``.

    Llobet et al. (2022) expected violations; identical to
    ``slidingRP.metrics.computeViol`` and ``matlab/computeViol.m``.
    """
    from scipy import stats

    Nc = spike_count * contamination_prop
    Nb = spike_count * (1 - contamination_prop)
    expected = 2 * ref_dur / rec_dur * Nc * (Nb + (Nc - 1) / 2)
    return 1 - stats.poisson.cdf(obs_viol, expected), expected


def slidingRP_from_acg(nACG, spike_count, rec_dur, cont_thresh=10.0,
                       conf_thresh=90.0, rp_reject=RP_REJECT):
    """Sliding RP metric from a precomputed ACG (fast analytical path).

    Mirrors ``slidingRP.metrics.slidingRP`` step for step, but takes the ACG
    (300 bins at 1/30000 s) rather than spike times, so the table can be built
    once and re-analysed many times.

    Returns a dict with max_conf, min_cont, rp_min_val, n_viol_short, passes,
    tau_first_pass and tau_pass0.
    """
    from scipy import stats

    nACG = np.asarray(nACG, dtype=np.float64)
    obs_viol = np.cumsum(nACG)

    # Honour the rp_reject argument rather than the module default, so the
    # tau_min sensitivity sweep actually varies tau_min.
    test_times = RP_CENTERS > rp_reject
    if not np.any(test_times):
        return dict(max_conf=0.0, min_cont=np.nan, rp_min_val=np.nan,
                    n_viol_short=int(np.sum(nACG[:int(np.argmax(RP_CENTERS > 0.002)) + 1])),
                    passes=False, tau_first_pass=np.nan,
                    tau_pass0=tau_pass0(spike_count, rec_dur, cont_thresh, conf_thresh))

    conf_at_thresh = 100 * computeViol(obs_viol, spike_count, REF_DUR,
                                       cont_thresh / 100, rec_dur)[0]
    max_conf = float(np.max(conf_at_thresh[test_times]))
    passes = bool(max_conf >= conf_thresh)

    # Shortest tau_r above tau_min that reaches the confidence threshold.
    ok = test_times & (conf_at_thresh >= conf_thresh)
    tau_first_pass = float(RP_CENTERS[np.argmax(ok)]) if np.any(ok) else np.nan

    # Minimum confirmable contamination, analytical (matches
    # compute_min_contamination / computeMinContamination.m).
    lam = stats.chi2.ppf(conf_thresh / 100, 2 * (obs_viol + 1)) / 2
    disc = (spike_count - 0.5) ** 2 - lam * rec_dur / REF_DUR
    cmin = np.full(REF_DUR.shape, np.nan)
    good = disc >= 0
    cmin[good] = ((spike_count - 0.5) - np.sqrt(disc[good])) / spike_count * 100
    cmin_test = cmin[test_times]
    if np.all(np.isnan(cmin_test)):
        min_cont, rp_min_val = np.nan, np.nan
    else:
        i = int(np.nanargmin(cmin_test))
        min_cont = float(cmin_test[i])
        rp_min_val = float(RP_CENTERS[test_times][i])
        if min_cont > 35:
            min_cont, rp_min_val = np.nan, np.nan

    n_viol_short = int(np.sum(nACG[:int(np.argmax(RP_CENTERS > 0.002)) + 1]))

    return dict(max_conf=max_conf, min_cont=min_cont, rp_min_val=rp_min_val,
                n_viol_short=n_viol_short, passes=passes,
                tau_first_pass=tau_first_pass,
                tau_pass0=tau_pass0(spike_count, rec_dur, cont_thresh, conf_thresh))


def tau_pass0(spike_count, rec_dur, cont_thresh=10.0, conf_thresh=90.0):
    """Shortest violation-free window that would pass, in seconds.

    A unit with zero observed violations up to tau passes iff
    ``1 - exp(-Ve(tau)) >= gamma``; inverting the Llobet Ve gives

        tau_pass0 = -ln(1 - gamma) * D / (2 * Nc * (Nb + (Nc - 1)/2))

    Units whose tau_pass0 exceeds the tested window (10 ms) cannot pass
    whatever their ACG looks like. See 06_underpowered_outcome/plan.md.
    """
    C = cont_thresh / 100.0
    Nc = spike_count * C
    Nb = spike_count * (1 - C)
    denom = 2 * Nc * (Nb + (Nc - 1) / 2)
    if denom <= 0:
        return np.inf
    return -np.log(1 - conf_thresh / 100.0) * rec_dur / denom


def min_passing_fr(rec_dur, tau, cont_thresh=10.0, conf_thresh=90.0):
    """Minimum firing rate (spikes/s) for a violation-free unit to pass at tau.

    Exact inverse of the Llobet expected-violation count: solves
    ``Ve(tau) = -ln(1 - gamma)`` for N, a quadratic in N, then FR = N / D.
    """
    C = cont_thresh / 100.0
    k = -np.log(1 - conf_thresh / 100.0) * rec_dur / (2 * tau)
    # Ve = 2*tau/D * C*N * ((1-C)*N + (C*N - 1)/2)  ->  a*N^2 + b*N - k = 0
    a = C * (1 - C) + C * C / 2
    b = -C / 2
    N = (-b + np.sqrt(b * b + 4 * a * k)) / (2 * a)
    return N / rec_dur


def hill_llobet_from_acg(nACG, spike_count, rec_dur, rp_dur, cont_thresh=10.0):
    """Hill-Llobet fixed-RP point estimate from a precomputed ACG.

    Mirrors ``matlab/RPmetric_Classic.m`` (metricType 'Llobet'), including its
    inclusive bin selection ``sum(nACG(1:find(rp > RPdur, 1)))``.
    """
    idx = int(np.argmax(RP_CENTERS > rp_dur))
    obs_viol = float(np.sum(np.asarray(nACG)[:idx + 1]))
    Nc = spike_count * cont_thresh / 100
    Nb = spike_count * (1 - cont_thresh / 100)
    expected = 2 * rp_dur / rec_dur * Nc * (Nb + (Nc - 1) / 2)
    with np.errstate(invalid="ignore"):
        est = 1 - np.sqrt(1 - obs_viol * rec_dur / (spike_count ** 2 * rp_dur))
    return bool(obs_viol <= expected), float(est), obs_viol


# --- table schema -----------------------------------------------------------

TABLE_COLUMNS = [
    "dataset", "species", "animal", "session", "insertion", "cluster_id",
    "acronym", "beryl", "cosmos", "sorter_label", "n_spikes", "rec_dur_s",
    "firing_rate", "acg_source",
]
