"""
Sliding Refractory Period (RP) Quality Metric
==============================================

Author: Ryan A. Ressmeyer (with strong contribution from Claude Code)

This module implements a statistical quality metric for evaluating whether spike-sorted
neural units have contaminated refractory periods. A contaminated refractory period
indicates that spikes from other neurons or noise are being mis-attributed to the unit. 
These methods were adapted from the "slidingRP" metric developed by Nick Steinmetz's lab.

Algorithm
---------

**Step 1 — Spike Train Model**

The recorded spike train is modelled as a mixture of two independent processes:

- Base neuron spikes: N_b = (1 - C) * N_s spikes, which obey a true refractory period
  and therefore produce no close spike pairs among themselves.
- Contaminating spikes: N_c = C * N_s spikes, which are Poisson with no refractory
  period, where C is the contamination proportion (0 = clean, 1 = fully contaminated).

N_s is the total spike count and F_r is the total firing rate (Hz), so the recording
duration is approximately D = N_s / F_r.

**Step 2 — Observed Refractory Period Violations**

For each tested refractory period duration τ_r, the observed violation count V_o(τ_r)
is the cumulative sum of the one-sided autocorrelogram (ACG) from ref_acg_t_start up
to τ_r:

    V_o(τ_r) = Σ_{i=ref_acg_t_start}^{τ_r} n_ACG(i)

The ACG is one-sided (only positive lags), so each close spike pair is counted once.
ref_acg_t_start defaults to 0.25 ms to avoid bias from Kilosort's near-duplicate spike
removal (which suppresses spikes within ~0.25 ms of each other).  The effective tested
window is the adjusted refractory period τ_adj = τ_r − ref_acg_t_start.

**Step 3 — Expected Violations Under a Given Contamination Level**

For a one-sided ACG, two classes of spike pairs within [0, τ_adj] produce violations:

  - Contam–base pairs (either time ordering): 2 * N_c * N_b * τ_adj / D
  - Contam–contam pairs (unordered):         N_c * (N_c − 1) * τ_adj / D

Base–base pairs contribute zero violations because the base neuron has its own
refractory period.  Combining and substituting N_c = C*N_s, N_b = (1−C)*N_s:

    V_e = [2 * N_c * N_b + N_c * (N_c − 1)] * τ_adj / D
        ≈ C * (2 − C) * F_r * τ_adj * N_s

The factor (2 − C) reduces to ≈ 2 for small C and captures the contam–contam
self-interactions. This matches the formula in Llobet et al. (2022) and Steinmetz et al.

**Step 4 — Poisson Likelihood**

V_o is modelled as Poisson(V_e).  The likelihood of observing V_o or fewer violations
given contamination C is:

    P(X ≤ V_o | λ = V_e) = poisson.cdf(V_o, V_e)

If the true contamination were C, and we observe very few violations, this probability
will be small — the observations are unlikely under the hypothesis "contamination = C".

**Step 5 — Minimum Rejected Contamination**

For each cluster and each τ_r, find C_min: the smallest contamination for which the
Poisson likelihood falls below (1 − confidence), i.e. we can reject the hypothesis that
contamination is as large as C_min at the given confidence level.

Two equivalent methods are provided:

a) Binary search (compute_min_contam_props): bisects [0, max_contam_prop] to find C
   where poisson.cdf(V_o, V_e(C)) = 1 − confidence.

b) Analytical solution (compute_min_contam_props_analytical): uses the identity
       poisson.cdf(r, λ) = 1 − chi2.cdf(2λ, 2*(r+1))
   to find the critical Poisson rate directly:
       λ_crit = chi2.ppf(confidence, 2*(V_o+1)) / 2
   then inverts the quadratic  C * (2 − C) = k  (where k = λ_crit / (F_r * τ_adj * N_s))
   via the closed-form root:
       C = 1 − sqrt(1 − k)

Functions
---------
calc_ccgs                           : Efficient pairwise cross-correlogram computation.
calc_local_firing_rate              : Fast local (ACG-weighted) firing rate via searchsorted.
calc_rp_violations                  : Fast cumulative violation counts via k-th order ISIs.
refractory_violation_likelihood     : Poisson CDF likelihood for a contamination level.
compute_min_contam_props            : Binary search for minimum rejected contamination.
compute_min_contam_props_analytical : Analytical (exact) solution for the same quantity.
compute_rvl_tensor                  : Full likelihood tensor over a contam × RP grid.
plot_min_contam_prop                : Plot the minimum contamination curve vs. RP.
plot_rvl                            : Plot the full RVL likelihood heatmap.
"""
import numpy as np
from tqdm import tqdm
from scipy.stats import poisson, chi2
import matplotlib.pyplot as plt
from typing import Optional

def ensure_ndarray(x, dtype=None):
    """
    Ensures that the input is a numpy.ndarray. If it is a tensor, it is converted to a numpy array.

    Parameters:
    ----------
    x : numpy.ndarray, torch.Tensor, int, float, list, or tuple
        The input array or tensor.

    Returns:
    -------
    numpy.ndarray
        The input converted to a numpy array.
    """
    if isinstance(x, int) or isinstance(x, float):
        x = [x]
    if isinstance(x, list) or isinstance(x, tuple):
        x = np.array(x)
    if dtype is not None:
        x = x.astype(dtype)
    return x


def calc_ccgs(
    spike_times: np.ndarray,
    bin_edges: np.ndarray,
    spike_clusters: Optional[np.ndarray] = None,
    uids: Optional[np.ndarray] = None,
    progress: bool = False,
) -> np.ndarray:
    """
    Compute all pairwise cross-correlograms among a set of spike clusters.

    This function uses a highly efficient, vectorized algorithm that iterates
    by "shift" in the spike array rather than by individual spike pairs. It
    skips counting spike pairs with a time lag of zero, meaning the central
    bin of an autocorrelogram will not equal the total spike count.

    Originally inspired by phylib correlogram function.
    https://github.com/cortex-lab/phylib/blob/master/phylib/stats/ccg.py

    Parameters
    ----------
    spike_times : np.ndarray (n_spikes,)
        A 1D array of spike times in seconds. For best performance, these
        should be pre-sorted, though the function will sort them if needed.
    bin_edges : np.ndarray (n_bins + 1,)
        A 1D array defining the edges of the time lag bins, in seconds.
        Example: `np.linspace(-0.05, 0.05, 101)` creates 100 bins for a
        +/- 50 ms window.
    spike_clusters : np.ndarray (n_spikes,), optional
        A 1D integer array mapping each spike to a cluster ID. If None,
        all spikes are treated as belonging to a single cluster.
    uids : np.ndarray (n_clusters,), optional
        An ordered 1D array of the unique cluster IDs to include in the
        output. The order of this array determines the axes of the output
        matrix. Spikes belonging to clusters not in `uids` will be ignored.
        If None, all unique clusters from `spike_clusters` will be used,
        sorted in ascending order.
    progress : bool, default=False
        If True, display a `tqdm` progress bar during computation.

    Returns
    -------
    correlograms : np.ndarray
        A 3D array of shape `(n_clusters, n_clusters, n_bins)` containing
        integer counts. `correlograms[i, j, k]` is the number of spikes from
        cluster `uids[j]` that occurred in time bin `k` relative to a spike
        from cluster `uids[i]`.
    """
    # --- 1. Input Validation and Preparation ---
    spike_times = np.asarray(spike_times, dtype=np.float64)
    bin_edges = np.asarray(bin_edges, dtype=np.float64)
    assert spike_times.ndim == 1, "spike_times must be a 1D array."
    assert bin_edges.ndim == 1, "bin_edges must be a 1D array."
    assert np.all(np.diff(bin_edges) > 0), "Bin edges must be monotonically increasing."

    if spike_clusters is None:
        spike_clusters = np.zeros(len(spike_times), dtype=np.int32)
    else:
        spike_clusters = np.asarray(spike_clusters, dtype=np.int32)

    assert len(spike_times) == len(spike_clusters), \
        "Spike times and spike clusters must have the same length."

    # The algorithm requires sorted spike times
    if not np.all(np.diff(spike_times) >= 0):
        print("Spike times are not sorted, sorting...")
        sort_inds = np.argsort(spike_times)
        spike_times = spike_times[sort_inds]
        spike_clusters = spike_clusters[sort_inds]

    # --- 2. Cluster ID Filtering and Mapping ---
    if uids is None:
        uids = np.unique(spike_clusters)
        needs_reordering = False
    else:
        uids = np.asarray(uids, dtype=np.int32)
        needs_reordering = True

    # Filter spikes to only include those in the `uids` list.
    membership_mask = np.isin(spike_clusters, uids)
    spike_times = spike_times[membership_mask]
    spike_clusters = spike_clusters[membership_mask]

    n_spikes = len(spike_times)
    n_clusters = len(uids)

    # Map raw cluster IDs to 0-based indices using the fast, sorted method
    if needs_reordering:
        sort_indices = np.argsort(uids)
        sorted_uids = uids[sort_indices]
        unsort_indices = np.argsort(sort_indices)
    else:
        sorted_uids = uids  # Already sorted and unique

    spike_inds = np.searchsorted(sorted_uids, spike_clusters)

    # --- 3. Binning Setup ---
    n_bins = len(bin_edges) - 1
    ccgs = np.zeros((n_clusters, n_clusters, n_bins), dtype=np.int32)
    min_bin, max_bin = bin_edges[0], bin_edges[-1]

    # Optimization: Use a faster calculation for linearly spaced bins
    bin_diffs = np.diff(bin_edges)
    if np.allclose(bin_diffs, bin_diffs[0]):
        bin_width = bin_diffs[0]
        def digitize(x):
            bins = np.asarray((x - min_bin) / bin_width, dtype=np.int32)
            return np.clip(bins, 0, n_bins - 1)
    else:
        def digitize(x):
            return np.digitize(x, bin_edges) - 1

    # --- 4. Main Correlogram Loop ---
    shift = 1
    pos_mask = np.ones(n_spikes, dtype=bool)
    neg_mask = np.ones(n_spikes, dtype=bool)
    pbar = tqdm(total=1.0, desc="Calculating CCGs", disable=not progress)

    while True:
        # --- Positive Lags (Ref: i, Target: i + shift) ---
        pos_mask[-shift:] = False
        active_pos_indices = np.where(pos_mask[:-shift])[0]
        has_pos = len(active_pos_indices) > 0

        if has_pos:
            pos_dts = spike_times[active_pos_indices + shift] - spike_times[active_pos_indices]
            valid_pos = (min_bin < pos_dts) & (pos_dts < max_bin)
            if np.any(valid_pos):
                valid_indices = active_pos_indices[valid_pos]
                pos_i, pos_j = spike_inds[valid_indices], spike_inds[valid_indices + shift]
                pos_bins = digitize(pos_dts[valid_pos])
                ravel_inds = np.ravel_multi_index((pos_i, pos_j, pos_bins), ccgs.shape)
                counts = np.bincount(ravel_inds, minlength=ccgs.size)
                ccgs += counts.reshape(ccgs.shape)
            pos_mask[active_pos_indices] = pos_dts < max_bin

        # --- Negative Lags (Ref: i + shift, Target: i) ---
        neg_mask[:shift] = False
        active_neg_indices = np.where(neg_mask[shift:])[0]
        has_neg = len(active_neg_indices) > 0

        if has_neg:
            neg_dts = spike_times[active_neg_indices] - spike_times[active_neg_indices + shift]
            valid_neg = (min_bin < neg_dts) & (neg_dts < max_bin)
            if np.any(valid_neg):
                valid_indices = active_neg_indices[valid_neg]
                neg_i, neg_j = spike_inds[valid_indices + shift], spike_inds[valid_indices]
                neg_bins = digitize(neg_dts[valid_neg])
                ravel_inds = np.ravel_multi_index((neg_i, neg_j, neg_bins), ccgs.shape)
                counts = np.bincount(ravel_inds, minlength=ccgs.size)
                ccgs += counts.reshape(ccgs.shape)
            neg_mask[active_neg_indices + shift] = neg_dts > min_bin

        # --- Loop Control and Progress Update ---
        if not has_pos and not has_neg:
            pbar.n = 1.0
            pbar.set_description("Calculating CCGs: Done")
            break

        if progress:
            progress_val = 1.0 - (np.sum(pos_mask) + np.sum(neg_mask)) / (2 * n_spikes)
            pbar.n = np.round(progress_val, 3)
            pbar.set_description(f"Calculating CCGs: Shift {shift}")

        shift += 1

    pbar.close()

    # Reorder axes to match original `uids` order if necessary
    if needs_reordering:
        ccgs = ccgs[unsort_indices][:, unsort_indices]

    return ccgs

def refractory_violation_likelihood(
        n_violations, 
        contam_prop,
        refractory_period,
        firing_rate, 
        n_spikes, 
        ):
    '''
    Calculate the likelihood of an observed number of refractory period violations under a poisson 
    model of refractory violations. likelihood = P(X <= N_v | R_c, T_ref, F_r, N_s), where X is a 
    poisson random variable with rate R_c * T_ref * F_r * N_s and N_v is the observed number of
    refractory period violations in the cluster. R_c is a specified contamination rate,
    T_ref is the refractory period, F_r is the firing rate of the cluster, and N_s is the number of
    spikes in the cluster.

    Parameters
    ----------
    n_violations : array_like
        the observed number of violations
    contam_prop : array_like
        the contamination proportion to test (as a proportion of the firing rate)
    refractory_period : array_like
        the refractory period in seconds
    firing_rate : array_like
        the firing rate of the cluster in Hz
    n_spikes : array_like
        the number of spikes in the cluster

    Returns
    -------
    likelihood : float
        the likelihood of the observing the number of violations or less if the cluster was contaminated at the rate specified

    '''
    # rate of contaminated spikes per second
    contamination_firing_rate = firing_rate * contam_prop

    # expected number of violations in the one-sided ACG:
    # = (2*N_c*N_b + N_c*(N_c-1)) * tau / D = C*(2-C) * F_r * tau * N_s
    expected_violations = contamination_firing_rate * (2 - contam_prop) * refractory_period * n_spikes

    # likelihood of observing the number of violations or less
    likelihood = poisson.cdf(n_violations, expected_violations)

    return likelihood

def calc_local_firing_rate(spike_times, fr_est_dur=1.0):
    """
    Compute the local (autocorrelation-weighted) firing rate of a spike train.

    Produces the same result as
    ``calc_ccgs(spike_times, [0, fr_est_dur]).squeeze() / fr_est_dur / n_spikes``
    but entirely via vectorised NumPy, without the CCG histogram overhead.

    For each spike i, count the number of spikes j > i with
    ``spike_times[j] - spike_times[i] < fr_est_dur``.  Summing those counts
    and dividing by ``n_spikes * fr_est_dur`` gives the same ACG-density
    estimate, which correctly tracks the *local* firing rate:  a neuron active
    for only 10 min of a 2-hour recording returns its true ~20 Hz rate, not the
    global time-averaged ~1.6 Hz.

    Parameters
    ----------
    spike_times : np.ndarray (n_spikes,)
        Sorted spike times in seconds.
    fr_est_dur : float
        Half-window duration in seconds (default 1.0).

    Returns
    -------
    firing_rate : float
        Local firing rate in Hz.
    """
    n_spikes = len(spike_times)
    if n_spikes == 0:
        return 0.0
    window_ends = np.searchsorted(spike_times, spike_times + fr_est_dur, side='right')
    spikes_in_window = window_ends - np.arange(n_spikes) - 1
    return np.sum(spikes_in_window) / (n_spikes * fr_est_dur)


def calc_rp_violations(spike_times, refractory_periods, ref_acg_t_start):
    """
    Compute cumulative refractory period violation counts using k-th order ISIs.

    Produces the same result as ``np.cumsum(calc_ccgs(st, np.r_[ref_acg_t_start,
    refractory_periods]).squeeze())``, but without the 3-D histogram overhead of
    calc_ccgs.  The ACG at lag τ equals the sum of all k-th order ISIs falling in
    [0, τ]:

        ACG(τ) = Σ_{k=1}^{∞} ISI_k(τ)

    For a 10 ms window, k rarely exceeds 2–3, so the outer loop terminates
    quickly.  Cumulative counts are then read off in O(n log n) via np.searchsorted
    rather than by binning.

    Parameters
    ----------
    spike_times : np.ndarray (n_spikes,)
        Sorted spike times in seconds.
    refractory_periods : np.ndarray (n_refractory_periods,)
        Refractory period upper bounds in seconds (monotonically increasing).
    ref_acg_t_start : float
        Lower bound of the violation window in seconds.

    Returns
    -------
    n_violations : np.ndarray (n_refractory_periods,) int
        Cumulative count of spike pairs with lag in [ref_acg_t_start, τ_r]
        for each τ_r.
    """
    max_tau = refractory_periods[-1]
    all_dts = []

    shift = 1
    while True:
        dts = spike_times[shift:] - spike_times[:-shift]
        valid_dts = dts[dts <= max_tau]
        if len(valid_dts) == 0:
            break
        all_dts.append(valid_dts)
        shift += 1

    if not all_dts:
        return np.zeros(len(refractory_periods), dtype=np.intp)

    all_dts = np.concatenate(all_dts)
    valid_dts = all_dts[all_dts >= ref_acg_t_start]
    valid_dts.sort()
    return np.searchsorted(valid_dts, refractory_periods, side='right')


def compute_min_contam_props(spike_times, spike_clusters=None, cids=None,
                       refractory_periods=np.exp(np.linspace(np.log(0.5e-3), np.log(10e-3), 100)),
                       max_contam_prop=1,
                       fr_est_dur = 1,
                       confidence = 0.95,
                       ref_acg_t_start = .25e-3,
                       progress = False,
                       tol = 10e-4,
                       max_iter = 100):
    '''
    Compute the minimum contamination proportion that can be rejected for each cluster
    via binary search over a Poisson model of refractory period violations.

    For each cluster and each tested refractory period τ_r, finds the smallest C such
    that poisson.cdf(V_o, V_e(C)) < (1 - confidence), i.e. the hypothesis that
    contamination equals C can be rejected at the given confidence level.

    Parameters
    ----------
    spike_times : array-like (n_spikes,)
        Spike times in seconds.
    spike_clusters : array-like (n_spikes,)
        Cluster IDs for each spike. If None, all spikes are assumed to be in the same cluster.
    cids : array-like (n_clusters,)
        Cluster IDs to test. Results returned in order of cids. If None, all clusters are tested.
    refractory_periods : array-like (n_refractory_periods,)
        Refractory periods to test in seconds.
    max_contam_prop : float
        Maximum contamination proportion to test (as a proportion of the firing rate).
    fr_est_dur : float
        Duration of the firing rate estimation window in seconds.
    confidence : float
        Confidence level for the test (e.g. 0.95 means 95% confidence). The rejection
        threshold on the Poisson likelihood is (1 - confidence).
    ref_acg_t_start : float
        Start time for the refractory period autocorrelogram in seconds.
        (necessary because Kilosort removes "duplicate" spikes within a .25 ms window)
    progress : bool
        Show a progress bar.
    tol : float
        Convergence tolerance for the binary search. Stops early when the search interval
        width falls below this value.
    max_iter : int
        Maximum number of binary search iterations.

    Returns
    -------
    min_contam_props : array (n_clusters, n_refractory_periods)
        Minimum contamination proportion that can be rejected under a Poisson model of
        refractory violations.
    firing_rates : array (n_clusters,)
        Firing rates for each cluster.
    '''
    spike_times = ensure_ndarray(spike_times).squeeze()

    if spike_clusters is None:
        spike_clusters = np.zeros(len(spike_times), dtype=np.int32)
    spike_clusters = ensure_ndarray(spike_clusters, dtype=np.int32).squeeze()
    assert spike_clusters.ndim == 1
    assert len(spike_times) == len(spike_clusters), "Spike times and spike clusters must have the same length."

    if cids is not None:
        cids = ensure_ndarray(cids, dtype=np.int32)
        cids_check = np.unique(spike_clusters)
        assert np.all(np.in1d(cids, cids_check)), "Some clusters are not in spike_clusters."
    else:
        cids = np.unique(spike_clusters)

    assert np.all(refractory_periods > 0), "Refractory periods must be positive."
    assert np.all(np.diff(refractory_periods) > 0), "Refractory periods must be monotonic."
    assert max_contam_prop > 0, "Contamination test proportions must be positive."


    firing_rates = np.zeros(len(cids))
    min_contam_props = np.ones((len(cids), len(refractory_periods))) * max_contam_prop
    for iC in tqdm(range(len(cids)), disable=not progress, desc="Calculating contamination"):
        cid = cids[iC]
        st_clu = spike_times[spike_clusters == cid]
        n_spikes = len(st_clu)
        firing_rate = calc_local_firing_rate(st_clu, fr_est_dur)
        firing_rates[iC] = firing_rate
        n_violations = calc_rp_violations(st_clu, refractory_periods, ref_acg_t_start)
        adj_ref_periods = refractory_periods - ref_acg_t_start

        # Vectorized binary search across all refractory periods simultaneously
        left = np.zeros(len(refractory_periods))
        right = np.full(len(refractory_periods), float(max_contam_prop))
        mid = np.empty(len(refractory_periods))
        for _ in range(max_iter):
            mid = (left + right) / 2
            likelihood = refractory_violation_likelihood(
                n_violations, mid, adj_ref_periods, firing_rate, n_spikes)
            too_high = likelihood < (1 - confidence)
            right = np.where(too_high, mid, right)
            left = np.where(too_high, left, mid)
            if np.max(right - left) < tol:
                break
        min_contam_props[iC] = mid

    return min_contam_props, firing_rates

def compute_min_contam_props_analytical(spike_times, spike_clusters=None, cids=None,
                       refractory_periods=np.exp(np.linspace(np.log(0.5e-3), np.log(10e-3), 100)),
                       max_contam_prop=1,
                       fr_est_dur=1,
                       confidence=0.95,
                       ref_acg_t_start=.25e-3,
                       progress=False):
    '''
    Compute the minimum contamination proportion that can be rejected using the exact
    analytical solution.

    Derivation
    ----------
    We want the smallest contamination C such that the Poisson likelihood of the observed
    violations V_o falls below (1 − confidence), i.e. find λ_crit satisfying:

        P(Poisson(λ_crit) ≤ V_o) = 1 − confidence                             ... (1)

    **Step 1 — Relate the Poisson CDF to the chi-squared CDF.**

    The Poisson CDF can be written in terms of the regularised upper incomplete gamma
    function Q(a, x) = Γ(a, x) / Γ(a):

        P(Poisson(λ) ≤ r) = Q(r + 1, λ)                                       ... (2)

    The chi-squared distribution with ν degrees of freedom is Gamma(ν/2, 2), so its CDF
    is the regularised *lower* incomplete gamma function:

        chi2.cdf(x, ν) = γ(ν/2, x/2) / Γ(ν/2)                                ... (3)

    Because the upper and lower incomplete gamma functions sum to the complete gamma,
    Q(a, x) = 1 − γ(a, x)/Γ(a).  Substituting a = r+1 and x = λ into (2) and comparing
    with (3) at ν = 2*(r+1) and x = 2λ:

        P(Poisson(λ) ≤ r) = Q(r+1, λ)
                           = 1 − γ(r+1, λ) / Γ(r+1)
                           = 1 − chi2.cdf(2λ, df=2*(r+1))                     ... (4)

    **Step 2 — Solve for λ_crit analytically.**

    Substituting (4) into (1):

        1 − chi2.cdf(2λ_crit, df=2*(V_o+1)) = 1 − confidence
        chi2.cdf(2λ_crit, df=2*(V_o+1))     = confidence
        2λ_crit = chi2.ppf(confidence, df=2*(V_o+1))
        λ_crit  = chi2.ppf(confidence, df=2*(V_o+1)) / 2                      ... (5)

    **Step 3 — Invert the expected-violations formula to recover C.**

    From the module docstring, V_e = C*(2−C)*F_r*τ_adj*N_s.  Setting V_e = λ_crit and
    letting k = λ_crit / (F_r * τ_adj * N_s):

        C*(2 − C) = k
        C^2 − 2C + k = 0
        C = 1 − sqrt(1 − k)       (smaller root, valid for C ∈ [0, 1])        ... (6)

    When k > 1 (i.e. λ_crit implies contamination above 100%) the result is capped at
    max_contam_prop via np.maximum(1 − k, 0) before taking the square root.

    Parameters (same as compute_min_contam_props)
    ----------
    spike_times : array-like (n_spikes,)
    spike_clusters : array-like (n_spikes,)
    cids : array-like (n_clusters,)
    refractory_periods : array-like (n_refractory_periods,)
    max_contam_prop : float
    fr_est_dur : float
    confidence : float
        Confidence level for the test (e.g. 0.95 means 95% confidence). The rejection
        threshold on the Poisson likelihood is (1 - confidence).
    ref_acg_t_start : float
    progress : bool

    Returns
    -------
    min_contam_props : array (n_clusters, n_refractory_periods)
    firing_rates : array (n_clusters,)
    '''
    spike_times = ensure_ndarray(spike_times).squeeze()

    if spike_clusters is None:
        spike_clusters = np.zeros(len(spike_times), dtype=np.int32)
    spike_clusters = ensure_ndarray(spike_clusters, dtype=np.int32).squeeze()
    assert spike_clusters.ndim == 1
    assert len(spike_times) == len(spike_clusters)

    if cids is not None:
        cids = ensure_ndarray(cids, dtype=np.int32)
    else:
        cids = np.unique(spike_clusters)

    adj_ref_periods = refractory_periods - ref_acg_t_start

    firing_rates = np.zeros(len(cids))
    min_contam_props = np.ones((len(cids), len(refractory_periods))) * max_contam_prop
    for iC in tqdm(range(len(cids)), disable=not progress, desc="Calculating contamination (analytical)"):
        cid = cids[iC]
        st_clu = spike_times[spike_clusters == cid]
        n_spikes = len(st_clu)
        firing_rate = calc_local_firing_rate(st_clu, fr_est_dur)
        firing_rates[iC] = firing_rate
        n_violations = calc_rp_violations(st_clu, refractory_periods, ref_acg_t_start)

        # Analytical solution: invert C*(2-C) * F_r * tau * N_s = lambda_critical
        # lambda_critical = chi2.ppf(confidence, 2*(r+1)) / 2
        # C*(2-C) = k  =>  C^2 - 2C + k = 0  =>  C = 1 - sqrt(1 - k)
        lambda_critical = chi2.ppf(confidence, df=2 * (n_violations + 1)) / 2
        k = lambda_critical / (firing_rate * adj_ref_periods * n_spikes)
        contam = 1 - np.sqrt(np.maximum(1 - k, 0))
        min_contam_props[iC] = np.minimum(contam, max_contam_prop)

    return min_contam_props, firing_rates


def plot_min_contam_prop(spike_times, min_contam_props, refractory_periods,
                         n_bins = 50, max_contam_prop=1, acg_t_start = .25e-3, axs=None):
    '''
    Utility for plotting the minimum contamination proportion that can be rejected for each cluster in the dataset.
    
    Parameters
    ----------
    spike_times : array-like (n_spikes,)
        Spike times in seconds.
    min_contam_props : array (n_refractory_periods)
        Minimum contamination rate that can be rejected under a poisson model of refractory violations.
    refractory_periods : array-like (n_refractory_periods,)
        Refractory periods to test in seconds.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axs : list of matplotlib.axes.Axes (2,)
        The axes objects.

    '''

    isis = np.diff(spike_times) * 1000
    max_refrac = refractory_periods.max() * 1000
    min_isi = acg_t_start * 1000
    min_prop = min_contam_props.min()

    if axs is None:
        fig, axs = plt.subplots(1,1)
    else:
        fig = axs.get_figure()
    bins = np.linspace(min_isi, max_refrac, n_bins)
    axs.hist(isis, bins=bins, edgecolor='black', color='black', alpha=0.6)
    axs.set_xlim([min_isi, max_refrac])
    axs.set_ylabel('ISI count (spikes)')
    axs.set_xlabel('ISI / Refractory Period (ms)')
    axs2 = axs.twinx()
    axs2.plot(refractory_periods*1000, min_contam_props, color='red', linewidth=3.5)
    axs2.axhline(min_prop, color='red', linestyle='--', linewidth=2)
    yticks = np.concatenate([np.linspace(0, max_contam_prop, 6), [min_prop]])
    axs2.set_ylim([0, max_contam_prop])
    axs2.set_yticks(yticks)
    axs2.set_yticklabels(['0', '', '', '', '', '1', f'{min_prop:.4g}'])
    axs2.tick_params(axis='y', colors='red')
    axs2.set_ylabel('Minimum Rejected Contamination Proportion', color='red')

    return fig, axs

# Depricated code from Nick Steinmetz's lab (Sliding RP violations)
# https://github.com/SteinmetzLab/slidingRefractory/blob/1.0.0/python/slidingRP/metrics.py
def compute_rvl_tensor(spike_times, spike_clusters=None, cids=None,
                      refractory_periods=np.exp(np.linspace(np.log(0.5e-3), np.log(10e-3), 100)),
                      contamination_test_proportions=np.exp(np.linspace(np.log(5e-3), np.log(.35), 50)),
                      fr_est_dur = 1,
                      ref_acg_t_start = .25e-3, 
                      progress = False):
    '''
    Compute the likelihood of observing the number of refractory period violations or fewer for many clusters, refractory periods, and test comtamination rates.

    Parameters
    ----------
    spike_times : array-like (n_spikes,)
        Spike times in seconds.
    spike_clusters : array-like (n_spikes,)
        The cluster ids for each spike. If None, all spikes are assumed to belong to a single cluster.
    cids : array-like (n_clusters,)
        The list of *all* unique clusters, in any order. That order will be used in the output array. If None, order the clusters by their appearance in `spike_clusters`.
    refractory_periods : array-like (n_refrac,)
        The refractory periods to test, in seconds.
    contamination_test_proportions : array-like (n_contam,)
        The contamination rates to test, as a proportion of the firing rate.
    fr_est_dur : float
        The duration in seconds over which to estimate the firing rate. 
    ref_acg_t_start : float. Default is .25e-3
        The start time in seconds for the refractory period autocorrelogram. 
        Necessary for Kilosort4, which removes duplicate spikes in a .25 ms window, which negatively biases the refractory likelihood estimates.

    Returns
    -------
    rvl_tensor : array
        A `(n_clusters, n_refrac, n_contam)` array with the likelihood of observing the number of refractory period violations or less if the cluster was contaminated at the rate specified.
    '''
    spike_times = ensure_ndarray(spike_times).squeeze()

    if spike_clusters is None:
        spike_clusters = np.zeros(len(spike_times), dtype=np.int32)
    spike_clusters = ensure_ndarray(spike_clusters, dtype=np.int32).squeeze()
    assert spike_clusters.ndim == 1
    assert len(spike_times) == len(spike_clusters), "Spike times and spike clusters must have the same length."

    if cids is not None:
        cids = ensure_ndarray(cids, dtype=np.int32)
        cids_check = np.unique(spike_clusters)
        assert np.all(np.in1d(cids, cids_check)), "Some clusters are not in spike_clusters."
    else:
        cids = np.unique(spike_clusters)

    rvl_tensor = np.ones((len(cids), len(contamination_test_proportions), len(refractory_periods)))

    iter = range(len(cids))
    if progress:
        iter = tqdm(iter, desc="Calculating RVL tensor", position=0, leave=True)
    for iC in iter:
        cid = cids[iC]
        cluster_spikes = spike_times[spike_clusters == cid]
        n_spikes = len(cluster_spikes)
        firing_rate = calc_local_firing_rate(cluster_spikes, fr_est_dur)
        refractory_violations = calc_rp_violations(cluster_spikes, refractory_periods, ref_acg_t_start)

        rvl_tensor[iC] = refractory_violation_likelihood(
                            refractory_violations[None,:], 
                            contamination_test_proportions[:,None],
                            refractory_periods[None,:] - ref_acg_t_start, 
                            firing_rate, 
                            n_spikes)

    return rvl_tensor

def plot_rvl(cluster_spikes, likelihoods, refractory_periods, contamination_test_proportions, likelihood_threshold=0.05):
    min_refrac, max_refrac = refractory_periods.min(), refractory_periods.max()
    min_contam, max_contam = contamination_test_proportions.min(), contamination_test_proportions.max()

    isis = np.diff(cluster_spikes)

    min_likelihood_per_contam = np.min(likelihoods, axis=1)
    min_likelihood_per_contam[min_likelihood_per_contam < likelihood_threshold] = np.inf
    lowest_contam_idx = np.argmin(min_likelihood_per_contam)
    lowest_contam = contamination_test_proportions[lowest_contam_idx]
    lowest_contam_likelihood = likelihoods[lowest_contam_idx]

    fig, axs = plt.subplots(3, 1, figsize=(5, 12), height_ratios=[1, 1.5, 1])
    axs[0].hist(isis * 1000, bins=np.arange(0, max_refrac*1000, .33))
    axs[0].set_title(f'ISI distribution')
    axs[0].set_ylabel('Count')
    axs[0].set_xlabel('ISI (ms)')
    axs[0].set_xlim([0, max_refrac*1000])

    extent = [min_refrac*1000, max_refrac*1000, min_contam, max_contam]
    from matplotlib.image import NonUniformImage
    im = NonUniformImage(axs[1], extent=extent, interpolation='nearest', cmap='viridis')
    im.set_data(refractory_periods*1000, contamination_test_proportions, likelihoods)
    im.set_clim(0, 1)
    axs[1].add_image(im)
    axs[1].axhline(lowest_contam, color='red', linestyle='--')
    axs[1].set_xlim(extent[:2])
    axs[1].set_ylim(extent[2:])
    fig.colorbar(im, ax=axs[1], orientation='horizontal', label='Likelihood')
    axs[1].set_title(f'Likelihood of observed refractory period violations')
    axs[1].set_xlabel('Refractory period (ms)')
    axs[1].set_ylabel('Contamination rate')

    axs[2].semilogy(refractory_periods*1000, lowest_contam_likelihood)
    axs[2].axhline(likelihood_threshold, color='red', linestyle='--')
    axs[2].set_title(f'Highest contamination rate more than {likelihood_threshold*100:.1g}% likely: {lowest_contam*100:.3g}%')
    axs[2].set_xlabel('Refractory period (ms)')
    axs[2].set_ylabel('Likelihood')
    axs[2].set_xlim([min_refrac*1000, max_refrac*1000])

    plt.tight_layout()
    return fig, axs


if __name__ == "__main__":
    import time

    rng = np.random.default_rng(42)
    duration = 3600.0
    rp = 2.5e-3
    n_units = 500

    # Generate 100 units with varying firing rates and contamination fractions
    real_rates = rng.uniform(5.0, 40.0, size=n_units)
    real_rates[0] = 20.0  # Set unit 0 to a known firing rate for plotting
    contam_fracs = rng.uniform(0.01, 0.50, size=n_units)
    contam_fracs[0] = 0.10  # Set unit 0 to a known contamination fraction for plotting

    unit_spike_times = []  # per-unit spike times (for plotting)
    all_spike_times_list = []
    all_spike_clusters_list = []

    for uid in range(n_units):
        real_rate = real_rates[uid]
        contam_rate = real_rate * contam_fracs[uid]

        # Real spikes: uniform draw then remove RP violations
        real_spikes = np.sort(rng.uniform(0, duration, size=int(duration * real_rate)))
        rp_violations = np.diff(real_spikes, prepend=-np.inf) < rp
        real_spikes = real_spikes[~rp_violations]

        # Contamination spikes: pure Poisson (uniform), no refractory period
        contam_spikes = np.sort(rng.uniform(0, duration, size=int(contam_rate * duration)))

        unit_spikes = np.sort(np.concatenate([real_spikes, contam_spikes]))
        unit_spike_times.append(unit_spikes)
        all_spike_times_list.append(unit_spikes)
        all_spike_clusters_list.append(np.full(len(unit_spikes), uid, dtype=np.int32))

    spike_times = np.concatenate(all_spike_times_list)
    spike_clusters = np.concatenate(all_spike_clusters_list)
    sort_idx = np.argsort(spike_times, kind='stable')
    spike_times = spike_times[sort_idx]
    spike_clusters = spike_clusters[sort_idx]
    cids = np.arange(n_units, dtype=np.int32)

    print(f"Total spikes: {len(spike_times)} across {n_units} units")

    refractory_periods = np.exp(np.linspace(np.log(0.5e-3), np.log(10e-3), 100))

    confidence = 0.90

    # --- Analytical solution ---
    t0 = time.perf_counter()
    min_contam_analytical, firing_rates_analytical = compute_min_contam_props_analytical(
        spike_times, spike_clusters, cids=cids, confidence=confidence, progress=True)
    t_analytical = time.perf_counter() - t0
    print(f"Analytical:      {t_analytical:.3f} s")

    # --- Binary search approach ---
    t0 = time.perf_counter()
    min_contam_binary, firing_rates_binary = compute_min_contam_props(
        spike_times, spike_clusters, cids=cids, confidence=confidence, progress=True)
    t_binary = time.perf_counter() - t0
    print(f"Binary search:   {t_binary:.3f} s")

    # --- Full RVL tensor approach ---
    contamination_test_proportions = np.exp(np.linspace(np.log(5e-3), np.log(1.0), 500))
    t0 = time.perf_counter()
    rvl_tensor = compute_rvl_tensor(
        spike_times, spike_clusters, cids=cids,
        contamination_test_proportions=contamination_test_proportions,
        progress=True)
    t_tensor = time.perf_counter() - t0
    print(f"Full RVL tensor: {t_tensor:.3f} s")
    print(f"Analytical vs binary speedup:     {t_binary / t_analytical:.1f}x")
    print(f"Analytical vs RVL tensor speedup: {t_tensor / t_analytical:.1f}x")

    # --- Profiling: exclude calc_ccgs time ---
    # Precompute all CCG data once, then time only the core computation of each method.
    fr_est_dur = 1
    ref_acg_t_start = .25e-3
    adj_ref_periods = refractory_periods - ref_acg_t_start
    max_contam_prop = 1
    tol = 10e-4
    max_iter = 100

    # --- Firing rate computation: calc_ccgs vs calc_local_firing_rate ---
    print("\n--- Firing rate computation: calc_ccgs vs calc_local_firing_rate ---")
    t0 = time.perf_counter()
    for iC in range(len(cids)):
        st_clu = spike_times[spike_clusters == cids[iC]]
        _ = calc_ccgs(st_clu, [0, fr_est_dur]).squeeze() / fr_est_dur / len(st_clu)
    t_ccgs_fr = time.perf_counter() - t0

    t0 = time.perf_counter()
    for iC in range(len(cids)):
        st_clu = spike_times[spike_clusters == cids[iC]]
        _ = calc_local_firing_rate(st_clu, fr_est_dur)
    t_local_fr = time.perf_counter() - t0

    print(f"calc_ccgs (FR):          {t_ccgs_fr:.3f} s")
    print(f"calc_local_firing_rate:  {t_local_fr:.3f} s")
    print(f"FR speedup:              {t_ccgs_fr / t_local_fr:.1f}x")

    # --- Violations computation: calc_ccgs vs calc_rp_violations ---
    print("\n--- Violations computation: calc_ccgs vs calc_rp_violations ---")
    t0 = time.perf_counter()
    for iC in range(len(cids)):
        st_clu = spike_times[spike_clusters == cids[iC]]
        _ = np.cumsum(calc_ccgs(st_clu, np.r_[ref_acg_t_start, refractory_periods]).squeeze())
    t_ccgs_viol = time.perf_counter() - t0

    t0 = time.perf_counter()
    for iC in range(len(cids)):
        st_clu = spike_times[spike_clusters == cids[iC]]
        _ = calc_rp_violations(st_clu, refractory_periods, ref_acg_t_start)
    t_rp_viol = time.perf_counter() - t0

    print(f"calc_ccgs + cumsum:  {t_ccgs_viol:.3f} s")
    print(f"calc_rp_violations:  {t_rp_viol:.3f} s")
    print(f"Violations speedup:  {t_ccgs_viol / t_rp_viol:.1f}x")

    print("\n--- Precomputing data for core timing (FR via calc_ccgs, violations via calc_rp_violations) ---")
    t0 = time.perf_counter()
    precomputed = []  # list of (n_spikes, firing_rate, n_violations) per cluster
    for iC in tqdm(range(len(cids)), desc="Precomputing"):
        cid = cids[iC]
        st_clu = spike_times[spike_clusters == cid]
        n_spikes = len(st_clu)
        firing_rate = calc_local_firing_rate(st_clu, fr_est_dur)
        n_violations = calc_rp_violations(st_clu, refractory_periods, ref_acg_t_start)
        precomputed.append((n_spikes, firing_rate, n_violations))
    t_precompute = time.perf_counter() - t0
    print(f"Precomputation: {t_precompute:.3f} s")

    # Core of analytical method (no CCGs)
    t0 = time.perf_counter()
    min_contam_analytical_core = np.ones((len(cids), len(refractory_periods))) * max_contam_prop
    for iC, (n_spikes, firing_rate, n_violations) in enumerate(precomputed):
        lambda_critical = chi2.ppf(confidence, df=2 * (n_violations + 1)) / 2
        k = lambda_critical / (firing_rate * adj_ref_periods * n_spikes)
        contam = 1 - np.sqrt(np.maximum(1 - k, 0))
        min_contam_analytical_core[iC] = np.minimum(contam, max_contam_prop)
    t_analytical_core = time.perf_counter() - t0

    # Core of binary search method (no CCGs)
    t0 = time.perf_counter()
    min_contam_binary_core = np.ones((len(cids), len(refractory_periods))) * max_contam_prop
    for iC, (n_spikes, firing_rate, n_violations) in enumerate(precomputed):
        left = np.zeros(len(refractory_periods))
        right = np.full(len(refractory_periods), float(max_contam_prop))
        mid = np.empty(len(refractory_periods))
        for _ in range(max_iter):
            mid = (left + right) / 2
            likelihood = refractory_violation_likelihood(
                n_violations, mid, adj_ref_periods, firing_rate, n_spikes)
            too_high = likelihood < (1 - confidence)
            right = np.where(too_high, mid, right)
            left = np.where(too_high, left, mid)
            if np.max(right - left) < tol:
                break
        min_contam_binary_core[iC] = mid
    t_binary_core = time.perf_counter() - t0

    # Core of RVL tensor method (no CCGs)
    t0 = time.perf_counter()
    rvl_tensor_core = np.ones((len(cids), len(contamination_test_proportions), len(refractory_periods)))
    for iC, (n_spikes, firing_rate, n_violations) in enumerate(precomputed):
        rvl_tensor_core[iC] = refractory_violation_likelihood(
            n_violations[None, :],
            contamination_test_proportions[:, None],
            refractory_periods[None, :] - ref_acg_t_start,
            firing_rate,
            n_spikes)
    t_tensor_core = time.perf_counter() - t0

    print(f"\n--- Core computation times (precompute excluded, {t_precompute:.2f} s shared) ---")
    print(f"Analytical core:       {t_analytical_core*1000:.2f} ms")
    print(f"Binary search core:    {t_binary_core*1000:.2f} ms  ({t_binary_core/t_analytical_core:.1f}x vs analytical)")
    print(f"RVL tensor core:       {t_tensor_core*1000:.2f} ms  ({t_tensor_core/t_analytical_core:.1f}x vs analytical)")
    print(f"\nPrecompute fraction of total: analytical {t_precompute/(t_precompute+t_analytical_core)*100:.1f}%"
          f", binary {t_precompute/(t_precompute+t_binary_core)*100:.1f}%"
          f", tensor {t_precompute/(t_precompute+t_tensor_core)*100:.1f}%")

    # Derive min_contam_props from rvl_tensor for comparison
    rvl = rvl_tensor[0]  # (n_contam, n_refrac)
    min_contam_from_tensor = np.full(len(refractory_periods), contamination_test_proportions[-1])
    for iR in range(len(refractory_periods)):
        below = np.where(rvl[:, iR] < (1 - confidence))[0]
        if len(below) > 0:
            min_contam_from_tensor[iR] = contamination_test_proportions[below[0]]

    # --- Plot: 2x2 grid comparison for unit 0 ---
    # Top row: RVL tensor (current version) — heatmap + min contam curve
    # Bottom row: binary search + analytical (legacy methods)
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    unit_label = (
        f"Unit 0  (real: {real_rates[0]:.1f} Hz, RP: {rp*1000:.1f} ms, "
        f"{contam_fracs[0]*100:.1f}% contam)  |  {n_units} units total\n"
        f"precompute (shared): {t_precompute:.2f} s   |   "
        f"Analytical solution: {t_analytical_core*1000:.1f} ms   |   "
        f"Binary search: {t_binary_core*1000:.1f} ms ({t_binary_core/t_analytical_core:.1f}\u00d7)   |   "
        f"RVL Tensor: {t_tensor_core*1000:.1f} ms ({t_tensor_core/t_analytical_core:.1f}\u00d7)\n"
        f"Binary: {t_binary:.2f} s ({t_binary/t_analytical:.1f}\u00d7 slower)   |   "
        f"RVL tensor: {t_tensor:.2f} s ({t_tensor/t_analytical:.1f}\u00d7 slower)"
    )
    fig.suptitle(unit_label, fontsize=9)

    # Top-left: RVL tensor heatmap
    im = axs[0, 0].pcolormesh(refractory_periods * 1000, contamination_test_proportions,
                               rvl_tensor[0], cmap='viridis', vmin=0, vmax=1)
    axs[0, 0].contour(refractory_periods * 1000, contamination_test_proportions,
                      rvl_tensor[0], levels=[1 - confidence], colors='red', linewidths=1.5)
    fig.colorbar(im, ax=axs[0, 0], label='Likelihood')
    axs[0, 0].set_xlabel('Refractory period (ms)')
    axs[0, 0].set_ylabel('Contamination proportion')
    axs[0, 0].set_title(f'RVL tensor heatmap\n(confidence={confidence} contour in red)')

    # Top-right: RVL tensor min contam curve
    plot_min_contam_prop(unit_spike_times[0], min_contam_from_tensor, refractory_periods,
                         max_contam_prop=1.0, axs=axs[0, 1])
    axs[0, 1].set_title(f"RVL tensor (min contam curve)\nprecompute {t_precompute:.2f} s + core {t_tensor_core*1000:.1f} ms = {t_tensor:.2f} s total")

    # Bottom-left: Binary search
    plot_min_contam_prop(unit_spike_times[0], min_contam_binary[0], refractory_periods,
                         max_contam_prop=1.0, axs=axs[1, 0])
    axs[1, 0].set_title(f"Binary search\nprecompute {t_precompute:.2f} s + core {t_binary_core*1000:.1f} ms = {t_binary:.2f} s total")

    # Bottom-right: Analytical
    plot_min_contam_prop(unit_spike_times[0], min_contam_analytical[0], refractory_periods,
                         max_contam_prop=1.0, axs=axs[1, 1])
    axs[1, 1].set_title(f"Analytical\nprecompute {t_precompute:.2f} s + core {t_analytical_core*1000:.1f} ms = {t_analytical:.2f} s total")

    plt.tight_layout()
    plt.savefig("analytical_vs_binary_vs_tensor.png", dpi=150)
    print("\nPlot saved to analytical_vs_binary_vs_tensor.png")
    plt.show()

