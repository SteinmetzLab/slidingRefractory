"""Estimate a unit's apparent refractory ("ACG recovery") duration from its ACG.

Python port of ``matlab/estimate_refractory_period.m``, which is the estimator
used for the macaque data in Roth et al. (2026) Fig 1 and, following the
J Neurophysiol revision, the single estimator applied to every dataset.

This is a *new* module: it does not change the behaviour of anything already in
the package. The older :func:`slidingRP.metrics.compute_rf` (a 2-parameter
sigmoid at 5% recovery, used for the mouse data in the submitted manuscript) is
retained for provenance; the two disagree by a median 0.15 ms on identical
ACGs, which is why the revision uses one estimator throughout.

The quantity returned is explicitly operational: the time at which a sigmoid
fitted to the rising phase of the ACG has recovered ``recovery_frac`` of the way
from its lower to its upper asymptote. Bursting, firing rate, residual
contamination, spike detection and sorter behaviour all shape the short-lag
ACG, so this is an *apparent* recovery time and not a biophysical absolute
refractory period.

Algorithm (identical to the MATLAB original)
--------------------------------------------
1. Median-filter the ACG with a 0.83 ms window (odd number of bins, zero-padded
   edges, matching MATLAB ``medfilt1``).
2. Trough: the *maximum* of the filtered ACG within 0 to 0.5 ms (the maximum,
   not the minimum, so that subsequent peak detection is forced outside this
   window).
3. Peak: all local maxima of the filtered ACG that exceed both 10% of its full
   range and the trough value; take the one closest to zero lag.
4. Fit a 4-parameter sigmoid between trough time and peak time.
5. RP estimate = the time at which the fit reaches ``recovery_frac`` of the way
   from its lower to its upper asymptote.
6. Floor: if the estimate falls before the first ACG bin with a nonzero count,
   replace it with that bin's time (recorded in ``floor_applied``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, signal, special


@dataclass
class RPFit:
    """Result of :func:`estimate_rp`."""

    rp_ms: float                 # estimated recovery time (ms), NaN on failure
    params: np.ndarray | None    # [ymin, ymax, k, x0] of the fitted sigmoid
    success: bool
    rmse: float
    rsquared: float
    floor_applied: bool          # estimate was clipped to the first nonzero bin
    reason: str                  # why the fit failed, '' if it succeeded


def _sigmoid(x, ymin, ymax, k, x0):
    # expit is the overflow-safe 1/(1+exp(-z)); the optimiser explores large k,
    # where the naive form warns and returns inf/nan.
    return ymin + (ymax - ymin) * special.expit(k * (x - x0))


def _sigmoid_jac(x, ymin, ymax, k, x0):
    """Analytic Jacobian of _sigmoid w.r.t. [ymin, ymax, k, x0].

    Supplying it rather than letting least_squares difference numerically is
    ~3x faster, which matters when the estimator runs over ~10^6 units.
    """
    s = special.expit(k * (x - x0))
    ds = s * (1 - s)
    rng_ = ymax - ymin
    return np.column_stack([1 - s, s, rng_ * ds * (x - x0), -rng_ * ds * k])


def estimate_rp(acg, bin_centers=None, bin_size=1 / 30000, recovery_frac=0.10):
    """Estimate the apparent refractory / ACG recovery duration.

    Parameters
    ----------
    acg : array_like
        Single-sided autocorrelogram counts, starting at zero lag.
    bin_centers : array_like, optional
        Bin times in seconds. Defaults to ``arange(len(acg)) * bin_size``,
        which is the convention used by the MATLAB caller.
    bin_size : float
        ACG bin width in seconds (default one sample at 30 kHz).
    recovery_frac : float
        Fraction of the fitted range at which the RP is read off (0.10 in the
        manuscript; 0.05 and 0.20 are reported as sensitivity analyses).

    Returns
    -------
    RPFit
    """
    acg = np.asarray(acg, dtype=np.float64)
    if bin_centers is None:
        bin_centers = np.arange(acg.size) * bin_size
    bin_centers = np.asarray(bin_centers, dtype=np.float64)
    fail = lambda why: RPFit(np.nan, None, False, np.nan, np.nan, False, why)  # noqa: E731

    if acg.size < 3 or not np.any(acg > 0):
        return fail("empty acg")

    # 1. median filter, odd window, zero-padded edges (MATLAB medfilt1 default)
    nfilt = int(round(0.83e-3 / bin_size))
    if nfilt % 2 == 0:
        nfilt += 1
    nfilt = max(nfilt, 1)
    filt = signal.medfilt(acg, kernel_size=nfilt)

    # 2. trough: max of the filtered ACG within 0 to 0.5 ms
    narrow = (bin_centers >= 0) & (bin_centers <= 0.5e-3)
    if not np.any(narrow):
        return fail("no bins in 0-0.5 ms")
    min_value = float(np.max(filt[narrow]))
    idx_narrow = np.flatnonzero(narrow)
    last = idx_narrow[filt[narrow] == min_value][-1]      # MATLAB uses 'last'
    min_time = float(bin_centers[last])

    # 3. first valid peak
    peak_locs, _ = signal.find_peaks(filt)
    if peak_locs.size == 0:
        return fail("no peaks")
    peak_values = filt[peak_locs]
    peak_times = bin_centers[peak_locs]
    threshold = 0.1 * (np.max(filt) - np.min(filt))
    valid = (peak_values > threshold) & (peak_values > min_value)
    if not np.any(valid):
        return fail("no valid peaks")
    vt, vv = peak_times[valid], peak_values[valid]
    closest = int(np.argmin(np.abs(vt)))
    max_time, max_value = float(vt[closest]), float(vv[closest])

    # 4. fit the sigmoid between trough and peak
    win = (bin_centers >= min_time) & (bin_centers <= max_time)
    if np.count_nonzero(win) < 3:
        return fail("too few points to fit")
    x_fit, y_fit = bin_centers[win], filt[win]

    lb = np.array([0.0, min_value, 0.0, min_time])
    ub = np.array([max_value * 2, max_value * 2, np.inf, max_time])
    if not np.all(ub > lb):
        return fail("degenerate fit bounds")

    # Multi-start. The objective is badly conditioned (the parameters span 1 to
    # ~5000) and a single start from the MATLAB initial guess can stall: on 200
    # macaque ACGs, MATLAB's lsqcurvefit stopped at a higher residual than this
    # search on 199/200 units (median 31% higher), shifting the estimate by up
    # to 0.09 ms on 11% of units. Starts vary the steepness and the midpoint;
    # the lowest-cost fit wins, which removes the dependence on the start.
    span = max(max_time - min_time, 1e-12)
    starts = [np.array([min_value, max_value, k / span, min_time + f * span])
              for k in (10.0, 1.0, 3.0, 30.0, 100.0) for f in (0.5,)]
    starts += [np.array([min_value, max_value, 10.0 / span, min_time + f * span])
               for f in (0.25, 0.75)]

    res = None
    for p0 in starts:
        p0 = np.clip(p0, lb + 1e-12, ub - 1e-12)
        try:
            r = optimize.least_squares(
                lambda p: _sigmoid(x_fit, *p) - y_fit, p0, bounds=(lb, ub),
                jac=lambda p: _sigmoid_jac(x_fit, *p),
                method="trf", max_nfev=1000, x_scale="jac")
        except Exception:  # noqa: BLE001
            continue
        if res is None or r.cost < res.cost:
            res = r
    if res is None:
        return fail("fit error")
    if not res.success:
        return fail("fit did not converge")

    ymin, ymax, k, x0 = res.x
    y_pred = _sigmoid(x_fit, *res.x)
    rmse = float(np.sqrt(np.mean(res.fun ** 2)))
    ss_res = float(np.sum((y_fit - y_pred) ** 2))
    ss_tot = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
    rsq = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # 5. read off the recovery time
    target = recovery_frac * (ymax - ymin) + ymin
    if target <= ymin or k <= 0:
        return fail("invalid sigmoid for RP readout")
    ratio = (ymax - ymin) / (target - ymin) - 1
    if ratio <= 0:
        return fail("invalid ratio")
    rp_ms = (x0 - np.log(ratio) / k) * 1000

    # 6. floor at the first nonzero ACG bin
    floor_applied = False
    first_nonzero = int(np.flatnonzero(acg > 0)[0])
    t_first = bin_centers[first_nonzero] * 1000
    if rp_ms < t_first:
        rp_ms = float(t_first)
        floor_applied = True

    return RPFit(float(rp_ms), res.x, True, rmse, float(rsq), floor_applied, "")


def estimate_rp_many(acgs, bin_size=1 / 30000, recovery_frac=0.10, bin_centers=None):
    """Vectorised convenience wrapper: returns (rp_ms, floor_applied, rsquared)."""
    acgs = np.atleast_2d(np.asarray(acgs))
    n = acgs.shape[0]
    rp = np.full(n, np.nan)
    floored = np.zeros(n, dtype=bool)
    rsq = np.full(n, np.nan)
    for i in range(n):
        f = estimate_rp(acgs[i], bin_centers, bin_size, recovery_frac)
        rp[i], floored[i], rsq[i] = f.rp_ms, f.floor_applied, f.rsquared
    return rp, floored, rsq
