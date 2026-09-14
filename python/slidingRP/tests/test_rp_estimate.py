"""Tests for the canonical RP / ACG-recovery estimator (slidingRP.rp_estimate).

Reference values in ``spec_checks/`` were produced by
``matlab/estimate_refractory_period.m`` on 200 macaque LGN ACGs.
"""
from pathlib import Path

import numpy as np
import pytest

from slidingRP.rp_estimate import estimate_rp, estimate_rp_many

BS = 1 / 30000
BIN_CENTERS = np.arange(300) * BS + BS / 2          # histdiff bin centres
SPEC = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
            r"01_fig1_rp_durations/spec_checks")


def synthetic_acg(rp_ms=2.0, baseline=100.0, ramp_ms=0.2, overshoot=0.25):
    """Refractory ACG: zero below rp, a fast ramp, then a small overshoot.

    The overshoot matters: the estimator (like the MATLAB original) locates the
    fit window using the first local *peak* of the ACG, so a perfectly
    monotonic ACG that plateaus has no peak and returns NaN. Real ACGs
    essentially always have one; see test_monotonic_acg_has_no_peak.
    """
    t = BIN_CENTERS * 1000
    acg = baseline / (1 + np.exp(-(t - rp_ms) / ramp_ms))
    acg = acg * (1 + overshoot * np.exp(-((t - rp_ms - 0.8) ** 2) / (2 * 0.5 ** 2)))
    acg[t < rp_ms - 3 * ramp_ms] = 0
    return np.round(acg)


@pytest.mark.parametrize("true_rp", [1.0, 1.5, 2.0, 3.0, 5.0])
def test_recovers_synthetic_rp(true_rp):
    """On a clean sigmoid ACG the estimate lands near the true onset."""
    fit = estimate_rp(synthetic_acg(true_rp), BIN_CENTERS, BS, 0.10)
    assert fit.success
    # the 10% point of a logistic with steepness 1/ramp is x0 - ln(9)*ramp
    expected = true_rp - np.log(9) * 0.2
    assert abs(fit.rp_ms - expected) < 0.1, (fit.rp_ms, expected)
    assert fit.rsquared > 0.99


def test_monotonic_acg_has_no_peak():
    """A plateauing ACG yields no local peak, so no estimate (as in MATLAB)."""
    fit = estimate_rp(synthetic_acg(2.0, overshoot=0.0), BIN_CENTERS, BS)
    assert not fit.success and np.isnan(fit.rp_ms)
    assert fit.reason == "no peaks"


def test_flat_acg_fails_cleanly():
    assert np.isnan(estimate_rp(np.zeros(300), BIN_CENTERS, BS).rp_ms)
    flat = estimate_rp(np.full(300, 5.0), BIN_CENTERS, BS)
    assert np.isnan(flat.rp_ms) or flat.floor_applied


def test_recovery_fraction_is_monotonic():
    """A larger recovery fraction must give a later time on the same fit."""
    acg = synthetic_acg(2.0)
    r = [estimate_rp(acg, BIN_CENTERS, BS, f).rp_ms for f in (0.05, 0.10, 0.20)]
    assert r[0] < r[1] < r[2]


def test_floor_at_first_nonzero_bin():
    """An ACG with an early stray count floors the estimate at that bin."""
    acg = synthetic_acg(3.0)
    acg[6] = 1                       # stray count at 0.217 ms
    fit = estimate_rp(acg, BIN_CENTERS, BS, 0.10)
    assert fit.rp_ms >= BIN_CENTERS[6] * 1000 - 1e-9


@pytest.mark.skipif(not (SPEC / "scratch_macaque_acg.csv").exists(),
                    reason="macaque spec-check data not present")
def test_matches_matlab_on_macaque_acgs():
    """Agreement with matlab/estimate_refractory_period.m on real ACGs.

    The two are not bit-for-bit: this implementation multi-starts the sigmoid
    fit and reaches a lower residual than MATLAB's single-start lsqcurvefit on
    199/200 units, so a handful of units differ by up to ~0.1 ms. The test
    pins the median agreement and the worst case.
    """
    acgs = np.loadtxt(SPEC / "scratch_macaque_acg.csv", delimiter=",")
    ref = np.loadtxt(SPEC / "scratch_macaque_rp_matlab.csv", delimiter=",")
    for frac, col in ((0.10, 2), (0.05, 3), (0.20, 4)):
        rp, _, _ = estimate_rp_many(acgs, BS, frac, BIN_CENTERS)
        both = ~np.isnan(rp) & ~np.isnan(ref[:, col])
        d = np.abs(rp[both] - ref[both, col])
        assert both.mean() > 0.95
        assert np.median(d) < 2e-3, (frac, np.median(d))
        assert np.max(d) < 0.25, (frac, np.max(d))
        assert np.mean(d < 0.01) > 0.80, (frac, np.mean(d < 0.01))
