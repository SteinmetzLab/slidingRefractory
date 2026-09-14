"""Statistical power of the Sliding RP test: how much data is enough?

A unit can fail the Sliding RP test for two quite different reasons: refractory
period violations were observed that are incompatible with acceptable
contamination, or too few spikes were recorded to establish acceptable
contamination at all. The metric's pass/fail decision does not distinguish
them, and deliberately so: what counts as a plausible refractory window is the
user's assumption, not the method's. But the distinction matters for
interpretation, and the quantities here make it explicit.

The existing ``nViolShort`` / ``n_spikes_below2`` output already separates the
two cases empirically (a failing unit with zero violations below 2 ms failed
for want of power, not because violations were seen). The functions here give
the analytical counterpart, so a user can tell in advance whether a given
firing rate and recording duration can support the conclusion at all.

Under the Llobet model the expected violation count up to tau for a unit with
N spikes recorded over D seconds at contamination C is

    Ve(tau) = 2 * tau * Nc * (Nb + (Nc - 1) / 2) / D,   Nc = C*N, Nb = (1-C)*N

so a unit with *zero* observed violations up to tau passes at confidence gamma
iff ``1 - exp(-Ve(tau)) >= gamma``. Both functions below invert that relation.
"""
from __future__ import annotations

import numpy as np


def tau_pass0(n_spikes, rec_dur, cont_thresh=10.0, conf_thresh=90.0):
    """Shortest violation-free window (s) that would let a unit pass.

    ``tau_pass0 = -ln(1 - gamma) * D / (2 * Nc * (Nb + (Nc - 1)/2))``

    A unit whose ``tau_pass0`` exceeds the longest tested refractory period
    (10 ms by default) cannot pass however clean its autocorrelogram is: there
    is not enough data. One whose ``tau_pass0`` is below a refractory duration
    the user considers plausible has enough data for the test to be meaningful.

    Parameters
    ----------
    n_spikes : int or array_like
        Total spike count of the unit.
    rec_dur : float or array_like
        Recording duration in seconds.
    cont_thresh : float
        Maximum acceptable contamination (%), the C_thresh of the metric.
    conf_thresh : float
        Required confidence (%), the gamma_thresh of the metric.

    Returns
    -------
    float or ndarray
        tau_pass0 in seconds; ``inf`` when the unit has no spikes.

    Examples
    --------
    At the defaults (10% contamination, 90% confidence) over a 1 h recording,
    a 0.2 spikes/s unit needs a violation-free window of ~90 ms, a 0.5
    spikes/s unit ~14 ms, and a 1.1 spikes/s unit ~3 ms.
    """
    n = np.asarray(n_spikes, dtype=np.float64)
    d = np.asarray(rec_dur, dtype=np.float64)
    c = cont_thresh / 100.0
    nc, nb = n * c, n * (1 - c)
    denom = 2 * nc * (nb + (nc - 1) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = -np.log(1 - conf_thresh / 100.0) * d / denom
    out = np.where(denom > 0, out, np.inf)
    return float(out) if np.ndim(out) == 0 else out


def min_passing_fr(rec_dur, tau, cont_thresh=10.0, conf_thresh=90.0):
    """Minimum firing rate (spikes/s) for a violation-free unit to pass at tau.

    Exact inverse of ``Ve(tau) = -ln(1 - gamma)``, which is a quadratic in the
    spike count N:

        [C(1-C) + C^2/2] N^2 - (C/2) N - ln(1/(1-gamma)) D / (2 tau) = 0

    Returns ``N / D``. This is the analytical form of manuscript Fig 4g, which
    was previously obtained by simulation.
    """
    d = np.asarray(rec_dur, dtype=np.float64)
    t = np.asarray(tau, dtype=np.float64)
    c = cont_thresh / 100.0
    k = -np.log(1 - conf_thresh / 100.0) * d / (2 * t)
    a = c * (1 - c) + c * c / 2
    b = -c / 2
    n = (-b + np.sqrt(b * b + 4 * a * k)) / (2 * a)
    out = n / d
    return float(out) if np.ndim(out) == 0 else out


def min_passing_spikes(tau, cont_thresh=10.0, conf_thresh=90.0):
    """Minimum spike count for a violation-free unit to pass at tau.

    Note this depends on ``tau/D`` only through ``Ve``, so unlike
    :func:`min_passing_fr` it is *not* a function of recording duration alone;
    it returns N for a nominal D = 1 s and must be scaled by D.
    """
    return min_passing_fr(1.0, tau, cont_thresh, conf_thresh)
