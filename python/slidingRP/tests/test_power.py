"""Tests for slidingRP.power (tau_pass0 / min_passing_fr) and the appended
``tau_pass0`` output of ``slidingRP``.

These quantities answer "is there enough data for this test to say anything?",
which the pass/fail decision deliberately does not distinguish from "this unit
has violations". They change no decision; they are diagnostics.
"""
import numpy as np
import pytest

from slidingRP import metrics
from slidingRP.power import min_passing_fr, tau_pass0
from slidingRP.simulations import genST


def _brute_force_tau_pass0(n_spikes, rec_dur, cont=10.0, conf=90.0):
    """Smallest tau on the metric's own grid at which zero violations passes."""
    bs = 1 / 30000
    rp = np.arange(0, 10 / 1000, bs) + bs / 2
    ref = rp + bs / 2
    c = 100 * metrics.computeViol(np.zeros_like(ref), n_spikes / rec_dur,
                                  n_spikes, ref, cont / 100, rec_dur)
    ok = np.flatnonzero(c >= conf)
    return ref[ok[0]] if ok.size else np.inf


@pytest.mark.parametrize("fr,dur", [(0.2, 3600), (0.5, 3600), (1.0, 3600),
                                    (2.0, 3600), (5.0, 7200), (1.0, 14400)])
def test_tau_pass0_matches_grid_search(fr, dur):
    n = int(fr * dur)
    analytic = tau_pass0(n, dur)
    grid = _brute_force_tau_pass0(n, dur)
    if np.isfinite(grid):
        # analytic is continuous, the grid search lands on the next bin edge
        assert 0 <= grid - analytic < 1 / 30000 + 1e-12, (analytic, grid)
    else:
        assert analytic > 10 / 1000


def test_tau_pass0_reference_values():
    """Values quoted in the manuscript/response (1 h recording, defaults)."""
    d = 3600
    assert tau_pass0(int(0.2 * d), d) * 1000 == pytest.approx(84.2, rel=0.02)
    assert tau_pass0(int(0.5 * d), d) * 1000 == pytest.approx(13.47, rel=0.02)
    assert tau_pass0(int(1.1 * d), d) * 1000 == pytest.approx(2.78, rel=0.02)
    # a unit that cannot pass at all: tau_pass0 beyond the 10 ms window
    assert tau_pass0(int(0.2 * d), d) > 0.010


def test_tau_pass0_scaling():
    """tau_pass0 falls as 1/(N^2/D): doubling the rate quarters it."""
    d = 3600
    a, b = tau_pass0(int(1.0 * d), d), tau_pass0(int(2.0 * d), d)
    assert a / b == pytest.approx(4.0, rel=0.02)
    # doubling duration at fixed rate halves it
    c = tau_pass0(int(1.0 * 2 * d), 2 * d)
    assert a / c == pytest.approx(2.0, rel=0.02)


def test_min_passing_fr_inverts_tau_pass0():
    for dur in (1800, 3600, 7200):
        for tau in (0.001, 0.002, 0.003, 0.005):
            fr = min_passing_fr(dur, tau)
            assert tau_pass0(fr * dur, dur) == pytest.approx(tau, rel=1e-6)


def test_min_passing_fr_matches_fig4g():
    """Fig 4g: ~1.1 spikes/s at 1 h, 3 ms RP, 90% confidence, 10% contamination."""
    assert min_passing_fr(3600, 0.003) == pytest.approx(1.1, rel=0.06)
    # higher confidence needs more spikes; shorter RP needs more spikes
    assert min_passing_fr(3600, 0.003, conf_thresh=95) > min_passing_fr(3600, 0.003)
    assert min_passing_fr(3600, 0.001) > min_passing_fr(3600, 0.003)


def test_min_passing_fr_agrees_with_simulation():
    """A clean unit just above FR_min passes; just below, it does not."""
    rng = np.random.default_rng(0)
    dur, rp_true = 3600.0, 0.003
    fr_min = min_passing_fr(dur, rp_true)
    for scale, expected in ((1.4, True), (0.6, False)):
        n_pass = 0
        for _ in range(20):
            st = genST(fr_min * scale, dur, rp_true)
            passes = metrics.slidingRP(st, params={'recDur': dur})[5]
            n_pass += bool(passes)
        if expected:
            assert n_pass >= 18, (scale, n_pass)
        else:
            assert n_pass <= 2, (scale, n_pass)


def test_slidingRP_returns_tau_pass0():
    """The appended 8th output matches the standalone function."""
    st = genST(5.0, 3600.0, 0.002)
    out = metrics.slidingRP(st, params={'recDur': 3600.0})
    assert len(out) == 8
    assert out[7] == pytest.approx(tau_pass0(st.size, 3600.0))


def test_slidingRP_all_has_tau_pass0_column():
    st1 = genST(5.0, 1800.0, 0.002)
    st2 = genST(0.3, 1800.0, 0.002)
    times = np.concatenate([st1, st2])
    clusters = np.concatenate([np.zeros(st1.size, int), np.ones(st2.size, int)])
    order = np.argsort(times)
    tbl = metrics.slidingRP_all(times[order], clusters[order],
                                params={'recDur': 1800.0})
    assert 'tau_pass0' in tbl
    assert len(tbl['tau_pass0']) == 2
    # the low-rate unit needs a far longer violation-free window
    assert tbl['tau_pass0'][1] > tbl['tau_pass0'][0]
