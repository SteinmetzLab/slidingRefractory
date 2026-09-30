"""Tests for the censor parameter (sorter duplicate-removal window).

censor = w means the sorter deleted every spike within w of an earlier spike of
the same unit, so only max(tau_r - w, 0) of each window can hold an observed
violation. censor = 0 must reproduce the published method exactly.
"""
import numpy as np
import pytest

from slidingRP import metrics
from slidingRP.power import min_passing_fr, tau_pass0
from slidingRP.simulations import RPmetric_Classic, genST

DUR = 3600.0


def contaminated_train(rate, cont, rp, seed, censor=0.0):
    np.random.seed(seed)
    base = genST(rate * (1 - cont), DUR, rp)
    other = genST(rate * cont, DUR, 0.0) if cont > 0 else np.empty(0)
    st = np.sort(np.concatenate([base, other]))
    if censor > 0:                           # keep a spike only if >= censor after the last kept one
        keep = np.ones(st.size, bool)
        last = -np.inf
        for i, t in enumerate(st):
            if t - last < censor:
                keep[i] = False
            else:
                last = t
        st = st[keep]
    return st


def test_default_is_censor_zero():
    st = contaminated_train(5.0, 0.10, 0.003, 1)
    a = metrics.slidingRP(st, params={"recDur": DUR})
    b = metrics.slidingRP(st, params={"recDur": DUR, "censor": 0.0})
    for x, y in zip(a, b):
        assert (np.isnan(x) and np.isnan(y)) or x == y


def test_censor_uses_the_observable_window():
    """The confidence at each window equals the uncensored formula evaluated at
    tau_r - w, and a window inside the censor carries no confidence."""
    w = 0.0005
    st = contaminated_train(8.0, 0.05, 0.002, 2, censor=w)
    cm, cont, rp, nacg, _ = metrics.computeMatrix(st, {"recDur": DUR, "censor": w})
    ref_dur = rp + (rp[1] - rp[0]) / 2
    obs = np.cumsum(nacg)
    i = int(np.flatnonzero(np.isclose(cont, 10))[0])
    expect = 100 * metrics.computeViol(obs, None, st.size, np.clip(ref_dur - w, 0, None),
                                       0.10, DUR)
    np.testing.assert_allclose(cm[i], expect, atol=1e-12)
    assert np.all(cm[:, ref_dur <= w] == 0)


def test_matrix_and_scalar_paths_agree_with_censor():
    w = 0.00025
    st = contaminated_train(5.0, 0.08, 0.003, 3, censor=w)
    p = {"recDur": DUR, "censor": w}
    max_conf = metrics.slidingRP(st, params=p)[0]
    cm, cont, rp, _, _ = metrics.computeMatrix(st, p)
    i = int(np.flatnonzero(np.isclose(cont, 10))[0])
    assert max_conf == pytest.approx(cm[i, rp > 0.0005].max(), abs=1e-9)


def test_power_functions_shift_by_the_censor():
    n, w = 18000, 0.00025
    assert tau_pass0(n, DUR, censor=w) == pytest.approx(tau_pass0(n, DUR) + w)
    assert min_passing_fr(DUR, 0.003, censor=w) == pytest.approx(min_passing_fr(DUR, 0.003 - w))
    assert np.isinf(min_passing_fr(DUR, 0.0002, censor=w))
    st = contaminated_train(5.0, 0.0, 0.003, 4, censor=w)
    out = metrics.slidingRP(st, params={"recDur": DUR, "censor": w})
    assert out[7] == pytest.approx(tau_pass0(st.size, DUR) + w)


def test_hill_llobet_censor_shortens_the_window():
    st = contaminated_train(5.0, 0.10, 0.003, 5, censor=0.0005)
    p = {"recDur": DUR, "RPdur": 0.003}
    pass0, est0 = RPmetric_Classic(st, p)
    passc, estc = RPmetric_Classic(st, dict(p, censor=0.0005))
    # the same observed count against a shorter expected window: a higher estimate
    assert estc > est0
    assert RPmetric_Classic(st, dict(p, censor=0.0)) == (pass0, est0)


def test_ignoring_a_censor_is_anti_conservative_and_the_correction_fixes_it():
    """Units at 12% contamination (above the 10% threshold) with a 0.5 ms
    censor: the published method accepts essentially all of them; accounting
    for the censor brings acceptance back near the uncensored rate (about 7% at
    this setting, 5 spikes/s, 3 ms, 1 h)."""
    w, n_sim = 0.0005, 150
    ignored = corrected = 0
    for k in range(n_sim):
        st = contaminated_train(5.0, 0.12, 0.003, 100 + k, censor=w)
        ignored += metrics.slidingRP(st, params={"recDur": DUR})[5]
        corrected += metrics.slidingRP(st, params={"recDur": DUR, "censor": w})[5]
    assert ignored / n_sim > 0.9, ignored
    assert corrected / n_sim < 0.25, corrected
