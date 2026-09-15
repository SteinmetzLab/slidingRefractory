"""Validation tests for the resubmission analysis code.

These guard the two re-expressions the analyses depend on:
  * ``slidingRP_from_acg`` must reproduce ``slidingRP.metrics.slidingRP``;
  * ``corrected_confidence_from_acg`` must reproduce
    ``computeMatrix(..., correction=True)``;
and the statistical properties of the new generators.

Run with:  pytest roth-et-al-2026/resubmission/test_resub.py -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import (BIN_SIZE, N_BINS, RP_CENTERS, computeACG_samples,  # noqa: E402
                       hill_llobet_from_acg, slidingRP_from_acg)
from resub_simulations import (corrected_confidence_from_acg, evaluate_arms,  # noqa: E402
                               gen_bursting, gen_graded_rp, gen_hard_rp,
                               gen_modulated_pair, gen_nonoverlapping,
                               make_train, poisson_test_fixed_tau)

from slidingRP import metrics  # noqa: E402

RNG = np.random.default_rng(20260913)


@pytest.mark.parametrize("rate,rp,cont", [(5.0, 0.002, 0.0), (5.0, 0.002, 0.10),
                                          (1.0, 0.003, 0.05), (20.0, 0.0015, 0.15)])
def test_slidingRP_from_acg_matches_package(rate, rp, cont):
    """The precomputed-ACG path must equal the package's spike-time path."""
    dur = 3600.0
    st, _ = make_train("standard", rate, cont, dur, rp, rng=RNG)
    ref = metrics.slidingRP(st, params={"recDur": dur})
    nACG = metrics.computeACG(st, BIN_SIZE, N_BINS)
    got = slidingRP_from_acg(nACG, st.size, dur)
    assert got["max_conf"] == pytest.approx(ref[0], rel=1e-10, abs=1e-10)
    assert got["passes"] == ref[5]
    assert got["n_viol_short"] == ref[3]
    assert got["tau_pass0"] == pytest.approx(ref[7], rel=1e-12)
    if np.isnan(ref[1]):
        assert np.isnan(got["min_cont"])
    else:
        assert got["min_cont"] == pytest.approx(ref[1], rel=1e-9)
        assert got["rp_min_val"] == pytest.approx(ref[2], rel=1e-9)


def test_corrected_confidence_matches_full_matrix():
    """One-row FWER correction == the full-grid correction at that row."""
    dur = 7200.0
    st, _ = make_train("standard", 5.0, 0.10, dur, 0.002, rng=RNG)
    st = st[:8000]
    conf_mat, cont, rp, _, _ = metrics.computeMatrix(
        st, {"recDur": dur, "correction": True, "cont": np.array([10.0])})
    ref = float(conf_mat[0, 0])
    got = corrected_confidence_from_acg(
        metrics.computeACG(st, BIN_SIZE, N_BINS), st.size, dur, 10.0)
    assert got == pytest.approx(ref, rel=1e-9, abs=1e-9)


def test_corrected_is_more_conservative_than_nominal():
    """The FWER correction can only lower confidence."""
    dur = 3600.0
    for cont in (0.0, 0.08, 0.10, 0.12):
        st, _ = make_train("standard", 5.0, cont, dur, 0.002, rng=RNG)
        nACG = metrics.computeACG(st, BIN_SIZE, N_BINS)
        nominal = slidingRP_from_acg(nACG, st.size, dur)["max_conf"]
        corrected = corrected_confidence_from_acg(nACG, st.size, dur)
        assert corrected <= nominal + 1e-9, (cont, nominal, corrected)


def test_computeACG_samples_matches_integer_ground_truth():
    """Integer-lag ACG equals a brute-force pairwise count."""
    s = np.sort(RNG.integers(0, 300000, 4000))
    got = computeACG_samples(s, N_BINS)
    d = s[:, None] - s[None, :]
    d = d[(d > 0) & (d < N_BINS)]
    want = np.bincount(d, minlength=N_BINS)[:N_BINS]
    assert np.array_equal(got, want)


def test_hill_llobet_matches_matlab_bin_convention():
    """obsViol uses MATLAB's inclusive sum(nACG(1:find(rp>RPdur,1)))."""
    st, _ = make_train("standard", 5.0, 0.10, 3600.0, 0.002, rng=RNG)
    nACG = metrics.computeACG(st, BIN_SIZE, N_BINS)
    _, _, obs = hill_llobet_from_acg(nACG, st.size, 3600.0, 0.002)
    idx = int(np.argmax(RP_CENTERS > 0.002))
    assert obs == float(np.sum(nACG[:idx + 1]))


def test_poisson_test_fixed_tau_is_calibrated():
    """A single fixed-tau test at the true RP has the nominal false-acceptance
    rate at the contamination threshold (no multiplicity)."""
    dur, rp, gamma = 7200.0, 0.003, 90.0
    n_pass = 0
    n_sim = 300
    for _ in range(n_sim):
        st, _ = make_train("standard", 10.0, 0.10, dur, rp, rng=RNG)
        nACG = metrics.computeACG(st, BIN_SIZE, N_BINS)
        c = poisson_test_fixed_tau(nACG, st.size, dur, rp)
        n_pass += c >= gamma
    far = n_pass / n_sim
    # nominal 10%; binomial SE at n=300 is 1.7 points, allow 3.5 SE
    assert 0.04 < far < 0.16, far


# --- generators -----------------------------------------------------------

def test_gen_hard_rp_rate_and_refractoriness():
    st = gen_hard_rp(10.0, 600.0, 0.003, RNG)
    assert st.size / 600.0 == pytest.approx(10.0, rel=0.05)
    assert np.min(np.diff(st)) >= 0.003 - 1e-12


def test_gen_graded_rp_is_absolute_then_graded():
    """Hard below `rp`, then suppressed but non-zero over the ramp.

    This is the shape Nick described: a unit hard to 2 ms and then graded to
    3 ms. An earlier version of this generator centred a logistic hazard *on*
    `rp`, which left 27% of baseline firing at zero lag and so had no absolute
    refractory period at all; that made graded recovery look catastrophic for
    the metric when it is in fact a non-issue. The assertions below are what
    distinguishes the two.
    """
    hard = gen_hard_rp(10.0, 900.0, 0.002, RNG)
    graded = gen_graded_rp(10.0, 900.0, 0.002, 0.001, RNG)
    assert np.sum(np.diff(hard) < 0.002) == 0
    assert np.sum(np.diff(graded) < 0.002) == 0           # absolute part
    isi = np.diff(graded)
    ramp = np.mean((isi >= 0.002) & (isi < 0.003))        # suppressed ramp
    after = np.mean((isi >= 0.003) & (isi < 0.004))       # full hazard
    assert 0 < ramp < after, (ramp, after)
    assert graded.size / 900.0 == pytest.approx(10.0, rel=0.15)


def test_gen_bursting_makes_short_isis():
    st = gen_bursting(10.0, 900.0, 0.002, p_burst=0.3, rng=RNG)
    isi = np.diff(st)
    assert np.mean((isi > 0.003) & (isi < 0.006)) > 0.1
    assert np.min(isi) >= 0.002 - 1e-12
    assert st.size / 900.0 == pytest.approx(10.0, rel=0.15)


@pytest.mark.parametrize("rho", [-1.0, 0.0, 1.0])
def test_gen_modulated_pair_sign_of_correlation(rho):
    b, c, r = gen_modulated_pair(4.5, 0.5, 1800.0, 0.002, rho=rho, rng=RNG)
    assert b.size / 1800.0 == pytest.approx(4.5, rel=0.2)
    if rho > 0:
        assert r > 0.05, r
    elif rho < 0:
        assert r < -0.05, r
    else:
        assert abs(r) < 0.1, r


@pytest.mark.parametrize("overlap", [0.0, 1.0])
def test_gen_nonoverlapping(overlap):
    b, c = gen_nonoverlapping(4.5, 0.5, 1800.0, 0.002, overlap=overlap, rng=RNG)
    edges = np.arange(0, 1801, 20.0)
    cb = np.histogram(b, edges)[0] > 0
    cc = np.histogram(c, edges)[0] > 0
    shared = np.mean(cc[cb]) if cb.any() else 0
    if overlap == 1.0:
        assert shared > 0.9, shared
    else:
        assert shared < 0.2, shared


def test_nonoverlapping_hides_contamination():
    """The failure mode the manuscript describes but never tested: with no
    temporal overlap a heavily contaminated unit passes."""
    dur, rp = 3600.0, 0.002
    n_pass_overlap, n_pass_disjoint = 0, 0
    for _ in range(40):
        st, _ = make_train("nonoverlapping", 5.0, 0.20, dur, rp, rng=RNG, overlap=1.0)
        n_pass_overlap += evaluate_arms(st, dur)["sliding_pass_90"]
        st, _ = make_train("nonoverlapping", 5.0, 0.20, dur, rp, rng=RNG, overlap=0.0)
        n_pass_disjoint += evaluate_arms(st, dur)["sliding_pass_90"]
    assert n_pass_overlap <= 4, n_pass_overlap          # 20% contamination: rejected
    assert n_pass_disjoint >= 30, n_pass_disjoint       # hidden by disjoint activity


def test_make_train_realised_contamination():
    st, info = make_train("standard", 10.0, 0.10, 1800.0, 0.002, rng=RNG)
    assert info["realised_cont"] == pytest.approx(0.10, abs=0.02)
    assert st.size / 1800.0 == pytest.approx(10.0, rel=0.1)


def test_evaluate_arms_is_self_consistent():
    dur = 3600.0
    st, _ = make_train("standard", 5.0, 0.05, dur, 0.002, rng=RNG)
    out = evaluate_arms(st, dur, true_rp=0.002, with_estimator=True,
                        with_correction=True)
    assert out["sliding_pass_90"] == (out["sliding_max_conf"] >= 90)
    # gamma sweep must be monotone: passing at 95 implies passing at 90
    assert not (out["sliding_pass_95"] and not out["sliding_pass_90"])
    assert out["corrected_conf"] <= out["sliding_max_conf"] + 1e-9
    assert np.isfinite(out["tau_pass0"])


def test_rp_reject_argument_is_honoured():
    """Regression: slidingRP_from_acg once ignored rp_reject and used the module
    default, which silently turned the tau_min sensitivity sweep into a no-op."""
    dur = 3600.0
    st, _ = make_train("standard", 5.0, 0.08, dur, 0.0015, rng=RNG)
    nACG = metrics.computeACG(st, BIN_SIZE, N_BINS)
    confs = [slidingRP_from_acg(nACG, st.size, dur, rp_reject=t)["max_conf"]
             for t in (0.00025, 0.0005, 0.001, 0.002, 0.004)]
    # raising tau_min can only remove candidate windows, so confidence is
    # non-increasing, and over this range it must actually change
    assert all(a >= b - 1e-9 for a, b in zip(confs, confs[1:])), confs
    assert confs[0] - confs[-1] > 1e-6, confs
    # and it must still agree with the package at the package's default
    assert slidingRP_from_acg(nACG, st.size, dur)["max_conf"] == pytest.approx(
        metrics.slidingRP(st, params={"recDur": dur})[0], rel=1e-10)


# --- the two candidate timepoints added for the revision --------------------

def test_tau_last_pass_matches_the_package_sweep():
    """The vectorised last-accepted-tau must agree with a per-unit reference.

    ``add_tau_last`` recomputes the whole confidence sweep in one matrix
    operation over many units at once. This checks that against the package's
    own single-unit path, both for the acceptance flag and for the first and
    last accepted windows.
    """
    from add_tau_last import tau_last_pass

    dur = 3600.0
    acgs, ns = [], []
    for rate, cont, rp in [(5.0, 0.0, 0.002), (5.0, 0.10, 0.002),
                           (1.0, 0.05, 0.003), (20.0, 0.15, 0.0015),
                           (12.0, 0.0, 0.0025)]:
        st, _ = make_train("standard", rate, cont, dur, rp, rng=RNG)
        acgs.append(metrics.computeACG(st, BIN_SIZE, N_BINS))
        ns.append(st.size)
    acgs = np.vstack(acgs)
    ns = np.array(ns, float)
    durs = np.full(ns.shape, dur)

    last, first, acc, nbins = tau_last_pass(acgs, ns, durs)
    for i in range(acgs.shape[0]):
        ref = slidingRP_from_acg(acgs[i], ns[i], dur)
        assert acc[i] == ref["passes"]
        if ref["passes"]:
            assert first[i] == pytest.approx(ref["tau_first_pass"], abs=1e-12)
            assert last[i] >= first[i]
            assert 1 <= nbins[i] <= RP_CENTERS.size
        else:
            assert np.isnan(first[i]) and np.isnan(last[i])
            assert nbins[i] == 0


def test_tau_last_pass_is_the_last_not_the_end_of_the_first_run():
    """Constructed case: an accepted window, a gap, then another accepted one.

    The accepted set of tau_r need not be contiguous. This pins the documented
    behaviour -- the last element, not the end of the first run -- so a future
    change to that has to be deliberate.
    """
    from add_tau_last import tau_last_pass

    dur, n = 3600.0, 40000.0
    # Clean out to 2 ms, then one bin holding 300 violations. Just after that
    # bin the expectation is far below 300 so the unit is rejected; by 10 ms the
    # expectation has grown past 800 and it is accepted again. The accepted set
    # is therefore two runs, and the answer must come from the second.
    acg = np.zeros((1, N_BINS))
    acg[0, RP_CENTERS.searchsorted(0.002)] = 300
    last, first, acc, nbins = tau_last_pass(acg, np.array([n]), np.array([dur]))
    assert acc[0]
    assert first[0] == pytest.approx(RP_CENTERS[RP_CENTERS > 0.0005][0])
    # the first run ends at the burst, well before the answer
    assert last[0] > 0.009
    # and the accepted set really is broken in two, not one long run
    assert nbins[0] < np.sum(RP_CENTERS > 0.0005)


def test_standardization_is_identity_without_a_firing_rate_effect():
    """If the value does not depend on firing rate, standardizing changes nothing."""
    from fr_standardize import reference_weights, standardized_median

    rng = np.random.default_rng(0)
    fr = np.exp(rng.normal(np.log(8), 0.7, 40000))
    v = rng.gamma(4, 0.6, fr.size)               # independent of fr
    w = reference_weights(fr)
    std, cov = standardized_median(fr, v, w)
    assert cov == pytest.approx(1.0)
    assert std == pytest.approx(np.median(v), abs=0.02)


def test_standardization_removes_a_known_firing_rate_effect():
    """Two groups with the same conditional means but different rate mixes.

    Raw medians differ only because the groups sample firing rate differently;
    standardization must bring them back together. It does not close the gap
    completely, and that residual is a real property of the method rather than a
    bug: within a bin the low-rate group still sits at the low edge and the
    high-rate group at the high edge, so a coarse grid leaves some of the effect
    behind. Here it removes 88% of a 0.69 ms gap.
    """
    from fr_standardize import reference_weights, standardized_median

    rng = np.random.default_rng(1)

    def group(n, log_mu):
        fr = np.exp(rng.normal(log_mu, 0.5, n))
        v = 4.0 - 0.5 * np.log(fr) + rng.normal(0, 0.3, n)   # same law in both
        return fr, v

    fr_a, v_a = group(30000, np.log(4))
    fr_b, v_b = group(30000, np.log(16))
    w = reference_weights(np.concatenate([fr_a, fr_b]))
    raw_gap = abs(np.median(v_a) - np.median(v_b))
    std_a = standardized_median(fr_a, v_a, w)[0]
    std_b = standardized_median(fr_b, v_b, w)[0]
    assert raw_gap > 0.4, raw_gap
    assert abs(std_a - std_b) < 0.2 * raw_gap
