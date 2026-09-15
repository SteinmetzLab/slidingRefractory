"""Simulation engine for the J Neurophysiol resubmission.

Three things live here:

1. **Generators** beyond the manuscript's hard-RP-Poisson base plus RP-free
   Poisson contamination (work package 05): a contaminating neuron with its own
   refractory period, correlated and non-stationary rates, non-overlapping
   ("place cell") activity, graded refractory recovery, and bursting.

2. **Arms** (work package 04): the several ways one can turn an ACG into a
   pass/fail decision, all evaluated on the *same* spike train so differences
   are not simulation noise.

3. **A runner** that sweeps parameters and returns a tidy table.

Nothing here changes the published metric. The Sliding RP arm calls the
package's own code path (via the precomputed-ACG re-expression in
``acg_table``, which is checked against ``slidingRP.metrics.slidingRP`` by
``test_acg_table.py``).
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import (BIN_SIZE, N_BINS, REF_DUR, RP_CENTERS, TEST_TIMES,  # noqa: E402
                       computeViol, hill_llobet_from_acg, slidingRP_from_acg)

from slidingRP.metrics import computeACG, _fwer_correct  # noqa: E402
from slidingRP.rp_estimate import estimate_rp  # noqa: E402


# --------------------------------------------------------------------------
# 1. Generators
# --------------------------------------------------------------------------

def gen_hard_rp(rate, duration, rp=0.0, rng=None):
    """Renewal process: ISI = rp + Exponential(mu), rate-corrected for the dead
    time so the realised rate is `rate`. Identical to the manuscript's genST."""
    rng = rng or np.random.default_rng()
    if rate <= 0:
        return np.empty(0)
    r_sim = rate / (1 - rp * rate)
    mu = 1 / r_sim
    n = max(int(np.ceil(rate * duration * 2)), 16)
    isi = rp + rng.exponential(mu, n)
    st = np.cumsum(isi)
    while st.size and st[-1] < duration:
        extra = rp + rng.exponential(mu, n)
        st = np.concatenate([st, st[-1] + np.cumsum(extra)])
    return st[st < duration]


@lru_cache(maxsize=256)
def _graded_survival(rate, rp, width, dt=2e-5):
    """ISI survival function for absolute-then-relative refractoriness.

    The hazard is exactly zero below `rp` (the absolute refractory period),
    then rises linearly to its asymptote over `width` (the relative refractory
    period), and is flat thereafter. The asymptote is calibrated so the
    realised mean rate equals `rate`. ``width = 0`` reduces to the hard-RP
    renewal process used elsewhere.

    Cached on the parameters, since the calibration is the expensive part and
    does not depend on the random draws.

    Returns (t, S) with S decreasing, ready for inverse-transform sampling.
    """
    t = np.arange(0, max(20 / rate, rp + 20 * max(width, 1e-4)), dt)
    shape = np.clip((t - rp) / width, 0.0, 1.0) if width > 0 else (t >= rp).astype(float)
    lam = 1.0
    for _ in range(60):
        S = np.exp(-np.cumsum(lam * shape) * dt)
        mean_isi = np.trapezoid(S, t)
        if mean_isi <= 0:
            break
        new_lam = lam * mean_isi * rate
        if abs(new_lam - lam) < 1e-9 * lam:
            lam = new_lam
            break
        lam = new_lam
    S = np.exp(-np.cumsum(lam * shape) * dt)
    return t, S


def gen_graded_rp(rate, duration, rp=0.002, width=0.001, rng=None):
    """Renewal process with an absolute refractory period then a graded recovery.

    Spiking is impossible for `rp` after each spike; the hazard then rises
    linearly to its asymptote over the next `width` seconds. This is the
    "hard to 2 ms, then graded to 3 ms" shape, not a soft suppression that
    still permits spikes at zero lag.
    """
    rng = rng or np.random.default_rng()
    if rate <= 0:
        return np.empty(0)
    t, S = _graded_survival(float(rate), float(rp), float(width))
    n = max(int(np.ceil(rate * duration * 1.5)), 16)
    st = np.empty(0)
    last = 0.0
    while last < duration:
        u = rng.random(n)
        isi = np.interp(u, S[::-1], t[::-1])       # S is decreasing
        new = last + np.cumsum(isi)
        st = np.concatenate([st, new])
        last = new[-1]
    return st[st < duration]


def _enforce_dead_time(st, rp):
    """Greedily drop spikes closer than `rp` to the previously kept spike.

    A single neuron cannot violate its own absolute refractory period, however
    its spikes were generated; without this a burst partner inserted after one
    spike can land arbitrarily close to the next.

    Implemented as repeated vectorised passes: each pass removes the first
    spike of every too-close pair, which is equivalent to the sequential greedy
    rule and converges in a couple of passes because violations are sparse.
    """
    if rp <= 0 or st.size < 2:
        return st
    st = np.sort(st)
    while True:
        gaps = np.diff(st)
        bad = np.flatnonzero(gaps < rp)
        if bad.size == 0:
            return st
        # drop the later spike of each violating pair, taking every other index
        # so that a run of close spikes is thinned rather than emptied
        drop = bad[::2] + 1
        st = np.delete(st, drop)


def gen_bursting(rate, duration, rp=0.002, p_burst=0.2, burst_isi=(0.003, 0.006),
                 rng=None):
    """Hard-RP renewal process where each spike may be followed by a burst
    partner at a short ISI, producing the short-latency ACG peak seen in many
    real units.

    The merged train is passed through the neuron's own dead time, so the
    absolute refractory period is still respected (a burst partner that would
    land within `rp` of the next regular spike removes that spike).
    """
    rng = rng or np.random.default_rng()
    if rate <= 0:
        return np.empty(0)
    # generate slightly hot, since dead-time enforcement removes a few spikes
    base = gen_hard_rp(rate / (1 + p_burst), duration, rp, rng)
    if base.size == 0:
        return base
    take = rng.random(base.size) < p_burst
    extra = base[take] + rng.uniform(*burst_isi, size=int(take.sum()))
    st = _enforce_dead_time(np.concatenate([base, extra]), rp)
    return st[st < duration]


def _modulation(duration, tau_s, rng, dt=0.05):
    """Zero-mean, unit-variance low-pass Gaussian modulation signal.

    An AR(1) process with correlation time tau_s, generated with lfilter rather
    than a Python loop: at dt = 50 ms and a 2 h recording the loop version ran
    144,000 iterations per simulated train and dominated the whole sweep.
    """
    from scipy import signal as _sig

    n = int(np.ceil(duration / dt)) + 1
    x = rng.standard_normal(n)
    a = np.exp(-dt / tau_s)
    y = _sig.lfilter([np.sqrt(1 - a * a)], [1.0, -a], x)
    return y, dt


def gen_modulated_pair(base_rate, cont_rate, duration, rp, rho=0.0, amp=0.8,
                       tau_s=2.0, rng=None):
    """Base neuron and contaminant driven by a shared slow rate modulation.

    ``r_b(t) = R_b (1 + amp*s(t))`` and ``r_c(t) = R_c (1 + rho*amp*s(t))``.
    ``rho = 0`` is the manuscript's independent case, ``rho = +1`` perfectly
    positively correlated, ``rho = -1`` anti-correlated. Rates are clipped at
    zero, then the realised total rates are renormalised to the targets.

    Implemented by thinning: generate at the maximum rate, keep each spike with
    probability r(t)/r_max. The base train is generated with its hard RP before
    thinning, which slightly lengthens the effective dead time; that is the
    intended behaviour (a modulated neuron still cannot fire within its RP).

    Returns (base_spikes, contaminant_spikes, realised_rate_correlation).
    """
    rng = rng or np.random.default_rng()
    s, dt = _modulation(duration, tau_s, rng)

    def thin(rate, rp_, weight):
        if rate <= 0:
            return np.empty(0)
        prof = np.clip(1 + weight * amp * s, 0, None)
        prof = prof / prof.mean()
        rmax = prof.max()
        st = gen_hard_rp(rate * rmax, duration, rp_, rng)
        if st.size == 0:
            return st
        p = prof[np.minimum((st / dt).astype(int), prof.size - 1)] / rmax
        return st[rng.random(st.size) < p]

    base = thin(base_rate, rp, 1.0)
    cont = thin(cont_rate, 0.0, rho)
    # realised correlation of the two rate profiles in 100 ms bins
    edges = np.arange(0, duration + 0.1, 0.1)
    cb = np.histogram(base, edges)[0].astype(float)
    cc = np.histogram(cont, edges)[0].astype(float)
    r = np.corrcoef(cb, cc)[0, 1] if cb.std() > 0 and cc.std() > 0 else np.nan
    return base, cont, r


def gen_nonoverlapping(base_rate, cont_rate, duration, rp, overlap=0.0,
                       block_s=20.0, active_frac=0.5, rng=None):
    """Base and contaminant active in partially overlapping epochs.

    The base neuron is active in a random ``active_frac`` of blocks. The
    contaminant is active in a set of blocks the same size, chosen so that a
    fraction ``overlap`` of them coincide with the base's active blocks.
    ``overlap = 1`` reproduces the stationary case (both always on together);
    ``overlap = 0`` is the hippocampal place-cell failure mode the manuscript
    describes but does not test. Total spike counts are held at the targets so
    the contamination fraction is exactly as specified.
    """
    rng = rng or np.random.default_rng()
    n_blocks = max(int(np.ceil(duration / block_s)), 2)
    n_active = max(int(round(active_frac * n_blocks)), 1)
    base_blocks = rng.choice(n_blocks, n_active, replace=False)
    n_shared = int(round(overlap * n_active))
    shared = rng.choice(base_blocks, min(n_shared, n_active), replace=False)
    others = np.setdiff1d(np.arange(n_blocks), base_blocks)
    n_extra = n_active - shared.size
    extra = rng.choice(others, min(n_extra, others.size), replace=False) \
        if n_extra > 0 and others.size else np.empty(0, dtype=int)
    cont_blocks = np.concatenate([shared, extra]).astype(int)

    def in_blocks(rate, rp_, blocks):
        """Generate at an inflated rate, keep only spikes inside the blocks."""
        if rate <= 0 or blocks.size == 0:
            return np.empty(0)
        frac = blocks.size / n_blocks
        st = gen_hard_rp(rate / frac, duration, rp_, rng)
        blk = np.minimum((st / block_s).astype(int), n_blocks - 1)
        return st[np.isin(blk, blocks)]

    return in_blocks(base_rate, rp, base_blocks), in_blocks(cont_rate, 0.0, cont_blocks)


def make_train(model, total_rate, cont_prop, duration, rp, rng=None, **kw):
    """Build one contaminated spike train under the named model.

    Returns (spike_times, info dict). The contamination proportion is the
    fraction of the *total* train contributed by the contaminating source, as
    in the manuscript.
    """
    rng = rng or np.random.default_rng()
    base_rate = (1 - cont_prop) * total_rate
    cont_rate = cont_prop * total_rate
    info = {}

    if model == "standard":                       # manuscript's model
        b = gen_hard_rp(base_rate, duration, rp, rng)
        c = gen_hard_rp(cont_rate, duration, 0.0, rng)
    elif model == "single_neuron_contaminant":    # 05-A
        b = gen_hard_rp(base_rate, duration, rp, rng)
        c = gen_hard_rp(cont_rate, duration, kw.get("cont_rp", 0.0015), rng)
    elif model == "modulated":                    # 05-B1
        b, c, r = gen_modulated_pair(base_rate, cont_rate, duration, rp,
                                     rho=kw.get("rho", 0.0),
                                     amp=kw.get("amp", 0.8),
                                     tau_s=kw.get("tau_s", 2.0), rng=rng)
        info["realised_rate_corr"] = r
    elif model == "nonoverlapping":               # 05-B2
        b, c = gen_nonoverlapping(base_rate, cont_rate, duration, rp,
                                  overlap=kw.get("overlap", 0.0),
                                  block_s=kw.get("block_s", 20.0), rng=rng)
    elif model == "graded":                       # 05-C1
        b = gen_graded_rp(base_rate, duration, rp, kw.get("width", 0.001), rng)
        c = gen_hard_rp(cont_rate, duration, 0.0, rng)
    elif model == "bursting":                     # 05-C2
        b = gen_bursting(base_rate, duration, rp, kw.get("p_burst", 0.2), rng=rng)
        c = gen_hard_rp(cont_rate, duration, 0.0, rng)
    else:
        raise ValueError(f"unknown model {model!r}")

    st = np.sort(np.concatenate([b, c]))
    info["n_base"], info["n_cont"] = b.size, c.size
    info["realised_cont"] = c.size / max(st.size, 1)
    return st, info


# --------------------------------------------------------------------------
# 2. Arms
# --------------------------------------------------------------------------

def corrected_confidence_from_acg(nACG, n_spikes, rec_dur, cont_thresh=10.0,
                                  rp_reject=0.0005):
    """FWER-corrected confidence at a single contamination level.

    Calls the package's own ``_fwer_correct`` with a one-row matrix, so the
    result is identical to ``computeMatrix(..., correction=True)`` at that row
    but ~70x cheaper (the full grid tests 70 contamination levels and the
    pass/fail decision needs only C_thresh).
    """
    nACG = np.asarray(nACG, dtype=np.float64)
    obs = np.cumsum(nACG)
    exp_viol = computeViol(obs, n_spikes, REF_DUR, cont_thresh / 100, rec_dur)[1]
    nominal = 100 * (1 - stats.poisson.cdf(obs, exp_viol))
    corrected = _fwer_correct(nominal[np.newaxis, :], exp_viol[np.newaxis, :],
                              obs, RP_CENTERS, rp_reject)
    return float(corrected[0, 0])


def poisson_test_fixed_tau(nACG, n_spikes, rec_dur, tau, cont_thresh=10.0):
    """The Sliding RP statistic at ONE pre-specified tau_r (no sliding).

    Isolates the statistical treatment (Poisson confidence rather than a point
    estimate) from the benefit of searching over tau_r.
    """
    if not np.isfinite(tau) or tau <= 0:
        return np.nan
    idx = int(np.argmin(np.abs(RP_CENTERS - tau)))
    obs = float(np.sum(np.asarray(nACG)[:idx + 1]))
    conf, _ = computeViol(obs, n_spikes, REF_DUR[idx], cont_thresh / 100, rec_dur)
    return 100 * float(conf)


def evaluate_arms(st, rec_dur, true_rp=None, cont_thresh=10.0,
                  gammas=(50, 60, 70, 75, 80, 85, 90, 95, 99),
                  with_correction=False, with_estimator=False,
                  hl_rps=(0.002, 0.003)):
    """Evaluate every decision rule on one spike train.

    Returns a flat dict. Confidences are returned rather than pass/fail
    wherever possible, so a whole gamma sweep costs nothing extra.
    """
    n = st.size
    out = {"n_spikes": n, "rec_dur": rec_dur, "firing_rate": n / rec_dur}
    if n < 2:
        out.update({"sliding_max_conf": 0.0, "sliding_min_cont": np.nan,
                    "tau_Cmin": np.nan, "n_viol_short": 0})
        return out

    nACG = computeACG(st, BIN_SIZE, N_BINS)
    r = slidingRP_from_acg(nACG, n, rec_dur, cont_thresh=cont_thresh,
                           conf_thresh=90.0)
    out["sliding_max_conf"] = r["max_conf"]
    out["sliding_min_cont"] = r["min_cont"]
    out["tau_Cmin"] = r["rp_min_val"]
    out["n_viol_short"] = r["n_viol_short"]
    out["tau_pass0"] = r["tau_pass0"]
    for g in gammas:
        out[f"sliding_pass_{g}"] = r["max_conf"] >= g

    for rp_dur in hl_rps:
        p, est, obs = hill_llobet_from_acg(nACG, n, rec_dur, rp_dur, cont_thresh)
        tag = f"{rp_dur*1000:g}".replace(".", "p")
        out[f"hl{tag}_pass"] = p
        out[f"hl{tag}_est"] = est

    if true_rp is not None:
        p, est, _ = hill_llobet_from_acg(nACG, n, rec_dur, true_rp, cont_thresh)
        out["hl_oracle_pass"] = p
        out["hl_oracle_est"] = est
        c = poisson_test_fixed_tau(nACG, n, rec_dur, true_rp, cont_thresh)
        out["pt_oracle_conf"] = c
        for g in gammas:
            out[f"pt_oracle_pass_{g}"] = c >= g

    if with_estimator:
        fit = estimate_rp(nACG, RP_CENTERS, BIN_SIZE, 0.10)
        out["rp_est_ms"] = fit.rp_ms
        if np.isfinite(fit.rp_ms) and fit.rp_ms > 0:
            tau_e = fit.rp_ms / 1000
            p, est, _ = hill_llobet_from_acg(nACG, n, rec_dur, tau_e, cont_thresh)
            out["hl_est_pass"] = p
            out["hl_est_est"] = est
            c = poisson_test_fixed_tau(nACG, n, rec_dur, tau_e, cont_thresh)
            out["pt_est_conf"] = c
            for g in gammas:
                out[f"pt_est_pass_{g}"] = c >= g
        else:
            # estimation failed: a user with no RP estimate has no principled
            # fallback, so the arm counts as a rejection (the alternative,
            # falling back to 2 ms, is reported separately)
            out["hl_est_pass"] = False
            out["hl_est_est"] = np.nan
            out["pt_est_conf"] = np.nan
            for g in gammas:
                out[f"pt_est_pass_{g}"] = False

    if with_correction:
        cc = corrected_confidence_from_acg(nACG, n, rec_dur, cont_thresh)
        out["corrected_conf"] = cc
        for g in gammas:
            out[f"corrected_pass_{g}"] = cc >= g

    return out
