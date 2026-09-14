"""Run the resubmission simulation sweeps.

    python run_sims.py calibration      # 03: realised false-acceptance rates
    python run_sims.py searchspace      # 03: dependence on tau_max / bin width
    python run_sims.py corrected        # 03: FWER-corrected test across true RP
    python run_sims.py decomposition    # 04: fixed / oracle / estimated / sliding
    python run_sims.py mismatch         # 05: departures from the generative model
    python run_sims.py all

Each writes a tidy parquet to ``D:/temp/slidingRP_resub/sims/<name>.pqt`` with
one row per condition and pass fractions at every confidence threshold, plus a
Clopper-Pearson 95% interval on the headline rate.
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from resub_simulations import evaluate_arms, make_train  # noqa: E402

OUT = Path(r"D:/temp/slidingRP_resub/sims")
sweep = None  # bound in __main__ to _sweep_impl with the chosen n_jobs
GAMMAS = (50, 60, 70, 75, 80, 85, 90, 95, 99)
CONT_GRID = np.array([0, 2, 4, 6, 7, 8, 8.5, 9, 9.5, 10,
                      10.5, 11, 11.5, 12, 13, 14, 16, 18, 20]) / 100


def _one_condition(cond, n_sim, seed, arms_kw):
    """Run n_sim spike trains for one parameter combination; aggregate."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_sim):
        st, info = make_train(cond["model"], cond["total_rate"], cond["cont_prop"],
                              cond["rec_dur"], cond["rp_dur"], rng=rng,
                              **cond.get("model_kw", {}))
        out = evaluate_arms(st, cond["rec_dur"],
                            true_rp=cond["rp_dur"] if arms_kw.get("oracle") else None,
                            with_estimator=arms_kw.get("estimator", False),
                            with_correction=arms_kw.get("correction", False),
                            gammas=GAMMAS)
        out.update({k: v for k, v in info.items() if k.startswith("realised")})
        rows.append(out)
    df = pd.DataFrame(rows)
    agg = {k: v for k, v in cond.items() if k != "model_kw"}
    agg.update({f"{k}": str(v) for k, v in cond.get("model_kw", {}).items()})
    agg["n_sim"] = n_sim
    for c in df.columns:
        if df[c].dtype == bool:
            agg[c] = float(df[c].mean())
        elif np.issubdtype(df[c].dtype, np.number):
            agg[c + "_mean"] = float(np.nanmean(df[c]))
            agg[c + "_median"] = float(np.nanmedian(df[c]))
    # Clopper-Pearson interval on the headline Sliding RP pass rate at 90%
    k = int(df["sliding_pass_90"].sum()) if "sliding_pass_90" in df else 0
    lo, hi = stats.beta.ppf(0.025, k, n_sim - k + 1) if k else 0.0, \
        stats.beta.ppf(0.975, k + 1, n_sim - k) if k < n_sim else 1.0
    agg["sliding_pass_90_lo"] = float(np.nan_to_num(lo))
    agg["sliding_pass_90_hi"] = float(np.nan_to_num(hi, nan=1.0))
    return agg


def _sweep_impl(name, conditions, n_sim, arms_kw, n_jobs=11, seed0=20260913):
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"[{name}] {len(conditions)} conditions x {n_sim} sims "
          f"= {len(conditions)*n_sim:,} spike trains", flush=True)
    res = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_one_condition)(c, n_sim, seed0 + i, arms_kw)
        for i, c in enumerate(conditions))
    df = pd.DataFrame(res)
    df.to_parquet(OUT / f"{name}.pqt")
    print(f"[{name}] done in {(time.time()-t0)/60:.1f} min -> {OUT/f'{name}.pqt'}",
          flush=True)
    return df


def grid(**kw):
    """Cartesian product of keyword lists -> list of condition dicts."""
    keys = list(kw)
    return [dict(zip(keys, vals)) for vals in itertools.product(*(kw[k] for k in keys))]


# --------------------------------------------------------------------------

def run_calibration(n_sim=4000):
    """03 deliverable 1: realised false-acceptance rate vs nominal.

    Contamination is simulated exactly at the threshold (the definition point
    of the nominal rate) and at 0.8x the threshold (true acceptance / power).
    The confidence sweep is free: max_conf is computed once per train.
    """
    conds = grid(model=["standard"], total_rate=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0],
                 rp_dur=[0.001, 0.0015, 0.002, 0.003, 0.005],
                 rec_dur=[1800.0, 3600.0, 7200.0, 14400.0],
                 cont_prop=[0.10, 0.08])
    return sweep("calibration", conds, n_sim, dict(oracle=False))


def run_searchspace(n_sim=4000):
    """03 deliverable 2: does the realised rate depend on tau_max / resolution?

    The tested window and bin width are properties of the *analysis*, not the
    data, so the same simulated trains are re-analysed at each setting. That is
    done in the analysis notebook from the stored ACGs; here we vary them
    through the metric by simulating at several true RPs and recording the full
    confidence trace, so the sweep is a re-analysis rather than a re-simulation.
    """
    conds = grid(model=["standard"], total_rate=[5.0],
                 rp_dur=[0.0015, 0.003, 0.010],
                 rec_dur=[7200.0], cont_prop=[0.10, 0.08])
    return sweep("searchspace_trains", conds, n_sim, dict(oracle=False))


def run_corrected(n_sim=1500):
    """03 deliverable 3: is the full-window null least favourable?

    If the coupling argument holds, the FWER-corrected test should be exactly
    calibrated when the true RP spans the whole tested window (10 ms) and
    strictly conservative for shorter true RPs. Firing rates are capped at
    10 spikes/s: the first-passage DP's state space grows with the observed
    violation count, and at 20 spikes/s a single evaluation takes ~8 s.
    """
    conds = grid(model=["standard"], total_rate=[1.0, 2.0, 5.0, 10.0],
                 rp_dur=[0.001, 0.0015, 0.002, 0.003, 0.005, 0.010],
                 rec_dur=[3600.0, 7200.0], cont_prop=[0.10, 0.08])
    return sweep("corrected", conds, n_sim, dict(oracle=False, correction=True))


def run_decomposition(n_sim=600):
    """04: separate 'avoiding a misspecified RP' from 'sliding' and from
    'treating the count statistically'."""
    conds = grid(model=["standard"], total_rate=[0.5, 1.0, 2.0, 5.0, 10.0],
                 rp_dur=[0.0015, 0.002, 0.003, 0.005],
                 rec_dur=[7200.0], cont_prop=list(CONT_GRID))
    return sweep("decomposition", conds, n_sim,
                 dict(oracle=True, estimator=True))


def run_decomposition_realistic(n_sim=400):
    """04, second grid: the same arms on ACG shapes where RP estimation is hard."""
    conds = (grid(model=["graded"], total_rate=[1.0, 5.0], rp_dur=[0.002, 0.003],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                  model_kw=[{"width": 0.0005}, {"width": 0.001}, {"width": 0.002}])
             + grid(model=["bursting"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                    rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                    model_kw=[{"p_burst": 0.1}, {"p_burst": 0.3}]))
    return sweep("decomposition_realistic", conds, n_sim,
                 dict(oracle=True, estimator=True))


def run_mismatch(n_sim=1000):
    """05: departures from the generative model."""
    conds = []
    # A: contamination from a single neuron with its own refractory period
    conds += grid(model=["single_neuron_contaminant"], total_rate=[1.0, 5.0],
                  rp_dur=[0.0015, 0.003], rec_dur=[7200.0],
                  cont_prop=list(CONT_GRID),
                  model_kw=[{"cont_rp": 0.0015}, {"cont_rp": 0.0025}])
    # B1: correlated / non-stationary rates
    conds += grid(model=["modulated"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                  model_kw=[{"rho": r, "tau_s": t}
                            for r in (-1.0, -0.5, 0.0, 0.5, 1.0) for t in (2.0, 30.0)])
    # B2: non-overlapping activity (the place-cell failure mode)
    conds += grid(model=["nonoverlapping"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                  model_kw=[{"overlap": o} for o in (0.0, 0.25, 0.5, 0.75, 1.0)])
    # C: graded recovery and bursting
    conds += grid(model=["graded"], total_rate=[1.0, 5.0], rp_dur=[0.002, 0.003],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                  model_kw=[{"width": w} for w in (0.0005, 0.001, 0.002)])
    conds += grid(model=["bursting"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID),
                  model_kw=[{"p_burst": p} for p in (0.1, 0.3)])
    # reference: the manuscript's own model, same grid
    conds += grid(model=["standard"], total_rate=[1.0, 5.0], rp_dur=[0.002, 0.003],
                  rec_dur=[7200.0], cont_prop=list(CONT_GRID))
    return sweep("mismatch", conds, n_sim, dict(oracle=False))


JOBS = {
    "calibration": run_calibration,
    "searchspace": run_searchspace,
    "corrected": run_corrected,
    "decomposition": run_decomposition,
    "decomposition_realistic": run_decomposition_realistic,
    "mismatch": run_mismatch,
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("job", choices=list(JOBS) + ["all"])
    ap.add_argument("--n-sim", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=11)
    args = ap.parse_args()
    todo = list(JOBS) if args.job == "all" else [args.job]
    import functools
    for j in todo:
        fn = JOBS[j]
        globals()["sweep"] = functools.partial(_sweep_impl, n_jobs=args.n_jobs)
        fn(n_sim=args.n_sim) if args.n_sim else fn()
