"""Second pass over the per-unit ACG tables: metrics and RP estimates.

Takes the ``<id>.npy`` / ``<id>.pqt`` pairs written by ``build_ibl.py``,
``build_allen.py`` and ``build_local.py`` and adds, for every unit:

* Sliding RP at the paper defaults: max confidence, minimum confirmable
  contamination, the tau_r at which it is reached, the shortest tau_r that
  reaches threshold, the <2 ms violation count, and tau_pass0;
* Hill-Llobet pass/fail and point estimate at 2 ms and 3 ms, and at the unit's
  own estimated RP;
* Sliding RP pass at several tau_min values (the tau_min sweep of package 07);
* the canonical ACG-recovery estimate at 5%, 10% and 20% recovery, with the
  floor flag and fit quality (package 01).

The RP estimator is the expensive part (~0.15 s/unit), so it is applied only to
units that clear a light pre-filter; everything else gets NaN. Metrics are
computed for every unit.

Usage
-----
    python enrich_table.py ibl --n-jobs 11
    python enrich_table.py all
"""
from __future__ import annotations

import argparse
import sys
import zlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import (BIN_SIZE, N_BINS, RP_CENTERS, hill_llobet_from_acg,  # noqa: E402
                       slidingRP_from_acg, tau_pass0)

from slidingRP.rp_estimate import estimate_rp  # noqa: E402

TABLES = Path(r"D:/temp/slidingRP_resub/acg_tables")
OUT = Path(r"D:/temp/slidingRP_resub/enriched")
TAU_MINS = (0.00025, 0.0005, 0.00075, 0.001, 0.0015)

# Pre-filter for the (expensive) RP estimator. Deliberately loose: the
# inclusion rules for Fig 1 are applied later, on the estimates.
MIN_SPIKES_FOR_RP = 500
MIN_FR_FOR_RP = 0.5


def enrich_one(npy_path, estimate=True, sens_frac=0.12, seed=20260913):
    """Enrich one file. The 5%/20% recovery-fraction sensitivity is computed on
    a deterministic random subsample (`sens_frac`) rather than every unit: the
    estimator dominates the cost and the sensitivity analysis only needs a
    sample. rp_ms_10 is computed for every eligible unit."""
    pqt = npy_path.with_suffix(".pqt")
    acgs = np.load(npy_path)
    tbl = pd.read_parquet(pqt)
    if len(tbl) != acgs.shape[0]:
        raise ValueError(f"{npy_path.name}: {len(tbl)} rows vs {acgs.shape[0]} ACGs")

    n_spikes = tbl["n_spikes"].to_numpy()
    rec_dur = tbl["rec_dur_s"].to_numpy()
    out = {k: np.full(len(tbl), np.nan) for k in (
        "max_conf", "min_cont", "tau_Cmin", "tau_first_pass", "n_viol_short",
        "tau_pass0", "hl2_est", "hl3_est", "hl_est_est",
        "rp_ms_10", "rp_ms_05", "rp_ms_20", "rp_r2")}
    for k in ("passes", "hl2_pass", "hl3_pass", "hl_est_pass", "rp_floor"):
        out[k] = np.zeros(len(tbl), dtype=bool)
    for t in TAU_MINS:
        out[f"pass_taumin_{t*1000:g}".replace(".", "p")] = np.zeros(len(tbl), bool)
    out["n_viol_1ms"] = np.zeros(len(tbl), dtype=np.int64)
    out["n_viol_3ms"] = np.zeros(len(tbl), dtype=np.int64)
    out["first_nonzero_bin"] = np.full(len(tbl), -1, dtype=np.int64)
    out["acg_0_0p5ms"] = np.zeros(len(tbl), dtype=np.int64)
    out["acg_0p5_1ms"] = np.zeros(len(tbl), dtype=np.int64)

    # zlib.crc32, not hash(): Python's string hash is randomised per process
    rng = np.random.default_rng(zlib.crc32(npy_path.stem.encode()) + seed)
    sens_pick = rng.random(len(tbl)) < sens_frac

    i1 = int(np.argmax(RP_CENTERS > 0.001)) + 1
    i3 = int(np.argmax(RP_CENTERS > 0.003)) + 1
    ih = int(np.argmax(RP_CENTERS > 0.0005)) + 1
    bin_centers = RP_CENTERS

    for i in range(len(tbl)):
        a = acgs[i].astype(np.float64)
        n, d = int(n_spikes[i]), float(rec_dur[i])
        if n < 2 or d <= 0:
            continue
        r = slidingRP_from_acg(a, n, d)
        out["max_conf"][i] = r["max_conf"]
        out["min_cont"][i] = r["min_cont"]
        out["tau_Cmin"][i] = r["rp_min_val"]
        out["tau_first_pass"][i] = r["tau_first_pass"]
        out["n_viol_short"][i] = r["n_viol_short"]
        out["tau_pass0"][i] = r["tau_pass0"]
        out["passes"][i] = r["passes"]
        for rp_dur, tag in ((0.002, "hl2"), (0.003, "hl3")):
            p, est, _ = hill_llobet_from_acg(a, n, d, rp_dur)
            out[f"{tag}_pass"][i] = p
            out[f"{tag}_est"][i] = est
        for t in TAU_MINS:
            rr = slidingRP_from_acg(a, n, d, rp_reject=t)
            out[f"pass_taumin_{t*1000:g}".replace(".", "p")][i] = rr["passes"]
        out["n_viol_1ms"][i] = int(a[:i1].sum())
        out["n_viol_3ms"][i] = int(a[:i3].sum())
        nz = np.flatnonzero(a > 0)
        out["first_nonzero_bin"][i] = nz[0] if nz.size else -1
        out["acg_0_0p5ms"][i] = int(a[:ih].sum())
        out["acg_0p5_1ms"][i] = int(a[ih:i1].sum())

        if estimate and n >= MIN_SPIKES_FOR_RP and n / d >= MIN_FR_FOR_RP:
            f10 = estimate_rp(a, bin_centers, BIN_SIZE, 0.10)
            out["rp_ms_10"][i] = f10.rp_ms
            out["rp_floor"][i] = f10.floor_applied
            out["rp_r2"][i] = f10.rsquared
            if sens_pick[i]:
                out["rp_ms_05"][i] = estimate_rp(a, bin_centers, BIN_SIZE, 0.05).rp_ms
                out["rp_ms_20"][i] = estimate_rp(a, bin_centers, BIN_SIZE, 0.20).rp_ms
            if np.isfinite(f10.rp_ms) and f10.rp_ms > 0:
                p, est, _ = hill_llobet_from_acg(a, n, d, f10.rp_ms / 1000)
                out["hl_est_pass"][i] = p
                out["hl_est_est"][i] = est

    for k, v in out.items():
        tbl[k] = v
    return tbl


def _worker(p, estimate, sens_frac=0.12):
    try:
        t = enrich_one(p, estimate, sens_frac)
        return p.stem, t, None
    except Exception as e:  # noqa: BLE001
        return p.stem, None, f"{type(e).__name__}: {e}"


def run(dataset, n_jobs=11, estimate=True, sens_frac=0.12):
    src = TABLES / dataset
    dst = OUT / dataset
    dst.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in dst.glob("*.pqt")}
    files = [p for p in sorted(src.glob("*.npy")) if p.stem not in done]
    print(f"[{dataset}] {len(files)} files to enrich ({len(done)} done)", flush=True)
    if not files:
        return
    t0 = time.time()
    for res in Parallel(n_jobs=n_jobs, verbose=5, return_as="generator_unordered")(
            delayed(_worker)(p, estimate, sens_frac) for p in files):
        stem, tbl, err = res
        if err:
            print(f"  {stem} FAILED {err}", flush=True)
        else:
            tbl.to_parquet(dst / f"{stem}.pqt")
    print(f"[{dataset}] done in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", choices=["ibl", "allen", "steinmetz", "macaque", "all"])
    ap.add_argument("--n-jobs", type=int, default=11)
    ap.add_argument("--no-estimate", action="store_true")
    ap.add_argument("--sens-frac", type=float, default=0.12,
                    help="fraction of units given the 5%/20% recovery sensitivity")
    a = ap.parse_args()
    sets = ["macaque", "steinmetz", "allen", "ibl"] if a.dataset == "all" else [a.dataset]
    for s in sets:
        run(s, a.n_jobs, not a.no_estimate, a.sens_frac)
