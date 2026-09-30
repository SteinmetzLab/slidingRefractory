"""Re-score the censored datasets (Allen, macaque) with their censor window (08).

Every metric column the enrichment produced is recomputed from the stored
autocorrelograms with the dataset's censor, and the original value is kept as
``<column>_nocensor``. Downstream analyses read the corrected columns, so every
figure and table for these two datasets reflects the censor from then on.
IBL and Steinmetz have no censor and are left alone.

The censor used:
  macaque  7 / 30000 s. Kilosort 4 (which Nick believes produced it) removes
           same-cluster spikes closer than int(duplicate_spike_ms * fs / 1000)
           samples; its default is 0.25 ms, i.e. 7 samples at 30 kHz, and the
           data show exactly that floor.
  allen    the measured equivalent window, from the pooled autocorrelogram
           (censor_realdata.survival_profile): a 0.1 ms hard floor plus a shadow
           that is complete by 0.37 ms. For every tested window (> 0.5 ms) a
           single window of that length gives the same observable length as
           the full measured profile, so the package's scalar parameter is exact
           here.

Before anything is written, the recomputation with censor 0 is checked against
the stored columns on every unit (it must reproduce them).

Idempotent: always recomputes from the ACGs and the *_nocensor originals.

Run:  python apply_censor.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import BIN_SIZE, N_BINS, REF_DUR, RP_CENTERS, RP_REJECT  # noqa: E402
from acg_table import hill_llobet_from_acg, slidingRP_from_acg  # noqa: E402
from censor_realdata import observable_length, survival_profile  # noqa: E402
from censor_windows import ACG_DIR, ENRICHED, good_mask  # noqa: E402

OUT = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/08_censoring_window")
TAU_MINS = (0.00025, 0.0005, 0.00075, 0.001, 0.0015)
METRIC_COLS = ["max_conf", "min_cont", "tau_Cmin", "tau_first_pass", "tau_pass0", "passes",
               "hl2_pass", "hl2_est", "hl3_pass", "hl3_est", "hl_est_pass", "hl_est_est",
               "tau_last_pass", "n_accepted_bins"] + [
    f"pass_taumin_{t*1000:g}".replace(".", "p") for t in TAU_MINS]


def measured_equivalent(ds):
    """Equivalent hard window (s) from the dataset's pooled ACG profile."""
    pooled = np.zeros(N_BINS)
    for f in sorted((ENRICHED / ds).glob("*.pqt")):
        npy = ACG_DIR / ds / (f.stem + ".npy")
        e = pd.read_parquet(f)
        a = np.load(npy).astype(float)
        pooled += a[good_mask(ds, e)].sum(0)
    s = survival_profile(pooled, integer_lag=ds in ("ibl", "macaque"))
    L = observable_length(s)
    return float(REF_DUR[15] - L[15]), s


def score(acg, n, d, rp_ms, censor):
    """Every metric column for one unit, with a censor."""
    out = dict.fromkeys(METRIC_COLS, np.nan)
    if n < 2 or d <= 0:          # as the enrichment stored them
        out.update(passes=False, hl2_pass=False, hl3_pass=False, hl_est_pass=False,
                   n_accepted_bins=0)
        for t in TAU_MINS:
            out[f"pass_taumin_{t*1000:g}".replace(".", "p")] = False
        return out
    r = slidingRP_from_acg(acg, n, d, censor=censor)
    out.update(max_conf=r["max_conf"], min_cont=r["min_cont"], tau_Cmin=r["rp_min_val"],
               tau_first_pass=r["tau_first_pass"], tau_pass0=r["tau_pass0"] + censor,
               passes=r["passes"])
    for rp_dur, tag in ((0.002, "hl2"), (0.003, "hl3")):
        p, est, _ = hill_llobet_from_acg(acg, n, d, rp_dur, censor=censor)
        out[f"{tag}_pass"], out[f"{tag}_est"] = p, est
    if np.isfinite(rp_ms) and rp_ms > 0:
        p, est, _ = hill_llobet_from_acg(acg, n, d, rp_ms / 1000, censor=censor)
        out["hl_est_pass"], out["hl_est_est"] = p, est
    else:
        out["hl_est_pass"] = False
    for t in TAU_MINS:
        out[f"pass_taumin_{t*1000:g}".replace(".", "p")] = \
            slidingRP_from_acg(acg, n, d, rp_reject=t, censor=censor)["passes"]
    # last accepted window, as in add_tau_last.py but with the observable length
    c = 0.10
    nc, nb = n * c, n * (1 - c)
    lam = 2 * np.clip(REF_DUR - censor, 0, None) / d * nc * (nb + (nc - 1) / 2)
    conf = 1 - stats.poisson.cdf(np.cumsum(acg), lam)
    ok = (RP_CENTERS > RP_REJECT) & (conf >= 0.9)
    out["n_accepted_bins"] = int(ok.sum())
    out["tau_last_pass"] = float(RP_CENTERS[np.flatnonzero(ok)[-1]]) if ok.any() else np.nan
    return out


def score_shard(f, npy, censor):
    e = pd.read_parquet(f)
    acg = np.load(npy).astype(np.float64)
    orig = {c: (e[f"{c}_nocensor"] if f"{c}_nocensor" in e else e[c]).to_numpy()
            for c in METRIC_COLS}
    rows0, rows = [], []
    for i in range(len(e)):
        n, d = int(e.n_spikes.iat[i]), float(e.rec_dur_s.iat[i])
        rp_ms = float(e.rp_ms_10.iat[i]) if "rp_ms_10" in e else np.nan
        rows0.append(score(acg[i], n, d, rp_ms, 0.0))
        rows.append(score(acg[i], n, d, rp_ms, censor))
    z, r = pd.DataFrame(rows0), pd.DataFrame(rows)
    bad = []
    for c in METRIC_COLS:
        a, b = z[c].to_numpy(), orig[c]
        if a.dtype == bool or b.dtype == bool:
            same = np.array_equal(a.astype(bool), b.astype(bool))
        else:
            a, b = a.astype(float), b.astype(float)
            same = np.allclose(a, b, rtol=1e-9, atol=1e-12, equal_nan=True)
        if not same:
            bad.append(c)
    for c in METRIC_COLS:
        e[f"{c}_nocensor"] = orig[c]
        e[c] = r[c].to_numpy()
    e["censor_s"] = censor
    return f, e, bad


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    w_allen, s_allen = measured_equivalent("allen")
    censors = {"macaque": 7 / 30000, "allen": w_allen}
    (OUT / "censor_by_dataset.json").write_text(json.dumps(
        {"censor_s": censors, "ibl": 0.0, "steinmetz": 0.0,
         "allen_survival_profile_0_to_0p5ms": [round(float(v), 4) for v in s_allen[:15]],
         "notes": {"macaque": "Kilosort 4 default duplicate_spike_ms = 0.25 ms -> 7 samples",
                   "allen": "measured equivalent of the pooled-ACG shadow (unverified pipeline)"}},
        indent=2))
    print("censor (ms):", {k: round(v * 1000, 4) for k, v in censors.items()}, flush=True)
    for ds, w in censors.items():
        files = sorted((ENRICHED / ds).glob("*.pqt"))
        jobs = [(f, ACG_DIR / ds / (f.stem + ".npy"), w) for f in files]
        res = Parallel(n_jobs=4, verbose=0)(delayed(score_shard)(*j) for j in jobs)
        allbad = sorted({c for _, _, b in res for c in b})
        if allbad:
            raise RuntimeError(f"{ds}: censor-0 recomputation does not reproduce {allbad}")
        n0 = n1 = n = 0
        for f, e, _ in res:
            e.to_parquet(f)
            n += len(e)
            n0 += int(e.passes_nocensor.astype(bool).sum())
            n1 += int(e.passes.astype(bool).sum())
        print(f"{ds}: {len(files)} shards, {n:,} units re-scored with censor {w*1000:.4f} ms; "
              f"accepted {n0:,} -> {n1:,} (all units, any label)", flush=True)


if __name__ == "__main__":
    main()
