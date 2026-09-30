"""Re-score every real unit with its dataset's censor window accounted for (08).

The published method assumes violations can land anywhere in [0, tau_r]. With a
sorter censor window, the lags that can be observed are only those where the
sorter keeps spikes, so the expected count should use the observable length

    L(tau_r) = integral_0^tau_r s(t) dt

where s(t) in [0, 1] is the fraction of spike pairs at lag t that survive the
sorter. A hard window w is s = 0 below w and 1 above, giving L = tau_r - w; a
soft shadow (Allen) is a ramp. s(t) is measured per dataset from the pooled
autocorrelogram of its sorter-accepted units, relative to its level between 0.4
and 0.8 ms, over lags below 0.5 ms (beyond that s = 1). See
08_censoring_window/censoring_math.md.

For each dataset this reports the acceptance rate of sorter-accepted units with
the published method, with a hard window at the measured floor, and with the
measured profile, and how many units change decision.

Run:  python censor_realdata.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import BIN_SIZE, N_BINS, REF_DUR, RP_CENTERS, RP_REJECT  # noqa: E402
from censor_windows import ACG_DIR, ENRICHED, good_mask  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/08_censoring_window")
DATASETS = ("ibl", "allen", "steinmetz", "macaque")
HARD_BINS = {"ibl": 0, "allen": 3, "steinmetz": 0, "macaque": 7}   # from censor_windows.py
PROFILE_UPTO = 15            # bins (0.5 ms); s = 1 beyond
PLATEAU = slice(12, 24)      # 0.4-0.8 ms
REGIONS = ("Isocortex", "HPF", "TH")
# IBL and macaque ACGs use exact integer lags, whose bin 0 (lag 0) is excluded
# by construction for every unit; that is not the sorter's doing
INTEGER_LAG = {"ibl", "macaque"}
# Where the profile is not a clean sorter signature: Steinmetz's pooled ACG
# rises gently across the whole 0-1.5 ms range, so a deficit relative to the
# 0.4-0.8 ms level there is at least partly biology. Report it as an upper bound.
PROFILE_IS_UPPER_BOUND = {"steinmetz"}


def survival_profile(pooled, integer_lag=False):
    """s_k: fraction of pairs surviving the sorter at lag bin k."""
    plateau = np.median(pooled[PLATEAU])
    s = np.ones(N_BINS)
    if plateau > 0:
        s[:PROFILE_UPTO] = np.clip(pooled[:PROFILE_UPTO] / plateau, 0, 1)
    if integer_lag and pooled[1:PROFILE_UPTO].min() > 0:
        s[0] = 1.0        # zero lag is excluded by construction, not censored
    return s


def observable_length(s):
    """L_k = sum_{j<=k} s_j * bin: the observable part of each window."""
    return np.cumsum(s) * BIN_SIZE


def accept(acg, n, dur, L, cont=10.0, conf=90.0, tau_min=RP_REJECT):
    """Vectorised Sliding RP acceptance with window lengths L (per bin)."""
    c = cont / 100
    nc, nb = n[:, None] * c, n[:, None] * (1 - c)
    expected = 2 * L[None, :] / dur[:, None] * nc * (nb + (nc - 1) / 2)
    obs = np.cumsum(acg, axis=1)
    conf_ = 1 - stats.poisson.cdf(obs, expected)
    test = RP_CENTERS > tau_min
    return (conf_[:, test] >= conf / 100).any(axis=1)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    L = ["Real units re-scored with their dataset's censor window", "=" * 56, ""]
    rows = []
    for ds in DATASETS:
        files = sorted((ENRICHED / ds).glob("*.pqt"))
        pooled = np.zeros(N_BINS)
        shards = []
        for f in files:
            npy = ACG_DIR / ds / (f.stem + ".npy")
            if not npy.exists():
                continue
            e = pd.read_parquet(f)
            acg = np.load(npy).astype(float)
            if acg.shape[0] != len(e):
                continue
            g = good_mask(ds, e)
            pooled += acg[g].sum(0)
            shards.append((e[g].reset_index(drop=True), acg[g]))
        s = survival_profile(pooled, integer_lag=ds in INTEGER_LAG)
        L_plain = REF_DUR.copy()
        L_hard = np.clip(REF_DUR - HARD_BINS[ds] * BIN_SIZE, 0, None)
        L_prof = observable_length(s)
        w_eff = (REF_DUR[PROFILE_UPTO] - L_prof[PROFILE_UPTO]) * 1000
        res = []
        for e, acg in shards:
            n = e.n_spikes.to_numpy(float)
            dur = e.rec_dur_s.to_numpy(float)
            a0 = accept(acg, n, dur, L_plain)
            # after apply_censor.py the stored decision is the censored one and
            # the published decision is kept as passes_nocensor
            published = (e.passes_nocensor if "passes_nocensor" in e else e.passes)
            assert np.array_equal(a0, published.to_numpy(bool)), "must reproduce published decisions"
            res.append(pd.DataFrame(dict(cosmos=e.cosmos.values, a0=a0,
                                         a_hard=accept(acg, n, dur, L_hard),
                                         a_prof=accept(acg, n, dur, L_prof))))
        r = pd.concat(res, ignore_index=True)
        L += [f"--- {ds} ---",
              f"sorter-accepted units {len(r):,}; hard floor {HARD_BINS[ds]} samples "
              f"({HARD_BINS[ds] / 30:.3f} ms); measured survival profile s(t) over 0-0.5 ms:",
              "  " + " ".join(f"{v:.2f}" for v in s[:PROFILE_UPTO]),
              f"  equivalent hard window (lost observable length at 0.5 ms): {w_eff:.3f} ms",
              f"acceptance: published {r.a0.mean():.3f}   hard window {r.a_hard.mean():.3f}   "
              f"measured profile {r.a_prof.mean():.3f}"
              + ("   (profile is an UPPER BOUND here: see PROFILE_IS_UPPER_BOUND)"
                 if ds in PROFILE_IS_UPPER_BOUND else ""),
              f"accepted by the published method but not with the profile: "
              f"{(r.a0 & ~r.a_prof).sum():,} units "
              f"({(r.a0 & ~r.a_prof).sum() / max(r.a0.sum(), 1):.1%} of accepted); "
              f"the reverse: {(~r.a0 & r.a_prof).sum():,}"]
        for reg in REGIONS:
            g = r[r.cosmos == reg]
            if len(g) < 50:
                continue
            L.append(f"  {reg:10s} n={len(g):7,}  published {g.a0.mean():.3f}  "
                     f"hard {g.a_hard.mean():.3f}  profile {g.a_prof.mean():.3f}")
            rows.append(dict(dataset=ds, region=reg, n=len(g), published=g.a0.mean(),
                             hard=g.a_hard.mean(), profile=g.a_prof.mean(), w_eff_ms=w_eff))
        L.append("")
    pd.DataFrame(rows).to_csv(OUTDIR / "censor_realdata.csv", index=False)
    (OUTDIR / "censor_realdata.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
