"""Excess short-lag coincidence between real neighbouring IBL neurons (05).

The bias that correlated firing introduces into Sliding RP is set by kappa, the
excess coincidence between a neuron and its contaminant:

    kappa(window) = CCG count in the window / count expected if independent - 1,

where the expectation is N_a N_b * (window length, both signs of lag) / D. A
contaminated unit's violations from neuron-contaminant pairs are scaled by
(1 + kappa) over the lags the metric tests, so a unit at true contamination C
behaves roughly like one at (1 + kappa) C (see correlated_rates.md). The
simulations report kappa for every correlated model, so this puts real pairs and
simulated pairs on the same axis.

Unlike a binned spike-count correlation, kappa is not diluted by counting noise
and does not depend on firing rate, which is why the count correlations in the
semi-synthetic figure (panel c) cannot be compared with the simulated rho.

Pairs: sorter-"good" units (IBL label 1) firing at least 2 spikes/s on the 8
cached brain-wide-map insertions, separated by at most 400 um in depth. The same
shared-spike screen as correlation_pairs.py is applied (CCG mass within
+/-0.5 ms between 0.5x and 2x its 5-10 ms level).

Run:  python pair_kappa.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from correlation_pairs import CACHE, corr_at  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/05_model_mismatch")
FS = 30000.0
WINDOWS_MS = {"k_0p5_3": (0.5, 3.0), "k_3_30": (3.0, 30.0), "k_30_300": (30.0, 300.0)}
MAX_PAIRS_PER_PROBE = 800
MIN_FR = 2.0


def _count(a, b, lo, hi):
    """Pairs with lag b - a in [lo, hi) or in (-hi, -lo], lags in samples."""
    return int((np.searchsorted(b, a + hi) - np.searchsorted(b, a + lo)
                + np.searchsorted(b, a - lo, side="right")
                - np.searchsorted(b, a - hi, side="right")).sum())


def ccg_central_ratio(a, b):
    """correlation_pairs.ccg_central_ratio, vectorised: CCG density within
    +/-0.5 ms relative to its density at 5-10 ms."""
    c = int(round(0.0005 * FS))
    central = int((np.searchsorted(b, a + c, side="right") - np.searchsorted(b, a - c)).sum())
    base = _count(a, b, int(0.005 * FS), int(0.010 * FS))
    if base == 0:
        return np.nan
    return (central / (2 * c)) / (base / (2 * (int(0.010 * FS) - int(0.005 * FS))))


def kappa(a, b, lo_ms, hi_ms, dur_s):
    """Excess coincidence of trains a, b (sample indices) at |lag| in [lo, hi)."""
    lo, hi = int(round(lo_ms * FS / 1000)), int(round(hi_ms * FS / 1000))
    n = (np.searchsorted(b, a + hi) - np.searchsorted(b, a + lo)
         + np.searchsorted(b, a - lo, side="right") - np.searchsorted(b, a - hi, side="right"))
    obs = float(n.sum())
    expected = a.size * b.size * 2 * (hi - lo) / (dur_s * FS)
    return obs / expected - 1.0 if expected > 0 else np.nan


def main(seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for fp in sorted(CACHE.glob("*.npz")):
        z = np.load(fp, allow_pickle=True)
        samp, clu = z["samples"], z["clusters"]
        dur = float(z["rec_dur"])
        cid_all, depth_all, label_all = z["cluster_id"], z["depths"], z["label"]
        order = np.argsort(clu, kind="stable")
        clu, samp = clu[order], samp[order]
        cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
        trains = {int(c): np.sort(samp[s:s + n]) for c, s, n in zip(cids, starts, counts)}
        info = {int(c): (float(depth_all[i]), float(label_all[i])) for i, c in enumerate(cid_all)}
        good = [c for c in trains if c in info and info[c][1] >= 1
                and trains[c].size / dur >= MIN_FR and np.isfinite(info[c][0])]
        cand = [(good[i], good[j]) for i in range(len(good)) for j in range(i + 1, len(good))
                if abs(info[good[i]][0] - info[good[j]][0]) <= 400]
        if not cand:
            continue
        take = rng.choice(len(cand), min(MAX_PAIRS_PER_PROBE, len(cand)), replace=False)
        for k in take:
            ca, cb = cand[k]
            a, b = trains[ca], trains[cb]
            row = dict(pid=fp.stem, sep_um=abs(info[ca][0] - info[cb][0]),
                       fr_a=a.size / dur, fr_b=b.size / dur,
                       ccg_ratio=ccg_central_ratio(a, b),
                       r_0p1=corr_at(a, b, dur, 0.1), r_1=corr_at(a, b, dur, 1.0))
            for key, (lo, hi) in WINDOWS_MS.items():
                row[key] = kappa(a, b, lo, hi, dur)
            rows.append(row)
        print(f"{fp.stem}: {len(good)} good units, {len(take)} pairs", flush=True)
    df = pd.DataFrame(rows)
    df.to_parquet(Path(r"D:/temp/slidingRP_resub/sims/pair_kappa.pqt"))

    clean = df.ccg_ratio.between(0.5, 2.0)
    sets = {"all pairs <= 400 um": np.ones(len(df), bool),
            "annulus 50-150 um, CCG-clean": df.sep_um.between(50, 150) & clean,
            "150-400 um, CCG-clean": df.sep_um.between(150, 400) & clean,
            "within 50 um, CCG-clean": (df.sep_um < 50) & clean}
    L = ["Excess short-lag coincidence (kappa) between real IBL neighbours", "=" * 64, "",
         f"{len(df):,} pairs of sorter-good units (>= {MIN_FR:g} spikes/s) on "
         f"{df.pid.nunique()} insertions.", "",
         "kappa = CCG count / count expected if independent - 1, over |lag| windows.",
         "0 means independent; +0.1 means 10% more coincidences than independence.", ""]
    L.append(f"{'selection':30s} {'n':>5}   " + "   ".join(
        f"{k:>22s}" for k in ("0.5-3 ms median [IQR]", "3-30 ms", "30-300 ms")))
    for name, m in sets.items():
        g = df[m]
        cells = []
        for key in WINDOWS_MS:
            v = g[key].dropna()
            cells.append(f"{v.median():+.3f} [{v.quantile(.25):+.3f},{v.quantile(.75):+.3f}]")
        L.append(f"{name:30s} {len(g):5d}   " + "   ".join(f"{c:>22s}" for c in cells))
    g = df[sets["annulus 50-150 um, CCG-clean"]]
    v = g.k_0p5_3.dropna()
    L += ["", "Annulus pairs, 0.5-3 ms: percentiles of kappa",
          "  " + "  ".join(f"p{q:g} {v.quantile(q / 100):+.3f}" for q in (5, 10, 25, 50, 75, 90, 95)),
          f"  fraction with kappa < 0 (fewer coincidences than independence, which hides "
          f"contamination): {np.mean(v < 0):.1%}",
          "", "How the count correlations relate (annulus pairs):",
          f"  median r at 100 ms {g.r_0p1.median():+.3f}, at 1 s {g.r_1.median():+.3f}",
          f"  Spearman(kappa 0.5-3 ms, r 100 ms) = "
          f"{g.k_0p5_3.corr(g.r_0p1, method='spearman'):+.2f}; "
          f"Spearman(kappa 30-300 ms, r 100 ms) = {g.k_30_300.corr(g.r_0p1, method='spearman'):+.2f}"]
    (OUTDIR / "pair_kappa.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
