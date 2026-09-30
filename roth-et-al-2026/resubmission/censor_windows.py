"""What censor window did each dataset's spike sorter impose? (08)

Sorters commonly remove spikes that fall within a short exclusion window of an
earlier spike of the same unit (duplicate removal). If that window is w, no unit
can show an inter-spike interval below w, and the Sliding RP expectation, which
assumes violations can land anywhere in [0, tau_r], is then too large by the
factor tau_r / (tau_r - w). See 08_censoring_window/censoring_math.md.

The stored autocorrelograms have one-sample (1/30000 s) resolution, so the
window is visible directly:

  * pooled over every unit of a dataset, the ACG should be (near) zero below w
    and jump at w;
  * unit by unit, the first nonzero bin is the unit's minimum ISI; w is the
    floor of that distribution, and units with counts below it are candidates
    for having been merged after the censoring was applied (merging two
    censored units does not censor pairs across them).

IBL and macaque ACGs were computed from integer sample indices, so bin k is a
lag of exactly k samples. Allen and Steinmetz ACGs were computed from spike
times in seconds, so bin k is a lag in [k, k+1) samples; Allen's times are
also aligned to a common clock, so lags are not exact sample multiples.

Run:  python censor_windows.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402

ACG_DIR = Path(r"D:/temp/slidingRP_resub/acg_tables")
ENRICHED = Path(r"D:/temp/slidingRP_resub/enriched")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/08_censoring_window")
FS = 30000.0
NB_SHOW = 45                       # 0-1.5 ms
DATASETS = ("ibl", "allen", "steinmetz", "macaque")
NAMES = plotstyle.DATASET_NAMES


def good_mask(ds, e):
    """Sorter-accepted units, as in load_enriched.rule_common."""
    if ds == "ibl":
        return e["label"].astype(float).values >= 1
    if ds == "allen":
        return e.get("quality", pd.Series([""] * len(e))).astype(str).str.contains("good").values
    if ds == "steinmetz":
        return e["phy_annotation"].astype(float).values >= 2
    return np.ones(len(e), bool)


def load_dataset(ds):
    pooled_all = np.zeros(300)
    pooled_good = np.zeros(300)
    rows = []
    for f in sorted((ENRICHED / ds).glob("*.pqt")):
        npy = ACG_DIR / ds / (f.stem + ".npy")
        if not npy.exists():
            continue
        e = pd.read_parquet(f)
        acg = np.load(npy).astype(float)
        if acg.shape[0] != len(e):
            continue
        g = good_mask(ds, e)
        pooled_all += acg.sum(0)
        pooled_good += acg[g].sum(0)
        nz = acg > 0
        first = np.where(nz.any(1), nz.argmax(1), -1)
        for i in range(len(e)):
            rows.append(dict(insertion=f.stem, good=bool(g[i]), first_bin=int(first[i]),
                             n_spikes=int(e.n_spikes.iat[i]),
                             below5=float(acg[i, :5].sum()), below8=float(acg[i, :8].sum()),
                             total_0_1ms=float(acg[i, :30].sum())))
    return pooled_all, pooled_good, pd.DataFrame(rows)


def floor_estimate(pooled, frac=0.02):
    """First bin at which the pooled ACG reaches `frac` of its 1-1.5 ms level."""
    ref = pooled[30:45].mean()
    if ref <= 0:
        return np.nan
    above = np.flatnonzero(pooled[:45] >= frac * ref)
    return int(above[0]) if above.size else np.nan


def main():
    plotstyle.apply()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)
    L = ["Censor windows by dataset", "=" * 26, "",
         "Bins are one sample (1/30000 s = 33.3 us). For IBL and macaque bin k is a lag",
         "of exactly k samples (bin 0, zero lag, is excluded by construction). For Allen",
         "and Steinmetz bin k is a lag in [k, k+1) samples.", ""]
    fig, axs = plt.subplots(1, len(DATASETS), figsize=(14, 3.4))
    summary = []
    for ax, ds in zip(axs, DATASETS):
        pa, pg, u = load_dataset(ds)
        if u.empty:
            continue
        x = np.arange(NB_SHOW) / FS * 1000
        ref_a, ref_g = pa[30:45].mean(), pg[30:45].mean()
        ax.step(x, pa[:NB_SHOW] / ref_a, where="post", color="0.55", lw=1.3, label="All units")
        ax.step(x, pg[:NB_SHOW] / ref_g, where="post", color="k", lw=1.6, label="Sorter-accepted")
        ax.set_yscale("symlog", linthresh=1e-3)
        ax.set_ylim(0, 2)
        ax.set_xlim(0, NB_SHOW / FS * 1000)
        ax.set_xlabel("Lag (ms)")
        ax.set_title(NAMES.get(ds, ds), loc="left")
        if ds == DATASETS[0]:
            ax.set_ylabel("Pooled ACG / its 1-1.5 ms level")
        leg = ax.legend(handlelength=0, handletextpad=0, fontsize=8)
        for t, h in zip(leg.get_texts(), leg.legend_handles):
            t.set_color(h.get_color()); h.set_visible(False)

        w_all, w_good = floor_estimate(pa), floor_estimate(pg)
        ug = u[u.good]
        fb = ug.first_bin[ug.first_bin >= 0]
        per_ins = ug[ug.first_bin >= 0].groupby("insertion").first_bin.min()
        L += [f"--- {NAMES.get(ds, ds)} ---",
              f"units: {len(u):,} ({u.good.sum():,} sorter-accepted), insertions: {u.insertion.nunique()}",
              "pooled ACG, first 12 bins (counts), sorter-accepted units:",
              "  " + " ".join(f"{int(v)}" for v in pg[:12]),
              "  relative to the 1-1.5 ms level: " + " ".join(f"{v / ref_g:.3f}" for v in pg[:12]),
              f"pooled floor (first bin >= 2% of the 1-1.5 ms level): all units {w_all}, "
              f"accepted {w_good} samples",
              "unit minimum ISI (first nonzero bin), sorter-accepted units:",
              "  " + ", ".join(f"bin {int(b)}: {int((fb == b).sum()):,}" for b in sorted(fb.unique())[:10]),
              f"  minimum across units {int(fb.min()) if len(fb) else 'n/a'}; median {fb.median():.0f}",
              f"per-insertion minimum ISI (bins): " + ", ".join(
                  f"{int(b)}: {int((per_ins == b).sum())}" for b in sorted(per_ins.unique())[:10]),
              ""]
        summary.append(dict(dataset=ds, pooled_floor_all=w_all, pooled_floor_good=w_good,
                            unit_min_bin=int(fb.min()) if len(fb) else np.nan,
                            frac_units_below_floor=float(np.mean(fb < (w_good or 0)))
                            if len(fb) and np.isfinite(w_good) else np.nan))
    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "censor_windows")
    pd.DataFrame(summary).to_csv(OUTDIR / "censor_windows.csv", index=False)
    (OUTDIR / "censor_windows.txt").write_text("\n".join(L))
    print("\n".join(L))
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
