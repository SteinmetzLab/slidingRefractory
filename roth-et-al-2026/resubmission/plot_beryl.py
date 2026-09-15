"""Beryl-level supplement: recovery time and acceptance across brain regions.

Approved as a supplement to consider (root plan, Q10). The Cosmos level used in
Fig 1 pools each region very coarsely -- "thalamus" spans LGd, LP, VPM, PO and
the rest. With the full IBL brain-wide map there are enough units to look one
level finer, at the Beryl parcellation, which is where a reader would go to ask
whether a particular nucleus is short or long.

Two panels: estimated recovery time per Beryl region (the Fig 1 quantity), and
the Sliding RP acceptance rate per region against Hill-Llobet at 3 ms, which
shows where a fixed threshold costs the most.

Run:  python plot_beryl.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402
from load_enriched import load_all, rule_common  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"01_fig1_rp_durations")
MIN_UNITS = 150
MIN_INSERTIONS = 5


def bootstrap_ci(x, n=2000, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    bs = np.median(x[rng.integers(0, x.size, (n, x.size))], axis=1)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main():
    plotstyle.apply()
    df = load_all(verbose=False)
    if df.empty:
        print("no enriched tables")
        return
    inc = df[rule_common(df, require_pass=True)]

    rows = []
    for reg, g in inc.groupby("beryl"):
        if not reg or reg in ("root", "void") or len(g) < MIN_UNITS:
            continue
        if g.insertion_key.nunique() < MIN_INSERTIONS:
            continue
        lo, hi = bootstrap_ci(g.rp_ms_10)
        acc = df[(df.beryl == reg) & ((df.sorter_label.isna()) | (df.sorter_label > 0))]
        rows.append(dict(beryl=reg, cosmos=g.cosmos.mode().iat[0], n=len(g),
                         n_ins=g.insertion_key.nunique(),
                         n_animals=g.animal_key.nunique(),
                         median=float(np.nanmedian(g.rp_ms_10)), lo=lo, hi=hi,
                         frac_below2=float(np.nanmean(g.rp_ms_10 < 2)),
                         fr=float(np.nanmedian(g.firing_rate)),
                         n_acc=len(acc),
                         sliding=float(acc.passes.mean()) if len(acc) else np.nan,
                         hl3=float(acc.hl3_pass.mean()) if len(acc) else np.nan))
    t = pd.DataFrame(rows).sort_values("median")
    t["rescued"] = t.sliding - t.hl3
    t.to_csv(OUTDIR / "beryl_table.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(11.5, max(4.5, 0.21 * len(t))),
                            gridspec_kw={"width_ratios": [1.25, 1]})
    y = np.arange(len(t))[::-1]
    cols = [plotstyle.REGION_COLORS.get(c, "0.5") for c in t.cosmos]

    ax = axs[0]
    for yy, r, c in zip(y, t.itertuples(), cols):
        ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=2.4, solid_capstyle="butt")
        ax.plot(r.median, yy, "o", color="w", mec=c, mew=1.3, ms=5)
    ax.axvline(2, color="0.7", ls=":", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.beryl}  ({r.cosmos}, n={r.n:,})" for r in t.itertuples()],
                       fontsize=6)
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_title("a  Median (95% CI) by Beryl region", loc="left")

    ax = axs[1]
    ax.barh(y + 0.2, t.sliding, 0.4, color="#1b9e77", label="Sliding RP")
    ax.barh(y - 0.2, t.hl3, 0.4, color="#7570b3", label="Hill-Llobet, 3 ms")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("Proportion of units accepted")
    ax.set_title("b  Acceptance rate", loc="left")
    ax.legend(fontsize=7, loc="lower right")

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "beryl")

    L = ["Beryl-level recovery times and acceptance rates", "=" * 52, "",
         f"Regions with at least {MIN_UNITS} units passing the common inclusion rule",
         f"and at least {MIN_INSERTIONS} insertions. {len(t)} regions qualify.", ""]
    L.append(f"{'region':10s} {'cosmos':10s} {'n':>7} {'ins':>5} {'animals':>8} "
             f"{'median':>8} {'[95% CI]':>18} {'<2ms':>7} {'FR':>7} "
             f"{'slidRP':>8} {'HL3ms':>7} {'rescued':>8}")
    for r in t.itertuples():
        L.append(f"{r.beryl:10s} {r.cosmos:10s} {r.n:7,} {r.n_ins:5d} {r.n_animals:8d} "
                 f"{r.median:8.3f} {f'[{r.lo:.3f}, {r.hi:.3f}]':>18} "
                 f"{r.frac_below2:7.1%} {r.fr:7.1f} {r.sliding:8.3f} {r.hl3:7.3f} "
                 f"{r.rescued:8.3f}")
    sh = t.nsmallest(5, "median")
    lo_ = t.nlargest(5, "median")
    L += ["", f"Shortest: " + ", ".join(f"{r.beryl} ({r.median:.2f} ms)" for r in sh.itertuples()),
          f"Longest : " + ", ".join(f"{r.beryl} ({r.median:.2f} ms)" for r in lo_.itertuples()),
          "",
          "Most rescued by Sliding RP relative to a fixed 3 ms threshold:",
          "  " + ", ".join(f"{r.beryl} (+{r.rescued:.2f})"
                           for r in t.nlargest(6, 'rescued').itertuples()),
          "",
          "Read with care: these are apparent ACG recovery times under one",
          "operational definition, not biophysical refractory periods, and the",
          "regions differ in firing rate, sorter yield and sample size."]
    (OUTDIR / "beryl_numbers.txt").write_text("\n".join(L))
    print("\n".join(L[:40]))


if __name__ == "__main__":
    main()
