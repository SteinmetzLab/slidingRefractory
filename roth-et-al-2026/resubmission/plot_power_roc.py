"""ROC curves: Sliding RP against the Poisson test at the true RP and at fixed
windows, compared at matched false acceptance (04). Reads power_roc.pqt.

For each firing rate and true RP, each arm's confidence is thresholded at every
value: true acceptance is the fraction of 8%-contaminated trains accepted, false
acceptance the fraction of 12%-contaminated trains accepted. The area under the
curve is the probability that a random 8% train gets a higher confidence than a
random 12% train (ties count half).

The key comparison: at the false acceptance Sliding RP actually has at gamma =
90, how much true acceptance does the oracle test reach?

Run:  python plot_power_roc.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402

SIMS = Path(r"D:/temp/slidingRP_resub/sims")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/04_hill_llobet_decomposition")
# drawing order: fixed windows underneath, the true-RP test and Sliding RP on top
# (where the true RP is 3 ms the fixed 3 ms curve is identical to the true-RP one)
ARMS = [("fixed_0.5ms", "Poisson test, fixed 0.5 ms", "0.72", "-", 1.3),
        ("fixed_1ms", "Poisson test, fixed 1 ms", "0.55", "-", 1.3),
        ("fixed_2ms", "Poisson test, fixed 2 ms", "0.38", "-", 1.3),
        ("fixed_3ms", "Poisson test, fixed 3 ms", "0.15", "-", 1.3),
        ("oracle", "Poisson test, true RP", "#1f78b4", "-", 2.0),
        ("sliding", "Sliding RP", "#1b9e77", "-", 2.2)]
GAMMA = 90.0


def roc(pos, neg):
    """(false acceptance, true acceptance) over every threshold, from the
    confidences of positives (8% trains) and negatives (12% trains)."""
    th = np.unique(np.concatenate([pos, neg, [np.inf]]))[::-1]
    tpr = np.array([np.mean(pos >= t) for t in th])
    fpr = np.array([np.mean(neg >= t) for t in th])
    return np.r_[0.0, fpr, 1.0], np.r_[0.0, tpr, 1.0]


def auc(pos, neg):
    """Mann-Whitney AUC with ties counting half."""
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method="average").to_numpy()
    rp = ranks[:pos.size].sum()
    return (rp - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size)


def tpr_at(pos, neg, target_fpr):
    """Best true acceptance with false acceptance <= target (no randomisation)."""
    th = np.unique(np.concatenate([pos, neg]))
    best = 0.0
    for t in th:
        if np.mean(neg >= t) <= target_fpr + 1e-12:
            best = max(best, np.mean(pos >= t))
    return best


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "power_roc.pqt")
    rates, rps = sorted(d.rate.unique()), sorted(d.rp.unique())
    fig, axs = plt.subplots(len(rps), len(rates), figsize=(3.1 * len(rates), 3.0 * len(rps) + 0.7),
                            sharex=True, sharey=True, squeeze=False)
    rows = []
    for i, rp in enumerate(rps):
        for j, fr in enumerate(rates):
            ax = axs[i, j]
            g = d[np.isclose(d.rp, rp) & np.isclose(d.rate, fr)]
            pos = {k: g[np.isclose(g.cont, 0.08)][k].to_numpy() for k, *_ in ARMS}
            neg = {k: g[np.isclose(g.cont, 0.12)][k].to_numpy() for k, *_ in ARMS}
            at10 = {k: g[np.isclose(g.cont, 0.10)][k].to_numpy() for k, *_ in ARMS}
            for k, lab, c, ls, lw in ARMS:
                x, y = roc(pos[k], neg[k])
                ax.plot(100 * x, 100 * y, ls, color=c, lw=lw,
                        label=lab if (i == 0 and j == 0) else None, drawstyle="steps-post")
            for k, mk in (("sliding", "o"), ("oracle", "s")):
                c = dict((a[0], a[2]) for a in ARMS)[k]
                ax.plot(100 * np.mean(neg[k] >= GAMMA), 100 * np.mean(pos[k] >= GAMMA), mk,
                        color=c, ms=7, mec="white", mew=1.2, zorder=5)
            ax.plot([0, 100], [0, 100], color="0.6", lw=0.8, ls=":", zorder=0)
            ax.set_title(f"{'abcdefghijkl'[i * len(rates) + j]}  {fr:g} spikes/s, RP {rp*1000:g} ms",
                         loc="left", fontsize=9)
            if j == 0:
                ax.set_ylabel("True acceptance at 8% (%)")
            if i == len(rps) - 1:
                ax.set_xlabel("False acceptance at 12% (%)")
            fa_s = np.mean(neg["sliding"] >= GAMMA)
            row = dict(rate=fr, rp_ms=rp * 1000,
                       sliding_fa=100 * fa_s, sliding_ta=100 * np.mean(pos["sliding"] >= GAMMA),
                       sliding_size10=100 * np.mean(at10["sliding"] >= GAMMA),
                       oracle_fa=100 * np.mean(neg["oracle"] >= GAMMA),
                       oracle_ta=100 * np.mean(pos["oracle"] >= GAMMA),
                       oracle_size10=100 * np.mean(at10["oracle"] >= GAMMA),
                       oracle_ta_at_sliding_fa=100 * tpr_at(pos["oracle"], neg["oracle"], fa_s))
            for k, *_ in ARMS:
                row[f"auc_{k}"] = auc(pos[k], neg[k])
            rows.append(row)
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=6,
               fontsize=8, frameon=False, handlelength=2.5, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.55 / (3.0 * len(rps) + 0.7)))
    plotstyle.save(fig, OUTDIR / "figures" / "power_roc")

    t = pd.DataFrame(rows)
    t.to_csv(OUTDIR / "power_roc.csv", index=False)
    L = ["Matched-false-acceptance comparison (ROC), 2 h, 2000 trains per point", "=" * 70, "",
         "At gamma = 90: Sliding RP's false acceptance (12%) and true acceptance (8%),",
         "the oracle's, and the oracle's true acceptance AT Sliding RP's false acceptance.",
         "Also the realised rate at exactly 10% contamination (the size).", "",
         f"{'rate':>5} {'RP':>4} | {'SRP FA':>7} {'SRP TA':>7} {'SRP @10%':>9} | {'orc FA':>7} "
         f"{'orc TA':>7} {'orc @10%':>9} | {'orc TA @ SRP FA':>16}"]
    for r in t.itertuples():
        L.append(f"{r.rate:5g} {r.rp_ms:4g} | {r.sliding_fa:6.1f}% {r.sliding_ta:6.1f}% "
                 f"{r.sliding_size10:8.1f}% | {r.oracle_fa:6.1f}% {r.oracle_ta:6.1f}% "
                 f"{r.oracle_size10:8.1f}% | {r.oracle_ta_at_sliding_fa:15.1f}%")
    L += ["", "Area under the ROC curve (0.5 = chance, 1 = perfect):", "",
          f"{'rate':>5} {'RP':>4} | " + " ".join(f"{k:>12s}" for k, *_ in ARMS)]
    for _, r in t.iterrows():
        L.append(f"{r['rate']:5g} {r['rp_ms']:4g} | " + " ".join(
            f"{r['auc_' + k]:12.3f}" for k, *_ in ARMS))
    (OUTDIR / "power_roc.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
