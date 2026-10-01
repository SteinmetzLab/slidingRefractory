"""Work package 04: what does the sliding actually buy?

Reviewer 2 major 2 asks us to separate the benefit of Sliding RP from the
simpler benefit of not using a badly misspecified refractory period. The arms,
all evaluated on the *same* simulated trains so differences are not simulation
noise:

  Hill-Llobet, 3 ms     the manuscript's comparison: a fixed conventional value
  Hill-Llobet, oracle   the true simulated refractory period (an upper bound on
                        what any RP-estimation scheme could achieve)
  Hill-Llobet, estimated  the manuscript's own ACG estimator, per unit
  Poisson test, oracle  the Sliding RP statistic at the true RP, without sliding
  Poisson test, estimated  the same at the estimated RP
  Sliding RP            the method as published

The two Poisson-test arms are the ones that isolate the question. Comparing
them with the Hill-Llobet arms at the same tau separates "treat the count
statistically" from "use the right tau"; comparing them with Sliding RP
separates "use the right tau" from "search over tau".

Run:  python plot_decomposition.py
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
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"04_hill_llobet_decomposition")
G = 90

ARMS = [
    ("hl3_pass", "Hill-Llobet, fixed 3 ms", "#7570b3", "-"),
    ("hl_oracle_pass", "Hill-Llobet, true RP", "#d95f02", "-"),
    ("hl_est_pass", "Hill-Llobet, estimated RP", "#d95f02", "--"),
    (f"pt_oracle_pass_{G}", "Poisson test, true RP", "#1f78b4", "-"),
    (f"pt_est_pass_{G}", "Poisson test, estimated RP", "#1f78b4", "--"),
    (f"sliding_pass_{G}", "Sliding RP", "#1b9e77", "-"),
]
MARK = {"-": "o", "--": "s"}
RATES = (0.5, 1.0, 2.0, 5.0)


def draw_arm(ax, x, y, c, ls, lab=None, ms=3.5):
    """One arm: solid line + filled circle, or dashed line + open square."""
    mk = MARK[ls]
    ax.plot(x, y, ls=ls, color=c, lw=1.8 if ls == "-" else 1.6,
            dashes=(5, 2.5) if ls == "--" else (None, None),
            marker=mk, ms=ms, mfc=c if mk == "o" else "white", mec=c, mew=1.1, label=lab)


FAR_C, TAR_C = 0.12, 0.08     # false / true acceptance points (Fig 4 convention)

def op(sub, col, cont):
    s = sub[np.isclose(sub.cont_prop, cont)]
    return s[col].mean() * 100 if len(s) and col in s else np.nan


def main():
    plotstyle.apply()
    p = SIMS / "decomposition.pqt"
    if not p.exists():
        print("decomposition.pqt not found; run:  python run_sims.py decomposition")
        return
    d = pd.read_parquet(p)
    p05 = SIMS / "decomposition_0p5.pqt"
    if p05.exists():
        d = pd.concat([d, pd.read_parquet(p05)], ignore_index=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)

    # --- full acceptance curves: every firing rate x true RP --------------
    rates = [r for r in RATES if r in set(d.total_rate)]
    fig, axs = plt.subplots(len(rates), 2, figsize=(10, 2.55 * len(rates) + 0.9),
                            sharex=True, sharey=True, squeeze=False)
    for i_r, fr in enumerate(rates):
        for j_c, rp in enumerate((0.0015, 0.003)):
            ax = axs[i_r, j_c]
            sub = d[np.isclose(d.rp_dur, rp) & np.isclose(d.total_rate, fr)]
            for col, lab, c, ls in ARMS:
                if col not in sub:
                    continue
                g = sub.groupby("cont_prop")[col].mean() * 100
                draw_arm(ax, g.index * 100, g.values, c, ls,
                         lab if (i_r == 0 and j_c == 0) else None)
            ax.axvline(10, color="0.6", ls=":", lw=1)
            ax.axvline(FAR_C * 100, color="#c2410c", ls=":", lw=0.8, alpha=0.6)
            ax.axvline(TAR_C * 100, color="#1b9e77", ls=":", lw=0.8, alpha=0.6)
            ax.set_title(f"{'abcdefgh'[2 * i_r + j_c]}  {fr:g} spikes/s, true RP "
                         f"{rp*1000:g} ms" + ("  (3 ms is too long)" if rp < 0.003 else ""),
                         loc="left")
            if j_c == 0:
                ax.set_ylabel("Units accepted (%)")
            if i_r == len(rates) - 1:
                ax.set_xlabel("True contamination (%)")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=3,
               fontsize=8.5, handlelength=3.2, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.75 / (2.55 * len(rates) + 0.9)))
    plotstyle.save(fig, OUTDIR / "figures" / "decomposition_curves")

    # --- summary panels -----------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
    axs = [None, None] + list(axs)

    # --- c: operating points by firing rate --------------------------------
    ax = axs[2]
    for col, lab, c, ls in ARMS:
        if col not in d:
            continue
        xs, ys = [], []
        for fr in sorted(d.total_rate.unique()):
            # one true RP (the Fig 4 setting), not pooled across 1.5, 3 and 5 ms
            sub = d[(d.total_rate == fr) & np.isclose(d.rp_dur, 0.003)]
            xs.append(op(sub, col, FAR_C))
            ys.append(op(sub, col, TAR_C))
        draw_arm(ax, xs, ys, c, ls, lab, ms=5)
        ax.annotate(f"{sorted(d.total_rate.unique())[0]:g}", (xs[0], ys[0]),
                    fontsize=5.5, color=c, xytext=(2, 2),
                    textcoords="offset points")
    ax.plot([0, 100], [0, 100], "k--", lw=1)
    ax.set_xlabel("False acceptance rate (%)")
    ax.set_ylabel("True acceptance rate (%)")
    ax.set_title("a  Operating points, true RP 3 ms", loc="left")
    ax.text(0.98, 0.04, "each arm traced over 0.5, 1, 2, 5 spikes/s\n(labelled at its lowest rate)",
            transform=ax.transAxes, ha="right", fontsize=7, color="0.35")
    ax.legend(fontsize=7, handlelength=3.2, loc="center right")

    # --- d: estimator error -------------------------------------------------
    ax = axs[3]
    if "rp_est_ms_median" in d:
        for rp, c in zip(sorted(d.rp_dur.unique()),
                         plotstyle.confidence_cmap(len(d.rp_dur.unique()))):
            sub = d[d.rp_dur == rp]
            g = sub.groupby("total_rate")["rp_est_ms_median"].median()
            ax.plot(g.index, g.values - rp * 1000, "o-", color=c, ms=4,
                    label=f"{rp*1000:g}")
        ax.axhline(0, color="0.6", ls=":", lw=1)
        ax.set_xscale("log")
        ax.set_xlabel("Firing rate (spikes/s)")
        ax.set_ylabel("Estimated minus true RP (ms)")
        ax.set_title("b  Estimator bias (fails at 0.5 spikes/s)", loc="left")
        plotstyle.plain_log_ticks(ax, ticks=(0.5, 1, 2, 5))
        ax.legend(title="True RP (ms)", fontsize=6.5)

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "decomposition")

    # --- numbers ------------------------------------------------------------
    L = ["Decomposing the Hill-Llobet comparison", "=" * 45, "",
         f"Sliding RP and the Poisson-test arms at gamma = {G}; 600 trains per",
         "contamination level; all arms on the same trains. 'False' is the",
         f"acceptance rate at {FAR_C*100:g}% contamination, 'true' at {TAR_C*100:g}% (the",
         "manuscript's Fig 4 convention).", ""]
    for rp in sorted(d.rp_dur.unique()):
        L += [f"True refractory period {rp*1000:g} ms:", ""]
        L.append(f"  {'arm':28s} " + "".join(
            f"{f'FR {fr:g}: false/true':>22}" for fr in sorted(d.total_rate.unique())))
        for col, lab, _, _ in ARMS:
            if col not in d:
                continue
            cells = []
            for fr in sorted(d.total_rate.unique()):
                sub = d[(d.rp_dur == rp) & (d.total_rate == fr)]
                cells.append(f"{op(sub, col, FAR_C):8.1f} /{op(sub, col, TAR_C):7.1f}")
            L.append(f"  {lab:28s} " + "".join(f"{c:>22}" for c in cells))
        L.append("")
    if "rp_est_ms_median" in d:
        L += ["Estimator accuracy (median estimated minus true RP, ms):", ""]
        for rp in sorted(d.rp_dur.unique()):
            row = []
            for fr in sorted(d.total_rate.unique()):
                sub = d[(d.rp_dur == rp) & (d.total_rate == fr)]
                row.append(f"FR {fr:g}: {sub.rp_est_ms_median.median() - rp*1000:+.3f}")
            L.append(f"  true RP {rp*1000:4.1f} ms   " + "   ".join(row))
    (OUTDIR / "decomposition_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
