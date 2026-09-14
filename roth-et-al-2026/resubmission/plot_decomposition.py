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
    ("hl3_pass", "Hill-Llobet, 3 ms", "#7570b3", "-"),
    ("hl_oracle_pass", "Hill-Llobet, oracle RP", "#d95f02", "-"),
    ("hl_est_pass", "Hill-Llobet, estimated RP", "#d95f02", "--"),
    (f"pt_oracle_pass_{G}", "Poisson test, oracle RP", "#1f78b4", "-"),
    (f"pt_est_pass_{G}", "Poisson test, estimated RP", "#1f78b4", "--"),
    (f"sliding_pass_{G}", "Sliding RP", "#1b9e77", "-"),
]


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
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)

    fig, axs = plt.subplots(1, 4, figsize=(14, 3.4))

    # --- a, b: acceptance curves where 3 ms is wrong, and where it is right
    for ax, rp, tag in ((axs[0], 0.0015, "a"), (axs[1], 0.003, "b")):
        sub = d[(d.rp_dur == rp) & (d.total_rate == 5.0)]
        for col, lab, c, ls in ARMS:
            if col not in sub:
                continue
            g = sub.groupby("cont_prop")[col].mean() * 100
            ax.plot(g.index * 100, g.values, ls, color=c, marker="o", ms=3,
                    label=lab)
        ax.axvline(10, color="0.6", ls=":", lw=1)
        ax.set_xlabel("True contamination (%)")
        ax.set_ylabel("Units accepted (%)")
        ax.set_title(f"{tag}  True RP {rp*1000:g} ms, 5 spikes/s", loc="left")
        if tag == "a":
            ax.legend(fontsize=6.5)

    # --- c: operating points by firing rate --------------------------------
    ax = axs[2]
    for col, lab, c, ls in ARMS:
        if col not in d:
            continue
        xs, ys = [], []
        for fr in sorted(d.total_rate.unique()):
            sub = d[d.total_rate == fr]
            xs.append(op(sub, col, 0.10))
            ys.append(op(sub, col, 0.08))
        ax.plot(xs, ys, ls, color=c, marker="o", ms=5, alpha=0.85, label=lab)
        ax.annotate(f"{sorted(d.total_rate.unique())[0]:g}", (xs[0], ys[0]),
                    fontsize=5.5, color=c, xytext=(2, 2),
                    textcoords="offset points")
    ax.plot([0, 100], [0, 100], "k--", lw=1)
    ax.set_xlabel("False acceptance rate (%)")
    ax.set_ylabel("True acceptance rate (%)")
    ax.set_title("c  Operating points (1, 2, 5 spikes/s)", loc="left")

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
        ax.set_title("d  Estimator bias", loc="left")
        ax.legend(title="True RP (ms)", fontsize=6.5)

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "decomposition")

    # --- numbers ------------------------------------------------------------
    L = ["Decomposing the Hill-Llobet comparison", "=" * 45, "",
         f"Sliding RP and the Poisson-test arms at gamma = {G}; 600 trains per",
         "contamination level; all arms on the same trains. 'False' is the",
         "acceptance rate at exactly 10% contamination, 'true' at 8%.", ""]
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
                cells.append(f"{op(sub, col, 0.10):8.1f} /{op(sub, col, 0.08):7.1f}")
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
