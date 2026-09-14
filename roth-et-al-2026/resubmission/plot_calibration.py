"""Figures for work package 03: calibration of the confidence parameter.

Panels
------
a  Realised false-acceptance rate vs nominal (1 - gamma), on log-log axes with
   the identity line, one trace per true refractory duration.
b  The same, split by firing rate: the discrepancy grows with statistical power.
c  Ratio realised/nominal as a function of gamma, showing that the gap is
   proportionally worst exactly where users operate (gamma = 90 to 99).
d  ROC-style operating points: true acceptance at 0.8*C_thresh against false
   acceptance at C_thresh, for the confidence family and for Hill-Llobet.
e  Corrected test: realised false-acceptance rate vs true refractory duration,
   testing whether the full-window null is least favorable (see
   03_confidence_calibration/least_favorable_null.md).
f  The power cost of the correction.

Run:  python plot_calibration.py
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
FIGS = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
            r"03_confidence_calibration/figures")
GAMMAS = np.array([50, 60, 70, 75, 80, 85, 90, 95, 99])


def _rates(df, col_tpl=lambda g: f"sliding_pass_{g}"):
    return np.array([df[col_tpl(g)].mean() * 100 for g in GAMMAS])


def main():
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)
    cal = pd.read_parquet(SIMS / "calibration.pqt")
    boundary = cal[cal.cont_prop == 0.10]
    power = cal[cal.cont_prop == 0.08]
    nominal = 100 - GAMMAS

    fig = plt.figure(figsize=(10.5, 6.6))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.36)

    # --- a: realised vs nominal, by true RP -------------------------------
    ax = fig.add_subplot(gs[0, 0])
    rps = sorted(boundary.rp_dur.unique())
    cols = plotstyle.confidence_cmap(len(rps))
    for rp, c in zip(rps, cols):
        sel = boundary[(boundary.rp_dur == rp) & (boundary.total_rate >= 2)
                       & (boundary.rec_dur >= 3600)]
        ax.plot(nominal, _rates(sel), "o-", color=c, ms=3.5,
                label=f"{rp*1000:g}")
    ax.plot([1, 50], [1, 50], "k--", lw=1, label="Nominal")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Nominal false acceptance rate (%)")
    ax.set_ylabel("Realized false acceptance rate (%)")
    ax.set_title("a  Calibration by refractory duration", loc="left")
    leg = ax.legend(title="True RP (ms)", loc="upper left", ncol=2)
    leg._legend_box.align = "left"

    # --- b: by firing rate -------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    frs = sorted(boundary.total_rate.unique())
    cols = plotstyle.confidence_cmap(len(frs))
    for fr, c in zip(frs, cols):
        sel = boundary[(boundary.total_rate == fr) & (boundary.rec_dur >= 3600)]
        ax.plot(nominal, _rates(sel), "o-", color=c, ms=3.5, label=f"{fr:g}")
    ax.plot([1, 50], [1, 50], "k--", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Nominal false acceptance rate (%)")
    ax.set_ylabel("Realized false acceptance rate (%)")
    ax.set_title("b  Calibration by firing rate", loc="left")
    leg = ax.legend(title="Firing rate (spikes/s)", loc="upper left", ncol=2)
    leg._legend_box.align = "left"

    # --- c: ratio vs gamma -------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    for rp, c in zip(rps, plotstyle.confidence_cmap(len(rps))):
        sel = boundary[(boundary.rp_dur == rp) & (boundary.total_rate >= 2)
                       & (boundary.rec_dur >= 3600)]
        ax.plot(GAMMAS, _rates(sel) / nominal, "o-", color=c, ms=3.5,
                label=f"{rp*1000:g}")
    ax.axhline(1, color="k", ls="--", lw=1)
    ax.set_xlabel("Confidence threshold (%)")
    ax.set_ylabel("Realized / nominal rate")
    ax.set_title("c  Inflation factor", loc="left")
    ax.legend(title="True RP (ms)", ncol=2)

    # --- d: ROC-style operating points -------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    b = boundary[(boundary.total_rate >= 2) & (boundary.rec_dur >= 3600)]
    p = power[(power.total_rate >= 2) & (power.rec_dur >= 3600)]
    fa, ta = _rates(b), _rates(p)
    cols = plotstyle.confidence_cmap(len(GAMMAS))
    ax.plot(fa, ta, "-", color="0.6", zorder=1)
    for i, g in enumerate(GAMMAS):
        ax.plot(fa[i], ta[i], "o", color=cols[i], ms=6, zorder=2,
                label=f"{g}" if g in (50, 90, 99) else None)
    ax.plot(b.hl2_pass.mean() * 100, p.hl2_pass.mean() * 100, "rx", ms=9,
            mew=2, label="Hill-Llobet, 2 ms")
    ax.plot(b.hl3_pass.mean() * 100, p.hl3_pass.mean() * 100, "r+", ms=11,
            mew=2, label="Hill-Llobet, 3 ms")
    ax.plot([0, 100], [0, 100], "k--", lw=1)
    ax.set_xlabel("False acceptance rate (%), 10% contamination")
    ax.set_ylabel("True acceptance rate (%), 8% contamination")
    ax.set_title("d  Operating points", loc="left")
    ax.legend(loc="lower right")

    # --- e, f: the corrected test ------------------------------------------
    corr_path = SIMS / "corrected.pqt"
    if corr_path.exists():
        cr = pd.read_parquet(corr_path)
        crb = cr[cr.cont_prop == 0.10]
        crp = cr[cr.cont_prop == 0.08]

        ax = fig.add_subplot(gs[1, 1])
        for g, c in zip((80, 90, 95), plotstyle.confidence_cmap(3)):
            y = [crb[crb.rp_dur == rp][f"corrected_pass_{g}"].mean() * 100
                 for rp in sorted(crb.rp_dur.unique())]
            ax.plot(np.array(sorted(crb.rp_dur.unique())) * 1000, y, "o-",
                    color=c, ms=4, label=f"{g}")
            ax.axhline(100 - g, color=c, ls=":", lw=1)
        ax.set_xlabel("True refractory duration (ms)")
        ax.set_ylabel("Realized false acceptance rate (%)")
        ax.set_title("e  FWER-corrected test", loc="left")
        ax.legend(title="Confidence (%)")

        ax = fig.add_subplot(gs[1, 2])
        rr = sorted(crb.rp_dur.unique())
        ax.plot(np.array(rr) * 1000,
                [crp[crp.rp_dur == rp]["sliding_pass_90"].mean() * 100 for rp in rr],
                "o-", color="0.3", ms=4, label="Standard")
        ax.plot(np.array(rr) * 1000,
                [crp[crp.rp_dur == rp]["corrected_pass_90"].mean() * 100 for rp in rr],
                "s-", color="#1b9e77", ms=4, label="Corrected")
        ax.set_xlabel("True refractory duration (ms)")
        ax.set_ylabel("True acceptance rate (%)")
        ax.set_title("f  Power cost of the correction", loc="left")
        ax.legend()
    else:
        for i, msg in ((1, "corrected.pqt not yet written"), (2, "")):
            ax = fig.add_subplot(gs[1, i])
            ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=8)
            ax.axis("off")

    plotstyle.save(fig, FIGS / "calibration")

    # --- a compact numeric summary for the response letter -----------------
    lines = ["Realized false-acceptance rate at C = C_thresh",
             "(FR >= 2 spikes/s, duration >= 1 h; n = 4000 per condition)", ""]
    lines.append(f"{'gamma':>6} {'nominal':>9} {'realized':>10} {'ratio':>7}")
    for g, nom, real in zip(GAMMAS, nominal, _rates(b)):
        lines.append(f"{g:>6} {nom:>8.0f}% {real:>9.1f}% {real/nom:>6.2f}x")
    m = boundary[(boundary.total_rate == 5) & (boundary.rec_dur == 7200)
                 & (boundary.rp_dur == 0.003)]
    lines += ["", "Manuscript condition (5 spikes/s, 2 h, RP 3 ms, gamma 90): "
              f"{m.sliding_pass_90.iloc[0]*100:.1f}% "
              f"(95% CI {m.sliding_pass_90_lo.iloc[0]*100:.1f}-"
              f"{m.sliding_pass_90_hi.iloc[0]*100:.1f}%)"]
    (FIGS.parent / "calibration_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
