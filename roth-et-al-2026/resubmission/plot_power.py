"""Figure and lookup table for work package 06: how much data is enough.

Replaces the simulated Fig 4g with the analytical curve and writes the
rate x duration lookup table the
reviewer asks for ("whether a 0.2-spike/s neuron recorded for one hour can
meaningfully be evaluated with this metric at all": it cannot).

Panels
------
a  Minimum firing rate for acceptance vs recording duration, one trace per
   confidence threshold: the analytical replacement for Fig 4g.
b  The same, by assumed refractory duration: short-RP units need more spikes,
   which is the reviewer's point about the populations that motivate the paper.
c  tau_pass0, the shortest violation-free window that would let a unit be
   accepted, against firing rate; the shaded band above 10 ms marks units that
   cannot be accepted whatever their autocorrelogram looks like.

The simulation check on panel a lives in its own figure,
plot_fig4g_validation.py, since it validates the closed form rather than
reporting a result.

Run:  python plot_power.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402

from slidingRP.power import min_passing_fr, tau_pass0  # noqa: E402

SIMS = Path(r"D:/temp/slidingRP_resub/sims")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"06_underpowered_outcome")


def main():
    plotstyle.apply()
    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    dur_h = np.linspace(0.25, 5, 120)

    # --- a: analytical Fig 4g ---------------------------------------------
    ax = axs[0]
    confs = [50, 60, 70, 80, 90, 95, 99]
    cols = plotstyle.confidence_cmap(len(confs))
    for g, c in zip(confs, cols):
        ax.plot(dur_h, min_passing_fr(dur_h * 3600, 0.003, conf_thresh=g),
                color=c, label=f"{g}")
    # simulation overlay: the firing rate at which uncontaminated simulated
    # units cross 50% acceptance (full validation in plot_fig4g_validation.py)
    val = SIMS / "fig4g_validation.pqt"
    if val.exists():
        v = pd.read_parquet(val)
        v = v[v.rp_dur == 0.003]
        first = True
        for rd, grp in v.groupby("rec_dur"):
            grp = grp.sort_values("total_rate")
            above = grp[grp.sliding_pass_90 >= 0.5]
            if len(above):
                ax.plot(rd / 3600, above.total_rate.iloc[0], "kv", ms=6,
                        zorder=5, label="Simulation" if first else None)
                first = False
    ax.set_xlabel("Recording duration (h)")
    ax.set_ylabel("Minimum firing rate for acceptance (spikes/s)")
    ax.set_title("a  Power at a 3 ms refractory period", loc="left")
    ax.legend(title="Confidence (%)", ncol=2)

    # --- b: by assumed refractory duration ---------------------------------
    ax = axs[1]
    taus = [0.001, 0.0015, 0.002, 0.003, 0.005]
    for t, c in zip(taus, plotstyle.confidence_cmap(len(taus))):
        ax.plot(dur_h, min_passing_fr(dur_h * 3600, t), color=c,
                label=f"{t*1000:g}")
    ax.set_xlabel("Recording duration (h)")
    ax.set_ylabel("Minimum firing rate for acceptance (spikes/s)")
    ax.set_title("b  Shorter refractory periods need more spikes", loc="left")
    ax.legend(title="Refractory period (ms)", ncol=2)

    # --- c: tau_pass0 vs firing rate ---------------------------------------
    ax = axs[2]
    fr = np.logspace(np.log10(0.1), np.log10(20), 200)
    for dh, c in zip([0.5, 1, 2, 4], plotstyle.confidence_cmap(4)):
        ax.plot(fr, tau_pass0(fr * dh * 3600, dh * 3600) * 1000, color=c,
                label=f"{dh:g}")
    ax.axhspan(10, 1e4, color="0.85", zorder=0)
    ax.text(0.13, 22, "Cannot be accepted: beyond the tested window",
            fontsize=7, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(0.01, 1000)
    ax.set_xlabel("Firing rate (spikes/s)")
    ax.set_ylabel("Shortest clean window for acceptance, tau_pass0 (ms)")
    ax.set_title("c  Is this unit evaluable at all?", loc="left")
    ax.legend(title="Duration (h)", ncol=2)

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "power")

    # --- lookup table ------------------------------------------------------
    rows = []
    for dh in (0.5, 1, 2, 4):
        for t in (0.001, 0.002, 0.003):
            for g in (80, 90, 95):
                rows.append(dict(duration_h=dh, tau_ms=t * 1000, confidence=g,
                                 min_fr=round(float(min_passing_fr(
                                     dh * 3600, t, conf_thresh=g)), 3)))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUTDIR / "min_passing_fr_table.csv", index=False)

    lines = ["Minimum firing rate (spikes/s) for a unit with no observed",
             "refractory period violations to be accepted, at a 10% contamination",
             "threshold. Below this rate the metric cannot establish acceptable",
             "contamination however clean the autocorrelogram is.", ""]
    piv = tab.pivot_table(index="duration_h", columns=["tau_ms", "confidence"],
                          values="min_fr")
    lines.append(piv.to_string())
    lines += ["", "tau_pass0: the shortest violation-free window that would let a",
              "unit be accepted (10% contamination, 90% confidence, 1 h recording):", ""]
    for f in (0.2, 0.3, 0.5, 0.8, 1.0, 1.1, 1.5, 2.0, 5.0):
        t = tau_pass0(f * 3600, 3600) * 1000
        tag = "  <- cannot be accepted (window is 10 ms)" if t > 10 else ""
        lines.append(f"  {f:4.1f} spikes/s -> {t:8.2f} ms{tag}")
    lines += ["",
              "The reviewer's example: a 0.2 spikes/s unit recorded for one hour",
              "would need a violation-free window of 84 ms, eight times the",
              "longest refractory period the algorithm tests. It cannot be",
              "evaluated by this metric, and no choice of parameters changes that."]
    (OUTDIR / "power_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
