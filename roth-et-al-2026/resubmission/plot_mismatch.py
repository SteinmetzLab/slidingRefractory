"""Figure for work package 05: departures from the assumed generative model.

Reviewer 2 major 3 and Reviewer 1's correlation minor. The manuscript's
simulations use a hard-refractory renewal base neuron plus refractory-free
Poisson contamination, which is the model the statistic assumes. These are the
departures.

Panels
------
a  Acceptance against true contamination for the shape departures: the
   manuscript's own model, a contaminating neuron with its own refractory
   period, bursting, and graded refractory recovery.
b  Correlated rates: false and true acceptance against the rate correlation
   rho between base neuron and contaminant.
c  Non-overlapping activity: the same against the fraction of epochs in which
   both are active. This is the hippocampal place-cell failure mode.
d  Graded recovery after a hard refractory period: the same against the width
   of the relative-refractory ramp. Acceptance is essentially unchanged, and
   marginally higher for wide ramps.
e  Everything at once, against the manuscript's model as reference.

Run:  python plot_mismatch.py
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
              r"05_model_mismatch")
GAMMA = 90


def far(sub):
    """Acceptance at exactly the contamination threshold: false acceptance."""
    s = sub[np.isclose(sub.cont_prop, 0.10)]
    return s[f"sliding_pass_{GAMMA}"].mean() * 100 if len(s) else np.nan


def tar(sub):
    """Acceptance at 0.8x the threshold: true acceptance."""
    s = sub[np.isclose(sub.cont_prop, 0.08)]
    return s[f"sliding_pass_{GAMMA}"].mean() * 100 if len(s) else np.nan


def curve(ax, sub, label, color, ls="-"):
    g = sub.groupby("cont_prop")[f"sliding_pass_{GAMMA}"].mean() * 100
    ax.plot(g.index * 100, g.values, ls, color=color, marker="o", ms=3,
            label=label)


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "mismatch.pqt")
    # The graded-recovery rows in mismatch.pqt were produced by the first
    # version of gen_graded_rp, whose hazard was nonzero at every lag (27% of
    # baseline at zero lag for a 2 ms width) so the simulated neurons had no
    # absolute refractory period at all. Replace them with the corrected model
    # (hard refractory period, then a linear ramp) from run_graded.py.
    gf = SIMS / "graded_fixed.pqt"
    if gf.exists():
        g = pd.read_parquet(gf)
        d = pd.concat([d[d.model != "graded"], g[g.model == "graded"]],
                      ignore_index=True)
    for c in ("rho", "tau_s", "overlap", "width", "p_burst", "cont_rp"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)
    ref = d[d.model == "standard"]

    fig = plt.figure(figsize=(12, 6.8))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)

    # --- a: shape departures ----------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    curve(ax, ref, "Manuscript model", "k")
    curve(ax, d[(d.model == "single_neuron_contaminant") & (d.cont_rp == 0.0025)],
          "Single-neuron contaminant", "#d95f02")
    curve(ax, d[(d.model == "bursting") & (d.p_burst == 0.3)], "Bursting", "#7570b3")
    curve(ax, d[(d.model == "graded") & (d.width == 0.002)],
          "Graded recovery (2 ms ramp)", "#1b9e77")
    ax.axvline(10, color="0.6", ls="--", lw=1)
    ax.set_xlabel("True contamination (%)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("a  Departures in spike-train shape", loc="left")
    ax.legend(fontsize=7)

    # --- b: correlated rates ----------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    m = d[d.model == "modulated"]
    rhos = sorted(m.rho.dropna().unique())
    ax.plot(rhos, [far(m[m.rho == r]) for r in rhos], "o-", color="#d95f02",
            ms=5, label="False acceptance (10% contam.)")
    ax.plot(rhos, [tar(m[m.rho == r]) for r in rhos], "s-", color="#1b9e77",
            ms=5, label="True acceptance (8% contam.)")
    ax.axhline(far(ref), color="#d95f02", ls=":", lw=1)
    ax.axhline(tar(ref), color="#1b9e77", ls=":", lw=1)
    ax.set_xlabel("Rate correlation, rho")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("b  Correlated firing rates", loc="left")
    ax.legend(fontsize=7)

    # --- c: non-overlapping activity --------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    n = d[d.model == "nonoverlapping"]
    ovs = sorted(n.overlap.dropna().unique())
    ax.plot(ovs, [far(n[n.overlap == o]) for o in ovs], "o-", color="#d95f02", ms=5)
    ax.plot(ovs, [tar(n[n.overlap == o]) for o in ovs], "s-", color="#1b9e77", ms=5)
    ax.axhline(far(ref), color="#d95f02", ls=":", lw=1)
    ax.set_xlabel("Fraction of epochs shared")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("c  Non-overlapping activity", loc="left")
    ax.annotate("place-cell\nfailure mode", xy=(0.0, far(n[n.overlap == 0.0])),
                xytext=(0.18, 72), fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=0.8))

    # --- d: graded recovery -----------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    g = d[d.model == "graded"]
    ws = sorted(g.width.dropna().unique())
    ax.plot([w * 1000 for w in ws], [far(g[g.width == w]) for w in ws], "o-",
            color="#d95f02", ms=5, label="False acceptance")
    ax.plot([w * 1000 for w in ws], [tar(g[g.width == w]) for w in ws], "s-",
            color="#1b9e77", ms=5, label="True acceptance")
    ax.axhline(far(ref), color="#d95f02", ls=":", lw=1)
    ax.axhline(tar(ref), color="#1b9e77", ls=":", lw=1)
    ax.set_xlabel("Relative-refractory ramp width (ms)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("d  Graded recovery after a hard 2 ms period", loc="left")
    ax.legend(fontsize=7)

    # --- e: summary --------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1:])
    rows = [("Manuscript model", ref)]
    rows += [(f"Contaminant RP {c*1000:g} ms",
              d[(d.model == "single_neuron_contaminant") & (d.cont_rp == c)])
             for c in sorted(d[d.model == "single_neuron_contaminant"].cont_rp.dropna().unique())]
    rows += [(f"Bursting p={p:g}", d[(d.model == "bursting") & (d.p_burst == p)])
             for p in sorted(d[d.model == "bursting"].p_burst.dropna().unique())]
    rows += [(f"Graded width {w*1000:g} ms", d[(d.model == "graded") & (d.width == w)])
             for w in ws]
    rows += [(f"rho = {r:+g}", d[(d.model == "modulated") & (d.rho == r)])
             for r in rhos]
    rows += [(f"Overlap {o:g}", d[(d.model == "nonoverlapping") & (d.overlap == o)])
             for o in ovs]
    y = np.arange(len(rows))[::-1]
    ax.barh(y + 0.19, [far(s) for _, s in rows], 0.36, color="#d95f02",
            label="False acceptance (10% contam.)")
    ax.barh(y - 0.19, [tar(s) for _, s in rows], 0.36, color="#1b9e77",
            label="True acceptance (8% contam.)")
    ax.axvline(far(ref), color="#d95f02", ls=":", lw=1)
    ax.axvline(tar(ref), color="#1b9e77", ls=":", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("Units accepted (%)")
    ax.set_title("e  All departures against the manuscript's model "
                 "(dotted lines)", loc="left")
    ax.legend(fontsize=7, loc="lower right")

    plotstyle.save(fig, OUTDIR / "figures" / "mismatch")

    # --- numbers -----------------------------------------------------------
    L = ["Departures from the assumed generative model", "=" * 50, "",
         f"Sliding RP at the defaults, gamma = {GAMMA}. 'False acceptance' is the",
         "acceptance rate for units simulated at exactly the 10% contamination",
         "threshold; 'true acceptance' is the rate at 8%. Each entry pools 1000",
         "trains per contamination level over firing rates 1 and 5 spikes/s.", "",
         f"REFERENCE, the manuscript's own model: false {far(ref):.1f}%, "
         f"true {tar(ref):.1f}%", ""]
    for name, sub in rows[1:]:
        L.append(f"  {name:28s} false {far(sub):5.1f}%   true {tar(sub):5.1f}%")
    L += ["",
          "Reading:",
          "",
          "1. A contaminating neuron with its own refractory period is mildly",
          "   anti-conservative (false acceptance 22.9-24.6% against 20.5%).",
          "   This is expected: the Llobet expected-violation count includes a",
          "   contaminant-contaminant term that a single neuron cannot produce,",
          "   so the test over-estimates the violations it should have seen.",
          "",
          "2. Bursting barely matters (17.6-18.2%), slightly conservative. The",
          "   short-latency ACG peak sits above the refractory window the metric",
          "   finds, so it does not corrupt the decision.",
          "",
          "3. Correlated rates matter a great deal and in the direction the",
          "   Discussion predicts. Negative correlation hides contamination",
          "   (false acceptance 62.7% at rho = -1); positive correlation makes",
          "   the test conservative (1.8% at rho = +1) but destroys power (true",
          "   acceptance falls to 4.5%). Note that rho is the underlying RATE",
          "   correlation; see correlation_bridge.txt for how that maps onto a",
          "   measurable spike-count correlation, which is much smaller.",
          "",
          "4. Non-overlapping activity is the most severe failure and it is the",
          "   one the manuscript names but never tested. With no shared epochs,",
          "   89.1% of units at the contamination threshold are accepted. The",
          "   transition is sharp: 63.3% at 25% overlap, 11.9% at 50%.",
          "",
          "5. Graded refractory recovery is a NON-issue, once it is modelled",
          "   correctly. With an absolute refractory period followed by a",
          "   linear ramp to baseline, acceptance is indistinguishable from the",
          "   hard-RP reference and if anything slightly higher: at 5 spikes/s,",
          "   true acceptance runs 77.6% (hard) to 77.0 / 78.1 / 80.0% for ramps",
          "   of 0.5 / 1 / 2 ms, with false acceptance flat at 26-29%. The extra",
          "   bins between the absolute RP and full recovery carry fewer",
          "   violations than baseline, so they give the sliding search a few",
          "   more windows that might work, each less likely to help than the",
          "   last -- exactly as one would predict.",
          "",
          "   An earlier version of this analysis reported that graded recovery",
          "   made the metric uninformative. That was an artifact: the first",
          "   generator used a logistic hazard centred on the refractory period,",
          "   which leaves 27% of baseline firing at zero lag and a minimum ISI",
          "   of zero. Those units genuinely violate at every lag, so rejecting",
          "   them was correct behaviour and said nothing about graded recovery."]
    (OUTDIR / "mismatch_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
