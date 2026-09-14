"""Figure for the semi-synthetic validation (work package 05, Reviewer 2 major 3).

Panels
------
a  Pass rate against injected contamination for Sliding RP and Hill-Llobet at
   2 and 3 ms. The leftmost point is the uncontaminated recipient, where the
   fixed-RP methods reject the majority of units that are clean by construction.
b  Minimum confirmable contamination against the known injected fraction.
c  The realized rate correlation between recipient and donor in these real
   recordings, for nearby and distant donors.
d  Pass rate at the contamination threshold as a function of that correlation:
   the assumption of uncorrelated rates matters, and in the direction the
   manuscript predicts.

Run:  python plot_semisynthetic.py
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


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "semisynthetic.pqt")
    r = d[d["mode"] == "random"].copy()
    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(1, 4, figsize=(13.5, 3.2))

    # --- a: pass rate vs injected fraction --------------------------------
    ax = axs[0]
    f = np.sort(r.f_injected.unique()) * 100
    for col, lab, c in (("passes", "Sliding RP", "#1b9e77"),
                        ("hl2_pass", "Hill-Llobet, 2 ms", "#d95f02"),
                        ("hl3_pass", "Hill-Llobet, 3 ms", "#7570b3")):
        y = [r[r.f_injected == v / 100][col].mean() for v in f]
        ax.plot(f, np.array(y) * 100, "o-", color=c, ms=4, label=lab)
    ax.axvline(10, color="0.6", ls="--", lw=1)
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Units passing (%)")
    ax.set_title("a  Semi-synthetic injection", loc="left")
    ax.legend()

    # --- b: estimate accuracy ---------------------------------------------
    ax = axs[1]
    base = r.recipient_base_cont.median()
    med = [r[r.f_injected == v / 100].min_cont.median() for v in f]
    q1 = [r[r.f_injected == v / 100].min_cont.quantile(.25) for v in f]
    q3 = [r[r.f_injected == v / 100].min_cont.quantile(.75) for v in f]
    ax.fill_between(f, q1, q3, color="#1b9e77", alpha=0.25, lw=0)
    ax.plot(f, med, "o-", color="#1b9e77", ms=4, label="Estimated (C_min)")
    ax.plot(f, f + base, "k--", lw=1,
            label=f"Injected + recipient baseline ({base:.1f}%)")
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Minimum confirmable contamination (%)")
    ax.set_title("b  Estimate against known ground truth", loc="left")
    ax.legend()

    # --- c: realised correlation in real data ------------------------------
    ax = axs[2]
    pairs = d.drop_duplicates(["pid", "recipient", "donor"])
    bins = np.linspace(-0.4, 0.5, 28)
    for k, c in (("near", "#d95f02"), ("far", "#1f78b4")):
        v = pairs[pairs.donor_kind == k].rate_corr.dropna()
        n, e = np.histogram(v, bins=bins)
        ax.stairs(n / n.sum(), e, color=c, lw=1.6,
                  label=f"{k} (median {v.median():+.3f})")
    ax.axvline(0, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Rate correlation, 100 ms bins")
    ax.set_ylabel("Proportion of pairs")
    ax.set_title("c  How correlated are real neighbours?", loc="left")
    ax.legend()

    # --- d: pass rate vs correlation at threshold --------------------------
    ax = axs[3]
    at = r[r.f_injected == 0.10].copy()
    edges = [-1, -0.02, 0.02, 0.08, 0.2, 1]
    at["bin"] = pd.cut(at.rate_corr, edges)
    xs, ys, ns = [], [], []
    for b, g in at.groupby("bin", observed=True):
        xs.append(g.rate_corr.median())
        ys.append(g.passes.mean() * 100)
        ns.append(len(g))
    ax.plot(xs, ys, "o-", color="#1b9e77", ms=5)
    for x, y, n in zip(xs, ys, ns):
        ax.annotate(f"n={n}", (x, y), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=6.5)
    ax.axvline(0, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Rate correlation, 100 ms bins")
    ax.set_ylabel("Units passing (%) at 10% injected")
    ax.set_title("d  Correlation shifts the decision", loc="left")

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "semisynthetic")

    # --- numbers -----------------------------------------------------------
    lines = ["Semi-synthetic contamination on real IBL recordings",
             "=" * 60, "",
             f"{d.recipient.nunique()} recipient units over {d.pid.nunique()} probe "
             f"insertions; {len(d):,} injected trains.",
             "Recipients were required to pass Sliding RP at a strict setting",
             "(5% contamination, 99% confidence) and fire at least 2 spikes/s,",
             "so they are clean by construction with a quantified baseline",
             f"(median C_min {base:.2f}%).", "",
             "Pass rate against injected contamination:", ""]
    lines.append(f"{'injected':>9} {'Sliding RP':>11} {'HL 2 ms':>9} {'HL 3 ms':>9} "
                 f"{'C_min median':>13}")
    for v in f:
        g = r[r.f_injected == v / 100]
        lines.append(f"{v:>8.0f}% {g.passes.mean()*100:>10.1f}% "
                     f"{g.hl2_pass.mean()*100:>8.1f}% {g.hl3_pass.mean()*100:>8.1f}% "
                     f"{g.min_cont.median():>12.2f}%")
    z = d[d.f_injected == 0]
    lines += ["",
              "The uncontaminated row is the striking one: on units that are clean",
              f"by construction, Sliding RP passes {z.passes.mean():.0%} while",
              f"Hill-Llobet rejects {1-z.hl2_pass.mean():.0%} of them at 2 ms and",
              f"{1-z.hl3_pass.mean():.0%} at 3 ms. Their median firing rate is",
              f"{z.fr_recipient.median():.1f} spikes/s, so this is not a low-power effect:",
              "it is the fixed-RP assumption rejecting neurons whose refractory",
              "period is shorter than assumed.", "",
              "Estimate accuracy: the minimum confirmable contamination tracks the",
              "known injected fraction plus the recipient's own baseline to within",
              "about one percentage point up to 12% injected, then under-estimates",
              "(16.2% estimated at 20% injected).", "",
              "Rate correlation between recipient and donor (100 ms bins):", ""]
    for k, g in pairs.groupby("donor_kind"):
        c = g.rate_corr
        lines.append(f"  {k:5s}: median {c.median():+.4f}, IQR "
                     f"[{c.quantile(.25):+.4f}, {c.quantile(.75):+.4f}], "
                     f"range [{c.min():+.3f}, {c.max():+.3f}], n = {len(c)}")
    lines += ["",
              "Effect of that correlation at 10% injected contamination:", ""]
    for b, g in at.groupby("bin", observed=True):
        lines.append(f"  r in {str(b):16s} n={len(g):4d}  passing "
                     f"{g.passes.mean()*100:5.1f}%  C_min median {g.min_cont.median():5.2f}%")
    lines += ["",
              "This is the assumption of uncorrelated rates being tested directly on",
              "real data, and it matters: over a correlation range of only about",
              "-0.1 to +0.2, which is entirely ordinary for neighbouring neurons,",
              "the pass rate for identically contaminated units runs from 69% down",
              "to 21%, and the contamination estimate from 8.4% to 13.9%. Negative",
              "correlation hides contamination and positive correlation makes the",
              "test conservative, exactly as the Discussion predicts, but the effect",
              "is larger than the text implies and should be stated quantitatively."]
    (OUTDIR / "semisynthetic_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
