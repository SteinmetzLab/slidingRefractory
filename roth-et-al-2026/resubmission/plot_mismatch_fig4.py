"""Model-mismatch figure, redone at the manuscript's Fig 4 setting (05).

The first version pooled 1 and 5 spikes/s and 2 and 3 ms refractory periods
into each curve, which is why its "manuscript model" curve was not sigmoidal,
and its correlated-rate model was mislabeled (rho scaled the contaminant's
modulation depth on a shared signal; it was not a correlation). Here every
curve is one condition at 5 spikes/s, 3 ms, 1 h, false acceptance is read at
12% contamination and true acceptance at 8% (the manuscript's Fig 4
convention), and the correlated model has a real correlation coefficient.

a  Shape departures, one curve each.
b  Correlated rates: acceptance against contamination for rho from -1 to +1.
c  False and true acceptance against kappa, the excess neuron-contaminant
   coincidence, for both correlated models, with the prediction from kappa
   alone.
d  Non-overlapping activity.
e  Graded recovery after a hard 3 ms period.
f  Where real neighbouring IBL neurons sit on the kappa axis.

Run:  python plot_mismatch_fig4.py
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
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/05_model_mismatch")
G = "sliding_pass_90"
FAR_C, TAR_C = 0.12, 0.08


def at(sub, c):
    s = sub[np.isclose(sub.cont_prop, c)]
    return 100 * s[G].mean() if len(s) else np.nan


def curve(ax, sub, color, label, lw=1.6, ls="-", ms=3):
    g = sub.sort_values("cont_prop")
    ax.plot(g.cont_prop * 100, g[G] * 100, ls, marker="o", ms=ms, lw=lw, color=color, label=label)


def colored_legend(ax, **kw):
    leg = ax.legend(handlelength=0, handletextpad=0, **kw)
    for t, h in zip(leg.get_texts(), leg.legend_handles):
        t.set_color(h.get_color())
        h.set_visible(False)
    return leg


def predicted(ref, c, kappa, kappa_cc):
    """Acceptance predicted from kappa alone: a unit at contamination c with
    excess coincidences behaves like an independent one at c_eff, where
    c_eff (1 - c_eff/2) = c (1 - c)(1 + kappa) + (c^2 / 2)(1 + kappa_cc)."""
    target = c * (1 - c) * (1 + kappa) + c * c / 2 * (1 + kappa_cc)
    ce = 1 - np.sqrt(max(0.0, 1 - 2 * target))
    g = ref.sort_values("cont_prop")
    return 100 * np.interp(ce, g.cont_prop.values, g[G].values)


def main():
    plotstyle.apply()
    ref = pd.read_parquet(SIMS / "fig4_reference.pqt")
    d = pd.read_parquet(SIMS / "fig4_mismatch.pqt")
    for c in ("rho", "tau_s", "overlap", "width", "p_burst", "cont_rp"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15, 8.2))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.3)

    # --- a -------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    curve(ax, ref, "k", "Manuscript model", lw=2.2)
    curve(ax, d[(d.model == "single_neuron_contaminant") & np.isclose(d.cont_rp, 0.0015)],
          "#d95f02", "One-neuron contaminant, 1.5 ms RP")
    curve(ax, d[(d.model == "single_neuron_contaminant") & np.isclose(d.cont_rp, 0.0025)],
          "#e6ab02", "One-neuron contaminant, 2.5 ms RP")
    curve(ax, d[(d.model == "bursting") & np.isclose(d.p_burst, 0.3)], "#7570b3", "Bursting (30%)")
    curve(ax, d[(d.model == "graded") & np.isclose(d.width, 0.002)], "#1b9e77",
          "Graded recovery (2 ms ramp)")
    ax.axvline(10, color="0.6", ls="--", lw=1)
    ax.set_xlabel("True contamination (%)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("a  Departures in spike-train shape", loc="left")
    colored_legend(ax, fontsize=7.5, loc="upper right")

    # --- b -------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    c = d[d.model == "correlated"]
    rhos = sorted(c.rho.dropna().unique())
    cmap = plt.cm.RdBu_r
    for r in rhos:
        col = "0.45" if r == 0 else cmap(0.5 + 0.47 * np.sign(r) * abs(r) ** 0.5)
        curve(ax, c[np.isclose(c.rho, r)], col, f"\u03c1 = {r:+g}" if r else "\u03c1 = 0",
              lw=1.3, ms=2)
    curve(ax, ref, "k", "Manuscript model", lw=1.6, ls="--", ms=0)
    ax.axvline(10, color="0.6", ls="--", lw=1)
    ax.set_xlabel("True contamination (%)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("b  Correlated rates (\u03c1 = rate correlation)", loc="left")
    colored_legend(ax, fontsize=6.5, loc="lower left", ncol=2)

    # --- c -------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    m = d[d.model == "modulated"]
    for sub, mk, name in ((c, "o", "Correlated (true \u03c1)"), (m, "s", "Shared drive, scaled depth")):
        rows = []
        for r, g in sub.groupby("rho"):
            k = g.realised_kappa_mean.mean()
            kcc = g.realised_kappa_cc_mean.mean()
            rows.append((k, kcc, at(g, FAR_C), at(g, TAR_C)))
        rows.sort()
        k, kcc, far, tar = map(np.array, zip(*rows))
        ax.plot(k, far, mk, ms=5.5, color="#c2410c", mfc="#c2410c" if mk == "o" else "white",
                ls="none")
        ax.plot(k, tar, mk, ms=5.5, color="#1b9e77", mfc="#1b9e77" if mk == "o" else "white",
                ls="none")
        ax.plot([], [], mk, ms=5.5, color="0.3", mfc="0.3" if mk == "o" else "white",
                ls="none", label=name)
        if mk == "o":
            kk = np.linspace(k.min(), k.max(), 60)
            kc = np.interp(kk, k, kcc)
            ax.plot(kk, [predicted(ref, FAR_C, a, b) for a, b in zip(kk, kc)], "-", color="#c2410c",
                    lw=1, alpha=0.7)
            ax.plot([], [], "-", color="0.3", lw=1, label="Predicted from \u03ba alone")
            ax.plot(kk, [predicted(ref, TAR_C, a, b) for a, b in zip(kk, kc)], "-", color="#1b9e77",
                    lw=1, alpha=0.7)
    ax.axhline(at(ref, FAR_C), color="#c2410c", ls=":", lw=1)
    ax.axhline(at(ref, TAR_C), color="#1b9e77", ls=":", lw=1)
    ax.axvline(0, color="0.7", lw=1)
    ax.set_xlabel("Excess neuron-contaminant coincidence, \u03ba")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("c  What matters is \u03ba", loc="left")
    ax.legend(fontsize=7, loc="upper right")
    ax.text(-0.2, 40, f"False acceptance\n({FAR_C*100:g}% contamination)", color="#c2410c", fontsize=8)
    ax.text(0.2, 45, f"True acceptance\n({TAR_C*100:g}% contamination)", color="#1b9e77", fontsize=8)

    # --- d -------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    n = d[d.model == "nonoverlapping"]
    ovs = sorted(n.overlap.dropna().unique())
    ax.plot(ovs, [at(n[np.isclose(n.overlap, o)], FAR_C) for o in ovs], "o-", color="#c2410c",
            ms=5, label=f"False acceptance ({FAR_C*100:g}%)")
    ax.plot(ovs, [at(n[np.isclose(n.overlap, o)], TAR_C) for o in ovs], "s-", color="#1b9e77",
            ms=5, label=f"True acceptance ({TAR_C*100:g}%)")
    ax.axhline(at(ref, FAR_C), color="#c2410c", ls=":", lw=1)
    ax.axhline(at(ref, TAR_C), color="#1b9e77", ls=":", lw=1)
    ax.set_xlabel("Fraction of epochs shared")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("d  Non-overlapping activity", loc="left")
    colored_legend(ax, fontsize=7.5)

    # --- e -------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    g = d[d.model == "graded"]
    widths = sorted(g.width.dropna().unique())
    xs = [0] + [w * 1000 for w in widths]
    ax.plot(xs, [at(ref, FAR_C)] + [at(g[np.isclose(g.width, w)], FAR_C) for w in widths], "o-",
            color="#c2410c", ms=5, label=f"False acceptance ({FAR_C*100:g}%)")
    ax.plot(xs, [at(ref, TAR_C)] + [at(g[np.isclose(g.width, w)], TAR_C) for w in widths], "s-",
            color="#1b9e77", ms=5, label=f"True acceptance ({TAR_C*100:g}%)")
    ax.set_xlabel("Relative-refractory ramp width (ms)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_title("e  Graded recovery after a hard 3 ms period", loc="left")
    colored_legend(ax, fontsize=7.5)

    # --- f -------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    pk = pd.read_parquet(SIMS / "pair_kappa.pqt")
    clean = pk.ccg_ratio.between(0.5, 2.0) & pk.sep_um.between(50, 150)
    v = pk.loc[clean, "k_0p5_3"].dropna()
    v2 = pk.loc[clean, "k_30_300"].dropna()
    bins = np.linspace(-1.0, 2.5, 71)
    for vals, col, name in ((v, "k", "Lags 0.5-3 ms"), (v2, "0.6", "Lags 30-300 ms")):
        h, e = np.histogram(np.clip(vals, bins[0], bins[-1]), bins=bins)
        ax.stairs(h / h.sum(), e, color=col, lw=1.6, label=f"{name} (median {vals.median():+.2f})")
    ks = c.groupby("rho").realised_kappa_mean.mean()
    for r, k in ks.items():
        if r in (-1.0, -0.5, 0.5, 1.0):
            ax.axvline(k, color=cmap(0.5 + 0.47 * np.sign(r) * abs(r) ** 0.5), lw=1, ls="--")
            ax.text(k, ax.get_ylim()[1] * 0.98, f"\u03c1={r:+g}", rotation=90, fontsize=6.5,
                    va="top", ha="right", color="0.3")
    ax.set_xlabel("Excess coincidence, \u03ba")
    ax.set_ylabel("Proportion of pairs")
    ax.set_title(f"f  Real IBL neighbours (50-150 \u00b5m, n={clean.sum()})", loc="left")
    colored_legend(ax, fontsize=7.5, loc="center right")

    (OUTDIR / "figures").mkdir(exist_ok=True)
    plotstyle.save(fig, OUTDIR / "figures" / "mismatch_fig4")

    # --- numbers -------------------------------------------------------------
    L = ["Model mismatch at the Fig 4 setting (5 spikes/s, 3 ms, 1 h)", "=" * 58, "",
         f"False acceptance at {FAR_C*100:g}%, true acceptance at {TAR_C*100:g}%, "
         "Sliding RP at 10% / 90%.", "",
         f"Manuscript model: FA {at(ref, FAR_C):.1f}%  TA {at(ref, TAR_C):.1f}%  "
         f"(4000 trains per point)", ""]

    def row(name, sub):
        L.append(f"  {name:44s} FA {at(sub, FAR_C):5.1f}%   TA {at(sub, TAR_C):5.1f}%")
    for cr in (0.0015, 0.0025):
        row(f"single-neuron contaminant, RP {cr*1000:g} ms",
            d[(d.model == "single_neuron_contaminant") & np.isclose(d.cont_rp, cr)])
    for p in (0.1, 0.3):
        row(f"bursting p={p:g}", d[(d.model == "bursting") & np.isclose(d.p_burst, p)])
    for w in widths:
        row(f"graded ramp {w*1000:g} ms", g[np.isclose(g.width, w)])
    for o in ovs:
        row(f"non-overlapping, shared fraction {o:g}", n[np.isclose(n.overlap, o)])
    L += ["", "Correlated rates (true rho), with the excess coincidence kappa and the",
          "spike-count correlations the same pairs would show:", "",
          f"  {'rho':>6} {'kappa':>7} {'kappa_cc':>9} {'r 100ms':>8} {'r 1s':>7}   "
          f"{'FA':>6} {'TA':>6}   {'FA pred':>7} {'TA pred':>7}"]
    for r, gg in c.groupby("rho"):
        k, kcc = gg.realised_kappa_mean.mean(), gg.realised_kappa_cc_mean.mean()
        L.append(f"  {r:+6.2f} {k:+7.3f} {kcc:9.3f} {gg.realised_r_0p1_mean.mean():+8.3f} "
                 f"{gg.realised_r_1_mean.mean():+7.3f}   {at(gg, FAR_C):5.1f}% {at(gg, TAR_C):5.1f}%"
                 f"   {predicted(ref, FAR_C, k, kcc):6.1f}% {predicted(ref, TAR_C, k, kcc):6.1f}%")
    L += ["", "Shared-drive model (the first version's 'rho'):"]
    for r, gg in m.groupby("rho"):
        L.append(f"  {r:+6.2f} kappa {gg.realised_kappa_mean.mean():+7.3f}   "
                 f"FA {at(gg, FAR_C):5.1f}%  TA {at(gg, TAR_C):5.1f}%")
    L += ["", f"Real IBL neighbours, 50-150 um, CCG-clean, n = {clean.sum()}:",
          f"  kappa at 0.5-3 ms: median {v.median():+.3f}, IQR [{v.quantile(.25):+.3f}, "
          f"{v.quantile(.75):+.3f}], 5-95% [{v.quantile(.05):+.3f}, {v.quantile(.95):+.3f}]",
          f"  kappa at 30-300 ms: median {v2.median():+.3f}, IQR [{v2.quantile(.25):+.3f}, "
          f"{v2.quantile(.75):+.3f}]"]
    (OUTDIR / "mismatch_fig4_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
