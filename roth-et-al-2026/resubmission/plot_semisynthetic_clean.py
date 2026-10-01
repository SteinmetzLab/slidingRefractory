"""Figure and numbers for the cleaned-recipient semi-synthetic validation (05).

Reads semisynthetic_clean.pqt (semisynthetic_clean.py). A pair's acceptance at
an injection level is the mean over its random draws; a curve is the mean over
pairs; 95% intervals resample recipients (each with all its pairs), since pairs
sharing a recipient are not independent.

Panels
------
a  Acceptance against injected contamination, Sliding RP and Hill-Llobet at
   2 ms, main-test pairs (solid) and dip pairs (dashed). Both start at 100%:
   recipients are cleaned to 2 ms.
b  Sliding RP's minimum confirmable contamination against the injected fraction.
c  Distribution of the excess short-lag coincidence kappa(0-2 ms) between
   recipient and donor, by donor distance.
d  Median kappa profile at 0.1 ms resolution, by donor distance.
e, f  Acceptance curves by kappa(0-2 ms) bin, for Sliding RP and Hill-Llobet.
g  Sliding RP split by the 100 ms spike-count correlation, for comparison with
   the earlier version of this analysis.
h  Why some units are still accepted at high contamination: each pair's
   kappa at 0-0.5 ms against 0.5-2 ms, colored by Sliding RP's acceptance at
   20% injected; dip pairs lie above the dashed line.

Run:  python plot_semisynthetic_clean.py
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
GREEN, ORANGE = "#1b9e77", "#d95f02"
KAPPA_EDGES = [-np.inf, -0.3, -0.1, 0.1, 0.3, 1.0, np.inf]
KAPPA_LABELS = ["< -0.3", "-0.3 to -0.1", "-0.1 to 0.1", "0.1 to 0.3", "0.3 to 1", "> 1"]
R_EDGES = [-np.inf, -0.02, 0.02, 0.08, 0.2, np.inf]
R_LABELS = ["< -0.02", "-0.02 to 0.02", "0.02 to 0.08", "0.08 to 0.2", "> 0.2"]
PAIR = ["pid", "recipient", "donor"]
# a pair has a short-lag dip when its coincidence at 0-0.5 ms is under half that
# at 0.5-2 ms (most likely overlapping waveforms the sorter did not resolve; see
# plot_dip_pairs.py). Dip pairs are kept out of the main test and shown apart.
DIP_RATIO = 0.5


def bin_colors(n):
    """Diverging blue (negative) to red (positive) with a dark gray center, so
    no bin is too pale to read."""
    six = ["#08519c", "#4292c6", "0.35", "#fb6a4a", "#cb181d", "#67000d"]
    five = ["#08519c", "#4292c6", "0.35", "#fb6a4a", "#a50f15"]
    return {6: six, 5: five}[n]


def pair_table(d):
    """One row per pair and injection level: acceptance averaged over draws."""
    keys = PAIR + ["f_injected"]
    agg = d.groupby(keys).agg(passes=("passes", "mean"), hl2_pass=("hl2_pass", "mean"),
                              min_cont=("min_cont", "median")).reset_index()
    first = d.drop_duplicates(PAIR).set_index(PAIR)
    cols = ["donor_kind", "depth_sep_um", "fr_recipient", "fr_donor", "r_100ms", "r_1s",
            "k_0_0p5", "k_0p5_1", "k_1_2", "k_0_2", "k_0p5_2", "k_3_30", "k_30_300"]
    out = agg.join(first[cols], on=PAIR)
    out["dip"] = (1 + out.k_0_0p5) / (1 + out.k_0p5_2) < DIP_RATIO
    return out


def curve(p, col, seed=0, n_boot=300):
    """Mean acceptance by injection level, with a recipient-level bootstrap CI."""
    f = np.sort(p.f_injected.unique())
    tab = p.pivot_table(index=["pid", "recipient", "donor"], columns="f_injected", values=col)
    tab = tab[f]
    # recipient = (insertion, cluster id); cluster ids repeat across insertions
    rec_code = pd.factorize(pd.MultiIndex.from_arrays(
        [tab.index.get_level_values("pid"), tab.index.get_level_values("recipient")]))[0]
    vals = tab.to_numpy()
    m = np.nanmean(vals, axis=0)
    # resampling recipients with all their pairs = weighting each recipient's
    # pair sums by its multinomial draw count
    n_rec = rec_code.max() + 1
    ok = ~np.isnan(vals)
    sums = np.stack([np.bincount(rec_code, np.where(ok[:, k], vals[:, k], 0), n_rec)
                     for k in range(vals.shape[1])], 1)
    cnts = np.stack([np.bincount(rec_code, ok[:, k], n_rec) for k in range(vals.shape[1])], 1)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        w = np.bincount(rng.integers(n_rec, size=n_rec), minlength=n_rec)
        boots.append((w @ sums) / (w @ cnts))
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    return f * 100, m * 100, lo * 100, hi * 100, len(tab)


def color_text_legend(ax, **kw):
    leg = ax.legend(handlelength=0, handletextpad=0, frameon=False, **kw)
    for text, h in zip(leg.get_texts(), leg.legend_handles):
        text.set_color(h.get_color() if hasattr(h, "get_color") else h.get_facecolor()[0])
        h.set_visible(False)
    return leg


def split_panel(ax, p, col, by, edges, labels, title):
    q = p.copy()
    q["bin"] = pd.cut(q[by], edges, labels=labels)
    cols = bin_colors(len(labels))
    out = []
    for lab, c in zip(labels, cols):
        g = q[q.bin == lab]
        if g.empty:
            continue
        x, y, lo, hi, n = curve(g, col)
        ax.fill_between(x, lo, hi, color=c, alpha=0.15, lw=0)
        ax.plot(x, y, color=c, lw=1.8, label=f"{lab} (n = {n})")
        out.append((lab, n, x, y))
    ax.axvline(10, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_ylim(-2, 102)
    ax.set_title(title, loc="left")
    return out


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "semisynthetic_clean.pqt")
    rc = pd.read_parquet(SIMS / "semisynthetic_clean_recipients.pqt")
    p_all = pair_table(d)
    p, p_dip = p_all[~p_all.dip], p_all[p_all.dip]
    pairs = p.drop_duplicates(PAIR)
    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14, 7.4))
    gs = fig.add_gridspec(2, 4)
    ax = fig.add_subplot(gs[0, 0])

    # --- a: overall curves ---------------------------------------------------
    n_dip = p_dip.drop_duplicates(PAIR).shape[0]
    for col, lab, c in (("passes", "Sliding RP", GREEN), ("hl2_pass", "Hill-Llobet, 2 ms", ORANGE)):
        x, y, lo, hi, n = curve(p, col)
        ax.fill_between(x, lo, hi, color=c, alpha=0.2, lw=0)
        ax.plot(x, y, color=c, lw=2, label=f"{lab}")
        if n_dip:
            x, y, lo, hi, _ = curve(p_dip, col)
            ax.plot(x, y, color=c, lw=1.4, ls=(0, (4, 2)), label=f"{lab}, dip pairs")
    ax.axvline(10, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Units accepted (%)")
    ax.set_ylim(-2, 102)
    ax.set_title(f"a  {n} pairs (dashed: {n_dip} dip pairs)", loc="left")
    ax.legend(loc="upper right", fontsize=7, handlelength=2.6, frameon=False)

    # --- b: estimate ------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    f = np.sort(p.f_injected.unique())
    med = [p[p.f_injected == v].min_cont.median() for v in f]
    q1 = [p[p.f_injected == v].min_cont.quantile(.25) for v in f]
    q3 = [p[p.f_injected == v].min_cont.quantile(.75) for v in f]
    ax.fill_between(f * 100, q1, q3, color=GREEN, alpha=0.25, lw=0)
    ax.plot(f * 100, med, color=GREEN, lw=2, label="Minimum confirmable contamination")
    ax.plot([0, 25], [0, 25], "k--", lw=1, label="Injected")
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Contamination (%)")
    ax.set_title("b  Sliding RP's estimate", loc="left")
    ax.legend(loc="upper left", fontsize=7.5, handlelength=2.2)

    # --- c: kappa distribution -------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    bins = np.linspace(-1, 2, 61)
    for k, c in (("near", ORANGE), ("mid", "0.45"), ("far", "#1f78b4")):
        v = pairs[pairs.donor_kind == k].k_0_2.dropna().clip(-1, 2)
        if v.empty:
            continue
        n_, e = np.histogram(v, bins=bins)
        ax.stairs(n_ / n_.sum(), e, color=c, lw=1.6,
                  label=f"{k}: median {v.median():+.2f} (n = {len(v)})")
    ax.axvline(0, color="0.6", ls=":", lw=1)
    ax.set_xlabel(r"Excess coincidence $\kappa$, lags 0-2 ms (clipped at 2)")
    ax.set_ylabel("Proportion of pairs")
    ax.set_title("c  How correlated at short lags?", loc="left")
    color_text_legend(ax, loc="upper right", fontsize=7.5)

    # --- g: diagnosis of the residual acceptance ---------------------------------
    ax = fig.add_subplot(gs[1, 3])
    top = p_all[np.isclose(p_all.f_injected, 0.20)]
    xx = np.linspace(-1, 1, 50)
    ax.plot(xx, (1 + xx) / DIP_RATIO - 1, color="0.3", lw=0.8, ls="--")
    sc = ax.scatter(top.k_0_0p5.clip(-1, 3), top.k_0p5_2.clip(-1, 3), c=top.passes * 100,
                    cmap="viridis", s=4, vmin=0, vmax=100, lw=0, rasterized=True)
    ax.axhline(0, color="0.6", ls=":", lw=1)
    ax.axvline(0, color="0.6", ls=":", lw=1)
    ax.set_xlabel(r"$\kappa$, lags 0-0.5 ms")
    ax.set_ylabel(r"$\kappa$, lags 0.5-2 ms")
    ax.set_xlim(-1.05, 3.05)
    ax.set_ylim(-1.05, 3.05)
    ax.set_title("h  Who is accepted at 20% injected?", loc="left")
    ax.text(-0.95, 2.9, "dip pairs\n(above dashed line)", fontsize=7, va="top", color="0.3")
    cb = fig.colorbar(sc, ax=ax, fraction=0.05)
    cb.set_label("Sliding RP acceptance (%)")

    # --- d, e, f: split curves ---------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    sk = split_panel(ax, p, "passes", "k_0_2", KAPPA_EDGES, KAPPA_LABELS,
                     r"e  Sliding RP, by $\kappa$(0-2 ms)")
    color_text_legend(ax, loc="upper right", fontsize=7)
    ax = fig.add_subplot(gs[1, 1])
    hk = split_panel(ax, p, "hl2_pass", "k_0_2", KAPPA_EDGES, KAPPA_LABELS,
                     r"f  Hill-Llobet 2 ms, by $\kappa$(0-2 ms)")
    color_text_legend(ax, loc="upper right", fontsize=7)
    ax = fig.add_subplot(gs[1, 2])
    sr = split_panel(ax, p, "passes", "r_100ms", R_EDGES, R_LABELS,
                     "g  Sliding RP, by count correlation (100 ms)")
    color_text_legend(ax, loc="upper right", fontsize=7)

    # --- h: mean short-lag coincidence profile, by donor distance -----------------
    ax = fig.add_subplot(gs[0, 3])
    prof = d.drop_duplicates(PAIR).merge(p_all.drop_duplicates(PAIR)[PAIR + ["dip"]], on=PAIR)
    prof = prof[~prof.dip]
    lags = (np.arange(30) + 0.5) * 0.1
    for kind, c in (("near", ORANGE), ("mid", "0.45"), ("far", "#1f78b4")):
        g = prof[prof.donor_kind == kind]
        if g.empty:
            continue
        m = np.stack(g.kappa_profile.to_numpy())
        ax.plot(lags, np.nanmedian(m, axis=0), color=c, lw=1.8, label=f"{kind} (median)")
        ax.fill_between(lags, *np.nanpercentile(m, [25, 75], axis=0), color=c, alpha=0.15, lw=0)
    ax.axhline(0, color="0.6", ls=":", lw=1)
    ax.axvline(0.5, color="0.6", ls="--", lw=0.8)
    ax.text(0.52, ax.get_ylim()[1], r"$\tau_{\min}$", va="top", fontsize=8, color="0.4")
    ax.set_xlabel("Lag between recipient and donor spikes (ms)")
    ax.set_ylabel(r"Excess coincidence $\kappa$")
    ax.set_title("d  Short-lag coincidence profile (no dip pairs)", loc="left")
    color_text_legend(ax, loc="lower right", fontsize=7.5)

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "semisynthetic_clean")

    # --- numbers --------------------------------------------------------------------
    used = rc[rc.frac_deleted <= 0.10]
    L = ["Semi-synthetic contamination, cleaned recipients (semisynthetic_clean.py)", "=" * 72, "",
         f"{d.pid.nunique()} IBL insertions; {used.shape[0]} recipients (of {rc.shape[0]} firing "
         f">= 2 spikes/s; {(rc.frac_deleted > 0.10).sum()} excluded for losing > 10% of spikes);",
         f"{pairs.groupby(['pid', 'recipient']).ngroups} recipients in the main test, "
         f"{len(pairs)} recipient-donor pairs in the main test, after setting aside "
         f"{p_dip.drop_duplicates(PAIR).shape[0]} dip pairs "
         f"({', '.join(f'{k} {v}' for k, v in p_dip.drop_duplicates(PAIR).donor_kind.value_counts().items())}); main: ({', '.join(f'{k} {v}' for k, v in pairs.donor_kind.value_counts().items())});",
         f"{d.rep.nunique()} random draws per pair; {len(d):,} trains evaluated.",
         f"Spikes deleted to clean recipients to 2 ms: median {used.frac_deleted.median():.1%}, "
         f"95th percentile {used.frac_deleted.quantile(.95):.1%}.",
         f"Recipient firing rate: median {pairs.fr_recipient.median():.1f} spikes/s "
         f"(IQR {pairs.fr_recipient.quantile(.25):.1f}-{pairs.fr_recipient.quantile(.75):.1f}).", "",
         "Acceptance against injected contamination (mean over pairs, 95% CI over recipients):", "",
         f"{'injected':>9} {'Sliding RP':>22} {'HL 2 ms':>22} {'C_min median':>13}"]
    xs, ys, los, his, _ = curve(p, "passes")
    xh, yh, loh, hih, _ = curve(p, "hl2_pass")
    for i, v in enumerate(xs):
        L.append(f"{v:8.0f}% {ys[i]:7.1f}% [{los[i]:5.1f}, {his[i]:5.1f}] "
                 f"{yh[i]:7.1f}% [{loh[i]:5.1f}, {hih[i]:5.1f}] "
                 f"{p[np.isclose(p.f_injected, v / 100)].min_cont.median():12.2f}%")
    if len(p_dip):
        xd, yd, _, _, nd = curve(p_dip, "passes")
        xe, ye, _, _, _ = curve(p_dip, "hl2_pass")
        L += ["", f"Dip pairs (n = {nd}), set aside: acceptance at 0 / 5 / 10 / 15 / 20 / 25% injected", "",
              "  Sliding RP  " + " ".join(f"{yd[np.argmin(abs(xd - v))]:5.1f}" for v in (0, 5, 10, 16, 20, 25)),
              "  HL 2 ms     " + " ".join(f"{ye[np.argmin(abs(xe - v))]:5.1f}" for v in (0, 5, 10, 16, 20, 25))]
    L += ["", "Excess short-lag coincidence kappa in the main test, by donor distance (median [IQR]):", ""]
    for k, g in pairs.groupby("donor_kind"):
        L.append(f"  {k:5s} n={len(g):5d}  " + "  ".join(
            f"{w} {g[w].median():+.2f} [{g[w].quantile(.25):+.2f}, {g[w].quantile(.75):+.2f}]"
            for w in ("k_0_0p5", "k_0p5_2", "k_0_2", "k_3_30")))
    for name, rows in (("Sliding RP by kappa(0-2 ms)", sk), ("Hill-Llobet 2 ms by kappa(0-2 ms)", hk),
                       ("Sliding RP by 100 ms count correlation", sr)):
        L += ["", f"{name}: acceptance at 0 / 5 / 10 / 15 / 20 / 25% injected", ""]
        for lab, n, x, y in rows:
            at = [y[np.argmin(abs(x - v))] for v in (0, 5, 10, 16, 20, 25)]
            L.append(f"  {lab:>14s} n={n:5d}  " + " ".join(f"{a:5.1f}" for a in at))

    # residual acceptance at 20%: what explains it
    top = p_all[np.isclose(p_all.f_injected, 0.20)].copy()
    top["collision"] = top.dip
    top["slow_anti"] = top.k_3_30 < -0.3
    L += ["", "Residual Sliding RP acceptance at 20% injected, all pairs, by pair type:", ""]
    for lab, m in (("short-lag dip (dip pairs)", top.collision),
                   ("slow anticorrelation (kappa 3-30 ms < -0.3), no dip",
                    ~top.collision & top.slow_anti),
                   ("neither", ~top.collision & ~top.slow_anti)):
        g = top[m]
        L.append(f"  {lab:62s} {len(g) / len(top):6.1%} of pairs, accepted {g.passes.mean():6.1%}, "
                 f"share of all acceptances {g.passes.sum() / max(top.passes.sum(), 1e-12):6.1%}")
    (OUTDIR / "semisynthetic_clean_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
