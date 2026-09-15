"""Beryl-level supplement: recovery time and acceptance across brain regions.

Approved as a supplement to consider (root plan, Q10). The Cosmos level used in
Fig 1 pools each region very coarsely -- "thalamus" spans LGd, LP, VPM, PO and
the rest. With the full IBL brain-wide map there are enough units to look one
level finer, at the Beryl parcellation, which is where a reader would go to ask
whether a particular nucleus is short or long.

Three figures:

    beryl.pdf              the main supplement. Same box-and-whisker convention
                           as Fig 1 (median, bootstrapped CI, interquartile box,
                           5-95% whiskers) plus the acceptance rates.
    beryl_definitions.pdf  the same regions under all three candidate
                           timepoints, to show how much the picture depends on
                           which one you call the recovery time.
    beryl_fr.pdf           raw against firing-rate-standardized medians, which
                           is the quantitative version of "is this just a
                           firing-rate difference?".

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
from fr_standardize import (FR_EDGES, reference_weights, standardized_ci,  # noqa: E402
                            standardized_median)
from load_enriched import RP_COLUMNS, load_all, rule_common  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"01_fig1_rp_durations")
MIN_UNITS = 150
MIN_INSERTIONS = 5
MAIN_COL = "rp_ms_10"


def region_table(inc, df, weights):
    """One row per qualifying Beryl region, with everything the figures need."""
    rows = []
    for reg, g in inc.groupby("beryl"):
        if not reg or reg in ("root", "void") or len(g) < MIN_UNITS:
            continue
        if g.insertion_key.nunique() < MIN_INSERTIONS:
            continue
        v = g[MAIN_COL].astype(float)
        lo, hi = plotstyle.bootstrap_ci(v)
        std, cov = standardized_median(g.firing_rate, v, weights)
        slo, shi = standardized_ci(g.firing_rate, v, weights, n_boot=300)
        acc = df[(df.beryl == reg)
                 & (df.sorter_label.isna() | (df.sorter_label > 0))]
        row = dict(beryl=reg, cosmos=g.cosmos.mode().iat[0], n=len(g),
                   n_sess=g.session_key.nunique(),
                   n_ins=g.insertion_key.nunique(),
                   n_animals=g.animal_key.nunique(),
                   median=float(np.nanmedian(v)), lo=lo, hi=hi,
                   q25=float(np.nanpercentile(v, 25)),
                   q75=float(np.nanpercentile(v, 75)),
                   frac_below2=float(np.nanmean(v < 2)),
                   fr=float(np.nanmedian(g.firing_rate)),
                   std_median=std, std_lo=slo, std_hi=shi, fr_coverage=cov,
                   n_acc=len(acc),
                   sliding=float(acc.passes.mean()) if len(acc) else np.nan,
                   hl3=float(acc.hl3_pass.mean()) if len(acc) else np.nan)
        for col in RP_COLUMNS:
            row[f"med_{col}"] = float(np.nanmedian(g[col].astype(float)))
        rows.append(row)
    t = pd.DataFrame(rows).sort_values("median").reset_index(drop=True)
    t["rescued"] = t.sliding - t.hl3
    return t


def row_labels(t):
    return [f"{r.beryl}  ({r.cosmos}, n={r.n:,}, {r.n_ins} ins.)"
            for r in t.itertuples()]


def fig_main(inc, t):
    """Panel a, the Fig 1 box-and-whisker convention; panel b, acceptance."""
    fig, axs = plt.subplots(1, 2, figsize=(12, max(4.5, 0.235 * len(t))),
                            gridspec_kw={"width_ratios": [1.3, 1]})
    y = np.arange(len(t))[::-1]
    cols = [plotstyle.REGION_COLORS.get(c, "0.45") for c in t.cosmos]

    ax = axs[0]
    for yy, r, c in zip(y, t.itertuples(), cols):
        v = inc.loc[inc.beryl == r.beryl, MAIN_COL].astype(float)
        plotstyle.box_row(ax, v, yy, c, height=0.55)
    ax.axvline(2, color="0.7", ls=":", lw=1, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(row_labels(t), fontsize=6)
    ax.set_ylim(-1, len(t))
    ax.set_xlim(0.8, 8)
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_title("a  Median (95% CI), interquartile box, 5-95% whiskers",
                 loc="left")

    ax = axs[1]
    ax.barh(y + 0.2, t.sliding, 0.4, color="#1b9e77", label="Sliding RP")
    ax.barh(y - 0.2, t.hl3, 0.4, color="#7570b3", label="Hill-Llobet, 3 ms")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_ylim(-1, len(t))
    ax.set_xlabel("Proportion of units accepted")
    ax.set_title("b  Acceptance rate", loc="left")
    ax.legend(fontsize=7, loc="lower right")

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "beryl")


def fig_definitions(inc, t):
    """The same regions under all three candidate timepoints."""
    cols = list(RP_COLUMNS)
    fig, axs = plt.subplots(1, len(cols), figsize=(13, max(4.5, 0.235 * len(t))),
                            sharey=True)
    y = np.arange(len(t))[::-1]
    colors = [plotstyle.REGION_COLORS.get(c, "0.45") for c in t.cosmos]
    for k, (ax, col) in enumerate(zip(axs, cols)):
        for yy, r, c in zip(y, t.itertuples(), colors):
            v = inc.loc[inc.beryl == r.beryl, col].astype(float)
            plotstyle.box_row(ax, v, yy, c, height=0.55)
        ax.axvline(2, color="0.7", ls=":", lw=1, zorder=0)
        ax.set_xlabel(RP_COLUMNS[col].axis)
        ax.set_title(f"{'abc'[k]}  {RP_COLUMNS[col].math}", loc="left")
        # each definition lives on its own scale, so give each panel its own
        # limits; the rows are in a common order, which is what is being read
        allv = inc[col].astype(float)
        ax.set_xlim(max(0.3, np.nanpercentile(allv, 0.5)),
                    min(10.2, np.nanpercentile(allv, 99.5)))
        ax.set_ylim(-1, len(t))
    axs[0].set_yticks(y)
    axs[0].set_yticklabels(row_labels(t), fontsize=6)
    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "beryl_definitions")


def fig_firing_rate(t):
    """Raw against firing-rate-standardized medians."""
    fig, axs = plt.subplots(1, 2, figsize=(10, max(4.5, 0.235 * len(t))),
                            gridspec_kw={"width_ratios": [1.35, 1]})
    y = np.arange(len(t))[::-1]
    colors = [plotstyle.REGION_COLORS.get(c, "0.45") for c in t.cosmos]

    ax = axs[0]
    for yy, r, c in zip(y, t.itertuples(), colors):
        ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=2.2, alpha=0.5,
                solid_capstyle="butt")
        ax.plot(r.median, yy, "o", color="w", mec=c, mew=1.2, ms=5)
        if np.isfinite(r.std_median):
            ax.plot(r.std_median, yy, "|", color=c, ms=9, mew=1.8)
            ax.plot([r.median, r.std_median], [yy, yy], color=c, lw=0.7,
                    alpha=0.6)
    ax.axvline(2, color="0.7", ls=":", lw=1, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.beryl}  ({r.cosmos}, {r.fr:.0f} spikes/s"
         + ("" if r.fr_coverage > 0.8 else ", partial") + ")"
         for r in t.itertuples()], fontsize=6)
    ax.set_ylim(-1, len(t))
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_title("a  Circle, observed median; tick, firing-rate standardized",
                 loc="left")

    ax = axs[1]
    for r, c in zip(t.itertuples(), colors):
        if not np.isfinite(r.std_median):
            continue
        full = r.fr_coverage > 0.8
        ax.plot(r.fr, r.std_median - r.median, "o", ms=4.5,
                color=c if full else "w", mec=c, mew=1.1)
    ax.axhline(0, color="0.7", lw=1)
    ax.set_xscale("log")
    plotstyle.plain_log_ticks(ax)
    ax.set_xlabel("Median firing rate of the region (spikes/s)")
    ax.set_ylabel("Standardized minus observed median (ms)")
    ax.set_title("b  Size of the adjustment", loc="left")
    ax.text(0.02, 0.04, "Open: standardized over part of the\nreference "
            "firing-rate range only", transform=ax.transAxes, fontsize=6,
            va="bottom", color="0.35")
    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "beryl_fr")


def write_numbers(t, weights):
    L = ["Beryl-level recovery times and acceptance rates", "=" * 52, "",
         f"Regions with at least {MIN_UNITS} units passing the common inclusion",
         f"rule and at least {MIN_INSERTIONS} insertions. {len(t)} regions",
         "qualify (the same 80 either way: no region has 5 insertions but",
         "fewer than 5 sessions, or the reverse).",
         "",
         "Columns: n units / n sessions / n insertions / n animals; median",
         "estimated ACG recovery time with bootstrapped 95% CI; the fraction",
         "below 2 ms; median firing rate; the firing-rate-standardized median",
         "and the share of the reference firing-rate range it covers; and the",
         "Sliding RP and Hill-Llobet (3 ms) acceptance rates.", ""]
    L.append(f"{'region':10s} {'cosmos':10s} {'n':>7} {'sess':>5} {'ins':>5} "
             f"{'anim':>5} {'median':>8} {'[95% CI]':>18} {'<2ms':>7} "
             f"{'FR':>7} {'std':>7} {'cov':>5} {'slidRP':>8} {'HL3ms':>7} "
             f"{'rescued':>8}")
    for r in t.itertuples():
        L.append(f"{r.beryl:10s} {r.cosmos:10s} {r.n:7,} {r.n_sess:5d} "
                 f"{r.n_ins:5d} {r.n_animals:5d} {r.median:8.3f} "
                 f"{f'[{r.lo:.3f}, {r.hi:.3f}]':>18} {r.frac_below2:7.1%} "
                 f"{r.fr:7.1f} {r.std_median:7.3f} {r.fr_coverage:5.2f} "
                 f"{r.sliding:8.3f} {r.hl3:7.3f} {r.rescued:8.3f}")

    sh, lo_ = t.nsmallest(5, "median"), t.nlargest(5, "median")
    L += ["",
          "Shortest: " + ", ".join(f"{r.beryl} ({r.median:.2f} ms)"
                                   for r in sh.itertuples()),
          "Longest : " + ", ".join(f"{r.beryl} ({r.median:.2f} ms)"
                                   for r in lo_.itertuples()),
          "",
          "Most rescued by Sliding RP relative to a fixed 3 ms threshold:",
          "  " + ", ".join(f"{r.beryl} (+{r.rescued:.2f})"
                           for r in t.nlargest(6, "rescued").itertuples()),
          "", "Spread within each Cosmos parent (median of region medians):"]
    for c, g in t.groupby("cosmos"):
        L.append(f"  {c:10s} {len(g):3d} regions  "
                 f"{g['median'].median():.3f} ms  "
                 f"range {g['median'].min():.3f}-{g['median'].max():.3f}  "
                 f"mean rescue {g.rescued.mean():+.3f}")

    d = t.std_median - t["median"]
    L += ["", "Firing-rate standardization", "-" * 27,
          f"Region median firing rates span {t.fr.min():.1f} to {t.fr.max():.1f} "
          f"spikes/s, so unlike the Cosmos level the regions are not matched on",
          "rate and the adjustment is not negligible.",
          f"  mean |adjustment| {d.abs().mean():.3f} ms, largest "
          f"{d.abs().max():.3f} ms",
          f"  rank correlation of raw and standardized medians: "
          f"{t['median'].corr(t.std_median, method='spearman'):.3f}",
          f"  {(t.fr_coverage <= 0.8).sum()} regions cover 80% or less of the "
          f"reference firing-rate range; read those as applying to the part of",
          "  the range they do cover.", "",
          "  largest adjustments:"]
    for r in t.reindex(d.abs().sort_values(ascending=False).index).head(8).itertuples():
        L.append(f"    {r.beryl:8s} {r.cosmos:10s} FR {r.fr:5.1f}  "
                 f"raw {r.median:.3f} -> std {r.std_median:.3f} "
                 f"({r.std_median - r.median:+.3f}, coverage {r.fr_coverage:.2f})")

    L += ["", "Median under each candidate timepoint (ms)", "-" * 42,
          "How much does the region ordering depend on which timepoint you",
          "call the recovery time? Spearman correlation of the 80 region",
          "medians, between definitions:"]
    cols = list(RP_COLUMNS)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = t[f"med_{cols[i]}"], t[f"med_{cols[j]}"]
            L.append(f"  {RP_COLUMNS[cols[i]].text:28s} vs "
                     f"{RP_COLUMNS[cols[j]].text:28s} "
                     f"rho = {a.corr(b, method='spearman'):+.3f}")
    L += ["",
          f"{'region':10s} " + " ".join(f"{RP_COLUMNS[c].text[:22]:>24s}"
                                        for c in RP_COLUMNS)]
    for r in t.itertuples():
        L.append(f"{r.beryl:10s} " + " ".join(
            f"{getattr(r, 'med_' + c):24.3f}" for c in RP_COLUMNS))

    L += ["", "Read with care: these are apparent quiet-window durations under",
          "three operational definitions, not biophysical refractory periods.",
          "The curve-fitting estimator behind the first column is imperfect,",
          "which is one reason the other two are shown beside it."]
    (OUTDIR / "beryl_numbers.txt").write_text("\n".join(L))
    return L


def main():
    plotstyle.apply()
    df = load_all(verbose=False)
    if df.empty:
        print("no enriched tables")
        return
    inc = df[rule_common(df, require_pass=True)].copy()
    weights = reference_weights(inc.firing_rate)
    t = region_table(inc, df, weights)
    t.to_csv(OUTDIR / "beryl_table.csv", index=False)

    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)
    fig_main(inc, t)
    fig_definitions(inc, t)
    fig_firing_rate(t)
    L = write_numbers(t, weights)
    print("\n".join(L[:20]))
    print(f"... full table in {OUTDIR / 'beryl_numbers.txt'}")


if __name__ == "__main__":
    main()
