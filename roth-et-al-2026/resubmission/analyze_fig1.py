"""Work package 01: the Fig 1 re-analysis.

One estimator applied identically to every dataset, one common inclusion rule
alongside the two original ones, hierarchical statistics that respect the
nesting of units in insertions in animals in datasets, and the four diagnostics
Reviewer 2 asked for.

Outputs (into 01_fig1_rp_durations/):
    figures/fig1.pdf            the redrawn Fig 1 (estimate + CI + distribution)
    figures/fig1_diagnostics.pdf  floor rate, RP vs firing rate, sensitivity
    fig1_numbers.txt            every number quoted in the text
    fig1_units.pqt              the per-unit table the figure is built from

Run:  python analyze_fig1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402
from fr_standardize import (FR_EDGES, by_bin_table, reference_weights,  # noqa: E402
                            standardized_ci, standardized_median)
from load_enriched import (RP_COLUMNS, load_all, rule_common,  # noqa: E402
                           rule_mouse)

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"01_fig1_rp_durations")
REGIONS = ["Isocortex", "HPF", "TH"]
OUT_TABLE = Path(r"D:/temp/slidingRP_resub/fig1_units.pqt")


def bootstrap_ci(x, fn=np.median, n=2000, seed=0):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    bs = fn(x[rng.integers(0, x.size, (n, x.size))], axis=1)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def group_stats(df, by=("dataset", "cosmos")):
    rows = []
    for key, g in df.groupby(list(by)):
        lo, hi = bootstrap_ci(g.rp_ms_10)
        rows.append(dict(zip(by, key if isinstance(key, tuple) else (key,))) | dict(
            n_units=len(g), n_sessions=g.session_key.nunique(),
            n_insertions=g.insertion_key.nunique(),
            n_animals=g.animal_key.nunique(),
            median=float(np.nanmedian(g.rp_ms_10)), ci_lo=lo, ci_hi=hi,
            q25=float(np.nanpercentile(g.rp_ms_10, 25)),
            q75=float(np.nanpercentile(g.rp_ms_10, 75)),
            frac_below_2ms=float(np.nanmean(g.rp_ms_10 < 2)),
            frac_floored=float(np.nanmean(g.rp_floor)),
            median_fr=float(np.nanmedian(g.firing_rate)))
        )
    return pd.DataFrame(rows)


def hierarchical(df, lines):
    """Mixed-effects models that respect units-in-insertions-in-animals."""
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        lines.append("statsmodels not available; hierarchical models skipped")
        return

    d = df[df.cosmos.isin(REGIONS)].copy()
    d["log_rp"] = np.log(d.rp_ms_10)
    d["log_fr"] = np.log(d.firing_rate)

    # --- mouse only: region effect, insertions nested in animals ----------
    m = d[(d.species == "mouse")].copy()
    if len(m) > 100:
        lines += ["", "Mouse, region effect (units nested in insertions in animals)",
                  "  log(RP) ~ C(cosmos) + (1 | animal_key/insertion_key)"]
        try:
            fit = smf.mixedlm("log_rp ~ C(cosmos, Treatment('Isocortex'))", m,
                              groups=m["animal_key"],
                              re_formula="1", vc_formula={"ins": "0 + C(insertion_key)"}
                              ).fit(method="lbfgs", maxiter=200)
            lines.append("  " + str(fit.summary().tables[1]).replace("\n", "\n  "))
        except Exception as e:  # noqa: BLE001
            lines.append(f"  model failed: {type(e).__name__}: {e}")
            try:
                fit = smf.mixedlm("log_rp ~ C(cosmos, Treatment('Isocortex'))", m,
                                  groups=m["animal_key"]).fit()
                lines.append("  (animal-level random intercept only)")
                lines.append("  " + str(fit.summary().tables[1]).replace("\n", "\n  "))
            except Exception as e2:  # noqa: BLE001
                lines.append(f"  fallback failed too: {e2}")

        lines += ["", "  same, with log firing rate as a covariate:"]
        try:
            fit2 = smf.mixedlm("log_rp ~ C(cosmos, Treatment('Isocortex')) + log_fr",
                               m, groups=m["animal_key"]).fit()
            lines.append("  " + str(fit2.summary().tables[1]).replace("\n", "\n  "))
        except Exception as e:  # noqa: BLE001
            lines.append(f"  failed: {e}")

    # --- species comparison on the shared regions -------------------------
    s = d[d.cosmos.isin(["Isocortex", "TH"])].copy()
    if s.species.nunique() > 1:
        lines += ["", "Species comparison (Isocortex and TH only), animal random effect",
                  "  log(RP) ~ C(cosmos) * C(species) + (1 | animal_key)",
                  f"  n = {len(s):,} units, {s.insertion_key.nunique()} insertions, "
                  f"{s.animal_key.nunique()} animals "
                  f"({s[s.species=='macaque'].animal_key.nunique()} macaque)"]
        try:
            fit = smf.mixedlm("log_rp ~ C(cosmos) * C(species)", s,
                              groups=s["animal_key"]).fit()
            lines.append("  " + str(fit.summary().tables[1]).replace("\n", "\n  "))
        except Exception as e:  # noqa: BLE001
            lines.append(f"  failed: {type(e).__name__}: {e}")

    # --- the same comparison at the animal level (the reviewer's alternative)
    lines += ["", "Animal-level medians (one value per animal), Isocortex vs TH:"]
    per_animal = (d[d.cosmos.isin(["Isocortex", "TH"])]
                  .groupby(["species", "dataset", "animal_key", "cosmos"])
                  .rp_ms_10.median().reset_index())
    for sp in sorted(per_animal.species.unique()):
        for reg in ["Isocortex", "TH"]:
            v = per_animal[(per_animal.species == sp) & (per_animal.cosmos == reg)].rp_ms_10
            if len(v):
                lines.append(f"  {sp:8s} {reg:10s} n_animals={len(v):3d} "
                             f"median of animal medians {np.median(v):.3f} ms "
                             f"[{np.percentile(v,25):.3f}, {np.percentile(v,75):.3f}]")
    # paired within-insertion contrast: the cleanest region comparison
    lines += ["", "Within-insertion paired contrast (insertions sampling both "
              "Isocortex and TH):"]
    piv = (d[d.cosmos.isin(["Isocortex", "TH"]) & (d.species == "mouse")]
           .groupby(["insertion_key", "cosmos"]).rp_ms_10.median().unstack())
    piv = piv.dropna()
    if len(piv) > 5:
        diff = piv["TH"] - piv["Isocortex"]
        from scipy import stats as sps
        w = sps.wilcoxon(diff)
        lines.append(f"  n = {len(piv)} insertions; median(TH - Isocortex) = "
                     f"{np.median(diff):+.3f} ms; Wilcoxon p = {w.pvalue:.3g}")
    else:
        lines.append(f"  only {len(piv)} insertions sample both regions")


def make_figure(df, stats, path):
    """Fig 1c redrawn: central estimate with CI *and* the distribution."""
    plotstyle.apply()
    groups = [(ds, reg) for ds in ["steinmetz", "ibl", "allen", "macaque"]
              for reg in REGIONS
              if len(df[(df.dataset == ds) & (df.cosmos == reg)]) >= 20]
    fig, axs = plt.subplots(1, 2, figsize=(10, max(4, 0.34 * len(groups))),
                            gridspec_kw={"width_ratios": [1, 1.15]})

    # left: IBL histograms by region (as the submitted Fig 1b)
    ax = axs[0]
    for reg in REGIONS:
        v = df[(df.dataset == "ibl") & (df.cosmos == reg)].rp_ms_10.dropna()
        if len(v) < 20:
            continue
        n, edges = np.histogram(v, bins=np.arange(0, 5.01, 0.2))
        ax.stairs(n / n.sum(), edges, color=plotstyle.REGION_COLORS[reg], lw=1.8,
                  label=f"{reg} (n={len(v):,})")
        ax.axvline(np.median(v), color=plotstyle.REGION_COLORS[reg], lw=1.4, ls="--")
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_ylabel("Proportion of neurons")
    ax.set_title("a  IBL brain-wide map (mouse)", loc="left")
    ax.set_xlim(0.8, 5)
    ax.legend()

    # right: median + CI, with the distribution behind it
    ax = axs[1]
    ypos = np.arange(len(groups))[::-1]
    labels = []
    for y, (ds, reg) in zip(ypos, groups):
        g = df[(df.dataset == ds) & (df.cosmos == reg)]
        v = g.rp_ms_10.dropna().values
        c = plotstyle.REGION_COLORS[reg]
        plotstyle.box_row(ax, v, y, c, height=0.44)
        labels.append(f"{plotstyle.DATASET_NAMES.get(ds, ds)} {reg}  "
                      f"n={len(v):,} / {g.session_key.nunique()} sess.")
    ax.axvline(2, color="0.7", lw=1, ls=":", zorder=0)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_xlim(0.8, 5)
    ax.set_title("b  Median (95% CI), interquartile box, 5-95% whiskers", loc="left")
    fig.tight_layout()
    plotstyle.save(fig, path)


def fig_definitions(df, path):
    """The Fig 1 comparison under all three candidate timepoints."""
    plotstyle.apply()
    groups = [(ds, reg) for ds in ["steinmetz", "ibl", "allen", "macaque"]
              for reg in REGIONS
              if len(df[(df.dataset == ds) & (df.cosmos == reg)]) >= 20]
    cols = list(RP_COLUMNS)
    fig, axs = plt.subplots(1, len(cols), figsize=(12, max(4, 0.36 * len(groups))),
                            sharey=True)
    ypos = np.arange(len(groups))[::-1]
    for k, (ax, col) in enumerate(zip(axs, cols)):
        for y, (ds, reg) in zip(ypos, groups):
            v = df[(df.dataset == ds) & (df.cosmos == reg)][col].astype(float)
            plotstyle.box_row(ax, v, y, plotstyle.REGION_COLORS[reg])
        ax.axvline(2, color="0.7", lw=1, ls=":", zorder=0)
        ax.set_xlabel(RP_COLUMNS[col].axis)
        ax.set_title(f"{'abc'[k]}  {RP_COLUMNS[col].math}", loc="left")
        allv = df[col].astype(float)
        ax.set_xlim(max(0.3, np.nanpercentile(allv, 0.5)),
                    min(10.2, np.nanpercentile(allv, 99.5)))
    axs[0].set_yticks(ypos)
    axs[0].set_yticklabels(
        [f"{plotstyle.DATASET_NAMES.get(ds, ds)} {reg}" for ds, reg in groups],
        fontsize=8)
    fig.tight_layout()
    plotstyle.save(fig, path)


def fig_firing_rate(df, path):
    """Is the region difference a firing-rate difference? Two ways to look."""
    plotstyle.apply()
    w = reference_weights(df.firing_rate)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2),
                            gridspec_kw={"width_ratios": [1, 1, 1.1]})

    # a: the firing-rate distributions themselves. If these coincide, no
    # amount of firing-rate dependence can confound the region comparison.
    ax = axs[0]
    bins = np.logspace(np.log10(2), np.log10(200), 40)
    for reg in REGIONS:
        v = df[df.cosmos == reg].firing_rate.dropna()
        n, e = np.histogram(v, bins=bins)
        ax.stairs(n / n.sum(), e, color=plotstyle.REGION_COLORS[reg], lw=1.6,
                  label=f"{reg} (median {np.median(v):.1f})")
    ax.set_xscale("log")
    ax.set_xlabel("Firing rate (spikes/s)")
    ax.set_ylabel("Proportion of neurons")
    ax.set_title("a  Firing-rate distributions by region", loc="left")
    plotstyle.plain_log_ticks(ax)
    ax.legend(fontsize=7)

    # b: the interaction. Median recovery time inside each firing-rate bin.
    ax = axs[1]
    centers = np.sqrt(FR_EDGES[:-1] * np.minimum(FR_EDGES[1:], 60))
    for reg in REGIONS:
        g = df[df.cosmos == reg]
        idx = np.digitize(g.firing_rate, FR_EDGES) - 1
        med = [np.nanmedian(g.rp_ms_10[idx == b]) if (idx == b).sum() >= 30
               else np.nan for b in range(len(FR_EDGES) - 1)]
        ax.plot(centers, med, "o-", color=plotstyle.REGION_COLORS[reg], ms=4,
                label=reg)
    ax.set_xscale("log")
    ax.set_xlabel("Firing rate (spikes/s)")
    ax.set_ylabel("Median estimated ACG recovery time (ms)")
    ax.set_title("b  Recovery time within firing-rate bins", loc="left")
    plotstyle.plain_log_ticks(ax)
    ax.legend(fontsize=7)

    # c: observed against standardized median, per dataset and region
    ax = axs[2]
    groups = [(ds, reg) for ds in ["steinmetz", "ibl", "allen", "macaque"]
              for reg in REGIONS
              if len(df[(df.dataset == ds) & (df.cosmos == reg)]) >= 20]
    ypos = np.arange(len(groups))[::-1]
    for y, (ds, reg) in zip(ypos, groups):
        g = df[(df.dataset == ds) & (df.cosmos == reg)]
        c = plotstyle.REGION_COLORS[reg]
        raw = float(np.nanmedian(g.rp_ms_10))
        std, cov = standardized_median(g.firing_rate, g.rp_ms_10.astype(float), w)
        lo, hi = plotstyle.bootstrap_ci(g.rp_ms_10)
        ax.plot([lo, hi], [y, y], color=c, lw=2.2, alpha=0.5,
                solid_capstyle="butt")
        ax.plot(raw, y, "o", color="w", mec=c, mew=1.3, ms=5.5)
        if np.isfinite(std):
            ax.plot(std, y, "|", color=c, ms=10, mew=2)
    ax.axvline(2, color="0.7", lw=1, ls=":", zorder=0)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{plotstyle.DATASET_NAMES.get(ds, ds)} {reg}"
                        for ds, reg in groups], fontsize=7)
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_title("c  Circle, observed median; tick, standardized", loc="left")
    fig.tight_layout()
    plotstyle.save(fig, path)


def firing_rate_section(df, lines):
    """Everything needed to answer 'is this just a firing-rate difference?'."""
    w = reference_weights(df.firing_rate)
    lines += ["", "=" * 60,
              "Firing rate: is the region difference a firing-rate difference?",
              "",
              "Reference firing-rate distribution used for standardization "
              "(all included units):",
              "  " + "  ".join(f"{a:g}-{b:g}:{x:.3f}" for a, b, x in
                               zip(FR_EDGES[:-1], FR_EDGES[1:], w)),
              "",
              "Firing rate by region (spikes/s):"]
    for reg in REGIONS:
        g = df[df.cosmos == reg]
        lines.append(f"  {reg:10s} median {g.firing_rate.median():6.2f}  "
                     f"q25 {g.firing_rate.quantile(.25):5.2f}  "
                     f"q75 {g.firing_rate.quantile(.75):6.2f}")
    lines += ["", "Observed and firing-rate-standardized medians (ms):"]
    for col in RP_COLUMNS:
        lines.append(f"  {RP_COLUMNS[col].text}:")
        for reg in REGIONS:
            g = df[df.cosmos == reg]
            raw = float(np.nanmedian(g[col].astype(float)))
            std, cov = standardized_median(g.firing_rate, g[col].astype(float), w)
            lo, hi = standardized_ci(g.firing_rate, g[col].astype(float), w,
                                     n_boot=300)
            lines.append(f"    {reg:10s} observed {raw:6.3f}  standardized "
                         f"{std:6.3f} [{lo:.3f}, {hi:.3f}]  "
                         f"change {std - raw:+.3f}  coverage {cov:.2f}")
    lines += ["",
              "Median recovery time inside each firing-rate bin (the check for",
              "an interaction; a region difference that survives inside every",
              "bin is not a firing-rate artifact):"]
    t = by_bin_table(df[df.cosmos.isin(REGIONS)], "cosmos", "rp_ms_10")
    piv = t.pivot(index="fr_bin", columns="cosmos", values="median")
    lines.append("  " + piv.round(3).to_string().replace("\n", "\n  "))


def definitions_section(df, lines):
    """The three candidate timepoints side by side."""
    from scipy import stats as sps
    lines += ["", "=" * 60,
              "Three candidate timepoints",
              "",
              "None of these is a refractory period. Each is a different",
              "operational answer to 'how long is this unit quiet for', and the",
              "point of listing all three is to show how much of Fig 1 depends",
              "on which one is used. All are computed on the same units.",
              ""]
    for col, spec in RP_COLUMNS.items():
        v = df[col].astype(float)
        lines.append(f"  {spec.text:28s} median {v.median():6.3f}  "
                     f"q05 {v.quantile(.05):6.3f}  q95 {v.quantile(.95):6.3f}  "
                     f"at the 0.5 ms floor {np.mean(v < 0.55):5.1%}  "
                     f"at the 10 ms ceiling {np.mean(v > 9.9):5.1%}")
    lines += ["", "Median by region (ms):"]
    for reg in REGIONS:
        g = df[df.cosmos == reg]
        lines.append(f"  {reg:10s} " + "  ".join(
            f"{RP_COLUMNS[c].text} {g[c].astype(float).median():6.3f}"
            for c in RP_COLUMNS))
    lines += ["", "Unit-level Spearman correlation between definitions:"]
    cols = list(RP_COLUMNS)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = df[cols[i]].astype(float), df[cols[j]].astype(float)
            m = a.notna() & b.notna()
            r = sps.spearmanr(a[m], b[m])
            lines.append(f"  {RP_COLUMNS[cols[i]].text:28s} vs "
                         f"{RP_COLUMNS[cols[j]].text:28s} rho = "
                         f"{r.statistic:+.3f}")
    lines += ["",
              "Caveats specific to each: tau_r at C_min piles up at the tau_min",
              "boundary (0.5 ms) for the fraction shown above, because for many",
              "units the tightest contamination bound comes from the shortest",
              "window; the last accepted tau_r is right-censored at the 10 ms",
              "edge of the tested range for the fraction shown above."]


def main():
    df = load_all()
    if df.empty:
        print("no enriched tables yet; run enrich_table.py first")
        return
    lines = ["Fig 1 re-analysis", "=" * 60, ""]
    lines.append(f"Loaded {len(df):,} units from {df.dataset.nunique()} datasets, "
                 f"{df.insertion_key.nunique()} insertions, {df.animal_key.nunique()} animals")
    lines.append(f"RP estimate available for {df.rp_ms_10.notna().sum():,} units "
                 f"({df.rp_ms_10.notna().mean():.1%} of all, "
                 f"{df.loc[df.firing_rate>=2, 'rp_ms_10'].notna().mean():.1%} of those >= 2 spikes/s)")

    # --- inclusion rules ---------------------------------------------------
    sel_common = rule_common(df, require_pass=False)
    sel_common_pass = rule_common(df, require_pass=True)
    sel_nolabel = rule_common(df, require_pass=True, require_sorter_good=False)
    sel_mouse_rule = rule_mouse(df)
    lines += ["", "Inclusion rules:",
              f"  common rule, sorter-good, no Sliding RP filter: {sel_common.sum():,} units",
              f"  common rule, sorter-good + Sliding RP accept  : {sel_common_pass.sum():,} units",
              f"  same but WITHOUT the sorter label             : {sel_nolabel.sum():,} units",
              f"  submitted mouse rule (>=2 spikes/s)           : {sel_mouse_rule.sum():,} units",
              "",
              "  Sorter label vs metric acceptance (median estimate, ms):"]
    for nm, sl in (("sorter-good + accepted", sel_common_pass),
                   ("accepted only, no label", sel_nolabel),
                   ("sorter-good only", sel_common)):
        for ds, g in df[sl].groupby("dataset"):
            if len(g) > 200:
                lines.append(f"    {nm:24s} {ds:10s} n={len(g):7,} "
                             f"median {g.rp_ms_10.median():.3f}")

    inc = df[sel_common_pass].copy()
    inc.to_parquet(OUT_TABLE)

    lines += ["", "Median estimated ACG recovery time by dataset and region",
              "(common inclusion rule, accepted by Sliding RP)", ""]
    st = group_stats(inc)
    st = st[st.cosmos.isin(REGIONS)].sort_values(["dataset", "cosmos"])
    lines.append(st.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # --- the four requested diagnostics -----------------------------------
    lines += ["", "Diagnostic 1: how often the first-nonzero-bin floor fires"]
    for (ds, reg), g in inc[inc.cosmos.isin(REGIONS)].groupby(["dataset", "cosmos"]):
        lines.append(f"  {ds:10s} {reg:10s} {g.rp_floor.mean():.1%} of {len(g):,} units")

    lines += ["", "Diagnostic 2: dependence on firing rate (Spearman, within region)"]
    from scipy import stats as sps
    for reg in REGIONS:
        g = inc[inc.cosmos == reg]
        if len(g) > 50:
            r = sps.spearmanr(g.firing_rate, g.rp_ms_10, nan_policy="omit")
            lines.append(f"  {reg:10s} rho = {r.statistic:+.3f} (p = {r.pvalue:.2g}, "
                         f"n = {len(g):,})")
    lines.append("  median RP by firing-rate octile (all regions pooled):")
    q = pd.qcut(inc.firing_rate, 8, duplicates="drop")
    for iv, g in inc.groupby(q, observed=True):
        lines.append(f"    {iv} n={len(g):,} median RP {np.nanmedian(g.rp_ms_10):.3f} ms")

    lines += ["", "Diagnostic 3: do the conclusions survive a common inclusion rule?"]
    for name, sel in (("submitted mouse rule", sel_mouse_rule),
                      ("common rule, no acceptance filter", sel_common),
                      ("common rule + acceptance filter", sel_common_pass)):
        g = df[sel]
        meds = g[g.cosmos.isin(REGIONS)].groupby("cosmos").rp_ms_10.median()
        lines.append(f"  {name:28s} " +
                     "  ".join(f"{k} {v:.3f}" for k, v in meds.items()))

    lines += ["", "Diagnostic 4: sensitivity to the recovery fraction"]
    for frac, col in ((0.05, "rp_ms_05"), (0.10, "rp_ms_10"), (0.20, "rp_ms_20")):
        meds = inc[inc.cosmos.isin(REGIONS)].groupby("cosmos")[col].median()
        order = "<".join(meds.sort_values().index)
        lines.append(f"  {frac:.0%} recovery: " +
                     "  ".join(f"{k} {v:.3f}" for k, v in meds.items()) +
                     f"   ordering {order}")

    # --- hierarchical statistics ------------------------------------------
    lines += ["", "=" * 60, "Hierarchical statistics"]
    hierarchical(inc, lines)

    # --- the two additions Nick asked for ---------------------------------
    definitions_section(inc[inc.cosmos.isin(REGIONS)], lines)
    firing_rate_section(inc[inc.cosmos.isin(REGIONS)], lines)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "fig1_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines[:60]))
    print(f"\n... full output in {OUTDIR/'fig1_numbers.txt'}")

    (OUTDIR / "figures").mkdir(exist_ok=True)
    make_figure(inc, st, OUTDIR / "figures" / "fig1")
    fig_definitions(inc[inc.cosmos.isin(REGIONS)],
                    OUTDIR / "figures" / "fig1_definitions")
    fig_firing_rate(inc[inc.cosmos.isin(REGIONS)],
                    OUTDIR / "figures" / "fig1_firing_rate")


if __name__ == "__main__":
    main()
