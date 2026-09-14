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
from load_enriched import load_all, rule_common, rule_mouse  # noqa: E402

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
            n_units=len(g), n_insertions=g.insertion_key.nunique(),
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
        v = df[(df.dataset == ds) & (df.cosmos == reg)].rp_ms_10.dropna().values
        c = plotstyle.REGION_COLORS[reg]
        # distribution: box from the quartiles with 5-95% whiskers
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        p5, p95 = np.percentile(v, [5, 95])
        ax.plot([p5, p95], [y, y], color=c, lw=0.8, alpha=0.6, zorder=1)
        ax.add_patch(plt.Rectangle((q1, y - 0.22), q3 - q1, 0.44, facecolor=c,
                                   alpha=0.25, edgecolor="none", zorder=2))
        lo, hi = bootstrap_ci(v)
        ax.plot([lo, hi], [y, y], color=c, lw=2.6, zorder=3, solid_capstyle="butt")
        ax.plot(med, y, "o", color="w", mec=c, mew=1.4, ms=5.5, zorder=4)
        labels.append(f"{ds} {reg}  n={len(v):,} / {df[(df.dataset==ds)&(df.cosmos==reg)].insertion_key.nunique()} ins")
    ax.axvline(2, color="0.7", lw=1, ls=":", zorder=0)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Estimated ACG recovery time (ms)")
    ax.set_xlim(0.8, 5)
    ax.set_title("b  Median (95% CI), interquartile box, 5-95% whiskers", loc="left")
    fig.tight_layout()
    plotstyle.save(fig, path)


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
    sel_mouse_rule = rule_mouse(df)
    lines += ["", "Inclusion rules:",
              f"  common rule (no Sliding RP filter): {sel_common.sum():,} units",
              f"  common rule + Sliding RP pass     : {sel_common_pass.sum():,} units",
              f"  submitted mouse rule (>=2 sp/s)   : {sel_mouse_rule.sum():,} units"]

    inc = df[sel_common_pass].copy()
    inc.to_parquet(OUT_TABLE)

    lines += ["", "Median estimated ACG recovery time by dataset and region",
              "(common inclusion rule, Sliding RP pass)", ""]
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
                      ("common rule, no pass filter", sel_common),
                      ("common rule + pass filter", sel_common_pass)):
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

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "fig1_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines[:60]))
    print(f"\n... full output in {OUTDIR/'fig1_numbers.txt'}")

    (OUTDIR / "figures").mkdir(exist_ok=True)
    make_figure(inc, st, OUTDIR / "figures" / "fig1")


if __name__ == "__main__":
    main()
