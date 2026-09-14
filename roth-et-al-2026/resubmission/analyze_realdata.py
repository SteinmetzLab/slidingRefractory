"""Work package 02: what the sliding procedure actually changes on real data.

Reviewer 2, Major comment 5, and Reviewer 1's selection-bias minor.

Outputs (into 02_realdata_pass_rates/):
    figures/realdata.pdf      pass rates, selected tau_r, and who gets excluded
    realdata_numbers.txt      the numbers for the text and response letter

Run:  python analyze_realdata.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402
from load_enriched import load_all  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"02_realdata_pass_rates")
REGIONS = ["Isocortex", "HPF", "TH"]
TAU_MIN_COLS = [("pass_taumin_0p25", 0.25), ("pass_taumin_0p5", 0.5),
                ("pass_taumin_0p75", 0.75), ("pass_taumin_1", 1.0),
                ("pass_taumin_1p5", 1.5)]


def main():
    plotstyle.apply()
    df = load_all()
    if df.empty:
        print("no enriched tables yet")
        return
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)
    lines = ["Real-data analysis: what Sliding RP changes", "=" * 60, ""]

    # sorter-accepted units are the main-text population
    acc = df[(df.sorter_label.isna()) | (df.sorter_label > 0)].copy()
    lines.append(f"All units {len(df):,}; sorter-accepted {len(acc):,}")

    # --- 1. pass rates by dataset x region --------------------------------
    lines += ["", "Acceptance rates by dataset and region (sorter-accepted units)", ""]
    rows = []
    for (ds, reg), g in acc[acc.cosmos.isin(REGIONS)].groupby(["dataset", "cosmos"]):
        if len(g) < 50:
            continue
        rows.append(dict(
            dataset=ds, region=reg, n=len(g),
            sliding=g.passes.mean(), hl2=g.hl2_pass.mean(), hl3=g.hl3_pass.mean(),
            rescued_vs_hl3=((g.passes) & (~g.hl3_pass)).mean(),
            lost_vs_hl3=((~g.passes) & (g.hl3_pass)).mean(),
            rescued_vs_hl2=((g.passes) & (~g.hl2_pass)).mean(),
            lost_vs_hl2=((~g.passes) & (g.hl2_pass)).mean(),
            median_fr=g.firing_rate.median()))
    pr = pd.DataFrame(rows)
    lines.append(pr.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    pr.to_csv(OUTDIR / "pass_rates.csv", index=False)

    # --- 2. distribution of the selected tau_r ----------------------------
    lines += ["", "Selected tau_r for accepted units (ms)", ""]
    p = acc[acc.passes]
    for (ds, reg), g in p[p.cosmos.isin(REGIONS)].groupby(["dataset", "cosmos"]):
        if len(g) < 50:
            continue
        t = g.tau_Cmin.dropna() * 1000
        tf = g.tau_first_pass.dropna() * 1000
        lines.append(
            f"  {ds:10s} {reg:10s} n={len(g):6,}  tau_Cmin median {t.median():.2f} "
            f"[{t.quantile(.25):.2f}, {t.quantile(.75):.2f}]  "
            f"at tau_min (<0.6 ms) {np.mean(t < 0.6):.1%}  "
            f"tau_first_pass median {tf.median():.2f}")
    lines += ["", "NOTE for interpretation: in simulation tau_Cmin systematically",
              "underestimates the true refractory period by about 30% (true RP",
              "1/2/3/5 ms -> median tau_Cmin 0.85/1.48/2.10/3.10 ms), so these",
              "values should not be read as refractory-period estimates."]

    # --- 3. sorter short-lag shadow and the tau_min sweep -----------------
    lines += ["", "Short-lag structure by dataset (sorter-dependent)", ""]
    for ds, g in acc.groupby("dataset"):
        fz = g.first_nonzero_bin.replace(-1, np.nan) / 30.0   # ms
        lines.append(f"  {ds:10s} first nonzero ACG bin: median {fz.median():.3f} ms, "
                     f"{np.mean(fz > 0.5):.1%} of units above 0.5 ms")
    lines += ["", "Acceptance rate as tau_min is raised (sorter-accepted units)", ""]
    hdr = f"  {'dataset':10s} " + " ".join(f"{t:>6.2f}ms" for _, t in TAU_MIN_COLS)
    lines.append(hdr)
    for ds, g in acc.groupby("dataset"):
        vals = [g[c].mean() if c in g else np.nan for c, _ in TAU_MIN_COLS]
        lines.append(f"  {ds:10s} " + " ".join(f"{v:>8.3f}" for v in vals))
    lines.append("  (a large drop from 0.5 to 1.0 ms would mean many acceptances rely on "
                 "the 0.5-1 ms window, where sorter duplicate-removal acts)")

    # --- 4. selection bias by firing rate (Reviewer 1) --------------------
    lines += ["", "Who gets excluded: rejections with vs without observed violations", ""]
    for (ds, reg), g in acc[acc.cosmos.isin(REGIONS)].groupby(["dataset", "cosmos"]):
        if len(g) < 50:
            continue
        fail = g[~g.passes]
        clean_fail = (fail.n_viol_short == 0).mean() if len(fail) else np.nan
        underpowered = (fail.tau_pass0 > 0.010).mean() if len(fail) else np.nan
        lines.append(
            f"  {ds:10s} {reg:10s} fail {len(fail):6,}/{len(g):6,} "
            f"({len(fail)/len(g):.1%}); of failures, {clean_fail:.1%} had no "
            f"violation below 2 ms and {underpowered:.1%} could not have been accepted "
            f"at any tau_r")

    # --- figure ------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))

    ax = axs[0, 0]
    if len(pr):
        x = np.arange(len(pr))
        w = 0.26
        ax.bar(x - w, pr.sliding, w, label="Sliding RP", color="#1b9e77")
        ax.bar(x, pr.hl2, w, label="Hill-Llobet, 2 ms", color="#d95f02")
        ax.bar(x + w, pr.hl3, w, label="Hill-Llobet, 3 ms", color="#7570b3")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r.dataset}\n{r.region}" for r in pr.itertuples()],
                           fontsize=6.5)
        ax.set_ylabel("Proportion of units accepted")
        ax.set_title("a  Acceptance rates by dataset and region", loc="left")
        ax.legend()

    ax = axs[0, 1]
    for reg in REGIONS:
        g = p[(p.cosmos == reg)]
        t = (g.tau_Cmin.dropna() * 1000).values
        if t.size < 50:
            continue
        n, e = np.histogram(t, bins=np.logspace(np.log10(0.4), np.log10(10), 40))
        ax.stairs(n / n.sum(), e, color=plotstyle.REGION_COLORS[reg], lw=1.6,
                  label=f"{reg} (n={t.size:,})")
    ax.axvline(0.5, color="0.5", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("Selected refractory duration, tau_Cmin (ms)")
    ax.set_ylabel("Proportion of accepted units")
    ax.set_title("b  Where the metric finds its best window", loc="left")
    ax.legend()

    ax = axs[1, 0]
    for ds, g in acc.groupby("dataset"):
        vals = [g[c].mean() if c in g else np.nan for c, _ in TAU_MIN_COLS]
        ax.plot([t for _, t in TAU_MIN_COLS], vals, "o-", ms=4, label=ds)
    ax.set_xlabel("tau_min (ms)")
    ax.set_ylabel("Proportion of units accepted")
    ax.set_title("c  Sensitivity to the short-lag exclusion", loc="left")
    ax.legend()

    ax = axs[1, 1]
    bins = np.logspace(np.log10(0.05), np.log10(50), 24)
    for reg in REGIONS:
        g = acc[acc.cosmos == reg]
        if len(g) < 100:
            continue
        idx = np.digitize(g.firing_rate, bins)
        xs, ys = [], []
        for i in range(1, len(bins)):
            s = g[idx == i]
            if len(s) >= 30:
                xs.append(np.sqrt(bins[i - 1] * bins[i]))
                ys.append(s.passes.mean())
        ax.plot(xs, ys, "o-", ms=3.5, color=plotstyle.REGION_COLORS[reg], label=reg)
    ax.set_xscale("log")
    ax.set_xlabel("Firing rate (spikes/s)")
    ax.set_ylabel("Proportion accepted")
    ax.set_title("d  Acceptance rate depends strongly on firing rate", loc="left")
    ax.legend()

    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "realdata")

    (OUTDIR / "realdata_numbers.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
