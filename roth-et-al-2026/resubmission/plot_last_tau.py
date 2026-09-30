"""The longest window at which each accepted unit is still accepted (02).

Companion to the tau_Cmin histogram (realdata panel b). tau_Cmin says where the
metric found its tightest bound; the last accepted tau_r says how far the clean
stretch of the autocorrelogram extends. A unit whose last accepted window sits
just above tau_min got through only on the very shortest windows, which is
where a sorter censor window (08) inflates the confidence, so the proportion of
such units is compared across datasets with and without a censor.

Population: sorter-"good" units (IBL label 1, Allen 'good', Steinmetz
phy_annotation >= 2, all macaque units) that Sliding RP accepts, in isocortex,
hippocampal formation and thalamus, i.e. the Fig 1 sorter criterion. The looser
population of realdata panel b (any positive IBL label) is reported in the
numbers file too: it contains a broad peak at 1.40-1.53 ms that comes almost
entirely from IBL label-0.33 units (failing two of IBL's three quality
criteria) and is absent from label-1 units.

Run:  python plot_last_tau.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402
from load_enriched import load_all  # noqa: E402

OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/02_realdata_pass_rates")
REGIONS = ["Isocortex", "HPF", "TH"]
DATASETS = ["ibl", "allen", "steinmetz", "macaque"]
CENSOR_MS = {"ibl": 0.0, "allen": 0.1, "steinmetz": 0.0, "macaque": 0.233}   # 08
SHORT_MS = 0.75          # "accepted only on the shortest windows": last tau <= this
CEIL_MS = 9.9            # right-censored at the 10 ms edge


def main():
    plotstyle.apply()
    df = load_all(verbose=False)
    loose = df[(df.sorter_label.isna()) | (df.sorter_label > 0)]
    acc = df[(df.sorter_label.isna()) | (df.sorter_label >= 1)]
    p = acc[acc.passes & acc.cosmos.isin(REGIONS)]
    pl = loose[loose.passes & loose.cosmos.isin(REGIONS)]
    bins = np.logspace(np.log10(0.45), np.log10(10.5), 46)

    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6), sharey=True)
    L = ["Longest accepted window for accepted units", "=" * 43, "",
         f"'Short only' = last accepted tau_r <= {SHORT_MS} ms; 'ceiling' = accepted",
         "right up to the 10 ms edge of the tested range.", ""]
    L.append(f"{'dataset':10s} {'region':10s} {'n':>7} {'median':>7} {'short only':>11} "
             f"{'ceiling':>8} {'tau_Cmin<0.6':>13} {'censor':>7}")
    for ax, ds in zip(axs, DATASETS):
        for reg in REGIONS:
            g = p[(p.dataset == ds) & (p.cosmos == reg)]
            t = g.rp_last_ms.dropna().values
            if t.size < 50:
                continue
            n, e = np.histogram(np.clip(t, bins[0], bins[-1] * 0.999), bins=bins)
            ax.stairs(n / n.sum(), e, color=plotstyle.REGION_COLORS[reg], lw=1.6,
                      label=f"{reg} (n={t.size:,})")
            L.append(f"{ds:10s} {reg:10s} {t.size:7,} {np.median(t):7.2f} "
                     f"{np.mean(t <= SHORT_MS):11.1%} {np.mean(t > CEIL_MS):8.1%} "
                     f"{np.mean(g.rp_Cmin_ms < 0.6):13.1%} {CENSOR_MS[ds]:6.3f}")
        ax.axvline(0.5, color="0.5", ls=":", lw=1)
        ax.set_xscale("log")
        ax.set_xticks([0.5, 1, 2, 5, 10])
        ax.set_xticklabels(["0.5", "1", "2", "5", "10"])
        ax.set_xlabel(r"Longest accepted $\tau_r$ (ms)")
        cz = CENSOR_MS[ds]
        ax.set_title(f"{plotstyle.DATASET_NAMES[ds]}  (censor "
                     f"{'none' if cz == 0 else f'{cz:g} ms'})", loc="left")
        leg = ax.legend(handlelength=0, handletextpad=0, fontsize=8)
        for tx, h in zip(leg.get_texts(), leg.legend_handles):
            tx.set_color(h.get_edgecolor() if hasattr(h, "get_edgecolor") else h.get_color())
            h.set_visible(False)
    axs[0].set_ylabel("Proportion of accepted units")
    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "last_accepted_tau")

    for name, pop in (("sorter-good units (figure)", p), ("any positive label (panel b)", pl)):
        L += ["", f"Pooled over regions, {name}:"]
        for ds in DATASETS:
            t = pop[pop.dataset == ds].rp_last_ms.dropna()
            L.append(f"  {plotstyle.DATASET_NAMES[ds]:10s} n={len(t):7,}  short only "
                     f"{np.mean(t <= SHORT_MS):6.1%}   median {t.median():.2f} ms   "
                     f"censor {CENSOR_MS[ds]:g} ms")
    (OUTDIR / "last_accepted_tau.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
