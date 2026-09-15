"""Shared plotting conventions for the resubmission figures.

Arial throughout, top and right spines hidden, editable text in vector output,
sentence-case labels with units, firing rates in spikes/s (not Hz).
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

FIGDIR = None  # set by the caller


def apply():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
    })


# Colour conventions carried over from the manuscript figures
REGION_COLORS = {"Isocortex": "#d95f02", "HPF": "#1f78b4", "TH": "#1b9e77"}
DATASET_MARKERS = {"ibl": "o", "allen": "s", "steinmetz": "^", "macaque": "D"}
#: how each dataset is spelled in a figure label
DATASET_NAMES = {"ibl": "IBL", "allen": "Allen", "steinmetz": "Steinmetz",
                 "macaque": "Macaque"}


def bootstrap_ci(x, fn=np.median, n=2000, seed=0):
    """Percentile bootstrap interval for a summary statistic."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    bs = fn(x[rng.integers(0, x.size, (n, x.size))], axis=1)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def box_row(ax, values, y, color, height=0.44, ci=True, seed=0):
    """One horizontal box-and-whisker row, the Fig 1 convention.

    Thin line, 5th to 95th percentile; shaded box, interquartile range; thick
    bar, bootstrapped 95% confidence interval of the median; open circle, the
    median. Returns (median, ci_lo, ci_hi) so the caller can tabulate the same
    numbers it just drew.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return np.nan, np.nan, np.nan
    p5, q1, med, q3, p95 = np.percentile(v, [5, 25, 50, 75, 95])
    ax.plot([p5, p95], [y, y], color=color, lw=0.8, alpha=0.6, zorder=1)
    ax.add_patch(plt.Rectangle((q1, y - height / 2), q3 - q1, height,
                               facecolor=color, alpha=0.25, edgecolor="none",
                               zorder=2))
    lo, hi = bootstrap_ci(v, seed=seed) if ci else (np.nan, np.nan)
    if np.isfinite(lo):
        ax.plot([lo, hi], [y, y], color=color, lw=2.6, zorder=3,
                solid_capstyle="butt")
    ax.plot(med, y, "o", color="w", mec=color, mew=1.4, ms=5.5, zorder=4)
    return float(med), lo, hi


def plain_log_ticks(ax, which="x", ticks=(2, 3, 5, 10, 20, 50, 100)):
    """Label a log axis with plain numbers instead of 3 x 10^0."""
    a = ax.xaxis if which == "x" else ax.yaxis
    lo, hi = (ax.get_xlim() if which == "x" else ax.get_ylim())
    t = [v for v in ticks if lo <= v <= hi]
    a.set_major_locator(mpl.ticker.FixedLocator(t))
    a.set_minor_locator(mpl.ticker.NullLocator())
    a.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))


def confidence_cmap(n):
    """Teal ramp used for the confidence-threshold families in Fig 4."""
    return plt.cm.viridis(np.linspace(0.15, 0.85, n))


def save(fig, path, also_png=True):
    fig.savefig(str(path) + ".pdf", bbox_inches="tight")
    if also_png:
        fig.savefig(str(path) + ".png", bbox_inches="tight", dpi=200)
    print("wrote", path, flush=True)
