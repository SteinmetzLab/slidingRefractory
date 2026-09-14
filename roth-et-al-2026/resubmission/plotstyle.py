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


def confidence_cmap(n):
    """Teal ramp used for the confidence-threshold families in Fig 4."""
    return plt.cm.viridis(np.linspace(0.15, 0.85, n))


def save(fig, path, also_png=True):
    fig.savefig(str(path) + ".pdf", bbox_inches="tight")
    if also_png:
        fig.savefig(str(path) + ".png", bbox_inches="tight", dpi=200)
    print("wrote", path, flush=True)
