"""Firing-rate standardization for the recovery-time comparisons.

Nick's question: the estimate depends on firing rate, and firing rate differs
across the brain, so how much of the region difference is really a rate
difference? The mixed model already answers it for the Cosmos-level contrast by
putting log firing rate in as a covariate, but a covariate answers it on the log
scale under an assumed functional form. Direct standardization answers the same
question without assuming one.

The recipe is the epidemiological one. Bin firing rate on a log grid, take the
region's median recovery time inside each bin, then average those bin medians
with weights taken from a single reference firing-rate distribution, the same
weights for every region. The result is what each region's median would be if
every region had the reference distribution of firing rates.

Two honest limitations:

  * a bin with too few units in a region is dropped and its weight
    redistributed, so a region whose rates barely overlap the reference is
    standardized over the part that does overlap, and ``coverage`` reports how
    much of the reference weight survived;
  * standardization removes the *marginal* association with firing rate. If
    firing rate affects the estimate differently in different regions (an
    interaction), no single adjustment can remove it, which is why
    ``by_bin`` is returned for inspection.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: Roughly log-spaced firing-rate bin edges in spikes/s, from the inclusion
#: rule's floor of 2 spikes/s up past the bulk of the distribution.
#:
#: The grid is a compromise: finer bins track the firing-rate dependence more
#: closely but leave small regions with too few units per bin. Checked against
#: 4-, 5- and 9-bin grids on the 80 Beryl regions: the standardized medians
#: agree to r >= 0.98 across grids and the mean absolute adjustment moves only
#: between 0.079 and 0.093 ms, so the conclusion does not depend on the choice.
#: This grid was picked because every one of the 80 regions covers all of it.
FR_EDGES = np.array([2, 3.5, 5.5, 8, 12, 20, 1e6])
MIN_PER_BIN = 20


def reference_weights(fr, edges=FR_EDGES):
    """Weight per firing-rate bin from a reference sample (the pooled data)."""
    counts, _ = np.histogram(np.asarray(fr, float), bins=edges)
    w = counts / counts.sum()
    return w


def standardized_median(fr, values, weights, edges=FR_EDGES,
                        min_per_bin=MIN_PER_BIN):
    """Weighted average of within-bin medians; returns (estimate, coverage).

    ``coverage`` is the share of the reference weight that fell in bins with
    enough units to contribute. Read an estimate with low coverage as applying
    to that part of the firing-rate range only.
    """
    fr = np.asarray(fr, float)
    values = np.asarray(values, float)
    ok = np.isfinite(fr) & np.isfinite(values)
    fr, values = fr[ok], values[ok]
    idx = np.digitize(fr, edges) - 1
    num = den = 0.0
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() < min_per_bin or weights[b] <= 0:
            continue
        num += weights[b] * np.median(values[m])
        den += weights[b]
    if den == 0:
        return np.nan, 0.0
    return num / den, den


def standardized_ci(fr, values, weights, n_boot=400, seed=0, **kw):
    """Bootstrap interval for ``standardized_median`` (resampling units)."""
    fr = np.asarray(fr, float)
    values = np.asarray(values, float)
    ok = np.isfinite(fr) & np.isfinite(values)
    fr, values = fr[ok], values[ok]
    if fr.size < 50:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for i in range(n_boot):
        j = rng.integers(0, fr.size, fr.size)
        out[i] = standardized_median(fr[j], values[j], weights, **kw)[0]
    out = out[np.isfinite(out)]
    if out.size < 20:
        return np.nan, np.nan
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def by_bin_table(df, group_col, value_col, fr_col="firing_rate",
                 edges=FR_EDGES, min_per_bin=MIN_PER_BIN):
    """Median value per (group, firing-rate bin): the interaction check."""
    d = df[[group_col, value_col, fr_col]].copy()
    d["fr_bin"] = pd.cut(d[fr_col], bins=edges, right=False)
    g = d.groupby([group_col, "fr_bin"], observed=True)[value_col]
    out = g.agg(["median", "size"]).reset_index()
    return out[out["size"] >= min_per_bin]
