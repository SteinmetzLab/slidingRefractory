"""Load the enriched per-unit tables into one harmonised DataFrame.

Harmonises the dataset-specific columns onto a common schema so Fig 1 (work
package 01) and the real-data pass-rate analysis (work package 02) can treat
all five datasets the same way:

    dataset species animal session insertion cluster_id acronym cosmos beryl
    n_spikes rec_dur_s firing_rate sorter_label acg_source
    + all the metric and RP-estimate columns added by enrich_table.py

``sorter_label`` is put on a common 0-1 scale where possible: IBL's ``label``
is already 0/0.33/0.67/1; Allen's ``quality`` becomes 1 for 'good' and 0 for
'noise'; Steinmetz's ``_phy_annotation`` >= 2 becomes 1; the macaque data has
no sorter label (NaN).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ENRICHED = Path(r"D:/temp/slidingRP_resub/enriched")
BWM_FIXTURE = None  # resolved lazily


def _bwm_map():
    """pid -> (subject, lab, date) from the brain-wide map release fixture."""
    global BWM_FIXTURE
    if BWM_FIXTURE is None:
        try:
            import brainwidemap
            p = Path(brainwidemap.__file__).parent / "fixtures" / "2023_12_bwm_release.csv"
            BWM_FIXTURE = pd.read_csv(p).set_index("pid")
        except Exception:  # noqa: BLE001
            BWM_FIXTURE = pd.DataFrame()
    return BWM_FIXTURE


def load(dataset, verbose=True):
    files = sorted((ENRICHED / dataset).glob("*.pqt"))
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)

    if dataset == "ibl":
        bwm = _bwm_map()
        df["dataset"] = "ibl"
        df["species"] = "mouse"
        df["insertion"] = df["pid"]
        df["session"] = df["eid"]
        if len(bwm):
            df["animal"] = df["pid"].map(bwm["subject"])
            df["lab"] = df["pid"].map(bwm["lab"])
        else:
            df["animal"] = df["eid"]
        df["firing_rate"] = df["fr_recomputed"]
        df["sorter_label"] = df.get("label", np.nan)
    elif dataset == "allen":
        df["dataset"] = "allen"
        df["species"] = "mouse"
        df["insertion"] = df["probe_id"].astype(str)
        df["animal"] = df.get("specimen_id", df["session"]).astype(str)
        df["session"] = df["session"].astype(str)
        df["firing_rate"] = df["fr_recomputed"]
        q = df.get("quality", pd.Series([""] * len(df)))
        df["sorter_label"] = (q.astype(str).str.contains("good")).astype(float)
        df["cluster_id"] = df.get("cluster_id", df.get("id"))
    elif dataset == "steinmetz":
        df["firing_rate"] = df["fr_recomputed"]
        df["sorter_label"] = (df["phy_annotation"].astype(float) >= 2).astype(float)
    elif dataset == "macaque":
        df["firing_rate"] = df["fr_recomputed"]
        df["sorter_label"] = np.nan

    keep_first = ["dataset", "species", "animal", "session", "insertion",
                  "cluster_id", "acronym", "cosmos", "beryl", "n_spikes",
                  "rec_dur_s", "firing_rate", "sorter_label", "acg_source"]
    for c in keep_first:
        if c not in df:
            df[c] = np.nan
    rest = [c for c in df.columns if c not in keep_first]
    df = df[keep_first + rest]
    if verbose:
        print(f"{dataset}: {len(df):,} units, {df.insertion.nunique()} insertions, "
              f"{df.animal.nunique()} animals", flush=True)
    return df


def load_all(datasets=("ibl", "allen", "steinmetz", "macaque"), verbose=True):
    parts = [load(d, verbose) for d in datasets]
    parts = [p for p in parts if len(p)]
    if not parts:
        return pd.DataFrame()
    df = pd.concat(parts, ignore_index=True)
    # a single animal key that is unique across datasets
    df["animal_key"] = df["dataset"].astype(str) + ":" + df["animal"].astype(str)
    df["insertion_key"] = df["dataset"].astype(str) + ":" + df["insertion"].astype(str)
    if verbose:
        print(f"TOTAL {len(df):,} units, {df.insertion_key.nunique()} insertions, "
              f"{df.animal_key.nunique()} animals", flush=True)
    return df


# --- inclusion rules (work package 01) -------------------------------------

def rule_mouse(df):
    """The submitted manuscript's mouse rule: >= 2 spikes/s and RP > 1 ms."""
    return (df.firing_rate >= 2) & (df.rp_ms_10 > 1) & df.rp_ms_10.notna()


def rule_macaque(df, filt_col="acg_0p5_1ms"):
    """The submitted manuscript's macaque rule: post-RP rate at least 5x the
    pre-RP rate, and RP > 1 ms. Approximated here from the stored ACG summary
    columns; the exact version is recomputed from the ACG when needed."""
    pre = df["acg_0_0p5ms"].astype(float) / 15.0          # 15 bins in 0-0.5 ms
    post = df["n_viol_3ms"].astype(float) / 90.0          # 90 bins in 0-3 ms
    return (post > 5 * pre) & (df.rp_ms_10 > 1) & df.rp_ms_10.notna()


def rule_common(df, min_fr=2.0, min_spikes=5000, min_r2=0.5, require_pass=True):
    """A single rule applied to every dataset (the reviewer's request).

    Deliberately does not condition on the shape of the ACG being measured,
    other than requiring the fit to have succeeded. The Sliding RP requirement
    is optional and reported both ways, since the submitted IBL sample was
    conditioned on an earlier version of the metric.
    """
    m = ((df.firing_rate >= min_fr) & (df.n_spikes >= min_spikes)
         & (df.rp_ms_10 > 1) & df.rp_ms_10.notna() & (df.rp_r2 >= min_r2))
    if require_pass:
        m &= df.passes
    return m
