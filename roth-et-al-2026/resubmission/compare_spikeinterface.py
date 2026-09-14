"""Work package 07: how SpikeInterface's sliding RP differs from this package.

State of SpikeInterface `main` as of 2026-09-13 (checked against the source on
GitHub):

* PR #4682 (merged 2026-07-21) corrected `_compute_violations` to the Llobet
  expected-violation formula and extended the contamination grid to include
  35%, so the *formula* now matches the manuscript. No further fix is needed
  there.
* Remaining differences:
  - `bin_size_ms = 0.25` by default. Note SpikeInterface converts this to
    samples with `max(int(0.25/1000*fs), 1)`, which at 30 kHz truncates 7.5 to
    **7 samples = 0.2333 ms**, not 0.25 ms. This package uses one sample
    (1/30000 s = 0.0333 ms).
  - the pass comparison is `conf > 0.9`, this package uses `>=`.
  - the only returned value is the minimum contamination confirmable at 90%
    confidence; there is no max confidence, no selected tau_r, no pass flag and
    no short-ISI violation count.
  - correlograms are computed in integer sample space, which is what this
    package now does for IBL and macaque data as well.

This script quantifies the one difference that can change a decision: the ACG
bin width. It re-bins the stored 1-sample ACGs into SpikeInterface's 7-sample
bins and recomputes the metric, so the comparison isolates binning from every
other difference.

Run:  python compare_spikeinterface.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import computeViol, slidingRP_from_acg  # noqa: E402

TABLES = Path(r"D:/temp/slidingRP_resub/acg_tables")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"07_parameter_sensitivity_docs")
FS = 30000.0
SI_BIN_SAMPLES = int(0.25 / 1000 * FS)      # = 7, SpikeInterface's truncation


def rebin(acg, k=SI_BIN_SAMPLES):
    """Sum an ACG in k-sample blocks (SpikeInterface's coarser binning)."""
    n = (acg.size // k) * k
    return acg[:n].reshape(-1, k).sum(1)


def metric_at_bin(acg, n_spikes, rec_dur, k, cont_thresh=10.0, conf_thresh=90.0,
                  rp_reject=0.0005, strict_gt=False):
    """Sliding RP on a coarsely binned ACG, everything else identical."""
    counts = rebin(np.asarray(acg, float), k)
    bin_s = k / FS
    centers = np.arange(counts.size) * bin_s + bin_s / 2
    ref_dur = centers + bin_s / 2
    obs = np.cumsum(counts)
    conf = 100 * computeViol(obs, n_spikes, ref_dur, cont_thresh / 100, rec_dur)[0]
    test = centers > rp_reject
    if not test.any():
        return 0.0, np.nan, False
    max_conf = float(np.max(conf[test]))
    passes = max_conf > conf_thresh if strict_gt else max_conf >= conf_thresh
    lam = stats.chi2.ppf(conf_thresh / 100, 2 * (obs + 1)) / 2
    disc = (n_spikes - 0.5) ** 2 - lam * rec_dur / ref_dur
    cmin = np.full(ref_dur.shape, np.nan)
    good = disc >= 0
    cmin[good] = ((n_spikes - 0.5) - np.sqrt(disc[good])) / n_spikes * 100
    ct = cmin[test]
    min_cont = float(np.nanmin(ct)) if not np.all(np.isnan(ct)) else np.nan
    if np.isfinite(min_cont) and min_cont > 35:
        min_cont = np.nan
    return max_conf, min_cont, bool(passes)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for dataset in ("ibl", "allen", "steinmetz", "macaque"):
        files = sorted((TABLES / dataset).glob("*.npy"))[:6]
        for f in files:
            acgs = np.load(f)
            tbl = pd.read_parquet(f.with_suffix(".pqt"))
            n = tbl["n_spikes"].to_numpy()
            d = tbl["rec_dur_s"].to_numpy()
            take = np.flatnonzero((n > 500) & (n / d > 0.5))
            if take.size > 400:
                take = np.random.default_rng(0).choice(take, 400, replace=False)
            for i in take:
                a = acgs[i].astype(float)
                ours = slidingRP_from_acg(a, int(n[i]), float(d[i]))
                mc, cmin, ps = metric_at_bin(a, int(n[i]), float(d[i]),
                                             SI_BIN_SAMPLES, strict_gt=True)
                rows.append(dict(dataset=dataset, file=f.stem, idx=int(i),
                                 n_spikes=int(n[i]), fr=float(n[i] / d[i]),
                                 ours_conf=ours["max_conf"], ours_cont=ours["min_cont"],
                                 ours_pass=ours["passes"],
                                 si_conf=mc, si_cont=cmin, si_pass=ps))
    df = pd.DataFrame(rows)
    if df.empty:
        print("no ACG tables found yet")
        return
    df.to_parquet(Path(r"D:/temp/slidingRP_resub/sims/si_comparison.pqt"))

    lines = ["SpikeInterface vs this package: the effect of the ACG bin width",
             "=" * 66, "",
             "SpikeInterface main (after PR #4682, merged 2026-07-21) uses the same",
             "Llobet expected-violation formula as this package. The remaining",
             "difference that can change a decision is the ACG bin width: it",
             "defaults to bin_size_ms = 0.25, which its sample conversion",
             "max(int(0.25/1000*fs), 1) truncates to 7 samples = 0.2333 ms at",
             "30 kHz. This package tests every sample (0.0333 ms).",
             "",
             f"Compared on {len(df):,} real units from "
             f"{df.dataset.nunique()} datasets ({df.file.nunique()} recordings),",
             "re-binning the same autocorrelograms so only the bin width differs.",
             ""]
    agree = (df.ours_pass == df.si_pass).mean()
    lines.append(f"Accept/reject agreement: {agree:.3%} "
                 f"({(df.ours_pass != df.si_pass).sum():,} of {len(df):,} differ)")
    lines.append(f"  ours accepts, SI rejects : {((df.ours_pass) & (~df.si_pass)).sum():,}")
    lines.append(f"  SI accepts, ours rejects : {((~df.ours_pass) & (df.si_pass)).sum():,}")
    both = df.ours_cont.notna() & df.si_cont.notna()
    if both.any():
        d_ = df.loc[both, "si_cont"] - df.loc[both, "ours_cont"]
        lines.append(f"Minimum contamination: median difference (SI - ours) "
                     f"{d_.median():+.3f} percentage points, "
                     f"IQR [{d_.quantile(.25):+.3f}, {d_.quantile(.75):+.3f}]")
    dc = df.si_conf - df.ours_conf
    lines.append(f"Max confidence: median difference {dc.median():+.2f} points, "
                 f"IQR [{dc.quantile(.25):+.2f}, {dc.quantile(.75):+.2f}]")
    lines += ["", "By dataset:", ""]
    for ds, g in df.groupby("dataset"):
        lines.append(f"  {ds:10s} n={len(g):5,}  agreement {(g.ours_pass==g.si_pass).mean():.3%}  "
                     f"ours accepts {g.ours_pass.mean():.3f}  SI accepts {g.si_pass.mean():.3f}")
    lines += ["", "Reading: a coarser bin cannot resolve short refractory windows,",
              "so it tends to lump genuine clean windows together with the bins",
              "just past the refractory period. The direction and size of the",
              "effect are given above; the Methods should state the bin width",
              "each implementation uses rather than implying they are identical."]
    (OUTDIR / "spikeinterface_comparison.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
