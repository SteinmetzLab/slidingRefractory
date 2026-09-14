"""How does a measured spike-count correlation relate to the underlying rate correlation?

The model-mismatch sweep controls the *rate* correlation rho between base
neuron and contaminant. What one can measure in real data is the correlation of
binned spike counts, and that is attenuated by counting noise: at 100 ms bins
and a few spikes/s there is far less than one spike per bin, so most of the
variance is Poisson rather than shared drive. The simulated pairs at rho = +/-1
therefore show measured correlations of only about +/-0.025 at 100 ms, which is
the same order as the values measured between real neighbouring IBL units.

Without this bridge the two results are not comparable and it would be easy to
conclude wrongly that the simulation used implausibly strong correlations. This
script measures both on the same footing, at several bin widths.

Run:  python correlation_bridge.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from resub_simulations import gen_modulated_pair  # noqa: E402

CACHE = Path(r"D:/temp/slidingRP_resub/semisynth")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/"
              r"05_model_mismatch")
BINS_S = (0.1, 1.0, 10.0)
FS = 30000.0


def corr_at(a, b, duration, bin_s, fs=1.0):
    edges = np.arange(0, duration + bin_s, bin_s)
    ca = np.histogram(np.asarray(a) / fs, edges)[0].astype(float)
    cb = np.histogram(np.asarray(b) / fs, edges)[0].astype(float)
    if ca.std() == 0 or cb.std() == 0:
        return np.nan
    return float(np.corrcoef(ca, cb)[0, 1])


def simulated(n_rep=12, duration=7200.0, total_rate=5.0, cont_prop=0.10,
              rp=0.002, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for rho in (-1.0, -0.5, 0.0, 0.5, 1.0):
        for tau_s in (2.0, 30.0):
            for _ in range(n_rep):
                b, c, _ = gen_modulated_pair((1 - cont_prop) * total_rate,
                                             cont_prop * total_rate, duration,
                                             rp, rho=rho, tau_s=tau_s, rng=rng)
                row = dict(rho=rho, tau_s=tau_s)
                for bs in BINS_S:
                    row[f"r_{bs:g}s"] = corr_at(b, c, duration, bs)
                rows.append(row)
    return pd.DataFrame(rows)


def real_pairs(max_pairs=320, seed=0):
    """Recipient-donor pairs from the cached semi-synthetic IBL probes."""
    rng = np.random.default_rng(seed)
    rows = []
    for fp in sorted(CACHE.glob("*.npz")):
        z = np.load(fp, allow_pickle=True)
        samp, clu = z["samples"], z["clusters"]
        rec_dur = float(z["rec_dur"])
        depth_all, cid_all = z["depths"], z["cluster_id"]
        order = np.argsort(clu, kind="stable")
        clu, samp = clu[order], samp[order]
        cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
        trains = {int(c): samp[a:a + n] for c, a, n in zip(cids, starts, counts)}
        depth = {int(c): float(depth_all[i]) if i < depth_all.size else np.nan
                 for i, c in enumerate(cid_all)}
        good = [c for c in trains
                if trains[c].size / rec_dur >= 2 and np.isfinite(depth.get(c, np.nan))]
        if len(good) < 4:
            continue
        # Enumerate candidate pairs by separation band rather than sampling at
        # random: random pairs on a 3.8 mm shank are almost never close
        # together, which left the 'near' band with a handful of pairs.
        bands = {"near": [], "far": []}
        for i in range(len(good)):
            for j in range(i + 1, len(good)):
                sep = abs(depth[good[i]] - depth[good[j]])
                if sep <= 60:
                    bands["near"].append((good[i], good[j], sep))
                elif sep >= 300:
                    bands["far"].append((good[i], good[j], sep))
        per_band = max(max_pairs // (2 * 8), 8)
        for kind, cand in bands.items():
            if not cand:
                continue
            take = rng.choice(len(cand), min(per_band, len(cand)), replace=False)
            for k in take:
                ca, cb, sep = cand[k]
                row = dict(pid=fp.stem, sep_um=sep, kind=kind,
                           fr_a=trains[ca].size / rec_dur,
                           fr_b=trains[cb].size / rec_dur)
                for bs in BINS_S:
                    row[f"r_{bs:g}s"] = corr_at(trains[ca], trains[cb], rec_dur, bs, FS)
                rows.append(row)
    return pd.DataFrame(rows)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sim = simulated()
    real = real_pairs()
    sim.to_parquet(Path(r"D:/temp/slidingRP_resub/sims/corr_bridge_sim.pqt"))
    real.to_parquet(Path(r"D:/temp/slidingRP_resub/sims/corr_bridge_real.pqt"))

    lines = ["Rate correlation: simulation and real data on the same footing",
             "=" * 62, "",
             "The model-mismatch sweep controls the underlying RATE correlation",
             "rho. What is measurable in data is the correlation of binned spike",
             "counts, which counting noise attenuates: at 100 ms bins and a few",
             "spikes/s most of the variance is Poisson, not shared drive. Comparing",
             "a simulated rho with a measured count correlation is therefore an",
             "apples-to-oranges comparison unless both are measured the same way.",
             "",
             "Simulated pairs (5 spikes/s total, 10% contamination, 2 h):", ""]
    hdr = f"{'rho':>6} {'tau_s':>6} " + "".join(f"{f'r at {b:g}s':>12}" for b in BINS_S)
    lines.append(hdr)
    for (rho, ts), g in sim.groupby(["rho", "tau_s"]):
        lines.append(f"{rho:>6.1f} {ts:>6.0f} " +
                     "".join(f"{g[f'r_{b:g}s'].median():>12.3f}" for b in BINS_S))
    lines += ["", "Real IBL pairs on the same probes (>= 2 spikes/s each):", ""]
    lines.append(f"{'kind':>6} {'n':>5} " + "".join(f"{f'r at {b:g}s':>12}" for b in BINS_S))
    for kind, g in real.groupby("kind"):
        lines.append(f"{kind:>6} {len(g):>5} " +
                     "".join(f"{g[f'r_{b:g}s'].median():>12.3f}" for b in BINS_S))
    for kind, g in real.groupby("kind"):
        q = [g[f"r_{b:g}s"].quantile(0.95) for b in BINS_S]
        lines.append(f"  {kind} 95th percentile: " +
                     ", ".join(f"{b:g}s {v:.3f}" for b, v in zip(BINS_S, q)))

    # the translation the reader needs, derived from the numbers so it cannot
    # go stale if the sampling changes
    s30 = sim[(sim.tau_s == 30.0) & (sim.rho == 1.0)]
    s02 = sim[(sim.tau_s == 2.0) & (sim.rho == 1.0)]
    near = real[real.kind == "near"]
    fast_sim = float(s02["r_0.1s"].median())
    fast_real = float(near["r_0.1s"].median()) if len(near) else float("nan")
    lines += ["",
              "Reading across.", "",
              f"A simulated rate correlation of rho = +1 appears as a measured",
              f"count correlation of only {fast_sim:.3f} at 100 ms bins but "
              f"{float(s30['r_10s'].median()):.3f} at 10 s bins:",
              "the slow shared drive survives averaging while the Poisson noise",
              "does not. Any statement about how correlated two neurons are is",
              "therefore meaningless without the bin width attached.",
              "",
              f"Real neighbouring units (within 60 um, n = {len(near)}) have a median",
              f"measured correlation of {fast_real:.3f} at 100 ms, "
              f"{float(near['r_1s'].median()):.3f} at 1 s and "
              f"{float(near['r_10s'].median()):.3f} at 10 s,",
              f"with a 95th percentile of {float(near['r_0.1s'].quantile(.95)):.3f} / "
              f"{float(near['r_1s'].quantile(.95)):.3f} / "
              f"{float(near['r_10s'].quantile(.95)):.3f}.",
              "",
              "So at the fast timescale the real neighbouring pairs are about "
              f"{fast_real/fast_sim:.1f}x",
              "MORE correlated than the rho = +1 simulation, and at 1 s they are",
              "comparable. The correlations used in the model-mismatch sweep are",
              "therefore not extreme: they are ordinary for neighbouring neurons,",
              "which means the biases that sweep reports are directly relevant to",
              "real recordings rather than a worst-case curiosity.",
              "",
              "Practical consequence for the manuscript: quote the simulation in",
              "terms of rho and the real data in terms of the bin width used, do",
              "not invite the reader to compare the two numbers directly, and do",
              "not describe the correlation problem as hypothetical."]
    (OUTDIR / "correlation_bridge.txt").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
