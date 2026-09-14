"""Semi-synthetic contamination on real IBL recordings.

Reviewer 2, Major comment 3, last paragraph:

    "An especially compelling validation, if feasible, would be semi-synthetic
    contamination: start with well-isolated recorded units and inject known
    fractions of spikes from another simultaneously recorded unit. This would
    provide realistic refractory structure, bursting, behavioral modulation, and
    correlation while preserving known contamination ground truth."

Procedure
---------
For each probe insertion:

1. **Recipients**: units that pass Sliding RP at a strict setting (5%
   contamination threshold, 99% confidence) and fire at least 2 spikes/s, so
   their own contamination is both small and quantified. Their minimum
   confirmable contamination at the strict setting is carried through as the
   baseline that adds to the injected fraction.
2. **Donors**: other units on the same probe. Two classes:
   - *near*: within `near_um` of the recipient's peak channel depth, which is
     where real misassignment happens and where firing is most correlated;
   - *far*: more than `far_um` away, an approximately independent contaminant.
3. **Injection**: a uniformly random subset of the donor's spikes is merged into
   the recipient's train so that the injected spikes are a fraction `f` of the
   combined total. Random thinning preserves the donor's refractory structure
   (no thinned train can have shorter ISIs than the donor did) along with its
   bursting and behavioral modulation. A contiguous-segment variant is also
   run, which additionally preserves the donor's local rate structure.
4. **Evaluation**: Sliding RP and Hill-Llobet at 2 and 3 ms on the merged train,
   plus the realized rate correlation between recipient and donor in 100 ms
   bins, so the result can be read against how correlated real neighbours are.

Usage
-----
    python semisynthetic.py --download --n-probes 8      # fetch and cache
    python semisynthetic.py --run                        # inject and evaluate
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import (BIN_SIZE, N_BINS, hill_llobet_from_acg,  # noqa: E402
                       slidingRP_from_acg)

CACHE = Path(r"D:/temp/slidingRP_resub/semisynth")
ONE_CACHE = Path(r"D:/temp/slidingRP_resub/ONE_semi")
F_GRID = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20])


def download(n_probes=8, seed=0):
    """Cache spike times, cluster ids and depths for a handful of BWM probes."""
    from one.api import ONE
    from brainwidemap import bwm_query
    from brainbox.io.one import SpikeSortingLoader

    CACHE.mkdir(parents=True, exist_ok=True)
    one = ONE(base_url="https://openalyx.internationalbrainlab.org",
              password="international", silent=True, cache_dir=str(ONE_CACHE))
    bwm = bwm_query(one)
    rng = np.random.default_rng(seed)
    pids = list(bwm.pid.unique())
    pick = [pids[i] for i in rng.choice(len(pids), n_probes * 2, replace=False)]

    saved = 0
    for pid in pick:
        if saved >= n_probes:
            break
        out = CACHE / f"{pid}.npz"
        if out.exists():
            saved += 1
            continue
        try:
            loader = SpikeSortingLoader(pid=pid, one=one)
            spikes, clusters, channels = loader.load_spike_sorting()
            clusters = loader.merge_clusters(spikes, clusters, channels)
            samp = np.asarray(one.load_dataset(loader.eid, "spikes.samples.npy",
                                               collection=loader.collection),
                              dtype=np.int64)
            np.savez_compressed(
                out,
                samples=samp, clusters=np.asarray(spikes.clusters, dtype=np.int32),
                rec_dur=float(np.max(spikes.times)),
                cluster_id=np.asarray(clusters["cluster_id"]),
                depths=np.asarray(clusters.get("depths", np.zeros(len(clusters["cluster_id"])))),
                acronym=np.asarray(clusters.get("acronym", np.array([""]))).astype(str),
                label=np.asarray(clusters.get("label", np.zeros(len(clusters["cluster_id"]))))),
            saved += 1
            print(f"cached {pid}: {samp.size:,} spikes", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"{pid} failed: {type(e).__name__}: {e}", flush=True)
    print(f"{saved} probes cached in {CACHE}", flush=True)


def _acg(samples):
    from acg_table import computeACG_samples
    return computeACG_samples(samples, N_BINS)


def _metrics(samples, rec_dur, cont_thresh=10.0, conf_thresh=90.0):
    a = _acg(samples)
    r = slidingRP_from_acg(a, samples.size, rec_dur, cont_thresh=cont_thresh,
                           conf_thresh=conf_thresh)
    hl2 = hill_llobet_from_acg(a, samples.size, rec_dur, 0.002, cont_thresh)
    hl3 = hill_llobet_from_acg(a, samples.size, rec_dur, 0.003, cont_thresh)
    return dict(max_conf=r["max_conf"], min_cont=r["min_cont"],
                tau_Cmin=r["rp_min_val"], n_viol_short=r["n_viol_short"],
                passes=r["passes"], tau_pass0=r["tau_pass0"],
                hl2_pass=hl2[0], hl2_est=hl2[1],
                hl3_pass=hl3[0], hl3_est=hl3[1])


def rate_correlation(a, b, rec_dur, bin_s=0.1, fs=30000.0):
    edges = np.arange(0, rec_dur + bin_s, bin_s)
    ca = np.histogram(a / fs, edges)[0].astype(float)
    cb = np.histogram(b / fs, edges)[0].astype(float)
    if ca.std() == 0 or cb.std() == 0:
        return np.nan
    return float(np.corrcoef(ca, cb)[0, 1])


def run(near_um=60.0, far_um=300.0, n_pairs_per_probe=25, seed=0,
        min_fr=2.0, strict_cont=5.0, strict_conf=99.0):
    rng = np.random.default_rng(seed)
    rows = []
    files = sorted(CACHE.glob("*.npz"))
    print(f"{len(files)} cached probes", flush=True)
    for fp in files:
        z = np.load(fp, allow_pickle=True)
        samp, clu = z["samples"], z["clusters"]
        rec_dur = float(z["rec_dur"])
        cid_all, depth_all = z["cluster_id"], z["depths"]
        acr_all = z["acronym"]
        order = np.argsort(clu, kind="stable")
        clu, samp = clu[order], samp[order]
        cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
        trains = {int(c): samp[a:a + n] for c, a, n in zip(cids, starts, counts)}
        depth = {int(c): float(depth_all[i]) if i < depth_all.size else np.nan
                 for i, c in enumerate(cid_all)}
        acr = {int(c): str(acr_all[i]) if i < acr_all.size else ""
               for i, c in enumerate(cid_all)}

        # candidate recipients: strictly clean and reasonably active
        recips = []
        for c, st in trains.items():
            fr = st.size / rec_dur
            if fr < min_fr or st.size < 1000:
                continue
            m = _metrics(st, rec_dur, cont_thresh=strict_cont, conf_thresh=strict_conf)
            if m["passes"]:
                recips.append((c, m["min_cont"]))
        if not recips:
            print(f"{fp.stem}: no strict-clean recipients", flush=True)
            continue

        rng.shuffle(recips)
        used = 0
        for c_r, base_cont in recips:
            if used >= n_pairs_per_probe:
                break
            st_r = trains[c_r]
            d_r = depth.get(c_r, np.nan)
            if not np.isfinite(d_r):
                continue
            near = [c for c in trains
                    if c != c_r and np.isfinite(depth.get(c, np.nan))
                    and abs(depth[c] - d_r) <= near_um and trains[c].size > 500]
            far = [c for c in trains
                   if c != c_r and np.isfinite(depth.get(c, np.nan))
                   and abs(depth[c] - d_r) >= far_um and trains[c].size > 500]
            for kind, pool in (("near", near), ("far", far)):
                if not pool:
                    continue
                c_d = int(pool[rng.integers(len(pool))])
                st_d = trains[c_d]
                corr = rate_correlation(st_r, st_d, rec_dur)
                for f in F_GRID:
                    for mode in ("random", "segment"):
                        if f == 0 and mode == "segment":
                            continue
                        n_inject = int(round(f * st_r.size / max(1 - f, 1e-9)))
                        n_inject = min(n_inject, st_d.size)
                        if n_inject == 0 and f > 0:
                            continue
                        if mode == "random":
                            take = rng.choice(st_d.size, n_inject, replace=False)
                            inj = st_d[np.sort(take)]
                        else:
                            if n_inject >= st_d.size:
                                inj = st_d
                            else:
                                s0 = int(rng.integers(0, st_d.size - n_inject))
                                inj = st_d[s0:s0 + n_inject]
                        merged = np.sort(np.concatenate([st_r, inj]))
                        m = _metrics(merged, rec_dur)
                        rows.append(dict(
                            pid=fp.stem, recipient=c_r, donor=c_d, donor_kind=kind,
                            mode=mode, f_injected=float(f),
                            n_recipient=int(st_r.size), n_injected=int(inj.size),
                            true_injected_frac=inj.size / merged.size,
                            recipient_base_cont=float(base_cont),
                            depth_sep_um=abs(depth[c_d] - d_r),
                            rate_corr=corr, rec_dur=rec_dur,
                            acronym_r=acr.get(c_r, ""), acronym_d=acr.get(c_d, ""),
                            fr_recipient=st_r.size / rec_dur,
                            fr_donor=st_d.size / rec_dur, **m))
                used += 1
        print(f"{fp.stem}: {used} recipient-donor sets, {len(rows)} rows so far",
              flush=True)

    df = pd.DataFrame(rows)
    out = Path(r"D:/temp/slidingRP_resub/sims/semisynthetic.pqt")
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    print(f"wrote {len(df)} rows to {out}", flush=True)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n-probes", type=int, default=8)
    ap.add_argument("--n-pairs", type=int, default=25)
    a = ap.parse_args()
    if a.download:
        download(a.n_probes)
    if a.run:
        t0 = time.time()
        run(n_pairs_per_probe=a.n_pairs)
        print(f"done in {(time.time()-t0)/60:.1f} min", flush=True)
