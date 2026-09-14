"""Build the per-unit ACG table for the IBL brain-wide map (public openalyx).

For every probe insertion in the 2023_12 BWM release (699 pids, 459 sessions,
139 mice), download the spike sorting, compute each unit's 0-10 ms
autocorrelogram at 1/30000 s from the integer ``spikes.samples``, store the
ACGs plus per-unit metadata, then delete the downloaded spike files.

Usage
-----
    python build_ibl.py                 # all pids, resumable
    python build_ibl.py --limit 3       # first 3 pids (smoke test)
    python build_ibl.py --workers 3     # parallel over pids
    python build_ibl.py --shard 0/3     # process pids 0, 3, 6, ... (for
                                        # running several processes by hand)

Outputs, per pid, in ``D:/temp/slidingRP_resub/acg_tables/ibl/``:
    <pid>.npy   int32 [n_units x 300] ACG counts
    <pid>.pqt   per-unit metadata (one row per unit, same order as the npy)

Notes
-----
* The ACG uses ``spikes.samples`` (integer sample index on the probe clock),
  not ``spikes.times`` (sync-corrected float seconds). These give different
  bin assignments for ~25% of spike pairs because of the seconds<->samples
  float round-trip; integer lags are exact. ``--acg-source times`` recomputes
  from float seconds for the comparison reported in the Methods.
* IBL's own pipeline metrics (``slidingRP_viol``, ``max_confidence``,
  ``min_contamination``, ``n_spikes_below2``) are carried through unchanged so
  the recomputation can be compared against them.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import N_BINS, BIN_SIZE, computeACG_samples  # noqa: E402

OUT_DIR = Path(r"D:/temp/slidingRP_resub/acg_tables/ibl")
ONE_CACHE = Path(r"D:/temp/slidingRP_resub/ONE")

# per-unit columns carried over from the IBL clusters table
CLUSTER_COLS = [
    "cluster_id", "acronym", "atlas_id", "label", "slidingRP_viol",
    "slidingRP_viol_forced", "max_confidence", "min_contamination",
    "n_spikes_below2", "firing_rate", "spike_count", "amp_median",
    "noise_cutoff", "presence_ratio", "contamination", "ks2_label",
    "x", "y", "z", "depths",
]


def group_by_cluster(clu):
    """Counting-sort order of spikes by cluster id; returns (order, cids, starts, ends).

    Faster than argsort for the 20M-element integer arrays here.
    """
    clu = np.asarray(clu)
    cids, inv, counts = np.unique(clu, return_inverse=True, return_counts=True)
    starts = np.zeros(cids.size, dtype=np.int64)
    np.cumsum(counts[:-1], out=starts[1:])
    order = np.argsort(inv, kind="stable")
    return order, cids, starts, starts + counts


def process_pid(pid, one, br, acg_source="samples", keep_files=False):
    from brainbox.io.one import SpikeSortingLoader
    from slidingRP.metrics import computeACG

    t0 = time.time()
    loader = SpikeSortingLoader(pid=pid, one=one)
    spikes, clusters, channels = loader.load_spike_sorting()
    clusters = loader.merge_clusters(spikes, clusters, channels)
    if spikes is None or "clusters" not in spikes or spikes.clusters is None:
        raise RuntimeError("no spikes")

    eid, pname = loader.eid, loader.pname
    rec_dur = float(np.max(spikes.times))

    if acg_source == "samples":
        lags = one.load_dataset(eid, "spikes.samples.npy", collection=loader.collection)
        lags = np.asarray(lags, dtype=np.int64)
    else:
        lags = np.asarray(spikes.times, dtype=np.float64)

    order, cids, starts, ends = group_by_cluster(spikes.clusters)
    lags = lags[order]

    acgs = np.zeros((cids.size, N_BINS), dtype=np.int32)
    for i, (a, b) in enumerate(zip(starts, ends)):
        if acg_source == "samples":
            acgs[i] = computeACG_samples(lags[a:b], N_BINS)
        else:
            acgs[i] = computeACG(lags[a:b], BIN_SIZE, N_BINS)

    meta = {}
    for c in CLUSTER_COLS:
        if c in clusters:
            v = np.asarray(clusters[c])
            meta[c] = v[np.searchsorted(np.asarray(clusters["cluster_id"]), cids)] \
                if c != "cluster_id" else cids
    tbl = pd.DataFrame(meta)
    tbl["n_spikes"] = (ends - starts).astype(np.int64)
    tbl["rec_dur_s"] = rec_dur
    tbl["fr_recomputed"] = tbl.n_spikes / rec_dur
    tbl["pid"] = pid
    tbl["eid"] = eid
    tbl["probe_name"] = pname
    tbl["collection"] = loader.collection
    tbl["acg_source"] = acg_source

    # Cosmos / Beryl remap from the Allen atlas id
    aid = np.asarray(tbl.get("atlas_id", np.full(len(tbl), np.nan)), dtype=float)
    ok = ~np.isnan(aid)
    for mapping in ("Cosmos", "Beryl"):
        out = np.array([""] * len(tbl), dtype=object)
        if ok.any():
            ids = br.remap(aid[ok].astype(int), source_map="Allen", target_map=mapping)
            out[ok] = br.id2acronym(ids)
        tbl[mapping.lower()] = out

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.save(OUT_DIR / f"{pid}.npy", acgs)
    tbl.to_parquet(OUT_DIR / f"{pid}.pqt")

    if not keep_files:
        ses_path = one.eid2path(eid)
        for sub in (Path(ses_path) / "alf",):
            for f in sub.rglob("spikes.*.npy"):
                try:
                    f.unlink()
                except OSError:
                    pass
    return len(tbl), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shard", type=str, default=None, help="i/n: take every n-th pid starting at i")
    ap.add_argument("--acg-source", default="samples", choices=["samples", "times"])
    ap.add_argument("--keep-files", action="store_true")
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--cache-dir", default=str(ONE_CACHE),
                    help="ONE cache; use a separate one per parallel shard")
    args = ap.parse_args()

    from one.api import ONE
    from brainwidemap import bwm_query
    from iblatlas.atlas import BrainRegions

    one = ONE(base_url="https://openalyx.internationalbrainlab.org",
              password="international", silent=True, cache_dir=args.cache_dir)
    br = BrainRegions()
    pids = list(bwm_query(one).pid.unique())
    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        pids = pids[i::n]
    if args.limit:
        pids = pids[: args.limit]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in OUT_DIR.glob("*.pqt")} if not args.redo else set()
    todo = [p for p in pids if p not in done]
    print(f"{len(pids)} pids assigned, {len(done)} already done, {len(todo)} to do", flush=True)

    n_units = 0
    t_start = time.time()
    for k, pid in enumerate(todo):
        try:
            nu, dt = process_pid(pid, one, br, args.acg_source, args.keep_files)
            n_units += nu
            el = time.time() - t_start
            rate = (k + 1) / el * 3600
            print(f"[{k+1}/{len(todo)}] {pid} {nu} units {dt:.0f}s "
                  f"({rate:.0f} pid/h, eta {(len(todo)-k-1)/max(rate,1e-9):.1f} h, "
                  f"{n_units:,} units so far)", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"[{k+1}/{len(todo)}] {pid} FAILED {type(e).__name__}: {e}", flush=True)
    print(f"done: {n_units:,} units in {(time.time()-t_start)/3600:.2f} h", flush=True)


if __name__ == "__main__":
    main()
