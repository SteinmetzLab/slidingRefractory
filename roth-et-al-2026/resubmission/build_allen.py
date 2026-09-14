"""Build the per-unit ACG table for the Allen Visual Coding Neuropixels dataset.

Streams each session NWB from the public Allen S3 bucket (no credentials),
computes every unit's 0-10 ms autocorrelogram at 1/30000 s, writes the ACGs and
per-unit metadata, then deletes the NWB (~2.5-3 GB each, 58 sessions).

Usage
-----
    python build_allen.py --limit 1     # smoke test
    python build_allen.py --shard 0/2   # parallel by hand
    python build_allen.py --local <nwb> # use an already-downloaded file

Outputs in ``D:/temp/slidingRP_resub/acg_tables/allen/``: ``<session>.npy`` /
``<session>.pqt`` as for the IBL builder.

Note: Allen NWBs store spike times as float seconds on the session master
clock; there is no integer sample column, so the ACG is computed from float
seconds (``acg_source = 'times'``). IBL and macaque use exact integer lags.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import N_BINS, BIN_SIZE  # noqa: E402

OUT_DIR = Path(r"D:/temp/slidingRP_resub/acg_tables/allen")
TMP_DIR = Path(r"D:/temp/slidingRP_resub/allen_tmp")
S3 = "https://allen-brain-observatory.s3-us-west-2.amazonaws.com/visual-coding-neuropixels/ecephys-cache"

UNIT_COLS = [
    "id", "cluster_id", "local_index", "peak_channel_id", "quality",
    "firing_rate", "isi_violations", "amplitude_cutoff", "presence_ratio",
    "snr", "isolation_distance", "d_prime", "nn_hit_rate", "nn_miss_rate",
    "l_ratio", "silhouette_score", "max_drift", "cumulative_drift",
    "waveform_duration", "waveform_halfwidth", "amplitude",
]


def download(url, dest, chunk=1 << 22):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for c in r.iter_content(chunk_size=chunk):
                f.write(c)
    return dest


def load_channels_table():
    """Cache-level channels.csv: the real CCF structure per channel.

    The NWB's own ``electrodes/location`` is coarse (most channels read
    'grey'); ``channels.csv`` from the ecephys cache carries
    ``ecephys_structure_acronym`` for all 123,224 channels.
    """
    path = TMP_DIR / "channels.csv"
    if not path.exists():
        download(f"{S3}/channels.csv", path)
    ch = pd.read_csv(path, usecols=["id", "ecephys_probe_id", "ecephys_structure_acronym",
                                    "anterior_posterior_ccf_coordinate",
                                    "dorsal_ventral_ccf_coordinate",
                                    "left_right_ccf_coordinate"])
    return ch.set_index("id")


def process_nwb(path, session_id, sessions_meta=None, channels_tbl=None):
    import h5py
    from slidingRP.metrics import computeACG

    with h5py.File(path, "r") as f:
        u = f["units"]
        st = u["spike_times"][:]
        idx = u["spike_times_index"][:]          # end offset of each unit
        starts = np.concatenate(([0], idx[:-1]))
        n_units = idx.size

        meta = {}
        for c in UNIT_COLS:
            if c in u:
                v = u[c][:]
                if v.dtype.kind == "S":
                    v = v.astype(str)
                if v.shape[0] == n_units:
                    meta[c] = v
        tbl = pd.DataFrame(meta)

        # electrode table -> region acronym per unit (via peak_channel_id)
        e = f["general/extracellular_ephys/electrodes"]
        eid_ = e["id"][:]
        loc = e["location"][:].astype(str)
        probe = e["probe_id"][:] if "probe_id" in e else np.zeros_like(eid_)
        order = np.argsort(eid_)
        pos = np.searchsorted(eid_[order], np.asarray(tbl["peak_channel_id"]))
        pos = np.clip(pos, 0, eid_.size - 1)
        sel = order[pos]
        ok = eid_[sel] == np.asarray(tbl["peak_channel_id"])
        tbl["nwb_location"] = np.where(ok, loc[sel], "")
        tbl["probe_id"] = np.where(ok, probe[sel], -1)

        # Real CCF structure from the cache-level channels table (the NWB's own
        # location field is coarse: most channels read 'grey').
        if channels_tbl is not None:
            j = channels_tbl.reindex(np.asarray(tbl["peak_channel_id"]))
            tbl["acronym"] = j["ecephys_structure_acronym"].fillna("").to_numpy()
            for src, dst in (("anterior_posterior_ccf_coordinate", "ccf_ap"),
                             ("dorsal_ventral_ccf_coordinate", "ccf_dv"),
                             ("left_right_ccf_coordinate", "ccf_lr")):
                tbl[dst] = j[src].to_numpy()
        else:
            tbl["acronym"] = tbl["nwb_location"]

        rec_dur = float(st.max())
        acgs = np.zeros((n_units, N_BINS), dtype=np.int32)
        n_spikes = np.zeros(n_units, dtype=np.int64)
        for i in range(n_units):
            s = st[starts[i]:idx[i]]
            n_spikes[i] = s.size
            if s.size > 1:
                acgs[i] = computeACG(s, BIN_SIZE, N_BINS)

    tbl["n_spikes"] = n_spikes
    tbl["rec_dur_s"] = rec_dur
    tbl["fr_recomputed"] = tbl.n_spikes / rec_dur
    tbl["session"] = session_id
    tbl["acg_source"] = "times"
    if sessions_meta is not None and session_id in sessions_meta.index:
        row = sessions_meta.loc[session_id]
        for c in ("specimen_id", "session_type", "sex", "genotype", "age_in_days"):
            if c in row:
                tbl[c] = row[c]
    return tbl, acgs


def remap_regions(tbl, br):
    """Allen acronym -> Cosmos / Beryl via iblatlas."""
    acr = np.asarray(tbl["acronym"], dtype=object)
    uniq = sorted({a for a in acr if a})
    cos, ber = {}, {}
    for a in uniq:
        try:
            aid = br.acronym2id(a)
            if len(aid) == 0:
                raise ValueError
            cos[a] = br.id2acronym(br.remap(aid[:1], source_map="Allen", target_map="Cosmos"))[0]
            ber[a] = br.id2acronym(br.remap(aid[:1], source_map="Allen", target_map="Beryl"))[0]
        except Exception:  # noqa: BLE001
            cos[a] = ber[a] = ""
    tbl["cosmos"] = [cos.get(a, "") for a in acr]
    tbl["beryl"] = [ber.get(a, "") for a in acr]
    return tbl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shard", type=str, default=None)
    ap.add_argument("--local", type=str, default=None)
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    from iblatlas.atlas import BrainRegions
    br = BrainRegions()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    sess_csv = TMP_DIR / "sessions.csv"
    if not sess_csv.exists():
        download(f"{S3}/sessions.csv", sess_csv)
    sessions = pd.read_csv(sess_csv).set_index("id")
    channels_tbl = load_channels_table()
    ids = list(sessions.index)
    if args.local:
        # session id comes from the filename, not the release order
        ids = [int(Path(args.local).stem.replace("session_", ""))]
    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        ids = ids[i::n]
    if args.limit:
        ids = ids[: args.limit]
    done = {int(p.stem) for p in OUT_DIR.glob("*.pqt")}
    todo = [s for s in ids if s not in done]
    print(f"{len(ids)} sessions assigned, {len(done)} done, {len(todo)} to do", flush=True)

    t_start = time.time()
    for k, sid in enumerate(todo):
        t0 = time.time()
        try:
            if args.local:
                nwb = Path(args.local)
            else:
                nwb = TMP_DIR / f"session_{sid}.nwb"
                if not nwb.exists():
                    download(f"{S3}/session_{sid}/session_{sid}.nwb", nwb)
            t_dl = time.time() - t0
            tbl, acgs = process_nwb(nwb, sid, sessions, channels_tbl)
            tbl = remap_regions(tbl, br)
            np.save(OUT_DIR / f"{sid}.npy", acgs)
            tbl.to_parquet(OUT_DIR / f"{sid}.pqt")
            if not args.keep and not args.local:
                nwb.unlink(missing_ok=True)
            el = time.time() - t_start
            print(f"[{k+1}/{len(todo)}] {sid}: {len(tbl)} units, dl {t_dl:.0f}s, "
                  f"total {time.time()-t0:.0f}s, eta {(len(todo)-k-1)*el/(k+1)/3600:.1f} h",
                  flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"[{k+1}/{len(todo)}] {sid} FAILED {type(e).__name__}: {e}", flush=True)
    print(f"done in {(time.time()-t_start)/3600:.2f} h", flush=True)


if __name__ == "__main__":
    main()
