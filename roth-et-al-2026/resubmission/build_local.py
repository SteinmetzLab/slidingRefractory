"""Build the per-unit ACG tables for the two local datasets.

    python build_local.py macaque      # D:/Horwitz  (7 recordings, 3 monkeys)
    python build_local.py steinmetz    # Steinmetz et al. 2019 (39 sessions, 10 mice)
    python build_local.py both

Macaque spike times are stored in integer samples at 30 kHz, so the ACG uses
exact integer lags (``acg_source = 'samples'``), as for IBL. Steinmetz 2019
provides float seconds only (``acg_source = 'times'``), as for Allen.

Macaque animal identity is inferred from the source folders recorded in
``script_fig1_hists.m`` and matched by unit count:
    LGN_1/2/3  = Osiris  (2024-04-24, 2024-05-07, 2024-05-14)
    V1_1, V1_2 = Q       (Q101525, Q122625)
    V1_3, V1_4 = W       (W040925, W123024)
This should be confirmed with the Horwitz lab.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import BIN_SIZE, N_BINS, computeACG_samples  # noqa: E402

OUT_DIR = Path(r"D:/temp/slidingRP_resub/acg_tables")
HORWITZ = Path(r"D:/Horwitz")
STEINMETZ = Path(r"D:/Dropbox/ucl/data/taskData")
FS = 30000.0

MACAQUE = {                       # name: (animal, region acronym, cosmos)
    "LGN_1": ("Osiris", "LGd", "TH"), "LGN_2": ("Osiris", "LGd", "TH"),
    "LGN_3": ("Osiris", "LGd", "TH"), "V1_1": ("Q", "VISp", "Isocortex"),
    "V1_2": ("Q", "VISp", "Isocortex"), "V1_3": ("W", "VISp", "Isocortex"),
    "V1_4": ("W", "VISp", "Isocortex"),
}


def build_macaque():
    out = OUT_DIR / "macaque"
    out.mkdir(parents=True, exist_ok=True)
    for name, (animal, acronym, cosmos) in MACAQUE.items():
        d = HORWITZ / name
        clu = np.load(d / "spike_clusters.npy").ravel()
        samp = np.load(d / "spike_times.npy").ravel().astype(np.int64)
        rec_dur = float(samp.max()) / FS
        order = np.argsort(clu, kind="stable")
        clu, samp = clu[order], samp[order]
        cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
        acgs = np.zeros((cids.size, N_BINS), dtype=np.int32)
        for i, (a, c) in enumerate(zip(starts, counts)):
            acgs[i] = computeACG_samples(samp[a:a + c], N_BINS)
        tbl = pd.DataFrame({
            "cluster_id": cids, "n_spikes": counts.astype(np.int64),
            "rec_dur_s": rec_dur, "fr_recomputed": counts / rec_dur,
            "dataset": "macaque", "species": "macaque", "animal": animal,
            "session": name, "insertion": name, "acronym": acronym,
            "cosmos": cosmos, "beryl": acronym, "acg_source": "samples",
        })
        np.save(out / f"{name}.npy", acgs)
        tbl.to_parquet(out / f"{name}.pqt")
        print(f"{name}: {len(tbl)} units, {rec_dur/60:.1f} min, animal {animal}",
              flush=True)


def _steinmetz_regions(session_dir, n_clusters):
    """Allen acronym per cluster from the peak channel's brain location."""
    peak = np.load(session_dir / "clusters.peakChannel.npy").ravel().astype(int)
    bl = pd.read_csv(session_dir / "channels.brainLocation.tsv", sep="\t")
    acr = bl["allen_ontology"].to_numpy().astype(str)
    idx = np.clip(peak - 1, 0, acr.size - 1)      # peakChannel is 1-based
    out = acr[idx]
    return out[:n_clusters] if out.size >= n_clusters else np.array([""] * n_clusters)


def build_steinmetz():
    out = OUT_DIR / "steinmetz"
    out.mkdir(parents=True, exist_ok=True)
    from slidingRP.metrics import computeACG
    from iblatlas.atlas import BrainRegions
    br = BrainRegions()
    sessions = sorted(p for p in STEINMETZ.iterdir() if p.is_dir())
    for d in sessions:
        try:
            st = np.load(d / "spikes.times.npy").ravel().astype(np.float64)
            clu = np.load(d / "spikes.clusters.npy").ravel().astype(np.int64)
            ann = np.load(d / "clusters._phy_annotation.npy").ravel()
        except FileNotFoundError as e:
            print(f"{d.name}: missing {e.filename}", flush=True)
            continue
        rec_dur = float(st.max())
        order = np.argsort(clu, kind="stable")
        clu, st = clu[order], st[order]
        cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
        acgs = np.zeros((cids.size, N_BINS), dtype=np.int32)
        for i, (a, c) in enumerate(zip(starts, counts)):
            if c > 1:
                acgs[i] = computeACG(st[a:a + c], BIN_SIZE, N_BINS)
        acr_all = _steinmetz_regions(d, ann.size)
        acr = np.array([acr_all[c] if c < acr_all.size else "" for c in cids])
        acr = np.array([a.split(" ")[0].strip() for a in acr])
        tbl = pd.DataFrame({
            "cluster_id": cids, "n_spikes": counts.astype(np.int64),
            "rec_dur_s": rec_dur, "fr_recomputed": counts / rec_dur,
            "dataset": "steinmetz", "species": "mouse",
            "animal": d.name.split("_")[0], "session": d.name, "insertion": d.name,
            "acronym": acr, "acg_source": "times",
            "phy_annotation": np.array([ann[c] if c < ann.size else -1 for c in cids]),
        })
        cos, ber = [], []
        for a in acr:
            try:
                aid = br.acronym2id(a)
                if len(aid) == 0:
                    raise ValueError
                cos.append(br.id2acronym(br.remap(aid[:1], source_map="Allen",
                                                  target_map="Cosmos"))[0])
                ber.append(br.id2acronym(br.remap(aid[:1], source_map="Allen",
                                                  target_map="Beryl"))[0])
            except Exception:  # noqa: BLE001
                cos.append("")
                ber.append("")
        tbl["cosmos"], tbl["beryl"] = cos, ber
        np.save(out / f"{d.name}.npy", acgs)
        tbl.to_parquet(out / f"{d.name}.pqt")
        print(f"{d.name}: {len(tbl)} units, {rec_dur/60:.1f} min, "
              f"cosmos {pd.Series(cos).value_counts().head(3).to_dict()}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["macaque", "steinmetz", "both"])
    a = ap.parse_args()
    if a.which in ("macaque", "both"):
        build_macaque()
    if a.which in ("steinmetz", "both"):
        build_steinmetz()
