"""Cache each semi-synthetic insertion's cluster template waveforms (05).

Used by plot_dip_pairs.py to tell a neuron split across two clusters (nearly
identical templates) from two different neurons whose overlapping spikes the
sorter failed to resolve (different templates on shared channels).

Saves <pid>_waveforms.npz next to the spike cache: clusters.waveforms
(n_clusters x samples x channels) and clusters.waveformsChannels.

Run:  D:/temp/slidingRP_resub/venv/Scripts/python fetch_waveforms.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

CACHE = Path(r"D:/temp/slidingRP_resub/semisynth")
ONE_CACHE = Path(r"D:/temp/slidingRP_resub/ONE_semi")


def main():
    from one.api import ONE
    from brainbox.io.one import SpikeSortingLoader
    one = ONE(base_url="https://openalyx.internationalbrainlab.org",
              password="international", silent=True, cache_dir=str(ONE_CACHE))
    for fp in sorted(CACHE.glob("*.npz")):
        if fp.stem.endswith("_waveforms"):
            continue
        out = CACHE / f"{fp.stem}_waveforms.npz"
        if out.exists():
            continue
        try:
            loader = SpikeSortingLoader(pid=fp.stem, one=one)
            # the spike cache used the most recent revision (no default revision is
            # set for these datasets), so take the most recent one here too, without
            # downloading the spikes again
            coll = f"alf/{loader.pname}/pykilosort"   # loader.collection is set only on load
            ds = one.list_datasets(loader.eid, collection=coll,
                                   filename="clusters.waveforms.npy")
            revs = sorted(part.strip("#") for d in ds for part in d.split("/")
                          if part.startswith("#") and part.endswith("#"))
            rev = revs[-1] if revs else None
            wf = one.load_dataset(loader.eid, "clusters.waveforms.npy",
                                  collection=coll, revision=rev)
            ch = one.load_dataset(loader.eid, "clusters.waveformsChannels.npy",
                                  collection=coll, revision=rev)
            np.savez_compressed(out, waveforms=np.asarray(wf, dtype=np.float32),
                                channels=np.asarray(ch))
            print(f"{fp.stem}: waveforms {np.shape(wf)}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"{fp.stem} failed: {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
