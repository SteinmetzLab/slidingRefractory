"""Semi-synthetic contamination on real IBL recordings, from cleaned recipients (05).

Second version of semisynthetic.py, changed after review (2026-10-01):

1. **Recipients are made clean rather than selected as clean.** Every unit firing
   at least ``MIN_FR`` spikes/s is a candidate. Its own violations are deleted:
   walking through the train, a spike is dropped if it falls within
   ``CLEAN_SAMPLES`` (60 samples = 2.0 ms, inclusive) of the last kept spike, as
   Kilosort's duplicate removal does. Afterwards the autocorrelogram is empty at
   every lag Hill-Llobet at 2 ms counts (its inclusive bin rule sums lags up to
   exactly 60 samples), so both Sliding RP and Hill-Llobet should accept at zero
   injection, apart from Sliding RP's power limit at low rates.

   The old version selected recipients by Sliding RP passing at a strict
   setting, which was circular (it guaranteed Sliding RP's 100% at zero
   injection) and left real sub-2 ms violations that Hill-Llobet then counted.

   "Clean" here means no violations, not no contamination: a contaminating spike
   that happens not to fall within 2 ms of another spike survives the cleaning.
   It then behaves exactly like one of the recipient's own spikes, which is what
   the test assumes of the base train, so it does not bias the comparison. The
   deleted fraction is recorded; recipients losing more than ``MAX_DELETED``
   are excluded as unlikely to be single units.

2. **Donors must have enough spikes** to reach the top of the injection grid
   (the old version silently capped the injection at the donor's spike count, so
   a third of its "20%" trains were far below 20%).

3. **Many more pairs, several draws per pair.** Each recipient is paired with up
   to ``N_NEAR`` donors within ``NEAR_UM`` of its depth, ``N_MID`` between
   ``NEAR_UM`` and ``FAR_UM``, and ``N_FAR`` beyond ``FAR_UM``. For each pair,
   ``N_REPS`` random orderings of the donor's spikes are drawn and the first
   ``n`` spikes of an ordering are injected at each level, so the injected sets
   are nested within a draw and each draw traces a monotone-input curve.

4. **Correlation measured where it acts.** Besides the spike-count correlation
   (100 ms and 1 s bins), each pair's excess short-lag coincidence

       kappa(lo, hi) = CCG count at |lag| in (lo, hi] / count expected if
                       independent - 1

   is recorded for several windows. Violations from recipient-donor pairs at lags
   up to tau scale with 1 + kappa(0, tau), so kappa(0, 2 ms) is the quantity that
   moves Hill-Llobet at 2 ms, and kappa at the short lags Sliding RP can choose
   (down to 0.5 ms) is what moves Sliding RP. See 05_model_mismatch/correlated_rates.md.

Usage
-----
    D:/temp/slidingRP_resub/venv/Scripts/python semisynthetic_clean.py --download 24   # top the cache up to 24 probes
    python semisynthetic_clean.py --run                    # inject and evaluate
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
from acg_table import (N_BINS, computeACG_samples, hill_llobet_from_acg,  # noqa: E402
                       slidingRP_from_acg)

CACHE = Path(r"D:/temp/slidingRP_resub/semisynth")
ONE_CACHE = Path(r"D:/temp/slidingRP_resub/ONE_semi")
OUT = Path(r"D:/temp/slidingRP_resub/sims/semisynthetic_clean.pqt")
FS = 30000.0

F_GRID = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 25]) / 100
CLEAN_SAMPLES = 60          # 2.0 ms at 30 kHz, inclusive
MIN_FR = 2.0                # recipient firing rate (spikes/s)
MIN_SPIKES = 1000
MAX_DELETED = 0.10          # recipients losing more than this are excluded
NEAR_UM, FAR_UM = 60.0, 300.0
N_NEAR, N_MID, N_FAR = 3, 2, 2
N_REPS = 4
KAPPA_WINDOWS_MS = {"k_0_0p5": (0.0, 0.5), "k_0p5_1": (0.5, 1.0), "k_1_2": (1.0, 2.0),
                    "k_0_2": (0.0, 2.0), "k_0p5_2": (0.5, 2.0),
                    "k_3_30": (3.0, 30.0), "k_30_300": (30.0, 300.0)}


# --- download ----------------------------------------------------------------

def download(target, seed=1):
    """Top the cache up to ``target`` brain-wide-map insertions.

    Needs the ONE/brainwidemap environment (D:/temp/slidingRP_resub/venv). Only
    spike samples, cluster ids, depths, acronyms and labels are kept; the ONE
    files for each session are deleted once its npz is written, so the ONE cache
    does not grow (about 1 GB per insertion otherwise).
    """
    from one.api import ONE
    from brainwidemap import bwm_query
    from brainbox.io.one import SpikeSortingLoader

    CACHE.mkdir(parents=True, exist_ok=True)
    one = ONE(base_url="https://openalyx.internationalbrainlab.org",
              password="international", silent=True, cache_dir=str(ONE_CACHE))
    bwm = bwm_query(one)
    have = {p.stem for p in CACHE.glob("*.npz")}
    pids = [p for p in bwm.pid.unique() if p not in have]
    rng = np.random.default_rng(seed)
    rng.shuffle(pids)
    for pid in pids:
        if len(have) >= target:
            break
        try:
            loader = SpikeSortingLoader(pid=pid, one=one)
            spikes, clusters, channels = loader.load_spike_sorting()
            clusters = loader.merge_clusters(spikes, clusters, channels)
            samp = np.asarray(one.load_dataset(loader.eid, "spikes.samples.npy",
                                               collection=loader.collection), dtype=np.int64)
            np.savez_compressed(
                CACHE / f"{pid}.npz",
                samples=samp, clusters=np.asarray(spikes.clusters, dtype=np.int32),
                rec_dur=float(np.max(spikes.times)),
                cluster_id=np.asarray(clusters["cluster_id"]),
                depths=np.asarray(clusters.get("depths", np.zeros(len(clusters["cluster_id"])))),
                acronym=np.asarray(clusters.get("acronym", np.array([""]))).astype(str),
                label=np.asarray(clusters.get("label", np.zeros(len(clusters["cluster_id"])))))
            have.add(pid)
            print(f"cached {pid}: {samp.size:,} spikes ({len(have)}/{target})", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"{pid} failed: {type(e).__name__}: {e}", flush=True)
        finally:
            try:   # remove this session's ONE files (only inside our own cache dir)
                sess = Path(one.eid2path(loader.eid))
                if ONE_CACHE in sess.parents and sess.exists():
                    shutil.rmtree(sess)
            except Exception:  # noqa: BLE001
                pass
    print(f"{len(have)} probes cached in {CACHE}", flush=True)


# --- helpers -----------------------------------------------------------------

def clean_train(s, dead=CLEAN_SAMPLES):
    """Drop every spike within ``dead`` samples (inclusive) of the last kept one."""
    s = np.sort(np.asarray(s, dtype=np.int64))
    keep = np.zeros(s.size, dtype=bool)
    last = None
    for i, t in enumerate(s.tolist()):
        if last is None or t - last > dead:
            keep[i] = True
            last = t
    return s[keep]


def kappa(a, b, lo_ms, hi_ms, dur_s):
    """Excess coincidence of sorted trains a, b (samples) at |lag| in (lo, hi]."""
    lo, hi = int(round(lo_ms * FS / 1000)), int(round(hi_ms * FS / 1000))
    n = (np.searchsorted(b, a + hi, side="right") - np.searchsorted(b, a + lo, side="right")
         + np.searchsorted(b, a - lo, side="left") - np.searchsorted(b, a - hi, side="left"))
    obs = float(n.sum())
    expected = a.size * b.size * 2 * (hi - lo) / (dur_s * FS)
    return obs / expected - 1.0 if expected > 0 else np.nan


def kappas(a, b, dur_s):
    """All KAPPA_WINDOWS_MS for sorted trains a, b, from cumulative CCG counts
    (pairs with 0 < |lag| <= h), which kappa() would give window by window."""
    edges = sorted({e for w in KAPPA_WINDOWS_MS.values() for e in w})
    cum = {}
    for e in edges:
        h = int(round(e * FS / 1000))
        if h == 0:
            cum[e] = 0.0
            continue
        n = (np.searchsorted(b, a + h, side="right") - np.searchsorted(b, a - h, side="left")).sum()
        n0 = (np.searchsorted(b, a, side="right") - np.searchsorted(b, a, side="left")).sum()
        cum[e] = float(n - n0)
    out = {}
    for k, (lo, hi) in KAPPA_WINDOWS_MS.items():
        nlo, nhi = int(round(lo * FS / 1000)), int(round(hi * FS / 1000))
        expected = a.size * b.size * 2 * (nhi - nlo) / (dur_s * FS)
        out[k] = (cum[hi] - cum[lo]) / expected - 1.0 if expected > 0 else np.nan
    return out


def kappa_profile(a, b, dur_s, step=3, top=90):
    """kappa in consecutive lag bins of ``step`` samples (0.1 ms) up to ``top``
    (3 ms): bin k covers |lag| in (k*step, (k+1)*step]."""
    _, _, lag = _pairs_within(a, b, top + 1)
    n = np.bincount((lag - 1) // step, minlength=top // step)[:top // step].astype(float)
    expected = a.size * b.size * 2 * step / (dur_s * FS)
    return (n / expected - 1.0).tolist() if expected > 0 else [np.nan] * (top // step)


def count_corr(a, b, dur_s, bin_s):
    edges = np.arange(0, dur_s + bin_s, bin_s)
    ca = np.histogram(a / FS, edges)[0].astype(float)
    cb = np.histogram(b / FS, edges)[0].astype(float)
    if ca.std() == 0 or cb.std() == 0:
        return np.nan
    return float(np.corrcoef(ca, cb)[0, 1])


def _pairs_within(x, y, n_bins=N_BINS):
    """All (i, j, lag) with 0 < |y[j] - x[i]| < n_bins, x and y sorted (samples)."""
    lo = np.searchsorted(x, y - (n_bins - 1), side="left")
    hi = np.searchsorted(x, y + (n_bins - 1), side="right")
    cnt = hi - lo
    j = np.repeat(np.arange(y.size), cnt)
    start = np.repeat(lo - np.r_[0, np.cumsum(cnt)[:-1]], cnt)
    i = start + np.arange(cnt.sum())
    lag = np.abs(y[j] - x[i])
    keep = (lag > 0) & (lag < n_bins)
    return i[keep], j[keep], lag[keep]


def _self_pairs(y, n_bins=N_BINS):
    """All (a, b, lag) with a < b and 0 < y[b] - y[a] < n_bins, y sorted."""
    A, B, L = [], [], []
    for shift in range(1, y.size):
        d = y[shift:] - y[:-shift]
        m = d < n_bins
        if not np.any(m):
            break
        m &= d > 0
        idx = np.nonzero(m)[0]
        A.append(idx); B.append(idx + shift); L.append(d[m])
    if not A:
        return (np.zeros(0, int),) * 3
    return np.concatenate(A), np.concatenate(B), np.concatenate(L)


def nested_acgs(recipient, donor, perm, n_list, n_bins=N_BINS):
    """ACGs of recipient + donor[perm[:n]] for every n in n_list, exactly.

    ACG(merged) = ACG(recipient) + CCG(recipient, injected) + ACG(injected), and
    donor spike j belongs to the injected set at level n iff its position in the
    permutation is below n, so every level comes from one pass over the pairs.
    Lag-0 coincidences are excluded, as in computeACG_samples.
    """
    acg_r = computeACG_samples(recipient, n_bins)
    rank = np.empty(donor.size, dtype=np.int64)
    rank[perm] = np.arange(donor.size)
    _, jc, lc = _pairs_within(recipient, donor, n_bins)
    a, b, ld = _self_pairs(donor, n_bins)
    rc, rd = rank[jc], np.maximum(rank[a], rank[b])
    return [acg_r + np.bincount(lc[rc < n], minlength=n_bins)[:n_bins]
            + np.bincount(ld[rd < n], minlength=n_bins)[:n_bins] for n in n_list]


def evaluate(train, dur, acg=None, n=None):
    a = computeACG_samples(train, N_BINS) if acg is None else acg
    n = train.size if n is None else n
    r = slidingRP_from_acg(a, n, dur)
    hl = hill_llobet_from_acg(a, n, dur, 0.002)
    return dict(passes=r["passes"], max_conf=r["max_conf"], min_cont=r["min_cont"],
                tau_Cmin=r["rp_min_val"], hl2_pass=hl[0], hl2_est=hl[1])


# --- the run -----------------------------------------------------------------

def _probe(fp, seed):
    rng = np.random.default_rng(seed)
    z = np.load(fp, allow_pickle=True)
    samp, clu = z["samples"], z["clusters"]
    dur = float(z["rec_dur"])
    order = np.argsort(clu, kind="stable")
    clu, samp = clu[order], samp[order]
    cids, starts, counts = np.unique(clu, return_index=True, return_counts=True)
    trains = {int(c): np.sort(samp[a:a + n]) for c, a, n in zip(cids, starts, counts)}
    cid_all = [int(c) for c in z["cluster_id"]]
    depth = dict(zip(cid_all, map(float, z["depths"])))
    acr = dict(zip(cid_all, map(str, z["acronym"])))
    label = dict(zip(cid_all, map(float, z["label"])))
    f_top = F_GRID.max()

    rows, recips = [], []
    for c, st in trains.items():
        if st.size < MIN_SPIKES or st.size / dur < MIN_FR or not np.isfinite(depth.get(c, np.nan)):
            continue
        cl = clean_train(st)
        recips.append(dict(pid=fp.stem, recipient=c, n_raw=st.size, n_clean=cl.size,
                           frac_deleted=1 - cl.size / st.size))
        if 1 - cl.size / st.size > MAX_DELETED:
            continue
        base = evaluate(cl, dur)
        need = int(np.ceil(f_top * cl.size / (1 - f_top)))
        donors = [d for d in trains if d != c and trains[d].size >= need
                  and np.isfinite(depth.get(d, np.nan))]
        sep = {d: abs(depth[d] - depth[c]) for d in donors}
        picks = []
        for kind, ok, k in (("near", lambda s: s <= NEAR_UM, N_NEAR),
                            ("mid", lambda s: NEAR_UM < s < FAR_UM, N_MID),
                            ("far", lambda s: s >= FAR_UM, N_FAR)):
            pool = [d for d in donors if ok(sep[d])]
            if pool:
                picks += [(kind, int(d)) for d in rng.choice(pool, min(k, len(pool)), replace=False)]
        for kind, d in picks:
            sd = trains[d]
            pair = dict(pid=fp.stem, recipient=c, donor=d, donor_kind=kind, depth_sep_um=sep[d],
                        n_recipient=cl.size, frac_deleted=1 - cl.size / st.size,
                        fr_recipient=cl.size / dur, fr_donor=sd.size / dur, rec_dur=dur,
                        acronym_r=acr.get(c, ""), acronym_d=acr.get(d, ""),
                        label_r=label.get(c, np.nan), label_d=label.get(d, np.nan),
                        r_100ms=count_corr(cl, sd, dur, 0.1), r_1s=count_corr(cl, sd, dur, 1.0))
            pair.update(kappas(cl, sd, dur))
            pair["kappa_profile"] = kappa_profile(cl, sd, dur)
            n_list = [int(round(f * cl.size / (1 - f))) for f in F_GRID]
            for rep in range(N_REPS):
                perm = rng.permutation(sd.size)
                acgs = nested_acgs(cl, sd, perm, n_list)
                for f, n_inj, acg in zip(F_GRID, n_list, acgs):
                    m = base if n_inj == 0 else evaluate(None, dur, acg, cl.size + n_inj)
                    rows.append(dict(pair, rep=rep, f_injected=float(f), n_injected=n_inj,
                                     true_injected_frac=n_inj / (cl.size + n_inj), **m))
    return rows, recips


def run(n_jobs=4, seed=0):
    from joblib import Parallel, delayed
    files = sorted(CACHE.glob("*.npz"))
    print(f"{len(files)} cached probes", flush=True)
    res = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_probe)(fp, seed + 7919 * i) for i, fp in enumerate(files))
    df = pd.DataFrame([r for rows, _ in res for r in rows])
    rc = pd.DataFrame([r for _, rec in res for r in rec])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT)
    rc.to_parquet(OUT.with_name("semisynthetic_clean_recipients.pqt"))
    print(f"wrote {len(df):,} rows ({df.groupby(['pid','recipient','donor']).ngroups} pairs, "
          f"{df.recipient.nunique()} recipients) to {OUT}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", type=int, default=0, help="top the cache up to N probes")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    if a.download:
        download(a.download)
    if a.run:
        t0 = time.time()
        run(a.jobs)
        print(f"done in {(time.time() - t0) / 60:.1f} min", flush=True)
