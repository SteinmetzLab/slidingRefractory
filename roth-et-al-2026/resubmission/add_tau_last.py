"""Add ``tau_last_pass`` to the enriched tables.

The enriched tables already carry two timepoints that the Sliding RP metric
produces on its own:

    tau_Cmin        the tau_r at which the minimum confirmable contamination
                    C_min is attained (the package's ``rp_min_val``)
    tau_first_pass  the shortest tau_r above tau_min at which the confidence at
                    C_thresh reaches gamma_thresh

Nick asked to also see the *last* such tau_r: the longest window over which the
unit still looks clean at the 10% / 90% operating point. As tau_r grows past a
unit's true refractory period real spikes start landing in the window, observed
violations climb faster than the linearly growing expectation, and the
confidence falls, so the last accepted tau_r is the natural Sliding-RP-native
answer to "how long does this unit stay quiet for".

Two properties to keep in mind when reading it:

  * it is right-censored at the edge of the tested window (10 ms), so a unit
    that is quiet throughout reports 9.98 ms rather than its true value;
  * the accepted set of tau_r need not be contiguous, and this takes the last
    element of it, not the end of the first run.

This is a read-only pass over the stored ACGs (about a minute for all four
datasets); nothing about the metric's behaviour changes.

Run:  python add_tau_last.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import BIN_SIZE, REF_DUR, RP_CENTERS, RP_REJECT  # noqa: E402

ACG_DIR = Path(r"D:/temp/slidingRP_resub/acg_tables")
ENRICHED = Path(r"D:/temp/slidingRP_resub/enriched")
CONT_THRESH = 10.0
CONF_THRESH = 90.0
CHUNK = 4000                      # units per poisson.cdf call


def tau_last_pass(acg, n_spikes, rec_dur, cont_thresh=CONT_THRESH,
                  conf_thresh=CONF_THRESH, rp_reject=RP_REJECT):
    """Vectorised over units: last accepted tau_r in seconds, NaN if none.

    Also returns the recomputed first accepted tau_r and acceptance flag, which
    are checked against the stored columns by the caller.
    """
    acg = np.asarray(acg, dtype=np.float64)
    n = np.asarray(n_spikes, dtype=np.float64)[:, None]
    d = np.asarray(rec_dur, dtype=np.float64)[:, None]
    c = cont_thresh / 100.0
    nc, nb = n * c, n * (1 - c)
    expected = 2 * REF_DUR[None, :] / d * nc * (nb + (nc - 1) / 2)
    obs = np.cumsum(acg, axis=1)
    conf = 1 - stats.poisson.cdf(obs, expected)

    ok = (RP_CENTERS > rp_reject)[None, :] & (conf >= conf_thresh / 100.0)
    any_ok = ok.any(axis=1)
    first = np.where(any_ok, RP_CENTERS[ok.argmax(axis=1)], np.nan)
    last_i = ok.shape[1] - 1 - ok[:, ::-1].argmax(axis=1)
    last = np.where(any_ok, RP_CENTERS[last_i], np.nan)
    n_accepted_bins = ok.sum(axis=1)
    return last, first, any_ok, n_accepted_bins


def process(dataset, verbose=True):
    files = sorted((ENRICHED / dataset).glob("*.pqt"))
    n_units = n_mismatch = 0
    t0 = time.time()
    for k, f in enumerate(files):
        npy = ACG_DIR / dataset / (f.stem + ".npy")
        if not npy.exists():
            print(f"  missing ACG for {f.stem}", flush=True)
            continue
        e = pd.read_parquet(f)
        acg = np.load(npy)
        if len(e) != acg.shape[0]:
            print(f"  row mismatch {f.stem}: {len(e)} vs {acg.shape[0]}", flush=True)
            continue

        last = np.full(len(e), np.nan)
        first = np.full(len(e), np.nan)
        acc = np.zeros(len(e), bool)
        nbins = np.zeros(len(e), int)
        ns = e.n_spikes.to_numpy(float)
        rd = e.rec_dur_s.to_numpy(float)
        for i in range(0, len(e), CHUNK):
            sl = slice(i, min(i + CHUNK, len(e)))
            last[sl], first[sl], acc[sl], nbins[sl] = tau_last_pass(
                acg[sl], ns[sl], rd[sl])

        # the recomputation must reproduce what is already stored
        stored_first = e.tau_first_pass.to_numpy(float)
        bad = ~(np.isclose(first, stored_first, equal_nan=True)
                & (acc == e.passes.to_numpy(bool)))
        n_mismatch += int(bad.sum())
        n_units += len(e)

        e["tau_last_pass"] = last
        e["n_accepted_bins"] = nbins
        e.to_parquet(f)
        if verbose and (k % 100 == 0 or k == len(files) - 1):
            print(f"  {dataset} {k + 1}/{len(files)} shards, {n_units:,} units, "
                  f"{time.time() - t0:.0f} s", flush=True)
    return n_units, n_mismatch


def main():
    total = bad = 0
    for ds in ("ibl", "allen", "steinmetz", "macaque"):
        if not (ENRICHED / ds).exists():
            continue
        n, m = process(ds)
        print(f"{ds}: {n:,} units, {m} disagreements with the stored columns",
              flush=True)
        total += n
        bad += m
    print(f"\ntotal {total:,} units, {bad} disagreements "
          f"(must be 0: the recomputation reproduces tau_first_pass and passes)")


if __name__ == "__main__":
    main()
