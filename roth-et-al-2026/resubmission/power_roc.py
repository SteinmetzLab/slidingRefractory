"""Does sliding buy power, or only move the operating point? (04)

At a fixed confidence threshold, Sliding RP accepts more clean units than the
Poisson test at the true refractory period, but it also accepts more
contaminated ones: the maximum over many windows is a multiple comparison.
Comparing the two at one threshold therefore mixes power with size. The fair
comparison is at matched false acceptance: trace each method's ROC curve (true
acceptance at 8% contamination against false acceptance at 12%) over every
confidence threshold, and compare the curves.

Theory says the oracle cannot be beaten under the simulated model: for lags up
to the true RP every violation comes from contamination, the cumulative count up
to the true RP is a sufficient statistic for the contamination rate, and the
one-sided Poisson test on it is uniformly most powerful (monotone likelihood
ratio). Sliding adds only shorter windows (no new information) and longer ones
(which include the neuron's own spikes). See sliding_power.md.

Arms, all on the same trains:
  sliding        max confidence over tau_r > 0.5 ms (the published statistic;
                 its FWER-corrected version is a monotone relabelling of the
                 same ordering, so it has the same ROC curve)
  oracle         Poisson test at the true RP
  fixed 0.5/1/2/3 ms   Poisson test at one window chosen in advance

Per-train confidences are kept, so any threshold can be applied afterwards.

Run:  python power_roc.py        (~5 min with 4 workers)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).parent))
from acg_table import BIN_SIZE, N_BINS, slidingRP_from_acg  # noqa: E402
from resub_simulations import make_train, poisson_test_fixed_tau  # noqa: E402

from slidingRP.metrics import computeACG  # noqa: E402

OUT = Path(r"D:/temp/slidingRP_resub/sims/power_roc.pqt")
RATES = (0.5, 1.0, 2.0, 5.0)
RPS = (0.0015, 0.003, 0.005)
CONTS = (0.08, 0.10, 0.12)
FIXED = (0.0005, 0.001, 0.002, 0.003)
DUR = 7200.0           # as in the decomposition sweep
N_SIM = 2000


def one(rate, rp, cont, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(N_SIM):
        st, _ = make_train("standard", rate, cont, DUR, rp, rng=rng)
        n = st.size
        acg = computeACG(st, BIN_SIZE, N_BINS)
        row = dict(rate=rate, rp=rp, cont=cont, n=n,
                   sliding=slidingRP_from_acg(acg, n, DUR)["max_conf"],
                   oracle=poisson_test_fixed_tau(acg, n, DUR, rp))
        for t in FIXED:
            row[f"fixed_{t*1000:g}ms"] = poisson_test_fixed_tau(acg, n, DUR, t)
        rows.append(row)
    return rows


def main():
    jobs = [(r, rp, c, 1_000_003 * i + 17) for i, (r, rp, c) in
            enumerate((r, rp, c) for r in RATES for rp in RPS for c in CONTS)]
    res = Parallel(n_jobs=4, verbose=5)(delayed(one)(*j) for j in jobs)
    df = pd.DataFrame([row for rows in res for row in rows])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT)
    print(f"wrote {OUT}: {len(df):,} trains")


if __name__ == "__main__":
    main()
