"""Graded refractory recovery, done properly (work package 05, case C1).

The first version of ``gen_graded_rp`` used a logistic hazard centred on the
refractory period, which leaves a substantial firing probability at every lag
(27% of baseline at zero lag for a 2 ms width) and a realised minimum ISI of
zero. The metric rejected those units, correctly, because they genuinely
violate at every lag; the apparent "collapse" was an artifact of the generator.

The intended model is an absolute refractory period followed by a relative one:
the hazard is exactly zero below ``rp`` and then rises linearly to its
asymptote over ``width``. Nick's prediction for that model is that acceptance
should if anything go *up* relative to a hard refractory period, because the
bins between ``rp`` and ``rp + width`` carry fewer violations than baseline and
so give the sliding search extra windows that might work, with each successive
bin less likely to help than the last.

This script tests exactly that: the hard-RP reference (width = 0) against three
ramp widths, at two firing rates.

Run:  python run_graded.py [--n-sim 800] [--n-jobs 2]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_sims as RS  # noqa: E402

CONT = [0, 0.04, 0.06, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14, 0.16, 0.20]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sim", type=int, default=800)
    ap.add_argument("--n-jobs", type=int, default=2)
    a = ap.parse_args()

    conds = RS.grid(model=["standard"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                    rec_dur=[7200.0], cont_prop=CONT)
    conds += RS.grid(model=["graded"], total_rate=[1.0, 5.0], rp_dur=[0.002],
                     rec_dur=[7200.0], cont_prop=CONT,
                     model_kw=[{"width": w} for w in (0.0005, 0.001, 0.002)])
    RS._sweep_impl("graded_fixed", conds, a.n_sim, dict(oracle=False),
                   n_jobs=a.n_jobs)


if __name__ == "__main__":
    main()
