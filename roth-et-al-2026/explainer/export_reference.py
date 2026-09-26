"""Write reference answers from the Python package for test_core.js.

For a set of spike trains (the Fig 2 neuron plus simulated units covering low
and high rates, short and long recordings, clean and heavily contaminated), save
the spike times as raw float64 and the package's answers as JSON:

  * the ACG (package computeACG) -- must match exactly;
  * slidingRP at the defaults (package) and at non-default thresholds
    (slidingRP_from_acg, itself tested equal to the package);
  * the full confidence matrix (package computeMatrix);
  * first and last accepted tau_r (add_tau_last);
  * Hill-Llobet at several fixed RPs;
  * a grid of Poisson tail probabilities and the lambda-for-confidence inverse;
  * acceptance proportions over 1000 simulated trains, which checks the
    JavaScript generator end to end against the Python one.

Run:  python export_reference.py      (then: node test_core.js)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "resubmission"))
from acg_table import hill_llobet_from_acg, slidingRP_from_acg  # noqa: E402
from add_tau_last import tau_last_pass  # noqa: E402
from resub_simulations import make_train  # noqa: E402

from slidingRP.metrics import computeACG, computeMatrix, slidingRP  # noqa: E402

OUT = Path(r"D:/temp/slidingRP_resub/explainer_ref")
BIN = 1 / 30000
NB = 300


def case_record(name, st, dur):
    st = np.sort(np.asarray(st, np.float64))
    st.tofile(OUT / f"{name}.f64")
    acg = computeACG(st, BIN, NB)
    n = st.size
    pk = slidingRP(st, params={"recDur": dur})
    alt = slidingRP_from_acg(acg, n, dur, cont_thresh=15.0, conf_thresh=80.0,
                             rp_reject=0.001)
    mat, cont, _, _, _ = computeMatrix(st, {"recDur": dur})
    last, first, acc, _ = tau_last_pass(acg[None, :], np.array([n]), np.array([dur]))
    hl = {}
    for rp_ms, c in ((2.0, 10.0), (3.0, 10.0), (1.5, 20.0), (0.2, 10.0)):
        p, est, obs = hill_llobet_from_acg(acg, n, dur, rp_ms / 1000, c)
        hl[f"{rp_ms}_{c}"] = dict(passes=bool(p), est=None if np.isnan(est) else est,
                                  obs=obs)

    def nn(x):
        return None if x is None or (isinstance(x, float) and np.isnan(x)) else float(x)

    return dict(
        name=name, n=int(n), dur=dur, acg=acg.astype(int).tolist(),
        default=dict(max_conf=float(pk[0]), min_cont=nn(pk[1]),
                     rp_min_val=nn(pk[2]), passes=bool(pk[5])),
        alt=dict(max_conf=alt["max_conf"], min_cont=nn(alt["min_cont"]),
                 rp_min_val=nn(alt["rp_min_val"]), passes=alt["passes"],
                 tau_first_pass=nn(alt["tau_first_pass"])),
        tau_first=nn(first[0]), tau_last=nn(last[0]),
        matrix=mat.tolist(), hl=hl)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cases = []
    st = np.load(HERE / "fig2_example_neuron_spike_times.npy")
    cases.append(case_record("fig2", st, 3600.0))

    rng = np.random.default_rng(20260926)
    specs = [("clean_5hz", 5.0, 0.0, 0.002, 3600.0),
             ("thresh_5hz", 5.0, 0.10, 0.002, 3600.0),
             ("low_rate", 0.8, 0.05, 0.003, 1800.0),
             ("high_rate", 40.0, 0.20, 0.0015, 3600.0),
             ("short_rp", 12.0, 0.04, 0.0008, 7200.0),
             ("dirty", 8.0, 0.35, 0.0025, 1200.0)]
    for name, rate, cont, rp, dur in specs:
        st, _ = make_train("standard", rate, cont, dur, rp, rng=rng)
        cases.append(case_record(name, st, dur))
        print(f"{name}: {st.size:,} spikes", flush=True)

    # Poisson tail grid, including large counts where naive methods break
    grid = []
    for k in (0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 100000, 1000000):
        for ratio in (0.5, 0.8, 0.95, 1.0, 1.05, 1.2, 2.0):
            lam = max(ratio * (k + 1), 1e-3)
            grid.append(dict(k=k, lam=lam, sf=float(stats.poisson.sf(k, lam)),
                             cdf=float(stats.poisson.cdf(k, lam))))
    inv = []
    for k in (0, 5, 75, 1000, 50000):
        for g in (0.5, 0.8, 0.9, 0.99):
            inv.append(dict(k=k, g=g, lam=float(stats.chi2.ppf(g, 2 * (k + 1)) / 2)))

    # end-to-end: acceptance over many simulated trains
    sim = dict(rate=5.0, cont=0.10, rp=0.002, dur=3600.0, n_trains=1000)
    n_srp = n_hl = 0
    for _ in range(sim["n_trains"]):
        st, _ = make_train("standard", sim["rate"], sim["cont"], sim["dur"], sim["rp"],
                           rng=rng)
        acg = computeACG(st, BIN, NB)
        n_srp += slidingRP_from_acg(acg, st.size, sim["dur"])["passes"]
        n_hl += hill_llobet_from_acg(acg, st.size, sim["dur"], 0.002, 10.0)[0]
    sim["p_sliding"] = n_srp / sim["n_trains"]
    sim["p_hl2"] = n_hl / sim["n_trains"]
    print(f"simulation: Sliding RP accepts {sim['p_sliding']:.3f}, "
          f"HL 2 ms accepts {sim['p_hl2']:.3f}")

    (OUT / "reference.json").write_text(json.dumps(
        dict(cases=cases, poisson=grid, inverse=inv, sim=sim)))
    print("wrote", OUT / "reference.json")


if __name__ == "__main__":
    main()
