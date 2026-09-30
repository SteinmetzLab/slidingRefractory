"""Python half of the Python/MATLAB parity check.

Writes parity_cases.mat: a set of spike trains, the option sets to run them
under, and every output of the Python implementation for each. check_parity.m
recomputes the same quantities in MATLAB and compares.

Trains: three clusters from the repository's unit-test IBL recording (the
regression clusters) plus simulated units spanning low and high rates, short
and long refractory periods, clean and contaminated, and one with a sorter
censor applied. Option sets: the defaults, two censor windows, and
non-default thresholds with a longer tau_min.

Run from the repository root:
    python parity/make_cases.py
then in MATLAB, with histdiff (cortex-lab/spikes) on the path:
    run parity/check_parity.m
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat

from slidingRP import metrics
from slidingRP.power import min_passing_fr, tau_pass0
from slidingRP.simulations import RPmetric_Classic, genST

HERE = Path(__file__).parent
ROOT = HERE.parent
REC = 3600.0


def censored(st, w):
    keep = np.ones(st.size, bool)
    last = -np.inf
    for i, t in enumerate(st):
        if t - last < w:
            keep[i] = False
        else:
            last = t
    return st[keep]


def simulated(rate, cont, rp, dur, seed, censor=0.0):
    np.random.seed(seed)
    base = genST(rate * (1 - cont), dur, rp)
    other = genST(rate * cont, dur, 0.0) if cont > 0 else np.empty(0)
    st = np.sort(np.concatenate([base, other]))
    return censored(st, censor) if censor > 0 else st


def main():
    trains, durs, names = [], [], []
    t = np.load(ROOT / "test-data" / "unit" / "spikes.times.npy")
    c = np.load(ROOT / "test-data" / "unit" / "spikes.clusters.npy")
    for cid in (274, 167, 275):
        st = np.sort(t[c == cid]).astype(np.float64)
        trains.append(st); durs.append(float(t.max())); names.append(f"ibl_cluster_{cid}")
    for name, args in [("clean_5hz_3ms", (5.0, 0.0, 0.003, REC, 1)),
                       ("cont10_5hz_3ms", (5.0, 0.10, 0.003, REC, 2)),
                       ("cont12_5hz_3ms", (5.0, 0.12, 0.003, REC, 3)),
                       ("low_0p8hz", (0.8, 0.05, 0.002, REC, 4)),
                       ("high_30hz_1p5ms", (30.0, 0.15, 0.0015, REC, 5)),
                       ("censored_0p25ms", (5.0, 0.12, 0.003, REC, 6, 0.00025))]:
        trains.append(simulated(*args)); durs.append(args[3]); names.append(name)

    opts = [dict(name="default", cont=10.0, conf=90.0, rpReject=0.0005, censor=0.0),
            dict(name="censor_0p25ms", cont=10.0, conf=90.0, rpReject=0.0005, censor=0.00025),
            dict(name="censor_0p5ms", cont=10.0, conf=90.0, rpReject=0.0005, censor=0.0005),
            dict(name="c15_g80_taumin1", cont=15.0, conf=80.0, rpReject=0.001, censor=0.0)]

    nT, nO = len(trains), len(opts)
    out = {k: np.full((nT, nO), np.nan) for k in
           ("max_conf", "min_cont", "rp_min_val", "n_below2", "passes", "tau_pass0")}
    mats = np.empty((nT, nO), dtype=object)
    for i, (st, dur) in enumerate(zip(trains, durs)):
        for j, o in enumerate(opts):
            p = {"recDur": dur, "censor": o["censor"]}
            r = metrics.slidingRP(st, params=p, conf_thresh=o["conf"], cont_thresh=o["cont"],
                                  rp_reject=o["rpReject"])
            out["max_conf"][i, j], out["min_cont"][i, j], out["rp_min_val"][i, j] = r[0], r[1], r[2]
            out["n_below2"][i, j], out["passes"][i, j], out["tau_pass0"][i, j] = r[3], float(r[5]), r[7]
            cm, _, _, _, _ = metrics.computeMatrix(st, dict(p, rpReject=o["rpReject"]))
            mats[i, j] = cm

    # the FWER-corrected path, on two trains, default settings (slow in both)
    corr_idx = [names.index("cont10_5hz_3ms"), names.index("censored_0p25ms")]
    corr = np.full((len(corr_idx), 2), np.nan)          # columns: censor 0, censor 0.25 ms
    for a, i in enumerate(corr_idx):
        for b, w in enumerate((0.0, 0.00025)):
            corr[a, b] = metrics.slidingRP(trains[i], params={"recDur": durs[i], "censor": w,
                                                              "correction": True})[0]

    # Hill-Llobet, spike-time path
    hl_cfg = [(m, rp, w) for m in ("Llobet", "Hill") for rp in (0.002, 0.003)
              for w in (0.0, 0.00025)]
    hl_pass = np.full((nT, len(hl_cfg)), np.nan)
    hl_est = np.full((nT, len(hl_cfg)), np.nan)
    for i, (st, dur) in enumerate(zip(trains, durs)):
        for k, (m, rp, w) in enumerate(hl_cfg):
            ps, est = RPmetric_Classic(st, {"metricType": m, "RPdur": rp, "recDur": dur,
                                            "censor": w, "contaminationThresh": 10})
            hl_pass[i, k], hl_est[i, k] = float(ps), est

    # power functions
    pg = np.array([(n, d, c, g, w) for n in (500, 18000, 200000) for d in (1800.0, 3600.0)
                   for c in (5.0, 10.0) for g in (80.0, 90.0) for w in (0.0, 0.00025)])
    tp0 = np.array([tau_pass0(n, d, c, g, censor=w) for n, d, c, g, w in pg])
    mfr = np.array([min_passing_fr(d, 0.002, c, g, censor=w) for n, d, c, g, w in pg])

    tr = np.empty((nT,), dtype=object)
    for i, st in enumerate(trains):
        tr[i] = st.reshape(-1, 1)
    savemat(HERE / "parity_cases.mat", {
        "trains": tr, "recDur": np.array(durs), "names": np.array(names, dtype=object),
        "opt_cont": np.array([o["cont"] for o in opts]),
        "opt_conf": np.array([o["conf"] for o in opts]),
        "opt_rpReject": np.array([o["rpReject"] for o in opts]),
        "opt_censor": np.array([o["censor"] for o in opts]),
        "opt_names": np.array([o["name"] for o in opts], dtype=object),
        **{f"py_{k}": v for k, v in out.items()}, "py_matrix": mats,
        "corr_idx": np.array(corr_idx) + 1, "py_corrected": corr,
        "hl_metric": np.array([m for m, _, _ in hl_cfg], dtype=object),
        "hl_rp": np.array([rp for _, rp, _ in hl_cfg]),
        "hl_censor": np.array([w for _, _, w in hl_cfg]),
        "py_hl_pass": hl_pass, "py_hl_est": hl_est,
        "pow_grid": pg, "py_pow_tau_pass0": tp0, "py_min_fr": mfr,
    }, do_compression=True)
    print(f"wrote {HERE / 'parity_cases.mat'}: {nT} trains x {nO} option sets, "
          f"{len(hl_cfg)} Hill-Llobet configurations, {len(pg)} power-function points")


if __name__ == "__main__":
    main()
