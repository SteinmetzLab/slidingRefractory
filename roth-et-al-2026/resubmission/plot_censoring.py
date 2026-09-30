"""Figure for 08: what a sorter censor window does to Sliding RP.

a  Acceptance against true contamination when the data have a censor window w
   and the published method ignores it (Fig 4 setting: 5 spikes/s, 3 ms, 1 h).
b  The same trains scored with the censor accounted for (tau_r -> tau_r - w).
c  False acceptance (12% contamination) and true acceptance (8%) against w,
   ignoring and accounting for the censor, for true RPs of 1.5 and 3 ms.
d  Real units re-scored with their dataset's measured censor.

Run:  python plot_censoring.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402

SIMS = Path(r"D:/temp/slidingRP_resub/sims")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/08_censoring_window")
FAR_C, TAR_C = 0.12, 0.08            # the manuscript's Fig 4 convention


def colored_legend(ax, **kw):
    leg = ax.legend(handlelength=0, handletextpad=0, **kw)
    for t, h in zip(leg.get_texts(), leg.legend_handles):
        t.set_color(h.get_color())
        h.set_visible(False)
    return leg


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "censoring.pqt")
    d["w_ms"] = pd.to_numeric(d.censor) * 1000
    ws = sorted(d.w_ms.unique())
    cols = {w: ("k" if w == 0 else plt.cm.Oranges(0.35 + 0.6 * i / (len(ws) - 2)))
            for i, w in enumerate([x for x in ws if x > 0])}
    cols[0.0] = "k"
    lab = lambda w: "No censor" if w == 0 else f"w = {w:.3g} ms"

    fig = plt.figure(figsize=(14.5, 7.2))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.3)
    ref = d[np.isclose(d.rp_dur, 0.003)]

    for k, (col, title) in enumerate((("srp_pass", "a  Censor ignored (published method)"),
                                      ("srp_pass_cc", "b  Censor accounted for"))):
        ax = fig.add_subplot(gs[0, k])
        for w in ws:
            g = ref[np.isclose(ref.w_ms, w)].sort_values("cont_prop")
            ax.plot(g.cont_prop * 100, g[col] * 100, "o-", ms=3, lw=1.6,
                    color=cols[w], label=lab(w))
        ax.axvline(10, color="0.6", ls="--", lw=1)
        ax.set_xlabel("True contamination (%)")
        ax.set_ylabel("Units accepted (%)")
        ax.set_ylim(-2, 102)
        ax.set_title(title, loc="left")
        if k == 0:
            colored_legend(ax, fontsize=8, loc="center right")

    ax = fig.add_subplot(gs[0, 2])
    for rp, mk in ((0.003, "o"), (0.0015, "s")):
        g = d[np.isclose(d.rp_dur, rp)]
        for col, ls, name in (("srp_pass", "-", "ignored"), ("srp_pass_cc", "--", "accounted for")):
            y = [100 * g[np.isclose(g.w_ms, w) & np.isclose(g.cont_prop, FAR_C)][col].mean() for w in ws]
            ax.plot(ws, y, ls, marker=mk, ms=4.5, color="#c2410c" if ls == "-" else "0.25",
                    label=f"{rp*1000:g} ms RP, censor {name}")
    ax.set_xlabel("Censor window in the data, w (ms)")
    ax.set_ylabel(f"False acceptance at {FAR_C*100:g}% (%)")
    ax.set_title("c  False acceptance against w", loc="left")
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 0])
    for rp, mk in ((0.003, "o"), (0.0015, "s")):
        g = d[np.isclose(d.rp_dur, rp)]
        for col, ls, name in (("srp_pass", "-", "ignored"), ("srp_pass_cc", "--", "accounted for")):
            y = [100 * g[np.isclose(g.w_ms, w) & np.isclose(g.cont_prop, TAR_C)][col].mean() for w in ws]
            ax.plot(ws, y, ls, marker=mk, ms=4.5, color="#1b9e77" if ls == "-" else "0.25",
                    label=f"{rp*1000:g} ms RP, censor {name}")
    ax.set_xlabel("Censor window in the data, w (ms)")
    ax.set_ylabel(f"True acceptance at {TAR_C*100:g}% (%)")
    ax.set_title("d  True acceptance against w", loc="left")
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 1])
    for rpk, col, name in (("hl3_pass", "#c2410c", "Hill-Llobet 3 ms, ignored"),
                           ("hl3_pass_cc", "0.25", "Hill-Llobet 3 ms, accounted for")):
        y = [100 * ref[np.isclose(ref.w_ms, w) & np.isclose(ref.cont_prop, FAR_C)][rpk].mean() for w in ws]
        ax.plot(ws, y, "o-" if "cc" not in rpk else "o--", ms=4.5, color=col, label=name)
    ax.set_xlabel("Censor window in the data, w (ms)")
    ax.set_ylabel(f"False acceptance at {FAR_C*100:g}% (%)")
    ax.set_title("e  The fixed-window method too (3 ms RP)", loc="left")
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 2])
    rd = pd.read_csv(OUTDIR / "censor_realdata.csv")
    rd = rd[rd.dataset.isin(["allen", "macaque"])]
    x = np.arange(len(rd))
    ax.bar(x - 0.2, rd.published, 0.38, color="0.65", label="Published method")
    ax.bar(x + 0.2, rd.profile, 0.38, color="#c2410c", label="Measured censor accounted for")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{plotstyle.DATASET_NAMES[r.dataset]}\n{r.region}" for r in rd.itertuples()],
                       fontsize=7.5)
    ax.set_ylabel("Proportion of sorter-accepted units accepted")
    ax.set_title("f  Censored datasets, re-scored", loc="left")
    ax.set_ylim(0, 0.85)
    ax.legend(fontsize=7.5, loc="upper left")

    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)
    plotstyle.save(fig, OUTDIR / "figures" / "censoring")

    L = ["Censor window simulations (5 spikes/s, 1 h, true RP 1.5 and 3 ms)", "=" * 66, "",
         f"{'RP':>5} {'w (ms)':>7} | {'ignored: FA 12%':>15} {'TA 8%':>6} | "
         f"{'accounted: FA 12%':>17} {'TA 8%':>6} | {'HL3 FA ign.':>11} {'HL3 FA acc.':>11}"]
    for rp in (0.0015, 0.003):
        g = d[np.isclose(d.rp_dur, rp)]
        for w in ws:
            h = g[np.isclose(g.w_ms, w)].set_index(np.round(g[np.isclose(g.w_ms, w)].cont_prop, 4))
            L.append(f"{rp*1000:4.1f}  {w:7.3f} | {100*h.loc[FAR_C,'srp_pass']:14.1f}% "
                     f"{100*h.loc[TAR_C,'srp_pass']:5.1f}% | {100*h.loc[FAR_C,'srp_pass_cc']:16.1f}% "
                     f"{100*h.loc[TAR_C,'srp_pass_cc']:5.1f}% | {100*h.loc[FAR_C,'hl3_pass']:10.1f}% "
                     f"{100*h.loc[FAR_C,'hl3_pass_cc']:10.1f}%")
    L += ["", "1000 trains per point; conditions differing only in w share trains, so the",
          "w = 0 row is the same trains before censoring."]
    (OUTDIR / "censoring_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
