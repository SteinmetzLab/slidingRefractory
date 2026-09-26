"""Animated explainer of the Sliding RP calculation: Fig 2, built step by step.

Renders an H.264 MP4 at 1920 x 1080, 30 fps, which PowerPoint plays directly.

The unit is the example neuron from the published Fig 2: a 10 spikes/s neuron
with a hard 2.5 ms refractory period plus a 1 spike/s Poisson contaminating
source (about 9% of all spikes), recorded for 1 h. Every number the movie shows
is computed from its spike times with the same formulas as slidingRP.metrics,
and the confidence matrix, maximum confidence and C_min are asserted against the
package before anything is drawn.

One convention to know when comparing with the paper: the published figure
printed 109 observed violations at 1.5 ms, while the current package gives 110
for the same spike train (75 and 157 at 1 and 2.2 ms agree). The figure was drawn
with an earlier ACG routine; the movie uses the current package.

A candidate refractory period tau_r here is the right edge of an ACG bin, which
is the window the package's expected-violation count uses, so the x position of
the confidence curve and the width of the shaded window are the same number.

Run:
    python make_movie.py                      # captioned version
    python make_movie.py --no-captions        # clean version to narrate over
    python make_movie.py --stills 40,300,1500 # PNG stills of chosen frames
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import imageio_ffmpeg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from scipy import stats  # noqa: E402

from slidingRP.metrics import computeACG, computeMatrix, slidingRP  # noqa: E402

HERE = Path(__file__).parent
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/explainer")

FPS = 30
SR = 30000.0
BIN = 1 / SR
NB = 300                                   # 0-10 ms, as in the package
RP_CENTERS = np.arange(0, 10 / 1000, BIN) + BIN / 2
REF_DUR = RP_CENTERS + BIN / 2             # window duration for bin k
TAU_MS = REF_DUR * 1000
TAU_MIN = 0.0005
REC_DUR = 3600.0
CONT = np.arange(0.5, 35.5, 0.5)           # package default grid (%)
GAMMA = 90.0
C_THRESH = 10.0
XMAX_MS = 5.0
NSHOW = int(round(XMAX_MS / 1000 / BIN))   # ACG bins drawn (0-5 ms)
TRUE_RP_MS = 2.5

# the paper's colors: three example windows, three contamination thresholds
EXAMPLES = [(1.0, "steelblue", "navy"),
            (1.5, "forestgreen", "darkgreen"),
            (2.2, "orchid", "darkorchid")]
SWEEP_FILL, SWEEP_LINE = "#f4a09a", "red"
C_COLORS = {15.0: "darkred", 10.0: "red", 7.5: "lightcoral"}


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------

def load_data():
    st = np.load(HERE / "fig2_example_neuron_spike_times.npy")
    n = st.size
    acg = computeACG(st, BIN, NB)
    obs = np.cumsum(acg)

    def expected(c_pct):
        c = c_pct / 100
        nc, nb = c * n, (1 - c) * n
        return 2 * REF_DUR / REC_DUR * nc * (nb + (nc - 1) / 2)

    lam = {c: expected(c) for c in CONT}
    mat = np.array([100 * (1 - stats.poisson.cdf(obs, lam[c])) for c in CONT])

    # the package must agree before we animate anything
    ref_mat, ref_cont, _, _, _ = computeMatrix(st, {"recDur": REC_DUR})
    assert np.allclose(ref_cont, CONT)
    assert np.allclose(mat, ref_mat, atol=1e-9), "confidence matrix disagrees"
    out = slidingRP(st, params={"recDur": REC_DUR})
    test = RP_CENTERS > TAU_MIN
    row = {c: mat[np.flatnonzero(CONT == c)[0]] for c in C_COLORS}
    assert abs(row[C_THRESH][test].max() - out[0]) < 1e-9, "max confidence disagrees"

    k_max = int(np.flatnonzero(test)[np.argmax(row[C_THRESH][test])])
    return dict(st=st, n=n, acg=acg, obs=obs, lam=lam, mat=mat, row=row,
                test=test, max_conf=float(out[0]), k_max=k_max,
                c_min=float(out[1]), tau_cmin_ms=float(out[2]) * 1000 + BIN * 500,
                accepted=bool(out[5]))


def k_of(tau_ms):
    """Bin whose right edge is tau_ms (-1 means an empty window)."""
    return int(np.clip(round(tau_ms / 1000 / BIN) - 1, -1, NB - 1))


def ease(p):
    return 0.5 - 0.5 * np.cos(np.pi * np.clip(p, 0, 1))


def seg(p, a, b):
    return float(np.clip((p - a) / (b - a), 0, 1))


# --------------------------------------------------------------------------
# the movie
# --------------------------------------------------------------------------

class Movie:
    def __init__(self, d, captions=True):
        self.d = d
        self.captions = captions
        self._style()
        self._build()
        self.scenes = self._timeline()
        self.starts = np.cumsum([0] + [s["n"] for s in self.scenes])
        self.n_frames = int(self.starts[-1])
        self.cur = -1
        self.cap_text, self.cap_frame = None, 0

    # ---- figure -----------------------------------------------------------
    @staticmethod
    def _style():
        plt.rcParams.update({
            "font.family": "Arial", "font.sans-serif": ["Arial", "DejaVu Sans"],
            "mathtext.fontset": "stix",
            "axes.spines.top": False, "axes.spines.right": False,
            "font.size": 15, "axes.labelsize": 17, "axes.titlesize": 17,
            "xtick.labelsize": 14, "ytick.labelsize": 14,
            "axes.linewidth": 1.2, "xtick.major.width": 1.2,
            "ytick.major.width": 1.2, "pdf.fonttype": 42, "ps.fonttype": 42,
        })

    def _build(self):
        d = self.d
        fig = plt.figure(figsize=(16, 9), dpi=120)
        fig.patch.set_facecolor("white")
        self.fig = fig
        y_top, h = (0.535, 0.305) if self.captions else (0.575, 0.325)
        y_bot = 0.085
        self.axA = fig.add_axes([0.07, y_top, 0.40, h])
        self.axB = fig.add_axes([0.07, y_bot, 0.40, h])
        self.axC = fig.add_axes([0.575, y_top, 0.335, h])
        self.axD = fig.add_axes([0.575, y_bot, 0.335, h])
        self.cax = fig.add_axes([0.925, y_bot, 0.011, h])

        # --- A: ACG with the sliding window ---
        ax = self.axA
        self.ymaxA = float(d["acg"][:NSHOW].max()) * 1.12
        self.bars = ax.bar(RP_CENTERS[:NSHOW] * 1000, np.zeros(NSHOW),
                           width=BIN * 1000, color="0.12", linewidth=0, zorder=2)
        ax.set_xlim(0, XMAX_MS)
        ax.set_ylim(0, self.ymaxA)
        ax.set_xlabel("Time from spike (ms)")
        ax.set_ylabel("Spike pairs per bin")
        self.win = Rectangle((0, 0), 0, self.ymaxA, facecolor="none",
                             edgecolor="none", alpha=0.35, zorder=1)
        ax.add_patch(self.win)
        self.win_line = ax.axvline(0, lw=2.4, color="none", zorder=3)
        self.tau_label = ax.text(0, self.ymaxA * 1.03, "", ha="center",
                                 va="bottom", fontsize=18, clip_on=False)
        self.count_label = ax.text(0.06, self.ymaxA * 0.86, "", ha="left",
                                   va="center", fontsize=20, zorder=4)

        # --- B: expected violations under the threshold ---
        ax = self.axB
        self.xB = np.arange(0, 601)
        (self.pmf_line,) = ax.plot([], [], lw=2.4)
        self.pmf_fill = None
        self.obs_line = ax.axvline(np.nan, ymax=0.76, color="0.1", lw=1.8, ls="--")
        ax.set_xlim(0, 260)
        ax.set_ylim(0, 0.05)
        ax.set_xlabel("Number of violations")
        ax.set_ylabel("Probability")
        ax.set_title(f"Violations expected if the unit were {C_THRESH:g}% contaminated",
                     loc="left", pad=10)
        self.mean_text = ax.text(0, 0, "", fontsize=15, va="center")
        self.obs_text = ax.text(0, 0, "", fontsize=15, va="bottom", ha="left",
                                color="0.1")
        self.prob_text = ax.text(0.99, 0.97, "", transform=ax.transAxes,
                                 ha="right", va="top", fontsize=16,
                                 linespacing=1.5)
        self.offscale = ax.text(0.99, 0.60, "", transform=ax.transAxes,
                                ha="right", va="center", fontsize=15, color="0.1")

        # --- C: confidence across candidate tau_r ---
        ax = self.axC
        ax.set_xlim(0, XMAX_MS)
        ax.set_ylim(0, 100)
        ax.set_xlabel(r"Tested refractory period $\tau_r$ (ms)")
        ax.set_ylabel("Confidence (%)")
        ax.set_title("Confidence across candidate refractory periods",
                     loc="left", pad=10)
        ax.axvspan(0, TAU_MIN * 1000, color="0.88", zorder=0, lw=0)
        ax.text(TAU_MIN * 500, 3, "Not tested", rotation=90, ha="center",
                va="bottom", fontsize=12, color="0.4")
        self.gamma_line = ax.axhline(GAMMA, color="k", lw=1.8, zorder=2)
        self.gamma_text = ax.text(XMAX_MS * 0.99, GAMMA + 1.5,
                                  f"Acceptance threshold ({GAMMA:g}%)",
                                  ha="right", va="bottom", fontsize=13)
        self.curves, self.c_labels = {}, {}
        for i, (c, col) in enumerate(C_COLORS.items()):
            lw = 3.0 if c == C_THRESH else 1.8
            (self.curves[c],) = ax.plot([], [], color=col, lw=lw, zorder=3)
            self.c_labels[c] = ax.text(3.35, 76 - 8.5 * i,
                                       rf"$C_{{\mathrm{{thresh}}}}$ = {c:g}%",
                                       color=col, fontsize=16, alpha=0)
        self.ex_marks = []
        for tau, fill, line in EXAMPLES:
            (m,) = ax.plot([], [], marker="X", ms=15, color=line, mec="white",
                           mew=1.2, ls="none", zorder=5)
            self.ex_marks.append(m)
        (self.star,) = ax.plot([], [], marker="*", ms=24, color="red",
                               mec="white", mew=1.4, ls="none", zorder=6)
        self.verdict = ax.text(3.35, 30, "", fontsize=16, va="center",
                               linespacing=1.4)

        # --- D: confidence matrix ---
        ax = self.axD
        x_edges = np.arange(NB + 1) * BIN * 1000
        y_edges = np.r_[CONT - 0.25, CONT[-1] + 0.25]
        self.mask_rows = np.arange(len(CONT))[:, None] * np.ones((1, NB))
        self.mesh = ax.pcolormesh(x_edges, y_edges,
                                  np.ma.masked_all(d["mat"].shape), cmap="viridis",
                                  vmin=0, vmax=100, shading="flat", rasterized=True)
        ax.set_xlim(0, XMAX_MS)
        ax.set_ylim(CONT[-1] + 0.25, 0)
        ax.set_xlabel(r"Tested refractory period $\tau_r$ (ms)")
        ax.set_ylabel("Contamination threshold (%)")
        ax.set_title("Confidence matrix", loc="left", pad=10)
        ax.add_patch(Rectangle((0, 0), TAU_MIN * 1000, CONT[-1] + 0.25,
                               facecolor="0.55", alpha=0.55, lw=0, zorder=3))
        cb = fig.colorbar(self.mesh, cax=self.cax)
        cb.set_label("Confidence (%)")
        cb.outline.set_visible(False)
        self.contour = ax.contour(TAU_MS - BIN * 500, CONT, d["mat"], levels=[GAMMA],
                                  colors="k", linewidths=2.4, zorder=4)
        self.contour.set_alpha(0)
        self.hlines = {c: ax.axhline(c, color=col, lw=2.2 if c == C_THRESH else 1.5,
                                     alpha=0, zorder=5)
                       for c, col in C_COLORS.items()}
        (self.cmin_mark,) = ax.plot([], [], marker="o", ms=13, color="white",
                                    mec="k", mew=2.2, ls="none", zorder=7)
        self.cmin_vline = ax.axvline(np.nan, color="white", ls="--", lw=1.6,
                                     zorder=6)
        self.cmin_text = ax.text(0.985, 0.04, "", transform=ax.transAxes,
                                 ha="right", va="bottom", fontsize=15, zorder=8,
                                 bbox=dict(boxstyle="round,pad=0.4", fc="white",
                                           ec="none", alpha=0.92))

        for a in (self.axB, self.axC, self.axD, self.cax):
            a.set_visible(False)

        self.cap = fig.text(0.07, 0.95, "", fontsize=20, va="center", ha="left",
                            color="0.1", linespacing=1.35)

    # ---- state setters ----------------------------------------------------
    def set_window(self, tau_ms, fill, line, show_count=True):
        k = k_of(tau_ms)
        tau = 0.0 if k < 0 else TAU_MS[k]
        n_obs = 0 if k < 0 else int(self.d["obs"][k])
        self.win.set_width(tau)
        self.win.set_facecolor(fill)
        self.win_line.set_xdata([tau, tau])
        self.win_line.set_color(line)
        self.tau_label.set_position((max(tau, 0.35), self.ymaxA * 1.03))
        self.tau_label.set_text(rf"$\tau_r$ = {tau:.2f} ms" if tau else "")
        self.tau_label.set_color(line)
        self.count_label.set_text(f"{n_obs:,}" if show_count and tau else "")
        self.count_label.set_color(line)
        return k

    def set_pmf(self, k, line, fill, reveal=1.0, shade=1.0, texts=1.0):
        d = self.d
        ax = self.axB
        if self.pmf_fill is not None:
            self.pmf_fill.remove()
            self.pmf_fill = None
        if k < 0:
            self.pmf_line.set_data([], [])
            for t in (self.mean_text, self.obs_text, self.prob_text, self.offscale):
                t.set_text("")
            self.obs_line.set_xdata([np.nan, np.nan])
            return
        lam = float(d["lam"][C_THRESH][k])
        o = int(d["obs"][k])
        x = self.xB
        y = stats.poisson.pmf(x, lam)
        m = max(2, int(len(x) * reveal))
        self.pmf_line.set_data(x[:m], y[:m])
        self.pmf_line.set_color(line)
        xl = ax.get_xlim()[1]
        if shade > 0:
            xs = x[x <= min(o, x[-1])]
            self.pmf_fill = ax.fill_between(xs, 0, y[:xs.size], color=fill,
                                            alpha=0.55 * shade, lw=0, zorder=1)
        p_le = 100 * stats.poisson.cdf(o, lam)
        conf = 100 - p_le
        vis = texts > 0
        peak = y.max()
        sd = np.sqrt(max(lam, 1e-9))
        mx = lam + 2.2 * sd + 0.02 * xl
        ha = "left"
        if mx > 0.72 * xl:
            mx, ha = lam - 2.2 * sd - 0.02 * xl, "right"
        self.mean_text.set_position((mx, 0.72 * peak))
        self.mean_text.set_ha(ha)
        self.mean_text.set_color(line)
        self.mean_text.set_text(f"Expected if {C_THRESH:g}% contaminated\n"
                                f"mean = {lam:.1f}" if reveal >= 1 else "")
        if shade > 0 and o <= xl:
            self.obs_line.set_xdata([o, o])
            self.obs_line.set_alpha(shade)
            self.obs_text.set_position((o + 0.012 * xl, 0.08 * ax.get_ylim()[1]))
            self.obs_text.set_text(f"Observed {o:,}")
            self.obs_text.set_alpha(shade)
            self.offscale.set_text("")
        else:
            self.obs_line.set_xdata([np.nan, np.nan])
            self.obs_text.set_text("")
            self.offscale.set_text(f"Observed {o:,} \u2192" if shade > 0 else "")
        if vis:
            self.prob_text.set_text(
                f"P(\u2264 {o:,} | {C_THRESH:g}% contamination) = {p_le:.1f}%\n"
                f"Confidence = 100% \u2212 {p_le:.1f}% = {conf:.1f}%")
            self.prob_text.set_alpha(texts)
        else:
            self.prob_text.set_text("")

    def trace_curve(self, c, upto_ms):
        m = TAU_MS <= upto_ms + 1e-12
        self.curves[c].set_data(TAU_MS[m], self.d["row"][c][m])

    def place_example(self, i, alpha=1.0):
        tau = EXAMPLES[i][0]
        k = k_of(tau)
        self.ex_marks[i].set_data([TAU_MS[k]], [self.d["row"][C_THRESH][k]])
        self.ex_marks[i].set_alpha(alpha)

    def example_numbers(self, i):
        k = k_of(EXAMPLES[i][0])
        o = int(self.d["obs"][k])
        lam = float(self.d["lam"][C_THRESH][k])
        return o, lam, float(self.d["row"][C_THRESH][k])

    # ---- timeline -----------------------------------------------------------
    def _timeline(self):
        d = self.d
        S = []

        def add(sec, step, enter=None, caption=None):
            S.append(dict(n=int(round(sec * FPS)), step=step, enter=enter,
                          caption=caption))

        ex_fill0, ex_line0 = EXAMPLES[0][1], EXAMPLES[0][2]

        # 1. the ACG
        def s_rise(p):
            h = d["acg"][:NSHOW] * ease(p)
            for b, v in zip(self.bars, h):
                b.set_height(v)
        add(2.5, s_rise, caption="The autocorrelogram (ACG) of one unit: how many "
            "pairs of spikes are separated by each time lag")
        add(1.5, lambda p: None, caption=lambda p: "The autocorrelogram (ACG) of one unit: "
            "how many pairs of spikes are separated by each time lag")

        # 2. a window, and the violations inside it
        def s_win1(p):
            self.set_window(EXAMPLES[0][0] * ease(p), ex_fill0, ex_line0)
        cap2 = ("Pick a candidate refractory period $\\tau_r$ and count the spikes "
                "that fall inside it: the observed violations")
        add(3.0, s_win1, caption=cap2)
        add(1.0, lambda p: None, caption=cap2)

        # 3. the expectation under the threshold
        o1, lam1, conf1 = self.example_numbers(0)

        def e_pmf():
            self.axB.set_visible(True)

        def s_pmf(p):
            self.set_pmf(k_of(EXAMPLES[0][0]), ex_line0, ex_fill0,
                         reveal=ease(p), shade=0, texts=0)
        add(3.0, s_pmf, enter=e_pmf,
            caption=f"If this unit were {C_THRESH:g}% contaminated, how many would "
            f"we expect?\nA Poisson count with mean {lam1:.1f}")

        # 4. how surprising is the observed count?
        def s_shade(p):
            self.set_pmf(k_of(EXAMPLES[0][0]), ex_line0, ex_fill0, reveal=1,
                         shade=ease(seg(p, 0, 0.4)), texts=ease(seg(p, 0.4, 0.7)))
        add(3.5, s_shade,
            caption=f"{o1} or fewer would happen {100 - conf1:.1f}% of the time at "
            f"{C_THRESH:g}% contamination,\nso the confidence that contamination is "
            f"below {C_THRESH:g}% is {conf1:.1f}%")

        # 5. onto the confidence curve
        def e_curve():
            self.axC.set_visible(True)

        def s_point1(p):
            self.place_example(0, alpha=ease(seg(p, 0.2, 0.6)))
        add(2.5, s_point1, enter=e_curve,
            caption=f"Plot that confidence against $\\tau_r$. The unit is accepted if "
            f"confidence reaches {GAMMA:g}% (black line)")

        # 6-7. two more windows
        def slide(i_from, i_to):
            t0, t1 = EXAMPLES[i_from][0], EXAMPLES[i_to][0]
            fill, line = EXAMPLES[i_to][1], EXAMPLES[i_to][2]

            def step(p):
                tau = t0 + (t1 - t0) * ease(seg(p, 0, 0.45))
                k = self.set_window(tau, fill, line)
                self.set_pmf(k, line, fill)
                self.place_example(i_to, alpha=ease(seg(p, 0.55, 0.75)))
            return step

        for i in (1, 2):
            o, lam, conf = self.example_numbers(i)
            add(4.0, slide(i - 1, i),
                caption=f"At $\\tau_r$ = {EXAMPLES[i][0]:g} ms: {o} observed, {lam:.1f} "
                f"expected, confidence {conf:.1f}%")

        # 8. the full sweep
        def e_sweep():
            for m in self.ex_marks:
                m.set_zorder(6)

        def s_sweep(p):
            xl = 260 + (600 - 260) * ease(seg(p, 0.0, 0.07))
            self.axB.set_xlim(0, xl)
            tau = XMAX_MS * seg(p, 0.03, 0.97)
            k = self.set_window(tau, SWEEP_FILL, SWEEP_LINE)
            if k >= 0:
                y = stats.poisson.pmf(self.xB, float(d["lam"][C_THRESH][k]))
                self.axB.set_ylim(0, max(0.05, 1.18 * y.max()))
            self.set_pmf(k, SWEEP_LINE, SWEEP_FILL)
            self.trace_curve(C_THRESH, tau)
            self.c_labels[C_THRESH].set_alpha(ease(seg(p, 0.05, 0.15)))

        def cap_sweep(p):
            tau = XMAX_MS * seg(p, 0.03, 0.97)
            if tau <= TRUE_RP_MS + 0.08:
                return ("Now slide $\\tau_r$ across every candidate duration, tracing "
                        "out the confidence curve")
            return (f"Past the unit's true refractory period ({TRUE_RP_MS:g} ms), "
                    "real spikes enter the window and confidence collapses")
        add(13.0, s_sweep, enter=e_sweep, caption=cap_sweep)

        # 9. the verdict
        tmax = TAU_MS[d["k_max"]]

        def s_verdict(p):
            a = ease(seg(p, 0, 0.35))
            self.star.set_data([tmax], [d["max_conf"]])
            self.star.set_alpha(a)
            word = "accepted" if d["accepted"] else "rejected"
            self.verdict.set_text(f"Maximum {d['max_conf']:.1f}%\n"
                                  f"at $\\tau_r$ \u2248 {tmax:.1f} ms: {word}")
            self.verdict.set_alpha(a)
        verdict_word = "Accepted" if d["accepted"] else "Rejected"
        add(4.5, s_verdict,
            caption=f"{verdict_word}: confidence passes {GAMMA:g}% at some $\\tau_r$, "
            "without assuming what the refractory period is")

        # 10. other thresholds
        test = d["test"]
        max15 = float(d["row"][15.0][test].max())
        max75 = float(d["row"][7.5][test].max())
        assert max15 >= GAMMA > max75, (max15, max75)

        def s_thresh(p):
            self.trace_curve(15.0, XMAX_MS * seg(p, 0.0, 0.42))
            self.c_labels[15.0].set_alpha(ease(seg(p, 0.05, 0.2)))
            self.trace_curve(7.5, XMAX_MS * seg(p, 0.5, 0.92))
            self.c_labels[7.5].set_alpha(ease(seg(p, 0.55, 0.7)))
        add(5.5, s_thresh,
            caption="Repeat for other contamination thresholds: below 15% is "
            "confirmed easily,\nbelow 7.5% never reaches the acceptance threshold")

        # 11. the whole matrix
        def e_matrix():
            self.axD.set_visible(True)
            self.cax.set_visible(True)

        def s_matrix(p):
            r = (len(CONT) + 1) * ease(seg(p, 0.0, 0.62))
            self.mesh.set_array(np.ma.masked_where(self.mask_rows >= r, d["mat"]))
            self.contour.set_alpha(ease(seg(p, 0.66, 0.82)))
            for c, h in self.hlines.items():
                h.set_alpha(0.95 * ease(seg(p, 0.84, 0.97)))
        add(7.5, s_matrix, enter=e_matrix,
            caption="Every threshold at once: the confidence matrix.\n"
            f"The black line is the {GAMMA:g}% confidence boundary")

        # 12. C_min
        def s_cmin(p):
            a = ease(seg(p, 0, 0.35))
            self.cmin_mark.set_data([d["tau_cmin_ms"]], [d["c_min"]])
            self.cmin_mark.set_alpha(a)
            self.cmin_vline.set_xdata([d["tau_cmin_ms"]] * 2)
            self.cmin_vline.set_alpha(a)
            self.cmin_text.set_text(
                f"Lowest contamination confirmed at {GAMMA:g}%:\n"
                f"$C_{{\\min}}$ = {d['c_min']:.1f}%  at  $\\tau_r$ \u2248 "
                f"{d['tau_cmin_ms']:.1f} ms")
            self.cmin_text.set_alpha(a)
            self.cmin_text.get_bbox_patch().set_alpha(0.92 * a)
        cap12 = (f"The lowest contamination the data can confirm at {GAMMA:g}% "
                 f"confidence: $C_{{\\min}}$ = {d['c_min']:.1f}%  "
                 "(the true value here is about 9%)")
        add(4.5, s_cmin, caption=cap12)
        add(3.0, lambda p: None, caption=cap12)
        return S

    # ---- frame driver -------------------------------------------------------
    def goto(self, f):
        i = int(np.searchsorted(self.starts, f, side="right") - 1)
        i = min(i, len(self.scenes) - 1)
        while self.cur < i:                      # finish any skipped scenes
            if self.cur >= 0:
                self.scenes[self.cur]["step"](1.0)
            self.cur += 1
            if self.scenes[self.cur]["enter"]:
                self.scenes[self.cur]["enter"]()
        sc = self.scenes[i]
        p = (f - self.starts[i]) / max(sc["n"] - 1, 1)
        sc["step"](p)
        if self.captions:
            c = sc["caption"](p) if callable(sc["caption"]) else sc["caption"]
            if c != self.cap_text:
                jumped = f != getattr(self, "last_f", -1) + 1
                self.cap_text, self.cap_frame = c, (f - 8 if jumped else f)
                self.cap.set_text(c or "")
            self.cap.set_alpha(min(1.0, (f - self.cap_frame + 1) / 8))
        self.last_f = f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-captions", action="store_true")
    ap.add_argument("--stills", default="", help="comma-separated frame numbers")
    args = ap.parse_args()

    d = load_data()
    print(f"{d['n']:,} spikes; max confidence {d['max_conf']:.2f}% at "
          f"{TAU_MS[d['k_max']]:.3f} ms; C_min {d['c_min']:.2f}% at "
          f"{d['tau_cmin_ms']:.3f} ms; accepted={d['accepted']}")
    for i, (tau, *_rest) in enumerate(EXAMPLES):
        k = k_of(tau)
        print(f"  tau {tau} ms: observed {int(d['obs'][k])}, expected "
              f"{d['lam'][C_THRESH][k]:.1f}, confidence {d['row'][C_THRESH][k]:.1f}%")

    mv = Movie(d, captions=not args.no_captions)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tag = "nocaptions" if args.no_captions else "captions"

    if args.stills:
        for f in sorted(int(x) for x in args.stills.split(",")):
            mv.goto(f)
            path = OUTDIR / f"_still_{tag}_{f:04d}.png"
            mv.fig.savefig(path, dpi=60)
            print("wrote", path)
        return

    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    writer = animation.FFMpegWriter(
        fps=FPS, codec="libx264", bitrate=-1,
        extra_args=["-pix_fmt", "yuv420p", "-profile:v", "high", "-crf", "18",
                    "-preset", "slow", "-threads", "2", "-movflags", "+faststart"])
    path = OUTDIR / f"sliding_rp_explainer_{tag}.mp4"
    with writer.saving(mv.fig, str(path), dpi=120):
        for f in range(mv.n_frames):
            mv.goto(f)
            writer.grab_frame()
            if f % 150 == 0:
                print(f"  frame {f}/{mv.n_frames}", flush=True)
    print(f"wrote {path} ({mv.n_frames / FPS:.1f} s)")


if __name__ == "__main__":
    main()
