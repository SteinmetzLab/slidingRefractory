"""One cross-correlogram per recipient-donor pair with a narrow short-lag dip (05).

A dip at |lag| < ~0.5 ms between nearby units most likely means the sorter
failed to resolve overlapping waveforms; a monosynaptic inhibitory connection
would instead give a trough displaced to one side of zero, starting after a
synaptic delay. These pairs hide injected contamination from Sliding RP's
shortest windows, so they are separated from the main semi-synthetic test and
inspected here one by one.

A pair is flagged as a dip pair when the coincidence at 0-0.5 ms is less than
half of that at 0.5-2 ms, (1 + kappa(0-0.5)) / (1 + kappa(0.5-2)) < DIP_RATIO.
The ratio, not kappa(0-0.5) alone, so that pairs whose firing is anticorrelated
at all lags (a genuine rate effect, which belongs in the main test) are not
flagged.

Output: 05_model_mismatch/figures/dip_pairs.pdf (multi-page), one panel per
pair: the signed CCG from -5 to +5 ms in 0.1 ms bins, as observed / expected if
independent (1 = no excess), with Sliding RP's acceptance at 10% and 20%
injected and the depth separation.

Run:  python plot_dip_pairs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).parent))
import plotstyle  # noqa: E402
from plot_semisynthetic_clean import PAIR, pair_table  # noqa: E402
from semisynthetic_clean import CACHE, FS, clean_train  # noqa: E402

SIMS = Path(r"D:/temp/slidingRP_resub/sims")
OUTDIR = Path(r"D:/Dropbox/papers/2026_SlidingRP/jNeurophysResubmission/05_model_mismatch")
DIP_RATIO = 0.5
HALF = 150            # 5 ms in samples
STEP = 3              # 0.1 ms


def dip_mask(p):
    return (1 + p.k_0_0p5) / (1 + p.k_0p5_2) < DIP_RATIO


def signed_ccg(a, b, dur):
    """CCG of donor b relative to recipient a, lags -5..5 ms, 0.1 ms bins,
    as observed / expected (lag-0 coincidences excluded)."""
    lo = np.searchsorted(b, a - HALF, side="left")
    hi = np.searchsorted(b, a + HALF, side="right")
    cnt = hi - lo
    i = np.repeat(np.arange(a.size), cnt)
    j = np.repeat(lo - np.r_[0, np.cumsum(cnt)[:-1]], cnt) + np.arange(cnt.sum())
    lag = b[j] - a[i]
    lag = lag[lag != 0]
    edges = np.arange(-HALF, HALF + STEP, STEP)
    n = np.histogram(lag, edges)[0].astype(float)
    # integer lags; bin [0, 3) has lost lag 0, so it spans 2 lags, not 3
    width = np.full(n.size, STEP, float)
    width[np.flatnonzero(edges[:-1] == 0)] = STEP - 1
    expected = a.size * b.size * width / (dur * FS)
    return (edges[:-1] + STEP / 2) / FS * 1000, n / expected


def shape_similarity(ccg, acg, x):
    """Correlation between a pair's CCG and the recipient's own (raw) ACG over
    0.1-5 ms, both folded to |lag|. Near 1 when the CCG reproduces the
    recipient's refractory trough and burst peaks, as it does when the two
    clusters are one neuron split in two."""
    m = np.abs(x) > 0.1
    fold = lambda v: np.array([np.mean(v[m][np.isclose(np.abs(x[m]), a)]) for a in np.unique(np.abs(x[m]).round(4))])
    a, b = fold(ccg), fold(acg)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan


def dip_shape(ccg, x):
    """Width and symmetry of the dip in a signed CCG (observed / expected).

    narrow: depth at |lag| 0.1-0.3 ms relative to 0.6-1.0 ms; a notch confined
            to the waveform-overlap zone gives a small value, a trough as wide
            as a refractory period gives a value near 1.
    asym:   (right - left) / (right + left) of the CCG at 0.1-1 ms; near 0 for a
            symmetric trough, near +/-1 when only one side is suppressed (as an
            inhibitory connection, or a burst split across two clusters, would give).
    """
    ax_ = np.abs(x)
    inner = np.mean(ccg[(ax_ > 0.1) & (ax_ < 0.3)])
    outer = np.mean(ccg[(ax_ > 0.6) & (ax_ < 1.0)])
    sel = (ax_ > 0.1) & (ax_ < 1.0)
    r, l = np.mean(ccg[sel & (x > 0)]), np.mean(ccg[sel & (x < 0)])
    return dict(inner=inner, outer=outer, narrow=inner / outer if outer > 0 else np.nan,
                asym=(r - l) / (r + l) if (r + l) > 0 else np.nan)


def load_templates(pid):
    """Cluster template waveforms (82 samples x 32 channels) and their channels,
    indexed like the npz cluster ids (fetch_waveforms.py)."""
    z = np.load(CACHE / f"{pid}_waveforms.npz")
    return z["waveforms"], z["channels"]


def template_similarity(wf, ch, a, b):
    """Cosine similarity of two clusters' templates over the union of their
    channels (zero where a template has no channel). A neuron split into two
    clusters gives nearly identical templates (near 1); two neurons whose
    footprints merely overlap give clearly lower values."""
    chans = np.union1d(ch[a], ch[b])
    va, vb = np.zeros((wf.shape[1], chans.size)), np.zeros((wf.shape[1], chans.size))
    va[:, np.searchsorted(chans, ch[a])] = wf[a]
    vb[:, np.searchsorted(chans, ch[b])] = wf[b]
    den = np.linalg.norm(va) * np.linalg.norm(vb)
    return float((va * vb).sum() / den) if den > 0 else np.nan


def main():
    plotstyle.apply()
    d = pd.read_parquet(SIMS / "semisynthetic_clean.pqt")
    p = pair_table(d)
    acc = p.pivot_table(index=PAIR, columns="f_injected", values="passes")
    pairs = p.drop_duplicates(PAIR).set_index(PAIR)
    dips = pairs[dip_mask(pairs)].reset_index().sort_values(["pid", "k_0_0p5"]).set_index(PAIR)
    print(f"{len(dips)} dip pairs of {len(pairs)} "
          f"({', '.join(f'{k} {v}' for k, v in dips.donor_kind.value_counts().items())})")

    per_page, ncol = 30, 6
    out = OUTDIR / "figures" / "dip_pairs.pdf"
    cache, sims = {}, []
    with PdfPages(out) as pdf:
        keys = list(dips.index)
        for start in range(0, len(keys), per_page):
            fig, axs = plt.subplots(per_page // ncol, ncol, figsize=(15, 13), sharex=True,
                                    squeeze=False)
            for ax, key in zip(axs.flat, keys[start:start + per_page]):
                pid, rcp, don = key
                if pid not in cache:
                    z = np.load(CACHE / f"{pid}.npz", allow_pickle=True)
                    cache.clear()
                    cache[pid] = (z["samples"], z["clusters"], float(z["rec_dur"]),
                                  *load_templates(pid))
                s, c, dur, wf, ch = cache[pid]
                tsim = template_similarity(wf, ch, rcp, don)
                raw = np.sort(s[c == rcp])
                a = clean_train(raw)
                b = np.sort(s[c == don])
                x, y = signed_ccg(a, b, dur)
                _, ya = signed_ccg(raw, raw, dur)     # recipient's own ACG, same scale
                sim = shape_similarity(y, ya, x)
                shp = dip_shape(y, x)
                sims.append(dict(pid=pid, recipient=rcp, donor=don, similarity=sim,
                                 template_sim=tsim, **shp,
                                 sep_um=dips.loc[key].depth_sep_um,
                                 acc10=acc.loc[key, 0.10], acc20=acc.loc[key, 0.20]))
                ax.stairs(y, np.r_[x - 0.05, x[-1] + 0.05], color="0.3", lw=0.7, fill=True,
                          alpha=0.55)
                ax.plot(x, ya, color="#7570b3", lw=0.9)
                ax.axhline(1, color="#d95f02", lw=0.8)
                for v in (-0.5, 0.5):
                    ax.axvline(v, color="#1b9e77", lw=0.8, ls=":")
                r = dips.loc[key]
                a10, a20 = acc.loc[key, 0.10], acc.loc[key, 0.20]
                ax.set_title(f"{pid[:8]} {rcp}/{don}, {r.depth_sep_um:.0f} µm, "
                             f"{r.fr_recipient:.1f}/{r.fr_donor:.1f} spikes/s\n"
                             f"Accepted {a10:.0%} at 10%, {a20:.0%} at 20%\n"
                             f"Template {tsim:.2f}, ACG similarity {sim:.2f}, asymmetry {shp['asym']:+.2f}",
                             fontsize=6.5)
                ax.tick_params(labelsize=7)
            for ax in axs.flat[len(keys[start:start + per_page]):]:
                ax.set_visible(False)
            for ax in axs[-1]:
                ax.set_xlabel("Donor lag re recipient (ms)", fontsize=8)
            for ax in axs[:, 0]:
                ax.set_ylabel("Observed / expected", fontsize=8)
            fig.suptitle(f"Recipient-donor pairs with a short-lag dip, (1 + kappa(0-0.5 ms)) / "
                         f"(1 + kappa(0.5-2 ms)) < {DIP_RATIO}; grouped by insertion, then sorted by kappa(0-0.5 ms). Gray: CCG; purple: the recipient's own ACG (raw, before cleaning), same scale. "
                         f"Page {start // per_page + 1}", fontsize=9)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)
    sim = pd.DataFrame(sims)
    sim.to_csv(OUTDIR / "dip_pairs_similarity.csv", index=False)
    print(f"wrote {out}")
    print(sim[["similarity", "template_sim", "narrow", "asym"]].describe().round(3))
    summary(pairs, sim, p)


def summary(pairs, sim, p):
    """One page: what the dip pairs are, against the other near pairs."""
    pairs = pairs.reset_index()
    pairs["is_dip"] = dip_mask(pairs)
    tsim = []
    for pid, g in pairs[pairs.donor_kind == "near"].groupby("pid"):
        wf, ch = load_templates(pid)
        tsim += [dict(pid=pid, recipient=r, donor=d, tsim=template_similarity(wf, ch, r, d),
                      is_dip=k) for r, d, k in zip(g.recipient, g.donor, g.is_dip)]
    t = pd.DataFrame(tsim)
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    ax = axs[0]
    bins = np.linspace(-1, 1, 41)
    for k, c, lab in ((False, "0.4", "other near pairs"), (True, "#d95f02", "dip pairs")):
        v = t[t.is_dip == k].tsim.dropna()
        n, e = np.histogram(v, bins)
        ax.stairs(n / n.sum(), e, color=c, lw=1.6, label=f"{lab} (n = {len(v)})")
    ax.set_xlabel("Template similarity (cosine)")
    ax.set_ylabel("Proportion of near pairs")
    ax.set_title("a  Templates of dip pairs", loc="left")
    ax.legend(frameon=False, fontsize=7.5)
    ax = axs[1]
    n, e = np.histogram(sim.asym.dropna(), np.linspace(-1, 1, 41))
    ax.stairs(n / n.sum(), e, color="#d95f02", lw=1.6)
    ax.set_xlabel("Dip asymmetry, (right - left) / (right + left)")
    ax.set_ylabel("Proportion of dip pairs")
    ax.set_title("b  Symmetric or one-sided?", loc="left")
    ax = axs[2]
    edges = np.array([0, 10, 20, 30, 40, 60, 100, 200, 300, 600, 1000, 4000])
    pairs["sep_bin"] = pd.cut(pairs.depth_sep_um, edges, include_lowest=True)
    g = pairs.groupby("sep_bin", observed=True).is_dip.agg(["mean", "size"])
    mid = [iv.mid for iv in g.index]
    ax.plot(mid, g["mean"] * 100, "o-", color="#d95f02", ms=4)
    ax.set_xscale("log")
    ax.set_xlabel("Depth separation of recipient and donor (µm)")
    ax.set_ylabel("Pairs with a dip (%)")
    ax.set_title("c  Dips are a near-neighbor effect", loc="left")
    ax = axs[3]
    q = p[dip_mask(p)].merge(sim[["pid", "recipient", "donor", "template_sim", "asym"]],
                             on=["pid", "recipient", "donor"])
    for lab, m, c in (("template > 0.5", q.template_sim > 0.5, "#7570b3"),
                      ("template -0.5 to 0.5", q.template_sim.between(-0.5, 0.5), "0.4"),
                      ("template < -0.5", q.template_sim < -0.5, "#1b9e77")):
        gg = q[m].groupby("f_injected").passes.mean() * 100
        ax.plot(gg.index * 100, gg.values, color=c, lw=1.8,
                label=f"{lab} (n = {q[m].drop_duplicates(PAIR).shape[0]})")
    ax.axvline(10, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Injected contamination (%)")
    ax.set_ylabel("Units accepted by Sliding RP (%)")
    ax.set_ylim(-2, 102)
    ax.set_title("d  Dip pairs, by template similarity", loc="left")
    leg = ax.legend(handlelength=0, handletextpad=0, frameon=False, fontsize=7.5)
    for text, h in zip(leg.get_texts(), leg.legend_handles):
        text.set_color(h.get_color())
        h.set_visible(False)
    fig.tight_layout()
    plotstyle.save(fig, OUTDIR / "figures" / "dip_pairs_summary")
    L = [f"Dip pairs: {len(sim)} of {len(pairs)} pairs.",
         f"Template similarity (near pairs): dip median {t[t.is_dip].tsim.median():.2f}, "
         f"other near median {t[~t.is_dip].tsim.median():.2f}; "
         f"dip pairs with similarity > 0.9: {(t[t.is_dip].tsim > 0.9).mean():.1%}, > 0.7: "
         f"{(t[t.is_dip].tsim > 0.7).mean():.1%}, < -0.5: {(t[t.is_dip].tsim < -0.5).mean():.1%} "
         f"(other near pairs < -0.5: {(t[~t.is_dip].tsim < -0.5).mean():.1%}).",
         f"One-sided dips (|asymmetry| > 0.5): {(sim.asym.abs() > 0.5).mean():.1%}.",
         f"Similarity of the CCG to the recipient's own ACG: median {sim.similarity.median():.2f}.",
         "Dip prevalence by depth separation (um): " + ", ".join(
             f"{iv.left:g}-{iv.right:g}: {r['mean']:.1%} (n={r['size']})" for iv, r in g.iterrows())]
    (OUTDIR / "dip_pairs_numbers.txt").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
