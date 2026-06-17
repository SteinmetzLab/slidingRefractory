# -*- coding: utf-8 -*-
"""
Spike-train simulations for validating the Sliding RP metric.

This is the Python twin of matlab/simulations + roth-et-al-2026/simulations:
- genST              : generate a simulated spike train (Poisson + hard RP).
- RPmetric_Classic   : Hill / Llobet point-estimate comparison metric.
- run_simulations    : sweep firing rate x RP x recording duration x
                       contamination x confidence, returning the same
                       pass-percentage columns as the MATLAB simDat table.
- plot_pass_vs_contamination : Fig 3 / S2 / Fig 4 style curves.

The Sliding RP metric itself (slidingRP) matches the authoritative MATLAB
implementation bit-for-bit on identical spike trains, so these simulations
reproduce the MATLAB figures (statistically) and serve as a cross-check.
"""
import numpy as np
from scipy import stats

from slidingRP.metrics import slidingRP, computeACG


def genST(rate, duration, rp=0, params=None):
    """Generate a simulated spike train (Poisson, optionally with a hard RP).

    Inter-spike intervals are rp + Exponential(mu), with mu rate-corrected so
    the realised total rate equals `rate` despite the dead time:
        rSim = rate / (1 - rp*rate);  mu = 1/rSim
    (matches matlab/simulations/genST.m and the manuscript Methods).

    Parameters
    ----------
    rate : float       desired total firing rate (spikes/s).
    duration : float   recording duration (s).
    rp : float         hard refractory period (s); default 0.
    params : dict       unused (kept for backward compatibility).

    Returns
    -------
    st : np.ndarray    spike times in [0, duration) seconds.
    """
    rSim = rate / (1 - (rp * rate))   # rate-corrected for the refractory period
    mu = 1 / rSim
    n = rate * duration
    isi = rp + np.random.exponential(mu, int(np.ceil(n * 2)))
    while np.sum(isi) < duration:
        isi = np.append(isi, np.random.exponential(mu))
    st = np.cumsum(isi)
    below = np.where(st < duration)[0]
    if below.size > 0:
        st = st[0:below[-1] + 1]   # keep the last in-window spike (matches MATLAB)
    else:
        st = np.array([])
    return st


def simulate_contaminated_train(base_rate, cont_rate, rec_dur, rp):
    """Base neuron (with RP) + contamination (no RP), concatenated and sorted."""
    st = genST(base_rate, rec_dur, rp)
    if cont_rate > 0:
        cont_st = genST(cont_rate, rec_dur, 0)
    else:
        cont_st = np.array([])
    return np.sort(np.concatenate((st, cont_st)))


def RPmetric_Classic(spikeTimes, params):
    """Hill / Llobet point-estimate contamination metric (prior art).

    Python twin of matlab/RPmetric_Classic.m. Unlike slidingRP, this uses a
    single fixed RP duration and a point estimate (not a statistical test).

    params keys: metricType ('Llobet' default | 'Hill'), contaminationThresh
    (%, default 10), recDur, RPdur (s, default 0.002), acgBinSize (default
    1/30000). For speed, a precomputed ACG may be passed via nACG + rp +
    spikeCount (as the MATLAB simulation does).

    Returns (passTest, estContam): passTest is obsViol <= expectedViol at RPdur;
    estContam is the inverted point estimate (proportion).
    """
    metricType = params.get('metricType', 'Llobet')
    contThresh = params.get('contaminationThresh', 10)
    RPdur = params.get('RPdur', 0.002)
    acgBinSize = params.get('acgBinSize', 1 / 30000)

    if 'nACG' in params:
        nACG = params['nACG']
        rp = params['rp']
        spikeCount = params['spikeCount']
        recDur = params.get('recDur', None)
        # match MATLAB sum(nACG(1:find(rp>RPdur,1))): inclusive of the first
        # bin whose centre exceeds RPdur (so nACG[:rpIdx+1] in 0-indexed Python).
        rpIdx = np.where(rp > RPdur)[0][0]
        obsViol = np.sum(nACG[:rpIdx + 1])
    else:
        spikeTimes = np.asarray(spikeTimes, dtype=np.float64)
        spikeCount = spikeTimes.size
        recDur = params.get('recDur', None)
        if recDur is None:
            recDur = float(np.max(spikeTimes)) if spikeCount else 0.0
        nbins = int(round(RPdur / acgBinSize))
        obsViol = int(np.sum(computeACG(spikeTimes, acgBinSize, nbins)))

    Nc = spikeCount * contThresh / 100
    Nb = spikeCount * (1 - contThresh / 100)

    # estContam is undefined (NaN) when observed violations exceed the formula's
    # valid range (very high contamination); passTest is unaffected.
    with np.errstate(invalid='ignore'):
        if metricType == 'Llobet':
            # contaminating spikes also violate with each other
            expectedViol = 2 * RPdur / recDur * Nc * (Nb + (Nc - 1) / 2)
            estContam = 1 - np.sqrt(1 - obsViol * recDur / (spikeCount ** 2 * RPdur))
        else:  # 'Hill': single other neuron, violations only with the base neuron
            expectedViol = 2 * RPdur / recDur * Nc * Nb
            estContam = 0.5 * (1 - np.sqrt(1 - 2 * obsViol * recDur / spikeCount ** 2 / RPdur))

    passTest = bool(obsViol <= expectedViol)
    return passTest, estContam


def default_sim_params():
    """Default sweep matching roth-et-al-2026/simulations/runSimulations.m."""
    return {
        'nSim': 1000,
        'RPdurs': np.array([1.5, 2, 3, 4, 5, 6]) / 1000,    # s
        'recDurs': np.array([0.5, 1, 2, 3]) * 3600,          # s
        'contProps': np.array([0, 2, 4, 6, 7, 8, 8.5, 9, 9.5, 10,
                               10.5, 11, 11.5, 12, 13, 14, 16, 18, 20]) / 100,
        'baseRates': np.array([0.5, 1, 2, 5, 10, 20]),       # spikes/s (total rate)
        'confThreshes': np.array([50, 60, 70, 80, 90]),
        'contThresh': 10,                                    # % acceptable contamination
    }


def run_simulations(params=None, progress=True, seed=None):
    """Sweep parameters and return pass-percentage arrays (Python twin of
    runSimulations.m). Returns a dict of flat arrays with the same columns as
    the MATLAB simDat table: total_rate, cont_prop, RP_dur, rec_dur, conf_level,
    passPct, passPctLlobet1_5/2/3, passPctHill1_5/2/3.
    """
    p = default_sim_params()
    if params:
        p.update(params)
    if seed is not None:
        np.random.seed(seed)

    nSim = p['nSim']
    contThresh = p['contThresh']
    cols = ['total_rate', 'cont_prop', 'RP_dur', 'rec_dur', 'conf_level',
            'passPct', 'passPctLlobet1_5', 'passPctLlobet2', 'passPctLlobet3',
            'passPctHill1_5', 'passPctHill2', 'passPctHill3']
    out = {c: [] for c in cols}

    rpDursClassic = [0.0015, 0.002, 0.003]

    total = (len(p['baseRates']) * len(p['contProps']) * len(p['RPdurs'])
             * len(p['recDurs']) * len(p['confThreshes']))
    done = 0
    for totalRate in p['baseRates']:
        for contProp in p['contProps']:
            baseRate = (1 - contProp) * totalRate
            contRate = contProp * totalRate
            for RPdur in p['RPdurs']:
                for recDur in p['recDurs']:
                    # generate all spike trains once per (rate, cont, RP, recDur);
                    # confidence threshold only changes the slidingRP decision.
                    passing = np.zeros((nSim, len(p['confThreshes'])))
                    llob = {d: np.zeros(nSim) for d in rpDursClassic}
                    hill = {d: np.zeros(nSim) for d in rpDursClassic}
                    for n in range(nSim):
                        combST = simulate_contaminated_train(baseRate, contRate, recDur, RPdur)
                        # slidingRP max-confidence (at the contamination threshold)
                        # does not depend on confThresh, so compute it once and
                        # threshold; >= matches the metric's pass rule.
                        max_conf = slidingRP(
                            combST, params={'recDur': recDur},
                            cont_thresh=contThresh)[0]
                        passing[n, :] = max_conf >= p['confThreshes']
                        # classic metrics share the 10 ms ACG (precomputed once)
                        nACG = computeACG(combST, 1 / 30000, 300)
                        rp = np.arange(0, 10 / 1000, 1 / 30000) + (1 / 30000) / 2
                        cp = {'recDur': recDur, 'contaminationThresh': contThresh,
                              'nACG': nACG, 'rp': rp, 'spikeCount': combST.size}
                        for d in rpDursClassic:
                            cp['metricType'] = 'Llobet'; cp['RPdur'] = d
                            llob[d][n] = RPmetric_Classic(None, cp)[0]
                            cp['metricType'] = 'Hill'
                            hill[d][n] = RPmetric_Classic(None, cp)[0]
                    for ci, conf in enumerate(p['confThreshes']):
                        out['total_rate'].append(totalRate)
                        out['cont_prop'].append(contProp)
                        out['RP_dur'].append(RPdur)
                        out['rec_dur'].append(recDur)
                        out['conf_level'].append(conf)
                        out['passPct'].append(passing[:, ci].mean() * 100)
                        out['passPctLlobet1_5'].append(llob[0.0015].mean() * 100)
                        out['passPctLlobet2'].append(llob[0.002].mean() * 100)
                        out['passPctLlobet3'].append(llob[0.003].mean() * 100)
                        out['passPctHill1_5'].append(hill[0.0015].mean() * 100)
                        out['passPctHill2'].append(hill[0.002].mean() * 100)
                        out['passPctHill3'].append(hill[0.003].mean() * 100)
                    done += len(p['confThreshes'])
                    if progress:
                        print(f'\r  simulations: {done}/{total}', end='', flush=True)
    if progress:
        print()
    return {c: np.array(v) for c, v in out.items()}


def plot_pass_vs_contamination(sim, ax=None, vary='total_rate',
                               fixed=None, metric='passPct'):
    """Fig 3 / S2 / Fig 4 style plot: pass percentage vs true contamination,
    one curve per value of `vary`, holding the other parameters at `fixed`.

    sim : dict from run_simulations.
    vary : column to make separate curves for (e.g. 'total_rate', 'RP_dur',
           'rec_dur', 'conf_level').
    fixed : dict of {column: value} to hold constant; sensible defaults used.
    metric : which pass column to plot (e.g. 'passPct', 'passPctLlobet3').
    """
    import matplotlib.pyplot as plt
    from statsmodels.stats.proportion import proportion_confint

    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    defaults = {'total_rate': 5, 'RP_dur': 0.003, 'rec_dur': 7200, 'conf_level': 90}
    fixed = {**defaults, **(fixed or {})}
    fixed.pop(vary, None)

    mask = np.ones(len(sim['cont_prop']), dtype=bool)
    for k, v in fixed.items():
        mask &= np.isclose(sim[k], v)

    vary_vals = np.unique(sim[vary][mask])
    # rough nSim recovery for CIs (pass% is a mean of nSim Bernoulli draws)
    for vv in vary_vals:
        m = mask & np.isclose(sim[vary], vv)
        order = np.argsort(sim['cont_prop'][m])
        x = sim['cont_prop'][m][order] * 100
        y = sim[metric][m][order]
        ax.plot(x, y, '.-', label=f'{vv:g}')
    ax.axvline(fixed.get('contThresh', 10), color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('True simulated contamination (%)')
    ax.set_ylabel('Percent pass (%)')
    ax.set_ylim(-5, 105)
    ax.legend(title=vary, frameon=False, fontsize=8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return ax
