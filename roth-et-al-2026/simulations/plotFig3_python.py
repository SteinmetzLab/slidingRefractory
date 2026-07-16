"""Python reproduction of Fig 3 (Sliding RP performance on simulated data).

Runs the Python simulation (slidingRP.simulations.run_simulations) and plots the
same panels as the MATLAB roth-et-al-2026/simulations/plotFig3AndS2.m:
  a) percent pass vs contamination, by firing rate (RP=3 ms, 2 h, conf=90)
  b) percent pass vs contamination, by RP duration (FR=5 sp/s, 2 h, conf=90)
  c) Hill-Llobet (3 ms) comparison, by firing rate

This is a cross-check that the Python and MATLAB metrics agree: on identical
spike trains the two implementations give identical pass decisions (verified
separately), so these curves match the MATLAB figure up to simulation noise.

Run:  python plotFig3_python.py   (nSim is reduced by default for speed)
"""
import numpy as np
import matplotlib.pyplot as plt
from slidingRP.simulations import run_simulations, plot_pass_vs_contamination


def make_fig3(nSim=200, savepath=None, seed=0):
    contProps = np.array([0, 2, 4, 6, 7, 8, 8.5, 9, 9.5, 10,
                          10.5, 11, 11.5, 12, 13, 14, 16, 18, 20]) / 100

    # Panel a/c: vary firing rate at RP = 3 ms
    simFR = run_simulations(dict(
        nSim=nSim, baseRates=np.array([0.5, 1, 2, 5, 10, 20]),
        RPdurs=np.array([0.003]), recDurs=np.array([7200.0]),
        confThreshes=np.array([90]), contProps=contProps, contThresh=10),
        progress=True, seed=seed)

    # Panel b: vary RP duration at FR = 5 sp/s
    simRP = run_simulations(dict(
        nSim=nSim, baseRates=np.array([5.0]),
        RPdurs=np.array([1.5, 2, 3, 4, 5, 6]) / 1000, recDurs=np.array([7200.0]),
        confThreshes=np.array([90]), contProps=contProps, contThresh=10),
        progress=True, seed=seed + 1)

    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    plot_pass_vs_contamination(simFR, ax=axs[0], vary='total_rate', metric='passPct')
    axs[0].set_title('Sliding RP by firing rate\n(RP=3 ms, 2 h, conf=90)')
    plot_pass_vs_contamination(simRP, ax=axs[1], vary='RP_dur', metric='passPct')
    axs[1].set_title('Sliding RP by RP duration\n(FR=5 sp/s, 2 h, conf=90)')
    plot_pass_vs_contamination(simFR, ax=axs[2], vary='total_rate', metric='passPctLlobet3')
    axs[2].set_title('Hill-Llobet (3 ms) by firing rate')
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath + '.png', dpi=150)
        fig.savefig(savepath + '.svg')
    return fig, (simFR, simRP)


if __name__ == "__main__":
    make_fig3(nSim=200, savepath="fig3_python")
    plt.show()
