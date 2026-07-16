# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`slidingRP` is a Python (and MATLAB) package for assessing spike-sorting quality using a sliding refractory period metric. It tests whether a neuron's spike train has contaminated refractory periods without assuming a fixed refractory period duration.

## Commands

### Run Python tests
```bash
pytest python/test_*
```

### Run a single Python test file
```bash
pytest python/slidingRP/tests/test_sliding_rp.py
```

### Run MATLAB tests
```matlab
results = runtests('matlab/tests/test_slidingRP.m')
```
Or by tag (e.g. only the dependency-free tests):
```matlab
results = runtests('matlab/tests/test_slidingRP.m', 'Tag', 'computeViol')
```
MATLAB tests require `histdiff` (cortex-lab/spikes) and `readNPY` (kwikteam/npy-matlab) on the MATLAB path for tests beyond the `computeViol` group.

### Install in development mode
```bash
pip install -e .
```

### Build and upload package to PyPI
```bash
rm -fR dist && rm -fR build
python setup.py sdist bdist_wheel
twine upload dist/*
```

## Architecture

### Python package layout (`python/slidingRP/`)

**Core modules:**
- `metrics.py` — The central module. Contains `slidingRP_all()` (batch), `slidingRP()` (single cluster), `computeMatrix()`, `computeViol()`, and `compute_rf()`.
- `simulations.py` — Spike train simulations used for validation and paper figures.
- `loadSaveData.py` — Data I/O utilities.

**Sub-packages:**
- `data_access/` — Lab-specific loaders (IBL via `one.api`, Steinmetz lab).
- `elts/` — Transforms from ACGs to refractory periods for low-temporal-resolution data.
- `tests/` — Unit tests; test data lives in `test-data/unit/` and `test-data/integration/`.

**MATLAB implementation** in `matlab/` mirrors the Python logic. `slidingRP_all.m` is the entry point; `slidingRP.m` handles single clusters.

**Paper code** in `roth-et-al-2026/` contains figure generation scripts tied to the manuscript.

### Core algorithm

**Motivation:** RP durations vary significantly across brain regions (thalamus shorter than cortex/hippocampus) and species (macaque shorter than mouse), with many neurons having RPs < 2 ms. Prior methods that assume a fixed RP of 2–3 ms incorrectly reject clean short-RP units and are statistically underpowered at low firing rates.

**Data flow:**
```
Spike times + Cluster IDs
  → ACG via computeACG() (histdiff-equivalent, bin size 1/30000 s, window 0–10 ms)
  → For each (τ_r, C) pair in the confidence matrix:
      V_e = 2τ_r × N_c × (N_b + (N_c−1)/2) / D   [Llobet formulation]
      V_o = cumsum(ACG up to τ_r)
      confidence γ = 1 − Poisson.cdf(V_o, V_e)
  → Pass if γ ≥ γ_thresh for any τ_r > τ_min at C = C_thresh
```

The Llobet formulation for V_e (expected violations) accounts for contaminating spikes also producing RP violations with each other, not just with base-neuron spikes.

**Algorithm outputs (per unit):**
- `max_confidence` (γ_max): highest confidence achieved at C_thresh across all τ_r
- `min_contamination` (C_min): lowest contamination confirmable at γ_thresh; NaN if > 35%
- `rp_min_val`: τ_r at which C_min is achieved (an estimate of the unit's true RP)
- `n_spikes_below2`: spike count with ISI < 2 ms (legacy diagnostic)
- `pass`: True if γ_max ≥ γ_thresh

**Key parameters** (defaults match 30 kHz Neuropixels data):
- `sampleRate`: 30000 Hz
- `binSizeCorr`: 1/30000 s
- `conf_thresh` (γ_thresh): 90 (%)
- `cont_thresh` (C_thresh): 10 (%)
- `τ_min`: 0.5 ms (bins below this are excluded from pass/fail to avoid noise)

**Known limitation on false acceptance rate:** The standard implementation does *not* apply a multiple-comparisons correction across the τ_r sweep. As a result, the empirical false acceptance rate at threshold contamination is ~30% (not 10%) for units near the threshold. A Markov-chain correction exists (see manuscript Methods) but is omitted from the standard code due to computational cost and because it requires knowing the true RP duration.

**Comparison to Hill-Llobet (prior method):** The Hill-Llobet method uses a fixed assumed RP duration and a point estimate of contamination (not a statistical test). It fails in two ways that Sliding RP avoids: (1) it always passes units with zero observed violations regardless of firing rate (no statistical power check), and (2) it incorrectly rejects all clean units whose true RP is shorter than the assumed duration.

### RP duration estimation (`compute_rf`)

Used to characterize RP durations for comparison across datasets (Fig. 1), not part of the pass/fail metric itself. Fits a 4-parameter sigmoid to the ACG and returns the time at which firing recovers to 10% of baseline. Steps: median-filter ACG → find first valid peak → fit sigmoid to ACG segment between trough and peak → extract RP_est.

### Dependencies

Required: `numpy`, `scipy`, `matplotlib`, `colorcet`, `statsmodels`

Optional (`[ibl]` extra): `ONE-api`, `ibllib`, `phylib` (IBL-specific data access in `data_access/`, `loadSaveData.py`, `elts/`). The core metric computes its ACG with the in-package `computeACG()` and no longer depends on `phylib`.
