# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - Unreleased

Major release aligning the Python implementation with the authoritative MATLAB
code and the 2026 manuscript. The expected-violations fix below **changes
computed confidence/contamination values**, so results differ from 1.x.

### Changed
- **Expected violations `Ve` corrected to the Llobet formulation**
  `Ve = 2·τr/D · Nc·(Nb + (Nc−1)/2)` (was `Nc·(Nb + Nc)`), accounting for
  contaminant–contaminant violations. Matches MATLAB and the manuscript.
- **ACG now computed in-package by `computeACG()`** (a `histdiff` equivalent),
  replacing `phylib.stats.correlograms`; reproduces the MATLAB ACG exactly.
- Contamination sweep is now `0.5:0.5:35` (70 levels); previously stopped at 34.5.
- Pass criterion standardised on `γ ≥ γ_thresh` (`>=`).
- `phylib` is no longer a core dependency (moved to the optional `[ibl]` extra).

### Added
- `compute_min_contamination()`: exact analytical minimum-contamination
  calculation, replacing the previous grid search.
- Optional family-wise (multiple-comparisons) correction across τr
  (`correction=True`, off by default).
- `forcePass` opt-in behaviour (default off).
- Parallel `slidingRP_all()` via `n_jobs`.
- Python simulation code reproducing the paper figures (Fig. 3 and
  cross-validation against MATLAB); MATLAB unit-test suite.
- PEP 621 `pyproject.toml` packaging and a GitHub Actions workflow publishing to
  PyPI on release via Trusted Publishing (OIDC).

### Fixed
- Hill formula corrected to `Nc·Nb`.
- `computeMatrix` half-bin offset in the observed-violations cumulative sum.
- `genST` (Python) now keeps the last in-window spike, matching MATLAB.

### Removed
- `setup.py` and `requirements.txt` (superseded by `pyproject.toml`).
- Stale IBL analysis scripts (`python/slidingRP/scripts/`), dead legacy code,
  and scratch artifacts.

## [1.0.0]

Initial tagged release.
