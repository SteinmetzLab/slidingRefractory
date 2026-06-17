# slidingRefractory — code review TODO

Review against `SlidingRPmanuscript_draft2026-03-06.pdf` (Methods: "The Sliding RP
metric", "The prior state of the art methods", "Spike train simulations").
Items are grouped by severity. Checkboxes are recommendations, not decisions —
add/override as you see fit.

Legend: 🔴 correctness / cross-implementation · 🟠 API / behavior · 🟢 cleanup.

## 🔴 High priority

> **DECISION (2026-06-13): Llobet is canonical.** Python `computeViol` has been
> changed to `2·refDur/D · Nc·(Nb + (Nc−1)/2)` (with `D = spikeCount/firingRate`)
> to match `matlab/computeViol.m` and the manuscript. **Still required:** (a) fix
> the test API (item 2) and run the metric under an env with `phylib`, then (b)
> regenerate the `EXPECTED` regression values in `tests/test_sliding_rp.py`. The
> code change was made but **not runnable/verified in the review environment**
> (phylib unavailable).

- [x] **1. Python and MATLAB disagree on the expected-violations formula `Ve`.**
  The manuscript and `matlab/computeViol.m` use the **Llobet** form
  `Ve = 2·τr/D · Nc·(Nb + (Nc−1)/2)`. But `python/.../metrics.py:computeViol`
  computes `expectedViol = contaminationRate · refDur · 2 · spikeCount`, which
  reduces to `Ve = 2·τr/D · Nc·Nt` (i.e. `Nc·(Nb + Nc)`, the "Hill 2-neuron"
  form — *not* Llobet). Verified numerically: at C=10% the Python `Ve` is ~5%
  larger, and near the pass/fail boundary the confidence differs materially
  (e.g. 0.59 vs 0.47 in an illustrative case). Consequences:
  - Python is systematically **more conservative** than MATLAB.
  - Python contradicts the Methods text, which states both implementations use
    Llobet.
  - The two reference implementations give different pass/fail decisions for
    units near threshold.
  Decide the canonical formula. If Llobet is canonical, change Python's
  `computeViol` to `2*refDur/recDur * Nc*(Nb + (Nc-1)/2)` and **regenerate the
  Python regression `EXPECTED` values** in `tests/test_sliding_rp.py` (they were
  produced by the current `Nc·Nt` code). If the IBL-production `Nc·Nt` form must
  stay for back-compat, then the manuscript text and `matlab/computeViol.m`
  should say so.

- [x] **2. The Python test suite (and `simulations.py`) call a stale API.**
  FIXED (2026-06-16): `slidingRP`/`slidingRP_all` now take an explicit
  `params=None` dict (MATLAB-like), so `slidingRP(st, params=params)` works.
  `tests/test_sliding_rp.py` updated + regenerated; `test_compute_rf.py` guarded
  with `importorskip`. Python suite now collects clean: **2 passed, 1 skipped**.
  (`simulations.py` still has the old API + the `confLlobet` NameError — covered
  by item 6; it's deprecated legacy.)

- [x] **3. `matlab/RPmetric_Classic.m` 'Hill' branch is internally inconsistent.**
  FIXED (2026-06-17): Hill `expectedViol` now uses `Nc*Nb` (was `Nc*(Nb+Nc)`),
  matching the manuscript and the estContam inversion. No test uses the 'Hill'
  path, and the paper figures plot only the Llobet comparison (`passPctLlobet*`),
  so figures are unaffected (the `passPctHill*` columns only change if
  runSimulations is re-run, and they are not plotted). Original finding below.

  ~~ORIGINAL:~~
  `passTest` uses `expectedViol = 2*RPdur/recDur * Nc.*(Nb + Nc)` (= `Nc·Nt`),
  but `estContam` in the same branch inverts `Ve = Nc·Nb` (the standard Hill
  quadratic), and the manuscript defines Hill as `Ve = 2τr·Nc·Nb/D`. The pass
  test should use `Nb` only: `... * expectedNc .* expectedNb`. As written, the
  'Hill' comparison is more lenient than textbook Hill. Check whether any
  published figure uses the 'Hill' (vs 'Llobet') comparison path before/after
  fixing.

---

## 🟠 Medium priority

- [x] **4. Python has no explicit recording-duration input.** FIXED (2026-06-16):
  Python `computeMatrix` now takes `params['recDur']` (default `max(spikeTimes)`,
  as MATLAB) and uses it directly in the Llobet `Ve`; `firingRate = n/recDur`.
  The old ACG-bin[1] rate estimate is gone.

- [x] **5. IBL-specific "force pass" hack is baked into the general
  `slidingRP`.** `metrics.py:225` hardcodes `(n_spikes_below2 == 0) and
  (firing_rate > 0.5)` → `pass_forced = True`, with a comment that it's tuned for
  IBL ~1 h recordings. This silently changes behavior for non-IBL data. Make it
  an opt-in parameter (off by default) and document it.

- [ ] **6. Build a Python simulation that matches the MATLAB one and reproduces
  the figures (DECISION 2026-06-17).** Rather than just deprecating the broken
  `simulations.py`, port `roth-et-al-2026/simulations/runSimulations.m` to Python
  (using the now-matching `slidingRP` + a Python `RPmetric_Classic` Hill/Llobet)
  and reproduce Fig 3 / S2 / Fig 4 from the Python output — an important extra
  cross-validation that the two sides agree. Plan: (a) clean Python `genST`
  (sample-consistent with MATLAB); (b) Python `RPmetric_Classic`; (c) a Python
  `runSimulations` writing the same table columns as `simDat.mat`; (d) Python
  figure scripts; (e) compare Python vs MATLAB pass-percentages on a shared
  parameter set. Replace the old broken `simulations.py`/`test_parfor.py` in the
  process.

- [x] **7. Pass/threshold comparison operator.** FIXED (2026-06-17): standardised
  on `>=` (per manuscript Outputs "at least equal to", and user decision). MATLAB
  `slidingRP` passTest and the corrected contamination find now use `>=`; Python
  fast-path pass uses `>=` (the corrected path's `pass_slidingRP_confmat` already
  did). Boundary-only effect; tests unaffected.

---

## 🟢 Low priority / cleanup

- [x] **8. `matlab/slidingRP_all.m:112` — code glued into a comment.** FIXED
  (2026-06-16): split `stCell = cell(numel(cids), 1);` onto its own line so the
  parfor preallocation is restored. Tests still 19/19.

- [x] **9. `genST` last-spike off-by-one between languages.** FIXED (2026-06-17):
  Python `simulations.py:genST` now keeps the last in-window spike (`...[-1] + 1`),
  matching MATLAB `st(1:find(st<duration,1,'last'))`.

- [x] **10. `compute_rf` flips a process-global warning filter.**
  `warnings.simplefilter("error", OptimizeWarning)` mutates global state as a
  side effect of calling the function. Wrap in a
  `with warnings.catch_warnings():` block instead.

- [x] **11. Resolve the `rp_min_val` / `timeOfLowestCont` status.** `metrics.py:193`
  has `# TODO this code could be deleted as rp_min_val is not informative`, but
  the manuscript (Outputs / Step 6) returns it as an RP-duration estimate. Decide
  whether it stays (and drop the comment) or goes (and update the manuscript).

- [x] **12. Corrected (Markov) variant untested and Python-less.** RESOLVED via
  item 25: the correction is now a `correction` flag on `computeMatrix`/`slidingRP`
  in BOTH languages (Python port verified bit-for-bit against MATLAB), and the
  MATLAB test suite has a `corrected` group. NOTE (still worth a manuscript
  cross-check): the boundary uses an adaptive `minPval` (min nominal p across
  valid bins) rather than the fixed `α = 1−γthresh` in the Methods eq for `c_k`;
  ported as-is from the original. Flag for the authors to confirm.

- [x] **13. metrics.py carries ~400 lines of commented-out legacy code.** The old
  `slidingRP_all`, `slidingRP`, `plotSlidingRP`, `fit_sigmoid`, `fitSigmoidACG`
  blocks are dead weight; git history preserves them. Remove for readability.

- [x] **14. Minor packaging:** `python/slidingRP/__init__.py` is empty — consider
  exposing `__version__` (setup.py is at 1.1.1) and the public API
  (`slidingRP`, `slidingRP_all`) for `from slidingRP import slidingRP`.

---

## 🔴 Discovered during testing (2026-06-13/14)

> **RESOLVED in code (2026-06-16):** Python now replicates MATLAB's `histdiff`
> binning exactly via a new `computeACG` (`floor((t_i−t_j)/binsize)`, verified
> bit-identical to histdiff on all 3 regression clusters, `sum|diff|=0`), and the
> phylib dependency was removed from `metrics.py`. Combined with the Llobet `Ve`
> (item 1) and explicit `recDur` (item 4), **Python now matches authoritative
> MATLAB bit-for-bit** on the regression clusters (e.g. cluster 275 confidence
> 27.4319% in both; was 47.8% under phylib). NB: both still bin float-seconds
> (the seconds↔samples quantization remains, identically, in both) — moving both
> to integer-sample binning is the future coordinated change (and the SI PR).

- [ ] **15. Python and MATLAB use different ACG engines → different violation
  counts.** MATLAB uses `histdiff` (cortex-lab/spikes); Python uses
  `phylib.stats.correlograms`. For regression cluster 275 they disagree:
  cumulative violations to 2 ms = **104 (histdiff) vs 99 (phylib)**, and the
  per-bin counts don't line up (histdiff first 8 bins `[0 0 0 0 0 0 0 1]` vs
  phylib `[0…0]`). (The apparent "2 ms boundary off-by-one" — idx 61 vs 60 — was
  a red herring: it's just MATLAB 1-indexing vs Python 0-indexing; both select
  the same bin. The `rp` bin-centre/edge conventions are identical.) Consequences,
  measured on the 3 regression clusters after the item-1 Llobet fix:

  | clu | old Python (Nc·Nt) | Python + Llobet | **MATLAB (authoritative)** |
  |-----|-----|-----|-----|
  | 167 maxConf | 0.80 | 0.75 | **0.53** |
  | 274 maxConf | 100.0 | 100.0 | **100.0** |
  | 275 maxConf | 53.9 | 47.8 | **27.4** |
  | 275 minCont | 15.0 | 16.0 | **19.0** |
  | 275 nBelow2 | 99 | 99 | **104** |
  | pass (167/274/275) | F/P/F | F/P/F | **F/P/F** |

  The Llobet fix moves Python toward MATLAB but a substantial gap remains because
  of the ACG-engine difference. **Pass/fail still agrees on these 3 clusters, but
  the confidence/contamination numbers do not, and near a real threshold the gap
  could flip decisions.** Implication: "make Python match the authoritative MATLAB"
  is *not* a one-liner — it requires reconciling the ACG computation and the
  bin-index conventions (`nViolShort`/2 ms boundary, `timeOfLowestCont`), then
  regenerating the Python `EXPECTED` values from MATLAB output. Needs a decision
  on how close Python must track MATLAB (bit-equivalent vs decision-equivalent).

  **Root cause (investigated 2026-06-16):** It is *not* a bin-edge/centre or
  coincidence-handling difference — it is how each engine maps seconds→bins:
  - **histdiff** bins the *exact float difference*: `bin = floor((t_i − t_j)·SR)`
    (faithful to the documented edges; counts each unordered pair once; exact
    coincidences `diff==0` excluded).
  - **phylib** *truncates each spike time to an integer sample first*
    (`(t·SR).astype(int64)`), then bins the *difference of truncated samples*:
    `bin = floor(t_i·SR) − floor(t_j·SR)` (same-sample pairs land in bin 0).
  - These differ by 0 or 1 per pair. The IBL `spikes.times` are **continuous
    sub-sample values** (drift-corrected/interpolated — fractional sample
    positions span [0,1), not grid-aligned), and *nearby* spikes share an almost
    identical sub-sample offset (slowly-varying drift correction). So near-pair
    differences land just below an integer (e.g. `7.9999…` samples): histdiff
    floors to 7, phylib (truncating both endpoints) recovers 8. **On cluster 275,
    essentially every sub-2 ms pair is binned one step higher by phylib.**
  - **Direction / effect:** phylib pushes short-lag violations to higher bins, so
    cumulative `obsViol(τr)` is *lower* at short τr → *higher* confidence →
    **phylib is less conservative** (275: 47.8 % vs MATLAB 27.4 %).
  - **Which is "correct"?** Neither, robustly — see RESOLUTION below. Each bins a
    different quantity (`floor((tᵢ−tⱼ)·SR)` vs `floor(tᵢ·SR)−floor(tⱼ·SR)`) and
    both are at the mercy of the seconds↔samples float round-trip; for
    sample-derived data histdiff is actually the *worse* approximation.

  **RESOLUTION (verified 2026-06-16):** The entire discrepancy is the
  **seconds↔samples float round-trip**, nothing else.
  - Fed *integer samples directly* (no `/SR`), **both engines are bit-exact**
    against the integer-lag ground truth (`sum|diff| = 0` for both phylib and
    histdiff). So neither has a binning bug per se.
  - Fed *seconds* `t = s/SR`, they diverge and **neither matches the true integer
    lags**: on a synthetic integer-sample train (first 80 bins), phylib was off by
    `sum|diff|=66`, histdiff by `406` — i.e. for sample-derived data **histdiff
    (floor of the float second-difference) is the *worse* of the two**, because
    `(sᵢ/SR − sⱼ/SR)·SR` underflows just below the integer and floors one bin low;
    phylib's truncate-both-then-subtract largely cancels that bias.
  - The user confirms the real data is recorded in **samples** then converted to
    seconds by `/SR`. So the true ACG is the integer-lag one, and the robust fix
    is unambiguous: **compute the ACG in integer-sample space** (pass sample
    indices / `round(t·SR)`, integer bin edges) in *both* languages. This makes
    MATLAB, Python, and SpikeInterface agree exactly and removes the float
    dependence entirely. (SpikeInterface already correlates in sample space — see
    item 17.) Note: `genST` (both languages) generates continuous-real spike times,
    not sample-quantized; regenerating sims in sample space would match the
    production convention (relates to todo #9).

- [ ] **16. Test-infrastructure status (verified in this env).** With `phylib`
  + `pytest` installed and `pip install -e .` done: **MATLAB 19/19 pass**
  (`histdiff`/`readNPY` from `githubDir`); **Python `test_sliding_rp.py` 2/2 fail**
  with `KeyError: 'sampleRate'` (item 2). Fixing item 2 unblocks the Python suite;
  the regression `EXPECTED` then needs regenerating (items 1 + 15).

---

## 🔴 Cross-package + PR review (2026-06-16)

- [ ] **17. SpikeInterface uses the `Nc·Nt` formula, NOT Llobet (verified from
  source).** `spikeinterface/metrics/quality/misc_metrics.py::_compute_violations`
  computes `expected_viol = firing_rate·C · ref_dur · 2 · spike_count` = `2·τr·Nc·Nt/D`
  — identical to the *old* Python `computeViol`, and its docstring says it was
  "adapted from slidingRefractory **1.0.0**". So the three implementations
  currently split: **MATLAB = Llobet**, **old Python = Nc·Nt**, **SpikeInterface =
  Nc·Nt**. Other SI differences worth knowing: default ACG `bin_size_ms = 0.25`
  (≈7 samples at 30 kHz, *not* 1/30000 s), `window_size_s = 1`, and it computes the
  correlogram in **integer-sample space** (`_compute_correlograms_numpy/numba` on
  `sample_index`) — i.e. SI already sidesteps the item-15 float issue. Decisions
  this forces: (a) if Llobet is canonical (it is, per item 1), a fix should be
  upstreamed to SpikeInterface, or the manuscript should note SI differs; (b) the
  manuscript states SI is a provided implementation — confirm what it actually
  ships and whether the 0.25 ms bin / sample-space ACG change results materially.

- [ ] **18. Review + adopt PR #6 ("Add ryan's sliding RP implementation",
  ryan-ressmeyer, open).** Adds `python/ryan_rp_methods.py` (1043 lines). Reviewed:
  it's high-quality, well-documented, vectorized, and **uses the Llobet-consistent
  expected-violations form** `Ve = C·(2−C)·F_r·τ·N_s` (= Llobet up to the
  negligible `Nc(Nc−1)` vs `Nc²` term), so it aligns with the authoritative MATLAB,
  not the `Nc·Nt` Python/SI form. Notable improvements worth pulling in:
  - **Analytical critical-contamination** (`compute_min_contam_props_analytical`):
    closed form via the Poisson↔χ² identity `λ_crit = chi2.ppf(conf, 2(V_o+1))/2`,
    then `C = 1 − sqrt(1 − k)`. Exact and much faster than the binary search or
    full tensor (paper-relevant: see item 19 speedups).
  - **Local firing-rate estimate** (`calc_local_firing_rate`, searchsorted-based):
    tracks the *local* rate so a unit active for 10 min of a 2 h recording isn't
    penalised — a genuine modelling improvement over `nACG[1]/n` (Python) / global
    `N/D`. NB: this *diverges from the manuscript's global-rate `D = N/F_r`
    assumption* — needs a conscious decision.
  - **`calc_rp_violations`** via k-th-order ISIs + searchsorted (fast, no histogram)
    and **`ref_acg_t_start = 0.25 ms`** default to avoid Kilosort duplicate-removal
    bias (a real-world correction the current code lacks).
  Caveats before merging: it's a standalone module (not integrated into the package
  API), has **no unit tests** (only a `__main__` benchmark/plot), adds a `tqdm`
  dependency, its CCG (`calc_ccgs`) bins float-second differences (same float
  caveat as item 15 unless moved to sample space), defaults to `confidence=0.95`
  (vs 0.90) and a log-spaced RP grid, and returns "min rejectable contamination"
  rather than pass/fail. **Action:** integrate into the package with tests, decide
  on the local-FR and `ref_acg_t_start` changes, and **port the chosen approach to
  MATLAB too** (user wants parity).

  **DECISION (2026-06-16, user):** adopt ONLY the **analytical critical-
  contamination** calculation (drop local-FR — intentionally excluded from MATLAB;
  drop the ISI-based counting — ACG is the mathematically correct route). Use the
  analytical C_min by default *whenever the full confidence matrix is not
  requested*, in both languages. VALIDATED (2026-06-16): exact-Llobet inversion
  `C = [(N−0.5) − sqrt((N−0.5)² − λ·D/τ)]/N` with `λ = chi2.ppf(conf, 2(V_o+1))/2`
  agrees with the grid matrix within grid resolution (cluster 275: analytical
  18.79% vs grid 19.0%; 274: 0.456% vs grid-floor 0.5%; 167: NaN both). NEXT:
  wire it into `slidingRP` (Python: always, since it never returns the matrix;
  MATLAB: when `nargout < confMatrix` / matrix not requested), then re-baseline
  both regression suites to the continuous C_min and cross-check Python==MATLAB.

  **DONE (2026-06-17):** Wired into `slidingRP` in both languages
  (`computeMinContamination.m` / `compute_min_contamination`); confidence/pass
  unchanged, `contamination`/`timeOfLowestCont` now analytical/continuous. The
  full matrix is built only on demand (MATLAB `nargout>=6`; Python via
  `computeMatrix`). MATLAB and Python agree bit-for-bit (274: 0.455893%, 275:
  18.791268%). Regression suites re-baselined; MATLAB 32 tests + Python 3 pass.
  Committed (689d05c).

---

## 🟠 Requested 2026-06-17

- [x] **25. Fold the correction into `slidingRP` as an option (off by default).**
  DONE (2026-06-17): removed `slidingRPCorrected.m` / `computeMatrixCorrected.m`;
  `computeMatrix`/`slidingRP` now take a `correction` flag (default off) in BOTH
  languages (the FWER Markov DP was ported to Python). Fig-S3 script + tests
  updated. Verified Python==MATLAB bit-for-bit with correction on (corrected
  confMatrix for cluster 275 matches to 3e-13; 275 corrected: conf 27.43%,
  cont 25%). Also aligned the Python `cont` grid to 0.5:0.5:35 (70 levels;
  Python previously stopped at 34.5, missing the 35% level). MATLAB 32 + Python
  4 tests pass.

## 🟢 Broader review dimensions (requested 2026-06-16)

- [ ] **19. Performance / parallelization.** (a) Adopt PR #6's analytical
  critical-contamination to replace the per-cell `computeViol` loop / full matrix
  where only pass/fail or `C_min` is needed (item 18). (b) Python `slidingRP_all`
  is a serial Python `for` over clusters — no parallelism (MATLAB has `parfor` via
  `slidingRP_all.m`); consider `joblib`/multiprocessing and/or vectorizing the ACG.
  (c) Profile: PR #6 notes ACG computation dominates, so the high-value win is a
  fast sample-space ACG, not micro-optimizing the confidence step.

- [~] **20. MATLAB ↔ Python usage parity.** PARTLY DONE (2026-06-16): Python now
  takes a `params` dict (like the MATLAB struct) with explicit `recDur`; the ACG
  and Llobet `Ve` match. STILL TODO: harmonise the return schema (MATLAB returns
  `[passTest, confidence, contamination, timeOfLowestCont, ...]`; Python returns
  `(max_conf, min_cont, rp_min_val, n_below2, fr, pass, pass_forced)`), the
  `pass_forced` IBL hack (item 5), and parameter naming, then document the mapping.

- [x] **21. Verify `roth-et-al-2026/` reproduces the paper figures.**
  **VERIFIED (2026-06-16)** against PDF figure renders (`_figcheck/`):
  - ✅ **Fig 3a–f + S2a–c** — `plotFig3AndS2.m` from committed `simDat.mat`
    (nSim=1000): curve shapes/orderings match (Sliding RP steepens with FR/RP/
    rec-dur; Llobet over-passes at low FR and over-rejects short RP).
  - ✅ **Fig 4a–f** — `plotFig4.m` from `simDat.mat` (teal confidence family +
    red Llobet pass-curves and ROC panels match). NB: Fig 4g is generated by
    `minimum_passing_FR.py`, and 4h,i by the IBL pipeline — not in `plotFig4.m`.
  - ✅ **Fig 1b,c (mouse)** — `plotFig1.m` cell 1 from the committed CSV: IBL
    histogram + median/CI ordering match (thalamus ~1.5–1.8 ms < cortex ~2.1–2.7
    < hippocampus ~2.5–3.2). Fig 1a (atlas slice + example ACGs) and the macaque
    panels need external data (`D:\Horwitz`, session files) — not reproducible here.
  - ✅ **Fig S3** — `runMultCompareCorrectionAndPlotS3.m` runs its own sims;
    a reduced `nSim=60` smoke run reproduced the qualitative result (corrected
    metric ≈10% pass at 10% contam for RP=10 ms; over-corrects to ~5% for RP=1.5
    ms). Full `nSim=1000` is slow (~2 h serial; the corrected/Markov metric
    dominates — that's why the script uses `parfor`).
  - [x] **Fig 2** — FIXED (2026-06-17): added a working Python `plotSlidingRP`
    to `metrics.py` (3 panels: ACG, confidence matrix + 90% iso-contour, and
    confidence traces; supports `plotXs`/`inputAxes`/`plotExtraContours`), and
    modernised `Fig2_metricExplanation.py` (`plotFig2` now uses `computeACG`/
    Llobet `Ve`, no phylib; imports `genST` explicitly; `os.path.join` paths).
    `runSaveFig2` runs end-to-end and produces a Fig-2-like figure (schematic +
    trace/matrix), verified visually.
  - Not locally reproducible (missing data, expected): Fig 1a/macaque, Fig 4g/h/i.

- [ ] **22. Documentation completeness.** Ensure each public function (both
  languages) has a complete docstring/help with params, returns, units, and a
  usage example; add/refresh a top-level README + a short "how to run the metric"
  example for Python and MATLAB; document defaults and the IBL-specific behaviors.

- [~] **23. Commenting clarity.** MATLAB library pass DONE (2026-06-16):
  `slidingRPCorrected.m` (added full help + fixed the function-name bug: was
  declared `slidingRP`), `computeMatrixCorrected.m` (added help header), `genST.m`
  (expanded help + labelled unreachable test cells), removed the `% TODO comment`
  placeholder. Core files (`computeViol`, `computeMatrix`, `slidingRP`,
  `slidingRP_all`, `RPmetric_Classic`, `plotSlidingRP`) were already well
  documented. STILL TODO: Python side — stale comments (the `# TODO this code
  could be deleted` at item 11, the half-bin "off slightly from matlab" in
  `computeMatrix`), explain the firing-rate-from-ACG trick and `pass_forced`
  rationale; and the research scripts in `roth-et-al-2026/` / `matlab/simulations`
  are only lightly commented.

- [ ] **24. Prune deprecated/useless Python.** Beyond item 13 (the ~400 commented
  lines in `metrics.py`) and item 6 (`simulations.py`): audit `python/plotSimsTest.py`,
  `python/simulationsFunctions.py`, `python/examples/`, and the many one-off
  `scripts/` (e.g. `explore_90conf_issues.py`, `scratch_ampFrACGFit.py`,
  `compare2msVersionIBLdata.py`) — move genuinely-needed analysis into
  `roth-et-al-2026/` and delete dead scratch. Clarify what is library vs paper vs
  throwaway.

---

## Things that look correct (verified, no action)

- MATLAB `computeViol.m` matches the manuscript Llobet `Ve` and `γ = 1 −
  poisscdf(Vo, Ve)`.
- Python and MATLAB pass-logic are equivalent in structure (min-contamination
  reachable at γthresh ⇔ max-confidence at the Cthresh row), given monotonicity
  of confidence in C — modulo the `>`/`>=` nit (item 7).
- The MATLAB simulation `genST.m` rate correction `rSim = rate/(1−rp·rate)` and
  ISI generation match the manuscript's `I = P + exprnd(µ)`, `µ = (1−PR)/R`.
- `cont = 0.5:0.5:35`, `rp` over 0–10 ms at 1/30000 s, and `rpReject = 0.5 ms`
  defaults match the Methods on both sides.
