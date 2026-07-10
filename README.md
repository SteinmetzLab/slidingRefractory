# slidingRefractory
Non-parametric code to perform a new test of whether neurons have contaminated refractory periods.

See Roth, Chapuis, Winter et al. (2026) for the full description of the method.

## Python

### Installation
Install using pip:
```commandline
pip install slidingRP
```

Install in development mode (clone the repository first):
```commandline
cd slidingRefractory
pip install -e .
```

### Minimal working example

```python
from pathlib import Path
import numpy as np
import pandas as pd
from slidingRP import metrics

# All clusters at once. Options go in a `params` dict; recDur (recording
# duration, s) is recommended — it defaults to max(spikeTimes).
spikes = pd.read_parquet(Path('test-data/integration') / 'spikes.pqt')
params = {'sampleRate': 30000, 'recDur': float(spikes.times.max())}
table = metrics.slidingRP_all(spikes.times.values, spikes.clusters.values,
                              params=params, n_jobs=1)   # n_jobs>1 parallelises

# One cluster
st = spikes.times.values[spikes.clusters.values == 0]
(max_conf, min_cont, rp_min_val, n_spikes_below2,
 firing_rate, passes, pass_forced) = metrics.slidingRP(st, params=params)
```

Per-cluster outputs (in both `slidingRP` and the `slidingRP_all` table):
- `max_confidence` / `max_conf` — max confidence (%) that contamination is below
  `cont_thresh` (default 10%).
- `value` / `passes` — pass/fail (`max_conf >= conf_thresh`, default 90%).
- `min_contamination` / `min_cont` — minimum confirmable contamination (%),
  computed analytically; NaN if > 35%.
- `rp_min_val` — tau_r (s) at which that minimum is reached (an RP estimate).
- `firing_rate` — spikes / recDur.

Options (via `params`): `recDur`, `sampleRate`, `binSizeCorr`, `correction`
(FWER multiple-comparisons correction, default off; Fig. S3), `forcePass`
(IBL-specific force-pass of low-rate units with zero violations, default off).
The Python results match the authoritative MATLAB implementation bit-for-bit on
identical spike trains. (The output *names/order* differ between languages — see
the `slidingRP` docstring for the mapping.)

### Run unit tests
```commandline
pytest python/slidingRP/tests/
```

### Build and upload package to PyPI
Packaging metadata lives in `pyproject.toml` (version, dependencies). Bump the
`version` there first, then:
```commandline
rm -fR dist build
python -m build
twine upload dist/*
```
(`python -m build` needs the `build` package: `pip install build twine`.)

---

## MATLAB

### Dependencies
The MATLAB implementation requires two external toolboxes on your MATLAB path:
- **[cortex-lab/spikes](https://github.com/cortex-lab/spikes)** — provides `histdiff` for ACG computation
- **[kwikteam/npy-matlab](https://github.com/kwikteam/npy-matlab)** — provides `readNPY` for loading test data

Clone both repositories and add them to your path, e.g.:
```matlab
addpath(genpath(fullfile(githubDir, 'spikes')))       % cortex-lab/spikes
addpath(genpath(fullfile(githubDir, 'npy-matlab')))   % kwikteam/npy-matlab
```

`githubDir` should be a function on your MATLAB path that returns the path to the folder containing your cloned repositories (including `spikes`, `npy-matlab`, and `slidingRefractory` itself). For example, create a file `githubDir.m` somewhere on your path:
```matlab
function d = githubDir()
d = 'C:/Users/you/github';  % folder that contains spikes/, npy-matlab/, slidingRefractory/
end
```

### Minimal working example

```matlab
% Load spike times and cluster IDs (times in seconds)
spikeTimes    = readNPY('spikes.times.npy');
spikeClusters = readNPY('spikes.clusters.npy');
spikeTimes    = double(spikeTimes) / 30000;  % convert samples to seconds if needed

% Run on all clusters with default parameters (contThresh=10%, confThresh=90%)
rpMetrics = slidingRP_all(spikeTimes, spikeClusters);

% Run on a single cluster
st = spikeTimes(spikeClusters == 42);
[passTest, confidence, contamination, timeOfLowestCont] = slidingRP(st);

% Visualise the result
plotSlidingRP(st);
```

Key output fields from `slidingRP_all`:
- `.confidence` — max confidence (%) that contamination is below threshold
- `.contamination` — minimum confirmable contamination (%); NaN if > 35%
- `.timeOfLowestCont` — estimated RP duration (seconds)

### Run unit tests

```matlab
% All tests (dependency-free tests run always; others skip if toolboxes absent)
results = runtests('matlab/tests/test_slidingRP.m')

% Only the dependency-free core tests
results = runtests('matlab/tests/test_slidingRP.m', 'Tag', 'computeViol')
```

The test suite is organised by tag:
- `computeViol` — formula and statistical logic; no external dependencies
- `computeMatrix` — output shape and range; requires `histdiff`
- `slidingRP` — synthetic spike train pass/fail; requires `histdiff`
- `slidingRP_all` — batch consistency; requires `histdiff`
- `regression` — real IBL test data; requires `histdiff` + `readNPY`

---

## Use of AI in building this repository

Claude Code (Anthropic) was used to:
- Refine and expand documentation across MATLAB functions and the README
- Search for and repair small bugs (e.g. an undefined variable in `computeMatrix.m`, an undefined variable in `plotFig3AndS2.m`)
- Produce the MATLAB unit test suite
- Create the simulation runner function and script (`runSimulations.m`, `script_runAndPlotSimulations.m`)

All code has been reviewed and verified by Nick Steinmetz, who takes responsibility for the content of the code.
