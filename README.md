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

TEST_DATA_PATH = Path('test-data/integration')
params = {'sampleRate': 30000, 'binSizeCorr': 1 / 30000}
spikes = pd.read_parquet(TEST_DATA_PATH / 'spikes.pqt')
table = metrics.slidingRP_all(spikes.times, spikes.clusters, **params)
```

### Run unit tests
```commandline
pytest python/test_*
```

### Upload package to PyPI
```commandline
rm -fR dist && rm -fR build
python setup.py sdist bdist_wheel
twine upload dist/*
```

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
