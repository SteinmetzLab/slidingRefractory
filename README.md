# slidingRefractory
Code to perform a new test of whether neurons have contaminated refractory periods, with a sliding window


## Python

### Installation
Install using pip:
```commandline
pip install slidingRP
```

Install using sources in development mode. First clone the repository.
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

# get the small test datasets from the github repository first
repo_path = "/home/ibladmin/Documents/PYTHON/int-brain-lab/slidingRefractory"
TEST_DATA_PATH = Path(repo_path).joinpath("test-data", "integration")

params = {'sampleRate': 30000, 'binSizeCorr': 1 / 30000}
spikes = pd.read_parquet(TEST_DATA_PATH.joinpath('spikes.pqt'))
table = metrics.slidingRP_all(spikes.times, spikes.clusters, **params)

assert np.allclose(pd.read_parquet(TEST_DATA_PATH.joinpath("rp_table.pqt")), pd.DataFrame(table), equal_nan=True)
```

### Contribute
#### Run unit tests
```commandline
 pytest python/test_*
```

#### Upload package
```commandline
rm -fR dist
rm -fR build
python setup.py sdist bdist_wheel
twine upload dist/*
```
