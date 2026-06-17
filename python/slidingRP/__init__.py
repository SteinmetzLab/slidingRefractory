"""slidingRP: a tuning-free spike-sorting quality metric from refractory-period
violations (Roth, Chapuis, Winter et al., 2026).

See slidingRP.metrics for the implementation. The MATLAB implementation in
matlab/ is the authoritative reference; this package matches it.
"""

__version__ = "1.1.1"

from slidingRP.metrics import (
    slidingRP,
    slidingRP_all,
    computeMatrix,
    computeViol,
    computeACG,
    compute_min_contamination,
    compute_rf,
)

__all__ = [
    "slidingRP",
    "slidingRP_all",
    "computeMatrix",
    "computeViol",
    "computeACG",
    "compute_min_contamination",
    "compute_rf",
    "__version__",
]
