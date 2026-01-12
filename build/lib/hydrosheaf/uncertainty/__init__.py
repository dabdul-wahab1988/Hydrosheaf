"""
Uncertainty quantification module for Hydrosheaf.

This module provides methods for quantifying uncertainty in transport parameters
and reaction extents through:
- Residual bootstrapping
- Bayesian MCMC inference
- Monte Carlo error propagation
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


@dataclass
class UncertaintyResult:
    """Container for uncertainty quantification results."""

    method: str  # "bootstrap", "bayesian", "monte_carlo"

    # Transport parameter uncertainty
    gamma_mean: Optional[float] = None
    gamma_std: Optional[float] = None
    gamma_ci_low: Optional[float] = None  # 2.5th percentile
    gamma_ci_high: Optional[float] = None  # 97.5th percentile

    f_mean: Optional[float] = None
    f_std: Optional[float] = None
    f_ci_low: Optional[float] = None
    f_ci_high: Optional[float] = None

    # Reaction extent uncertainty
    extents_mean: List[float] = field(default_factory=list)
    extents_std: List[float] = field(default_factory=list)
    extents_ci_low: List[float] = field(default_factory=list)
    extents_ci_high: List[float] = field(default_factory=list)

    # Posterior samples (optional, for Bayesian)
    gamma_samples: Optional[np.ndarray] = None
    extents_samples: Optional[np.ndarray] = None

    # Convergence diagnostics (Bayesian)
    r_hat: Optional[Dict[str, float]] = None
    ess: Optional[Dict[str, float]] = None

    # Bootstrap samples
    n_resamples: int = 0


__all__ = [
    "UncertaintyResult",
]
