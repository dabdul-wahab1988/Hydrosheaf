"""
Parameter and Observation definitions for PEST-style calibration.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class AdjustableParameter:
    """A parameter to be estimated."""

    name: str
    value: float
    lower_bound: float = -1.0e30
    upper_bound: float = 1.0e30
    log_transform: bool = False
    group: str = "default"
    prior_mean: Optional[float] = None
    prior_sigma: Optional[float] = None  # For regularization
    description: Optional[str] = None  # Human-readable description for documentation

    def to_internal(self, val: float) -> float:
        """Convert real value to optimization space (e.g. log)."""
        if self.log_transform:
            return np.log10(max(1e-10, val))
        return val

    def from_internal(self, internal_val: float) -> float:
        """Convert optimization space value to real value."""
        if self.log_transform:
            return 10.0**internal_val
        return internal_val


@dataclass
class Observation:
    """A target value to match."""

    name: str
    value: float
    weight: float = 1.0
    group: str = "default"
    time: Optional[float] = None  # For time-series plotting

    # Placeholder for the simulated value during run
    simulated: float = 0.0

    def residual(self) -> float:
        """Weighted residual (Simulated - Observed) * Weight."""
        return (self.simulated - self.value) * self.weight
