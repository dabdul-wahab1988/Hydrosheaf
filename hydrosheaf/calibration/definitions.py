"""
Parameter and Observation definitions for PEST-style calibration.
"""

from dataclasses import dataclass
from typing import Optional, Dict, List, Any
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
    fixed: bool = False
    tied_to: Optional[str] = None

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


@dataclass
class CalibrationResult:
    """Standard schema for single-point calibration results (e.g. GLM)."""
    success: bool
    phi: float
    optimal_parameters: Dict[str, float]
    residuals: Optional[Dict[str, float]] = None
    phi_by_group: Optional[Dict[str, float]] = None
    covariance_matrix_path: Optional[str] = None
    correlation_matrix_path: Optional[str] = None
    sensitivity_table: Optional[Dict[str, Any]] = None
    identifiability: Optional[Dict[str, float]] = None


@dataclass
class EnsembleResult:
    """Standard schema for ensemble-based calibration results (e.g. IES)."""
    success: bool
    prior_parameters: Optional[Dict[str, List[float]]] = None       # parameter name -> realization values
    posterior_parameters: Optional[Dict[str, List[float]]] = None   # parameter name -> realization values
    simulated_observations: Optional[Dict[str, List[float]]] = None # observation name -> realization values
    phi_history: Optional[List[float]] = None
    posterior_forecast_summaries: Optional[Dict[str, Dict[str, float]]] = None
    parcov_path: Optional[str] = None
    obscov_path: Optional[str] = None
    localizer_path: Optional[str] = None
    restart_from: Optional[str] = None


@dataclass
class SensitivityResult:
    """Standard schema for sensitivity analysis (e.g. SEN)."""
    success: bool
    first_order_sobol: Optional[Dict[str, float]] = None
    total_sobol: Optional[Dict[str, float]] = None
    morris_elementary_effects: Optional[Dict[str, List[float]]] = None
    ranked_importance: Optional[List[str]] = None
    recommended_calibratable_subset: Optional[List[str]] = None


@dataclass
class OptimizationResult:
    """Standard schema for management optimization (e.g. MOU / OPT)."""
    success: bool
    pareto_front: Optional[List[Dict[str, float]]] = None
    optimal_decision_variables: Optional[Dict[str, float]] = None
    constraints_status: Optional[Dict[str, float]] = None


@dataclass
class DataAssimilationResult:
    """Standard schema for data assimilation (e.g. DA)."""
    success: bool
    cycles: List[Dict[str, Any]]

