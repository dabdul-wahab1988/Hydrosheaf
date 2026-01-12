"""
Reactive transport integration module for Hydrosheaf.

This module provides forward validation by coupling Hydrosheaf inverse results
with reactive transport simulators:
- PHREEQC kinetics
- Forward-inverse consistency metrics
- Thermodynamic path validation
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class KineticParameters:
    """Rate law parameters for a single reaction."""

    reaction_name: str
    rate_constant: float  # mol/m²/s (or mol/L/s for homogeneous)
    surface_area: float  # m²/L (for mineral reactions)
    reaction_order: float = 1.0
    activation_energy: float = 0.0  # J/mol (for temperature dependence)

    # Optional: Arrhenius temperature correction
    # k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))
    reference_temp_k: float = 298.15


@dataclass
class ReactiveTransportResult:
    """Result of forward reactive transport simulation."""

    edge_id: str
    simulator: str  # "phreeqc_kinetic", "mt3dms", "pflotran"

    # Input from inverse model
    inverse_extents: List[float]
    inverse_residence_time_days: float

    # Forward simulation output
    forward_x_v: List[float]
    forward_si_trajectory: Dict[str, List[float]] = field(default_factory=dict)  # SI vs time for each mineral
    forward_time_steps: List[float] = field(default_factory=list)  # time points in days

    # Consistency metrics
    rmse: float = 0.0
    nse: float = 0.0
    pbias: float = 0.0
    thermodynamic_consistent: bool = False

    # Detailed comparison
    per_ion_error: List[float] = field(default_factory=list)
    per_ion_bias: List[float] = field(default_factory=list)


@dataclass
class ValidationSummary:
    """Network-level forward validation summary."""

    n_edges_validated: int
    n_edges_consistent: int
    mean_rmse: float
    mean_nse: float

    # Edges with poor consistency (for review)
    inconsistent_edges: List[str] = field(default_factory=list)
    edge_results: Dict[str, ReactiveTransportResult] = field(default_factory=dict)


__all__ = [
    "KineticParameters",
    "ReactiveTransportResult",
    "ValidationSummary",
]
