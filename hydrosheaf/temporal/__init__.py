"""
Temporal dynamics module for Hydrosheaf.

This module provides methods for handling time-series data along flow paths,
including:
- Time-series sample loading and interpolation
- Residence time estimation (cross-correlation, gradient, tracer decay)
- Temporal edge fitting with time-averaged parameters
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np


@dataclass
class TimeSeriesSample:
    """A chemical sample with timestamp."""

    sample_id: str
    node_id: str
    timestamp: datetime
    concentrations: List[float]  # length = len(ion_order), mmol/L

    # Optional metadata
    temperature_c: Optional[float] = None
    ph: Optional[float] = None
    ec_uS_cm: Optional[float] = None
    isotopes: Optional[Dict[str, float]] = None  # e.g., {"18O": -5.2, "2H": -35.0}


@dataclass
class TemporalNode:
    """A node with time-series data."""

    node_id: str
    samples: List[TimeSeriesSample] = field(default_factory=list)  # sorted by timestamp

    # Derived quantities (computed)
    mean_concentration: Optional[List[float]] = None
    std_concentration: Optional[List[float]] = None
    trend_coefficients: Optional[List[float]] = None  # linear trend per ion
    seasonal_amplitude: Optional[List[float]] = None


@dataclass
class TemporalEdgeResult:
    """Result of temporal edge fitting."""

    edge_id: str
    u: str
    v: str

    # Estimated residence time
    residence_time_days: float
    residence_time_method: str  # "gradient", "cross_correlation", "tracer_decay"
    transport_model: str  # "evap" or "mix"

    # Optional fields
    residence_time_uncertainty: Optional[float] = None

    # Time-averaged transport parameters
    gamma_mean: Optional[float] = None
    gamma_std: Optional[float] = None
    f_mean: Optional[float] = None
    f_std: Optional[float] = None

    # Time-series reaction extents
    reaction_extents_series: List[List[float]] = field(default_factory=list)
    reaction_extents_mean: List[float] = field(default_factory=list)
    reaction_extents_std: List[float] = field(default_factory=list)

    # Time points
    timestamps: List[datetime] = field(default_factory=list)

    # Fit quality
    total_residual_norm: float = 0.0
    per_time_residual: List[float] = field(default_factory=list)


__all__ = [
    "TimeSeriesSample",
    "TemporalNode",
    "TemporalEdgeResult",
]
