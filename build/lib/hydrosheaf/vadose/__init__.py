"""Vadose-zone (unsaturated) physics module.

This module provides open, 1D vertical-column unsaturated flow physics intended to
support Hydrosheaf's residence-time and connectivity inference between the vadose
zone (e.g., lysimeters) and groundwater.

The initial implementation focuses on:
  - 1D mixed-form Richards equation with van Genuchten-Mualem constitutive laws
  - Atmospheric top flux (rain/irrigation minus surface evaporation)
  - Root-zone sink term for transpiration (Feddes-style stress response)
  - Recharge time series at the bottom boundary (deep percolation)
  - Travel-time (TTD) summaries derived from advective pore velocities
"""

from .contracts import (
    VadoseForcingSample,
    VadoseLayer,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .io import (
    load_forcing_csv,
    load_links_csv,
    load_profile,
)
from .run import (
    VadoseEdgePrior,
    VadoseRunResult,
    build_vadose_edge_priors,
    run_vadose_profile,
)
from .calibrate import (
    ThetaObservation,
    VadoseCalibrationResult,
    calibrate_ks_and_kc,
    load_theta_observations_csv,
    scale_profile_ks,
)
from .nitrate import (
    BreakthroughPoint,
    BreakthroughSummary,
    NitrateLoadingSample,
    load_no3_loading_csv,
    predict_no3_breakthrough,
)

__all__ = [
    "VadoseForcingSample",
    "VadoseLayer",
    "VadoseLinksRow",
    "VadoseProfile",
    "VadoseRunConfig",
    "VadoseEdgePrior",
    "VadoseRunResult",
    "load_forcing_csv",
    "load_links_csv",
    "load_profile",
    "build_vadose_edge_priors",
    "run_vadose_profile",
    "ThetaObservation",
    "VadoseCalibrationResult",
    "calibrate_ks_and_kc",
    "load_theta_observations_csv",
    "scale_profile_ks",
    "BreakthroughPoint",
    "BreakthroughSummary",
    "NitrateLoadingSample",
    "load_no3_loading_csv",
    "predict_no3_breakthrough",
]
