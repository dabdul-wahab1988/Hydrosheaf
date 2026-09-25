"""
3D flow network module for Hydrosheaf.

This module extends Hydrosheaf from 2D horizontal flow networks to fully 3D
representation with vertical discretization, multi-aquifer systems, and
aquitard leakage.
"""

from .types_3d import (
    Node3D,
    Edge3D,
    LayeredAquiferSystem,
    Network3D,
)
from .build_3d import (
    attach_bedrock_elevation_raster,
    attach_configured_geophysics,
    attach_geophysical_observations,
    compute_geophysical_transit_time,
    geophysical_barrier_check,
    load_geophysical_observations,
)

__all__ = [
    "Node3D",
    "Edge3D",
    "LayeredAquiferSystem",
    "Network3D",
    "load_geophysical_observations",
    "attach_geophysical_observations",
    "attach_bedrock_elevation_raster",
    "attach_configured_geophysics",
    "geophysical_barrier_check",
    "compute_geophysical_transit_time",
]
