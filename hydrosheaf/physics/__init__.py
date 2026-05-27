"""Physics-based priors integration (optional).

This package provides adapters for importing external groundwater flow/transport
model outputs (e.g., MODFLOW/MODPATH) as priors over graph edges and travel times.
"""

from .priors import apply_physics_priors, load_physics_priors
from .modflow_head import (
    GridGeometry,
    build_grid_geometry_from_nodes,
    build_grid_geometry_from_params,
    cell_heads_to_sample_field,
    compute_head_gradient,
    map_gradient_to_nodes,
    parse_fhd,
    parse_hds,
    try_parse_cbc,
)

__all__ = [
    "apply_physics_priors",
    "build_grid_geometry_from_nodes",
    "build_grid_geometry_from_params",
    "cell_heads_to_sample_field",
    "compute_head_gradient",
    "GridGeometry",
    "load_physics_priors",
    "map_gradient_to_nodes",
    "parse_fhd",
    "parse_hds",
    "try_parse_cbc",
]
