"""Physics-based priors integration (optional).

This package provides adapters for importing external groundwater flow/transport
model outputs (e.g., MODFLOW/MODPATH) as priors over graph edges and travel times.
"""

from .priors import apply_physics_priors, load_physics_priors

__all__ = ["apply_physics_priors", "load_physics_priors"]

