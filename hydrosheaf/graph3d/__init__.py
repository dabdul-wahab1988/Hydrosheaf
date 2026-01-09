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

__all__ = [
    "Node3D",
    "Edge3D",
    "LayeredAquiferSystem",
    "Network3D",
]
