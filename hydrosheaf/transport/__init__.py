"""
Solute Transport Modeling Package for Hydrosheaf.

This package provides methods for simulating solute transport in the saturated zone,
including:
- 1D advection-dispersion transport via FloPy/MT3DMS
- Analytical ADE solutions (Green's functions)
- Coupling with vadose zone recharge
- First-order decay (denitrification)
"""

from .flopy_1d import (
    TransportResult,
    build_1d_transport_model,
    run_1d_transport,
    check_flopy_available,
)

from .coupling import (
    VadoseCouplingResult,
    couple_vadose_saturated,
)

from .analytical import (
    ade_1d_green,
    simulate_1d_ade
)

from .validation import TransportValidationProblem

__all__ = [
    "TransportResult",
    "build_1d_transport_model",
    "run_1d_transport",
    "check_flopy_available",
    "VadoseCouplingResult",
    "couple_vadose_saturated",
    "ade_1d_green",
    "simulate_1d_ade",
    "TransportValidationProblem",
]
