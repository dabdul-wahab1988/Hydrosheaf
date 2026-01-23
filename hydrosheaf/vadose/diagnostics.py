from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class FlowRegimeDiagnostic:
    advective_velocity_m_day: float
    tracer_velocity_m_day: float
    ratio_tracer_advective: float
    interpretation: str


def diagnose_flow_regime(
    profile_depth_m: float,
    advective_tt_days: float,
    tracer_tt_days: float
) -> FlowRegimeDiagnostic:
    """Compare advective (Darcy) and tracer-derived travel times to diagnose flow regime.
    
    Args:
        profile_depth_m: Thickness of the vadose zone.
        advective_tt_days: Mean travel time calculated from Richard's equation (q/theta).
        tracer_tt_days: Observed or calibrated travel time of peak tracer breakthrough.
    
    Returns:
        Diagnostic object with velocities and interpretation.
    """
    v_adv = float(profile_depth_m) / max(1e-3, float(advective_tt_days))
    v_tr = float(profile_depth_m) / max(1e-3, float(tracer_tt_days))
    ratio = v_tr / max(1e-9, v_adv)
    
    interp = "Matrix flow dominant."
    if ratio > 50.0:
        interp = "Extreme preferential flow (fracture/macropore dominant)."
    elif ratio > 10.0:
        interp = "Strong preferential flow; dual-domain model required."
    elif ratio > 2.0:
        interp = "Moderate heterogeneity; enhanced dispersion (high CV) recommended."
    
    return FlowRegimeDiagnostic(v_adv, v_tr, ratio, interp)
