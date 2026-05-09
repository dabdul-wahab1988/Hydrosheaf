"""Global Sensitivity Analysis for Hydrosheaf."""

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple
import copy
import numpy as np

from ..config import Config
from ..log import get_logger

logger = get_logger("uncertainty.sensitivity")

@dataclass
class SensitivityResult:
    parameter: str
    perturbation: float
    residual_change: float
    phase_stability: float  # Fraction of phases that remained selected
    phase_shifts: List[Tuple[str, str]] # (original, new) for any additions/removals

def analyze_sensitivity(
    solve_func,
    base_inputs: Dict[str, object],
    config: Config,
    parameters: List[str] = ["age_years", "cl"],
    scale: float = 0.1
) -> List[SensitivityResult]:
    """
    Perform a simple global sensitivity analysis by perturbing key inputs.
    solve_func: a callable that returns a solver result with residual and reaction extents.
    """
    if not config.sensitivity_analysis_enabled:
        return []

    results = []
    
    # Run baseline
    base_output = solve_func(**base_inputs)
    base_residual = getattr(base_output, "residual_norm", 0.0)
    base_phases = set()
    if hasattr(base_output, "reaction_fit"):
        # Assuming we can get labels from config elsewhere or passed in
        # For simplicity, we just use indices if labels aren't available
        extents = base_output.reaction_fit.extents
        base_phases = {i for i, e in enumerate(extents) if abs(e) > 1e-6}

    for param in parameters:
        # Perturb +scale
        perturbed_inputs = copy.deepcopy(base_inputs)
        if param in perturbed_inputs:
            val = perturbed_inputs[param]
            if isinstance(val, (int, float)):
                perturbed_inputs[param] = val * (1.0 + scale)
            
            p_output = solve_func(**perturbed_inputs)
            p_residual = getattr(p_output, "residual_norm", 0.0)
            p_phases = set()
            if hasattr(p_output, "reaction_fit"):
                extents = p_output.reaction_fit.extents
                p_phases = {i for i, e in enumerate(extents) if abs(e) > 1e-6}
            
            # Phase stability
            stability = 0.0
            if base_phases:
                common = base_phases.intersection(p_phases)
                stability = len(common) / len(base_phases)
            elif not p_phases:
                stability = 1.0
            
            results.append(SensitivityResult(
                parameter=param,
                perturbation=scale,
                residual_change=(p_residual - base_residual) / (base_residual + 1e-9),
                phase_stability=stability,
                phase_shifts=[] # TODO: Map indices to labels if possible
            ))

    return results
