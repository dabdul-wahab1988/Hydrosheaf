"""Global Sensitivity and Robustness Analysis for Hydrosheaf."""

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Any
import copy
import numpy as np

from ..config import Config
from ..log import get_logger

logger = get_logger("uncertainty.sensitivity")

@dataclass
class SensitivityResult:
    parameter: str
    method: str = "OAT"
    perturbation: float = 0.0
    residual_change: float = 0.0
    phase_stability: float = 0.0  # Fraction of phases that remained selected
    phase_shifts: List[Tuple[str, str]] = field(default_factory=list)

@dataclass
class RobustnessReport:
    """A comprehensive robustness report for a specific flow-path discovery."""
    base_minerals: List[str]
    n_trials: int
    # Probability of Inclusion (PI) for each mineral in the candidate library
    inclusion_probabilities: Dict[str, float]
    # Mean and Std of reaction extents across trials
    extent_means: Dict[str, float]
    extent_stds: Dict[str, float]
    # Overall confidence score (0.0 to 1.0)
    confidence_score: float

def analyze_sensitivity_mc(
    solve_func: Any,
    base_inputs: Dict[str, Any],
    config: Config,
    n_trials: int = 50,
    concentration_error_pct: float = 0.05,
    age_error_pct: float = 0.1,
) -> RobustnessReport:
    """
    Perform Monte Carlo Sensitivity Analysis (MCSA) by perturbing all inputs simultaneously.
    Returns a RobustnessReport with inclusion probabilities.
    """
    if not config.sensitivity_analysis_enabled:
        return RobustnessReport([], 0, {}, {}, {}, 0.0)

    # 1. Run Baseline
    base_output = solve_func(**base_inputs)
    base_labels = getattr(base_output, "z_labels", [])
    base_extents = getattr(base_output, "z_extents", [])
    base_minerals = [label for label, extent in zip(base_labels, base_extents) if abs(extent) > 1e-6]

    # 2. Run Trials
    trial_extents = [] # List of lists
    all_labels = base_labels

    for _ in range(n_trials):
        perturbed_inputs = copy.deepcopy(base_inputs)

        # Perturb concentrations (x_u, x_v if present)
        for key in ["x_u", "x_v"]:
            if key in perturbed_inputs:
                vec = np.array(perturbed_inputs[key])
                noise = np.random.normal(0, concentration_error_pct, size=vec.shape)
                perturbed_inputs[key] = (vec * (1.0 + noise)).tolist()

        # Perturb age
        if "residence_time_days" in perturbed_inputs and perturbed_inputs["residence_time_days"] is not None:
            tau = perturbed_inputs["residence_time_days"]
            noise = np.random.normal(0, age_error_pct)
            perturbed_inputs["residence_time_days"] = max(1.0, tau * (1.0 + noise))

        try:
            t_output = solve_func(**perturbed_inputs)
            trial_extents.append(getattr(t_output, "z_extents", []))
        except Exception:
            continue

    if not trial_extents:
        return RobustnessReport(base_minerals, 0, {}, {}, {}, 0.0)

    # 3. Aggregate Statistics
    extents_matrix = np.array(trial_extents)
    n_success = len(trial_extents)

    inclusion_probs = {}
    means = {}
    stds = {}

    for i, label in enumerate(all_labels):
        col = extents_matrix[:, i]
        is_included = np.abs(col) > 1e-6
        inclusion_probs[label] = float(np.sum(is_included) / n_success)
        means[label] = float(np.mean(col))
        stds[label] = float(np.std(col))

    # Confidence Score: Average PI of the minerals selected in the baseline
    conf = 0.0
    if base_minerals:
        conf = np.mean([inclusion_probs[m] for m in base_minerals])

    return RobustnessReport(
        base_minerals=base_minerals,
        n_trials=n_success,
        inclusion_probabilities=inclusion_probs,
        extent_means=means,
        extent_stds=stds,
        confidence_score=float(conf)
    )

def analyze_sensitivity_oat(
    solve_func: Any,
    base_inputs: Dict[str, Any],
    config: Config,
    parameters: List[str] = ["residence_time_days"],
    scale: float = 0.1
) -> List[SensitivityResult]:
    """Simple One-At-a-Time sensitivity analysis."""
    if not config.sensitivity_analysis_enabled:
        return []

    results = []
    base_output = solve_func(**base_inputs)
    base_residual = getattr(base_output, "anomaly_norm", 0.0)
    base_labels = getattr(base_output, "z_labels", [])
    base_extents = getattr(base_output, "z_extents", [])
    base_phases = {label for label, extent in zip(base_labels, base_extents) if abs(extent) > 1e-6}

    for param in parameters:
        if param not in base_inputs or base_inputs[param] is None:
            continue

        perturbed_inputs = copy.deepcopy(base_inputs)
        val = perturbed_inputs[param]
        if isinstance(val, (int, float)):
            perturbed_inputs[param] = val * (1.0 + scale)

            p_output = solve_func(**perturbed_inputs)
            p_residual = getattr(p_output, "anomaly_norm", 0.0)
            p_labels = getattr(p_output, "z_labels", [])
            p_extents = getattr(p_output, "z_extents", [])
            p_phases = {label for label, extent in zip(p_labels, p_extents) if abs(extent) > 1e-6}

            stability = 0.0
            if base_phases:
                common = base_phases.intersection(p_phases)
                stability = len(common) / len(base_phases)
            elif not p_phases:
                stability = 1.0

            results.append(SensitivityResult(
                parameter=param,
                method="OAT",
                perturbation=scale,
                residual_change=(p_residual - base_residual) / (base_residual + 1e-9),
                phase_stability=stability
            ))

    return results
