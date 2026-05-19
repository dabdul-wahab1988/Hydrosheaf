"""
M6 Robustness and Sensitivity Analysis.
Quantifies discovery stability across regularization paths and input variance.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
import copy

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.config import default_config
from hydrosheaf.models.reactions import build_reaction_dictionary, fit_reactions
from hydrosheaf.uncertainty.sensitivity import analyze_sensitivity_mc, analyze_sensitivity_oat
from hydrosheaf.inference.edge_fit import fit_edge

# Paths
BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_DIR / "results"
ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
PSI_ACTIVE_MINERALS = [
    "calcite_closed",
    "dolomite_closed",
    "albite",
    "anorthite",
    "k_feldspar",
    "biotite",
    "chlorite",
    "halite",
    "gypsum",
    "fluorite",
    "pyrite_oxidation_aerobic",
]


def configure_common(config, active_minerals=None):
    """Apply the shared M6 chemistry configuration without reusing mutated state."""
    config.ion_order = list(ION_ORDER)
    config.weights = [1.0] * len(ION_ORDER)
    config.conservative_weights = [1.0] * len(ION_ORDER)
    config.measured_ions = list(ION_ORDER)
    if active_minerals is not None:
        config.active_minerals = list(active_minerals)
    return config

def run_regularization_path(residual, matrix, labels, weights, site, config, penalty_scales):
    """Scenario 3.1: The L1 Discovery Plateau."""
    print(f"\n--- Running Regularization Path ({site}) ---")
    lambdas = np.logspace(-4, -1, 20)
    path_data = []

    for l in lambdas:
        fit = fit_reactions(
            residual, matrix, weights, 
            lambda_l1=l, 
            lambda_l2=config.lambda_l2, 
            penalty_scales=penalty_scales
        )
        row = {"site": site, "lambda": l, "residual_norm": fit.residual_norm, "aicc": fit.aicc}
        for label, extent in zip(labels, fit.extents):
            row[label] = extent
        path_data.append(row)
    return path_data

def run_system_state_aicc(residual, config, weights):
    """Scenario 3.2: AICc-Based System Resolution (Open vs Closed)."""
    print("\n--- Running AICc Model Selection (Open vs Closed) ---")

    open_config = configure_common(copy.deepcopy(config), ["calcite_open", "dolomite_open"])
    m_open, l_open, _, p_open = build_reaction_dictionary(open_config)
    fit_open = fit_reactions(residual, m_open, weights, lambda_l1=0.01, penalty_scales=p_open, lambda_l2=open_config.lambda_l2)

    closed_config = configure_common(copy.deepcopy(config), ["calcite_closed", "dolomite_closed"])
    m_closed, l_closed, _, p_closed = build_reaction_dictionary(closed_config)
    fit_closed = fit_reactions(residual, m_closed, weights, lambda_l1=0.01, penalty_scales=p_closed, lambda_l2=closed_config.lambda_l2)

    print(f"  Open System   - AICc: {fit_open.aicc:.2f}, BIC: {fit_open.bic:.2f}")
    print(f"  Closed System - AICc: {fit_closed.aicc:.2f}, BIC: {fit_closed.bic:.2f}")

def run_uncertainty_cascade(config, site):
    """Scenario 3.3 & 3.4: Monte Carlo Phase Stability Index (PSI)."""
    print(f"\n--- Running Monte Carlo Phase Stability Analysis ({site}) ---")
    psi_config = configure_common(copy.deepcopy(config), PSI_ACTIVE_MINERALS)
    psi_config.sensitivity_analysis_enabled = True

    # Differentiate sites slightly
    if site == "Manu":
        x_u = [1.0, 0.5, 2.0, 0.2, 2.5, 2.0, 0.5, 0.1, 0.0, 0.01]
        x_v = [1.0, 0.5, 2.4, 0.2, 2.8, 2.1, 0.5, 0.1, 0.0, 0.01]
    else:
        x_u = [0.8, 0.4, 1.5, 0.3, 2.0, 1.0, 0.2, 0.5, 0.1, 0.05]
        x_v = [1.2, 0.6, 2.0, 0.3, 2.5, 1.1, 0.3, 0.6, 0.1, 0.05]

    base_inputs = {
        "x_u": x_u,
        "x_v": x_v,
        "config": psi_config,
        "residence_time_days": 3650.0
    }

    report = analyze_sensitivity_mc(
        fit_edge,
        base_inputs,
        psi_config,
        n_trials=100,
        concentration_error_pct=0.05
    )

    psi_data = [{"site": site, "mineral": m, "psi": p, "mean": report.extent_means[m], "std": report.extent_stds[m]}
                for m, p in report.inclusion_probabilities.items()]
    return psi_data

def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    weights = [1.0] * len(ION_ORDER)

    config = configure_common(default_config(), PSI_ACTIVE_MINERALS)

    matrix, labels, _, penalty_scales = build_reaction_dictionary(config)

    all_path_data = []
    all_psi_data = []

    for site in ["Manu", "Talensi"]:
        if site == "Manu":
            residual = [0.5, 0.1, 0.2, 0.05, 1.0, 0.2, 0.1, 0.0, 0.0, 0.0]
        else:
            residual = [0.4, 0.2, 0.5, 0.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.0]
            
        path_data = run_regularization_path(residual, matrix, labels, weights, site, config, penalty_scales)
        all_path_data.extend(path_data)
        
        run_system_state_aicc(residual, config, weights)
        
        psi_data = run_uncertainty_cascade(config, site)
        all_psi_data.extend(psi_data)

    pd.DataFrame(all_path_data).to_csv(RESULT_DIR / "m6_regularization_path.csv", index=False)
    pd.DataFrame(all_psi_data).to_csv(RESULT_DIR / "m6_phase_stability_index.csv", index=False)
    print("\nM6 Robustness Analysis Complete.")

if __name__ == "__main__":
    main()
