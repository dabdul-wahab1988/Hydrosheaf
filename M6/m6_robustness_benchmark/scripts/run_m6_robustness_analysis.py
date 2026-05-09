"""
M6 Robustness and Sensitivity Analysis.
Quantifies discovery stability across regularization paths and input variance.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
import math

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

def run_regularization_path(residual, matrix, labels, weights):
    """Scenario 3.1: The L1 Discovery Plateau."""
    print("\n--- Running Regularization Path (lambda sweep) ---")
    lambdas = np.logspace(-4, -1, 20) # Focusing on the interesting range
    path_data = []

    for l in lambdas:
        fit = fit_reactions(residual, matrix, weights, lambda_l1=l)
        row = {"lambda": l, "residual_norm": fit.residual_norm, "aicc": fit.aicc}
        for label, extent in zip(labels, fit.extents):
            row[label] = extent
        path_data.append(row)

    df = pd.DataFrame(path_data)
    df.to_csv(RESULT_DIR / "m6_regularization_path.csv", index=False)

    # Identify the 'optimal' lambda based on AICc
    best_row = df.loc[df['aicc'].idxmin()]
    print(f"  Optimal Lambda (min AICc): {best_row['lambda']:.4f}")
    print(f"  Regularization path saved to {RESULT_DIR}")

def run_system_state_aicc(residual, config, weights):
    """Scenario 3.2: AICc-Based System Resolution (Open vs Closed)."""
    print("\n--- Running AICc Model Selection (Open vs Closed) ---")

    # Model 1: Open System
    config.active_minerals = ["calcite_open", "dolomite_open"]
    m_open, l_open, _, p_open = build_reaction_dictionary(config)
    fit_open = fit_reactions(residual, m_open, weights, lambda_l1=0.01, penalty_scales=p_open)

    # Model 2: Closed System
    config.active_minerals = ["calcite_closed", "dolomite_closed"]
    m_closed, l_closed, _, p_closed = build_reaction_dictionary(config)
    fit_closed = fit_reactions(residual, m_closed, weights, lambda_l1=0.01, penalty_scales=p_closed)

    print(f"  Open System   - AICc: {fit_open.aicc:.2f}, BIC: {fit_open.bic:.2f}")
    print(f"  Closed System - AICc: {fit_closed.aicc:.2f}, BIC: {fit_closed.bic:.2f}")

    winner = "Open" if fit_open.aicc < fit_closed.aicc else "Closed"
    print(f"  Statistically Selected System State: {winner}")

def run_uncertainty_cascade(config):
    """Scenario 3.3 & 3.4: Monte Carlo Phase Stability Index (PSI)."""
    print("\n--- Running Monte Carlo Phase Stability Analysis ---")
    config.sensitivity_analysis_enabled = True

    # Sample Case: Ghana-like evolution
    ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
    config.ion_order = ion_order
    config.weights = [1.0] * len(ion_order)
    config.conservative_weights = [1.0] * len(ion_order)

    x_u = [1.0, 0.5, 2.0, 0.2, 2.5, 2.0, 0.5, 0.1, 0.0, 0.01]
    # Add 0.3 mmol albite, 0.1 mmol halite
    x_v = [1.0, 0.5, 2.4, 0.2, 2.8, 2.1, 0.5, 0.1, 0.0, 0.01]

    base_inputs = {
        "x_u": x_u,
        "x_v": x_v,
        "config": config,
        "residence_time_days": 3650.0 # 10 years
    }

    report = analyze_sensitivity_mc(
        fit_edge,
        base_inputs,
        config,
        n_trials=100,
        concentration_error_pct=0.05
    )

    print(f"  Base Minerals Identified: {report.base_minerals}")
    print(f"  Overall Confidence Score: {report.confidence_score*100:.1f}%")

    # Print Inclusion Probabilities (PSI)
    print("\n  Phase Inclusion Probabilities (PSI):")
    for mineral, prob in report.inclusion_probabilities.items():
        if prob > 0:
            print(f"    - {mineral:15s}: {prob*100:5.1f}% (Mean Extent: {report.extent_means[mineral]:.4f})")

    # Save PSI Matrix
    psi_data = [{"mineral": m, "psi": p, "mean": report.extent_means[m], "std": report.extent_stds[m]}
                for m, p in report.inclusion_probabilities.items()]
    pd.DataFrame(psi_data).to_csv(RESULT_DIR / "m6_phase_stability_index.csv", index=False)

def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    # Setup a realistic test case (derived from Ghana proxies)
    ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
    # 0.5 mmol increase in Ca, 1.0 mmol increase in HCO3
    residual = [0.5, 0.1, 0.2, 0.05, 1.0, 0.2, 0.1, 0.0, 0.0, 0.0]
    weights = [1.0] * len(ion_order)

    config = default_config()
    config.ion_order = ion_order
    config.active_minerals = ["calcite", "dolomite", "albite", "halite", "pyrite_oxidation_aerobic"]
    config.measured_ions = ion_order

    matrix, labels, _, _ = build_reaction_dictionary(config)

    # 1. Lambda Path
    run_regularization_path(residual, matrix, labels, weights)

    # 2. AICc Selection
    run_system_state_aicc(residual, config, weights)

    # 3. Uncertainty Cascade (MC)
    run_uncertainty_cascade(config)

    print("\nM6 Robustness Analysis Complete.")

if __name__ == "__main__":
    main()
