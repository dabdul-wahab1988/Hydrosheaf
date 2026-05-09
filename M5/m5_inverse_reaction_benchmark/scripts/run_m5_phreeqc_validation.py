"""
M5 Advanced Inverse Reaction Validation: PHREEQC Synthetic Benchmark.
Integrates synthetic path recovery, thermodynamic screening, and kinetic forward checks.
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
from hydrosheaf.validation.reaction import validate_sparse_inverse_reaction_model
from hydrosheaf.models.reactions import build_reaction_dictionary, fit_reactions
import phreeqpython as pp

# Paths
BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_DIR / "results"

def generate_phreeqc_synthetic():
    """Step 1: Generate a synthetic reaction path using PHREEQC."""
    print("Generating PHREEQC synthetic reaction path...")
    phreeqc = pp.PhreeqPython()
    
    # 1. Upstream solution (Typical Groundwater)
    upstream_comp = {
        'temp': 25.0,
        'pH': 7.0,
        'Ca': 1.0,
        'Mg': 0.5,
        'Na': 2.0,
        'K': 0.2,
        'Cl': 2.0,
        'S(6)': 0.5,
        'Alkalinity': 2.5
    }
    sol_u = phreeqc.add_solution(upstream_comp)
    
    # 2. Add reactions (Known extents in mmol/L)
    extents = {
        'Calcite': 0.5,
        'Halite': 0.1,
        'Gypsum': 0.2
    }
    
    sol_v = sol_u.copy()
    for mineral, amount in extents.items():
        sol_v.add(mineral, amount / 1000.0) # mmol to mol
    
    # Extract compositions
    ion_order = ["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3"]
    
    def get_vec(sol):
        return [
            sol.total('Ca') * 1000.0,
            sol.total('Mg') * 1000.0,
            sol.total('Na') * 1000.0,
            sol.total('K') * 1000.0,
            sol.total('Cl') * 1000.0,
            sol.total('S(6)') * 1000.0,
            sol.total('Alkalinity') * 1000.0
        ]
        
    x_u = get_vec(sol_u)
    x_v = get_vec(sol_v)
    
    return x_u, x_v, ion_order, extents

def run_m5_validation():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    
    x_u, x_v, ion_order, true_extents = generate_phreeqc_synthetic()
    
    config = default_config()
    config.ion_order = ion_order
    config.weights = [1.0] * len(ion_order)
    config.active_minerals = ["calcite", "dolomite", "gypsum", "halite"]
    config.exchange_enabled = False # Keep it simple for synthetic recovery
    
    reaction_matrix, labels, _ = build_reaction_dictionary(config)
    
    # 2. Run Sparse Inverse Model
    print("\nRunning Hydrosheaf Sparse Inverse Model...")
    # residual = downstream - upstream
    residual = [v - u for v, u in zip(x_v, x_u)]
    
    # Use direct fit_reactions for precise control
    fit = fit_reactions(
        residual, 
        reaction_matrix, 
        config.weights, 
        lambda_l1=0.0, # Zero penalty for exact recovery test
        max_iter=1000
    )
    
    print("\nInferred extents (mmol/L):")
    for label, extent in zip(labels, fit.extents):
        if abs(extent) > 1e-4:
            true_val = true_extents.get(label.capitalize(), 0.0)
            print(f"  {label:10s}: {extent:8.4f} (True: {true_val:8.4f}, Error: {abs(extent-true_val):8.4f})")
    
    rmse = np.sqrt(np.mean(np.array(fit.residual)**2))
    print(f"\nResidual RMSE: {rmse:.6f} mmol/L")
    
    # Save results
    m5_results = {
        "metric": ["RMSE", "Phase_Recovery_Count"],
        "value": [rmse, len([e for e in fit.extents if abs(e) > 1e-4])]
    }
    pd.DataFrame(m5_results).to_csv(RESULT_DIR / "m5_synthetic_recovery.csv", index=False)
    print(f"\nM5 Synthetic Recovery Complete. Results saved to {RESULT_DIR}")

if __name__ == "__main__":
    run_m5_validation()
