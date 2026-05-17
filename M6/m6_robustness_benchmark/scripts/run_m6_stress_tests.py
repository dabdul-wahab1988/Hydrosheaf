"""
M6 Structural Robustness and Bias Stress Testing.
Addresses Flaws 1, 2, and 3 from the expert reviewer critique.
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

# Paths
BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_DIR / "results"

def run_leave_one_out_test(residual, config, weights):
    """
    Flaw 1: Structural Uncertainty (Leave-One-Out Mineral Test).
    Systematically removes one mineral from the library and observes AICc response.
    """
    print("\n--- Running Leave-One-Out (LOO) Structural Test ---")
    
    # Original baseline
    m_base, l_base, _, p_base = build_reaction_dictionary(config)
    fit_base = fit_reactions(residual, m_base, weights, lambda_l1=0.01, penalty_scales=p_base, lambda_l2=config.lambda_l2)
    base_aicc = fit_base.aicc
    base_minerals = {l_base[i] for i, e in enumerate(fit_base.extents) if abs(e) > 1e-6}
    
    loo_results = []
    active_minerals = list(config.active_minerals)
    
    for i, mineral in enumerate(active_minerals):
        # Create library with one mineral removed
        sub_minerals = [m for j, m in enumerate(active_minerals) if i != j]
        config_loo = copy.deepcopy(config)
        config_loo.active_minerals = sub_minerals
        
        m_loo, l_loo, _, p_loo = build_reaction_dictionary(config_loo)
        fit_loo = fit_reactions(residual, m_loo, weights, lambda_l1=0.01, penalty_scales=p_loo, lambda_l2=config.lambda_l2)
        
        # Calculate AICc delta
        delta_aicc = fit_loo.aicc - base_aicc
        
        # Is this mineral "Essential"? (If AICc increases significantly)
        is_essential = delta_aicc > 2.0
        
        loo_results.append({
            "removed_mineral": mineral,
            "delta_aicc": delta_aicc,
            "is_essential": is_essential,
            "was_in_baseline": mineral in base_minerals
        })
        
        print(f"  Removed {mineral:15s}: dAICc={delta_aicc:+6.2f} Essential? {'YES' if is_essential else 'NO'}")

    df = pd.DataFrame(loo_results)
    df.to_csv(RESULT_DIR / "m6_loo_structural_uncertainty.csv", index=False)

def run_regional_bias_test(config):
    """
    Flaw 2: Systematic Bias vs Random Noise.
    Applies a systematic shift to inputs and evaluates if the discovery flips.
    """
    print("\n--- Running Regional Bias Stress Test ---")
    ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
    config.ion_order = ion_order
    config.weights = [1.0] * len(ion_order)
    config.active_minerals = ["albite", "anorthite", "calcite", "dolomite", "pyrite_oxidation_aerobic"]
    
    # Baseline Case: Silicate weathering
    # 0.2 mmol increase in Na, 0.2 mmol increase in HCO3
    residual_base = [0.0, 0.0, 0.2, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    m, l, _, p = build_reaction_dictionary(config)
    fit_ref = fit_reactions(residual_base, m, config.weights, lambda_l1=0.01, penalty_scales=p, lambda_l2=config.lambda_l2)
    ref_minerals = {l[i] for i, e in enumerate(fit_ref.extents) if abs(e) > 1e-6}
    
    # Bias: Assume a systematic +0.1 mmol shift in ALL Calcium measurements due to regional hard water
    bias_residual = list(residual_base)
    bias_residual[ion_order.index("Ca")] += 0.1
    
    fit_bias = fit_reactions(bias_residual, m, config.weights, lambda_l1=0.01, penalty_scales=p, lambda_l2=config.lambda_l2)
    bias_minerals = {l[i] for i, e in enumerate(fit_bias.extents) if abs(e) > 1e-6}
    
    print(f"  Reference Minerals (No Bias): {ref_minerals}")
    print(f"  Discovery with +0.1 Ca Bias:  {bias_minerals}")
    
    stability = len(ref_minerals.intersection(bias_minerals)) / len(ref_minerals) if ref_minerals else 1.0
    print(f"  Discovery Stability under Bias: {stability*100:.1f}%")

def run_process_grouping_stability(config):
    """
    Flaw 3: L1-Jitter (Equifinality).
    Evaluates the stability of 'Groups' (e.g., Total Silicates) instead of just species.
    """
    print("\n--- Running Process Grouping Stability Analysis ---")
    # We define groups: Silicates (albite, anorthite), Carbonates (calcite, dolomite)
    groups = {
        "Silicates": ["albite", "anorthite", "k_feldspar"],
        "Carbonates": ["calcite", "dolomite", "calcite_open", "calcite_closed"]
    }
    
    # Use Monte Carlo to see if individual minerals flip but groups stay stable
    n_trials = 50
    group_hits = {"Silicates": 0, "Carbonates": 0}
    species_hits = {"albite": 0, "anorthite": 0}
    
    ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
    residual = [0.1, 0.0, 0.2, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0] # Mixed signal
    
    m, l, _, p = build_reaction_dictionary(config)
    
    for _ in range(n_trials):
        # Add random noise
        noisy_res = [r + np.random.normal(0, 0.02) for r in residual]
        fit = fit_reactions(noisy_res, m, [1.0]*10, lambda_l1=0.01, penalty_scales=p, lambda_l2=config.lambda_l2)
        active = {l[i] for i, e in enumerate(fit.extents) if abs(e) > 1e-6}
        
        # Check Species
        if "albite" in active: species_hits["albite"] += 1
        if "anorthite" in active: species_hits["anorthite"] += 1
        
        # Check Groups
        if any(s in active for s in groups["Silicates"]): group_hits["Silicates"] += 1
        if any(c in active for c in groups["Carbonates"]): group_hits["Carbonates"] += 1
        
    print(f"  Species Stability (Albite):    {species_hits['albite']/n_trials*100:5.1f}%")
    print(f"  Species Stability (Anorthite): {species_hits['anorthite']/n_trials*100:5.1f}%")
    print(f"  GROUP Stability (Silicates):   {group_hits['Silicates']/n_trials*100:5.1f}%")
    
    print("\n  INTERPRETATION: If GROUP stability > Species stability, the model suffers from ")
    print("  mineral equifinality but preserves physical process discovery.")

def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    config = default_config()
    config.ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
    config.active_minerals = ["albite", "anorthite", "calcite", "dolomite", "halite", "pyrite_oxidation_aerobic"]
    config.measured_ions = config.ion_order
    
    # Target Residual (Typical mixed granite signal)
    weights = [1.0] * 10
    residual = [0.3, 0.1, 0.2, 0.05, 0.8, 0.2, 0.1, 0.0, 0.0, 0.02]
    
    run_leave_one_out_test(residual, config, weights)
    run_regional_bias_test(config)
    run_process_grouping_stability(config)
    
    print("\nM6 Stress Testing Complete.")

if __name__ == "__main__":
    main()
