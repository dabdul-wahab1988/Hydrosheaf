"""
Technical Profiling for M2 Thesis Upgrades.
Quantifies performance gains and robustness of Logic Gates, Temporal Sheafs, and Sensitivity Analysis.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import time

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import default_config
from hydrosheaf.uncertainty.sensitivity import analyze_sensitivity

# Path to Northern Ghana data
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_ROOT = BENCHMARK_ROOT / "external" / "northern_ghana"
RESULT_DIR = EXTERNAL_ROOT / "results"

def load_pilot_samples():
    path = RESULT_DIR / "northern_ghana_samples.csv"
    if not path.exists():
        raise FileNotFoundError("Run run_e4c_northern_ghana_validation.py first to generate samples.")
    return pd.read_csv(path)

def load_pilot_edges():
    path = RESULT_DIR / "northern_ghana_edge_inputs.csv"
    if not path.exists():
        raise FileNotFoundError("Run run_e4c_northern_ghana_validation.py first to generate edge inputs.")
    return pd.read_csv(path)

def run_convergence_profile(samples, edges):
    """Test 2: Convergence Analysis (Logic Gates vs. Standard)"""
    print("\n--- Running Convergence Profile ---")
    
    edge_tuples = [(str(r['edge_id']), str(r['u']), str(r['v'])) for r in edges.to_dict("records")]
    sample_records = samples.to_dict("records")
    
    # Warmup
    config_warm = default_config()
    config_warm.phreeqc_enabled = False
    fit_network_pipeline(sample_records, edge_tuples, config_warm)

    # Mode 1: No Logic Gates
    config_standard = default_config()
    config_standard.use_thermodynamic_logic_gates = False
    config_standard.phreeqc_enabled = False
    
    # Mode 2: Logic Gates Enabled with Simulated SI
    config_gated = default_config()
    config_gated.use_thermodynamic_logic_gates = True
    config_gated.phreeqc_enabled = False
    config_gated.si_logic_gate_threshold = 0.5
    
    # Simulate SI mask: assume all wells are supersaturated with calcite to test pruning
    wells = set(edges['u']).union(set(edges['v']))
    sim_si_masks = {w: {"saturation_indices": {"calcite": 1.0, "dolomite": 1.0}} for w in wells}
    
    start = time.time()
    results_std, _ = fit_network_pipeline(sample_records, edge_tuples, config_standard)
    std_time = time.time() - start
    std_iters = np.mean([r.reaction_iterations for r in results_std if hasattr(r, 'reaction_iterations')])
    
    start = time.time()
    results_gated, _ = fit_network_pipeline(
        sample_records, edge_tuples, config_gated, 
        phreeqc_results=sim_si_masks # API passes this as pre_si_masks to directed_section
    )
    gated_time = time.time() - start
    gated_iters = np.mean([r.reaction_iterations for r in results_gated if hasattr(r, 'reaction_iterations')])
    
    print(f"Standard L1: {std_time:.4f}s, Avg Iters: {std_iters:.1f}")
    print(f"Gated L1:    {gated_time:.4f}s, Avg Iters: {gated_iters:.1f}")
    if std_iters > 0:
        print(f"Iteration Reduction: {((std_iters - gated_iters)/std_iters)*100:.1f}%")

def run_phase_stability_analysis(samples, edges):
    """Test 1: Phase Stability Index"""
    print("\n--- Running Phase Stability Analysis ---")
    config = default_config()
    config.sensitivity_analysis_enabled = True
    
    # Pick the top edge from Ghana
    top_edge = edges.iloc[0].to_dict()
    u_site = top_edge['u']
    v_site = top_edge['v']
    
    sample_u = samples[samples['well_id'] == u_site].iloc[0].to_dict()
    sample_v = samples[samples['well_id'] == v_site].iloc[0].to_dict()
    
    # Mock solve function for sensitivity wrapper
    from hydrosheaf.inference.edge_fit import fit_edge
    def solve_edge_wrapper(**kwargs):
        # Kwargs contains perturbed parameters like 'Cl'
        # We update the sample vector accordingly
        x_u_perturbed = list(kwargs.get('x_u', [0.0]*len(config.ion_order)))
        return fit_edge(x_u_perturbed, list(kwargs.get('x_v', [0.0]*len(config.ion_order))), config)

    # Initial vectors
    from hydrosheaf.data.schema import vector_from_sample
    x_u, _ = vector_from_sample(sample_u, config.ion_order, config.missing_policy, config.detection_limit_policy)
    x_v, _ = vector_from_sample(sample_v, config.ion_order, config.missing_policy, config.detection_limit_policy)
    
    sens_results = analyze_sensitivity(
        solve_edge_wrapper, 
        {"x_u": x_u, "x_v": x_v}, 
        config,
        parameters=["x_u"], # Perturb the upstream concentrations
        scale=0.1
    )
    
    for res in sens_results:
        print(f"Parameter: {res.parameter}")
        print(f"  Residual Change: {res.residual_change*100:.1f}%")
        print(f"  Phase Stability: {res.phase_stability*100:.1f}%")

def main():
    try:
        samples = load_pilot_samples()
        edges = load_pilot_edges()
        
        # Filter for a representative subset to speed up profiling
        samples_sub = samples.head(50)
        edges_sub = edges[edges['u'].isin(samples_sub['well_id']) & edges['v'].isin(samples_sub['well_id'])].head(20)
        
        if edges_sub.empty:
            # Fallback if head(50) doesn't have connections
            edges_sub = edges.head(20)
            connected_wells = set(edges_sub['u']).union(set(edges_sub['v']))
            samples_sub = samples[samples['well_id'].isin(connected_wells)]

        run_convergence_profile(samples_sub, edges_sub)
        run_phase_stability_analysis(samples_sub, edges_sub)
        
    except Exception as e:
        print(f"Profiling failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
