"""
Run M4 Topology Validation Benchmark using USGS MODPATH reference.
Quantifies Hydrosheaf's ability to reproduce advective connectivity.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
import json
import math

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.physics.modpath import (
    load_modpath_endpoint_records,
    load_modpath_pathline_points,
)
from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.validation import validate_independent_graph_against_modpath

# Paths
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
M2_MODPATH_DATA = REPO_ROOT / "M2" / "m2_benchmark" / "external" / "modpath" / "input" / "selected_output"
RESULT_DIR = BENCHMARK_ROOT / "results"

def load_modpath_reference():
    endpoint_file = M2_MODPATH_DATA / "base-MP.end"
    pathline_file = M2_MODPATH_DATA / "base-MP.path"
    
    if not endpoint_file.exists():
        raise FileNotFoundError(f"Reference MODPATH data not found at {endpoint_file}")
        
    endpoints = load_modpath_endpoint_records(str(endpoint_file))
    # Convert endpoints to a set of true directed edges (initial_cell -> final_cell)
    ref_edges = []
    nodes_info = {}
    for ep in endpoints:
        if ep.initial_cell is not None and ep.final_cell is not None and ep.initial_cell != ep.final_cell:
            u = f"cell_{int(ep.initial_cell)}"
            v = f"cell_{int(ep.final_cell)}"
            ref_edges.append((u, v))
            # Store coordinates for node proxy
            nodes_info[u] = {"id": u, "lat": ep.y0, "lon": ep.x0, "elevation": ep.z0}
            nodes_info[v] = {"id": v, "lat": ep.y, "lon": ep.x, "elevation": ep.z}
            
    return ref_edges, list(nodes_info.values())

def run_topology_benchmark():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Loading MODPATH reference data...")
    ref_edges, nodes = load_modpath_reference()
    print(f"Loaded {len(ref_edges)} reference flow paths over {len(nodes)} nodes.")
    
    scenarios_results = []
    
    # 1. Scenario: Spatial-Only Graph Inference (Section 2.6.1)
    print("\nRunning Scenario 2.6.1: Spatial-Only Inference...")
    samples = []
    for n in nodes:
        samples.append({
            "site_id": n["id"],
            "lat": n["lat"],
            "lon": n["lon"],
            "elevation": n.get("elevation", 0.0)
        })
        
    inferred_spatial_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=True)
    inferred_spatial = [(e.u, e.v) for e in inferred_spatial_obj]
    
    print(f"  Inferred {len(inferred_spatial)} spatial edges.")
    print(f"  Example reference edges: {ref_edges[:5]}")
    print(f"  Example inferred edges:  {inferred_spatial[:5]}")
    
    report_spatial = validate_independent_graph_against_modpath(inferred_spatial, ref_edges)
    m = report_spatial["metrics"]
    scenarios_results.append({
        "scenario": "2.6.1_spatial_only",
        "precision": m["precision"],
        "recall": m["recall"],
        "f1": m["f1"]
    })
    print(f"  F1: {m['f1']:.3f} (TP={m['tp']}, FP={m['fp']}, FN={m['fn']})")

    # 2. Scenario: Head-Gradient Constrained Inference (Section 2.6.2)
    print("\nRunning Scenario 2.6.2: Head-Gradient Constrained...")
    # By setting allow_uphill=False, we use the 'elevation' (proxy for head) to gate direction
    inferred_head_obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
    inferred_head = [(e.u, e.v) for e in inferred_head_obj]
    
    report_head = validate_independent_graph_against_modpath(inferred_head, ref_edges)
    m = report_head["metrics"]
    scenarios_results.append({
        "scenario": "2.6.2_head_constrained",
        "precision": m["precision"],
        "recall": m["recall"],
        "f1": m["f1"]
    })
    print(f"  F1: {m['f1']:.3f}")
    
    # 3. Scenario: Hydrostratigraphic Constrained Inference (Section 2.6.3)
    print("\nRunning Scenario 2.6.3: Hydrostratigraphic Constrained...")
    # Add fake 'aquifer_unit' based on depth to nodes
    for s in samples:
        # Simple binary split: shallow (< 20m) and deep (>= 20m)
        s["aquifer_unit"] = "shallow" if s["elevation"] > -20.0 else "deep"
        
    # Custom filter: only allow edges within same unit
    inferred_strat_obj = []
    # Re-run inference without filter first to get candidates
    candidates = infer_edges_from_coordinates(samples, max_neighbors=5, allow_uphill=False)
    for e in candidates:
        u_unit = next(s["aquifer_unit"] for s in samples if s["site_id"] == e.u)
        v_unit = next(s["aquifer_unit"] for s in samples if s["site_id"] == e.v)
        if u_unit == v_unit:
            inferred_strat_obj.append(e)
            
    inferred_strat = [(e.u, e.v) for e in inferred_strat_obj]
    
    report_strat = validate_independent_graph_against_modpath(inferred_strat, ref_edges)
    m = report_strat["metrics"]
    scenarios_results.append({
        "scenario": "2.6.3_hydrostratigraphic",
        "precision": m["precision"],
        "recall": m["recall"],
        "f1": m["f1"]
    })
    print(f"  F1: {m['f1']:.3f} (TP={m['tp']}, FP={m['fp']}, FN={m['fn']})")

    # Save results
    res_df = pd.DataFrame(scenarios_results)
    res_df.to_csv(RESULT_DIR / "m4_topology_benchmark_summary.csv", index=False)
    print(f"\nM4 Benchmark results saved to {RESULT_DIR}")

def run_sparsity_sensitivity(samples, ref_edges):
    """
    Evaluates how inference performance decays as the number of available nodes decreases.
    """
    print("\nRunning Scenario 2.6.4: Node Sparsity Sensitivity...")
    sensitivity_results = []
    
    # Total nodes available
    total_nodes = len(samples)
    node_fractions = [0.1, 0.25, 0.5, 0.75, 1.0]
    
    for frac in node_fractions:
        n_count = int(total_nodes * frac)
        if n_count < 5: continue
        
        # Run 5 random trials per fraction for robustness
        trials_f1 = []
        trials_recall = []
        for _ in range(5):
            # Sub-sample nodes
            sub_samples = pd.DataFrame(samples).sample(n_count).to_dict("records")
            sub_ids = {s["id"] for s in sub_samples}
            
            # Filter reference edges: both u and v must be in sub-sample
            sub_ref = [e for e in ref_edges if e[0] in sub_ids and e[1] in sub_ids]
            
            if not sub_ref:
                continue
                
            # Run Head-Gradient Constrained Inference (our standard M4 scenario)
            # Map 'id' to 'site_id' for consistency with build.py expectation
            for s in sub_samples: s["site_id"] = s["id"]
            
            inferred_obj = infer_edges_from_coordinates(sub_samples, max_neighbors=2, allow_uphill=False)
            inferred = [(e.u, e.v) for e in inferred_obj]
            
            report = validate_independent_graph_against_modpath(inferred, sub_ref)
            trials_f1.append(report["metrics"]["f1"])
            trials_recall.append(report["metrics"]["recall"])
            
        if trials_f1:
            sensitivity_results.append({
                "node_fraction": frac,
                "node_count": n_count,
                "mean_f1": np.mean(trials_f1),
                "std_f1": np.std(trials_f1),
                "mean_recall": np.mean(trials_recall)
            })
            print(f"  Frac {frac:.2f} ({n_count} nodes): Mean F1={np.mean(trials_f1):.3f}")

    res_df = pd.DataFrame(sensitivity_results)
    res_df.to_csv(RESULT_DIR / "m4_sparsity_sensitivity.csv", index=False)
    print(f"Sparsity sensitivity results saved to {RESULT_DIR}")

if __name__ == "__main__":
    try:
        # Load and run main benchmark
        RESULT_DIR.mkdir(parents=True, exist_ok=True)
        print("Loading MODPATH reference data...")
        ref_edges, nodes_info = load_modpath_reference()
        
        # Main scenarios
        run_topology_benchmark()
        
        # New Sparsity Analysis
        run_sparsity_sensitivity(nodes_info, ref_edges)
        
    except Exception as e:
        print(f"M4 Benchmark failed: {e}")
        import traceback
        traceback.print_exc()
