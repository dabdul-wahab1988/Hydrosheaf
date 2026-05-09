"""
Validation of Hydrosheaf on the 'talensi' Field Dataset.
Tests whether the framework reaches high accuracy (R2 ~0.999) using M2 thesis upgrades.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import math

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import default_config
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.inference.network_fit import infer_edges

def load_talensi_data():
    path = REPO_ROOT / "talensi.csv"
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")
    
    df = pd.read_csv(path)
    
    # Pre-process columns
    mapping = {
        "Code": "sample_id",
        "Town": "site_id",
        "Longitude": "lon",
        "Latitude": "lat",
        "Elevation": "elevation",
        "Temp": "temp_c"
    }
    df.rename(columns=mapping, inplace=True)
    
    # Ion conversion
    # Note: 'F' was missing in previous runs, checking header now.
    ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'Fe']
    for ion in ions:
        # Force numeric conversion, coercion errors to NaN
        df[ion] = pd.to_numeric(df[ion], errors='coerce')
        # Convert mg/L to mmol/L
        df[ion] = df[ion].apply(lambda x: mgL_to_mmolL(float(x), ion) if pd.notna(x) else 0.0)
    
    # Add dummy F if missing for Config consistency (if not in header)
    if 'F' not in df.columns:
        df['F'] = 0.0
    else:
        df['F'] = pd.to_numeric(df['F'], errors='coerce').fillna(0.0).apply(lambda x: mgL_to_mmolL(x, 'F'))
        
    # Coordinates mapping for inference
    df["easting"] = df["lon"]
    df["northing"] = df["lat"]
    df["screen_depth"] = 25.0 # Assumption for Talensi wells
    df["head_meas"] = df["elevation"]
    
    return df

def run_validation():
    df = load_talensi_data()
    print(f"Loaded {len(df)} samples from talensi.csv")
    
    # 1. Infer Topology (Elevation as Head Proxy)
    sample_dicts = df.to_dict("records")
    
    # 2. Configure Hydrosheaf with Thesis Upgrades
    config = default_config()
    config.ion_order = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe']
    config.weights = [1.0] * len(config.ion_order)
    config.active_minerals = [
        "calcite", "dolomite", "albite", "halite", 
        "fluorite", "k_feldspar", "biotite", "CO2(g)", 
        "pyrite_oxidation_aerobic",
        "SO4_input" 
    ]
    config.exchange_enabled = True
    config.edge_radius_km = 100.0
    config.edge_max_neighbors = 10
    config.latent_endmembers_enabled = True
    config.latent_endmembers_count = 3
    
    # Enable Upgrades
    config.use_thermodynamic_logic_gates = True
    config.use_temporal_sheaf_sections = True
    config.phreeqc_enabled = False
    
    # Technical Remediation: Honest Modeling
    config.honest_modeling = True
    config.measured_ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'Fe'] 

    print("Building geography-based topology...")
    inferred_edges_obj = infer_edges(
        sample_dicts, 
        method="probabilistic", 
        max_neighbors=10, 
        allow_uphill=True,
        head_key="head_meas",
        elevation_key="elevation",
        config=config
    )
    edge_tuples = [(str(e.edge_id), str(e.u), str(e.v)) for e in inferred_edges_obj]
    print(f"Inferred {len(edge_tuples)} flow-path edges using probabilistic logic.")
    
    # 3. Run Inversion
    print("Fitting reaction pathways...")
    results, extras = fit_network_pipeline(sample_dicts, edge_tuples, config)
    
    # 4. Analyze Accuracy
    r2_scores = [r.chemistry_r2 for r in results if not math.isnan(r.chemistry_r2)]
    objectives = [r.objective_score for r in results if not math.isnan(r.objective_score)]
    
    if not r2_scores:
        print("No valid fits produced.")
        return

    print("\n--- 'talensi' Dataset Accuracy Summary ---")
    print(f"Median Chemistry R2:  {np.median(r2_scores):.5f}")
    print(f"Max Chemistry R2:     {np.max(r2_scores):.5f}")
    print(f"Mean Objective Score: {np.mean(objectives):.5f}")
    
    # Identify Top Process Pathways
    print("\n--- Detailed Process Breakdown ---")
    
    # Real-to-Real
    real_fits = [r for r in results if not ("Virtual" in r.u or "Virtual" in r.v)]
    if real_fits:
        top_real = sorted(real_fits, key=lambda x: x.objective_score)[:3]
        print("\nTop 3 Real-to-Real Flow Paths:")
        for res in top_real:
            print(f"  * {res.edge_id}: R2={res.chemistry_r2:.4f}")
            for label, extent in zip(res.z_labels, res.z_extents):
                if abs(extent) > 0.01: print(f"    - {label}: {extent:+.3f} mmol/L")
    
    # Virtual-to-Real (Source Discovery)
    virt_fits = [r for r in results if "Virtual" in r.u and not "Virtual" in r.v]
    if virt_fits:
        top_virt = sorted(virt_fits, key=lambda x: x.objective_score)[:3]
        print("\nTop 3 Source-to-Well Flow Paths (Recharge Discovery):")
        for res in top_virt:
            print(f"  * {res.edge_id}: R2={res.chemistry_r2:.4f}")
            for label, extent in zip(res.z_labels, res.z_extents):
                if abs(extent) > 0.01: print(f"    - {label}: {extent:+.3f} mmol/L")

    # Final Verdict
    target_met = np.max(r2_scores) >= 0.999
    print(f"\nFinal Verdict: Accuracy target (0.999) reached? {'YES' if target_met else 'NO'}")

if __name__ == "__main__":
    try:
        run_validation()
    except Exception as e:
        print(f"Validation failed: {e}")
        import traceback
        traceback.print_exc()
