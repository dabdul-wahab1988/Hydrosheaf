"""
Improve Model Fit by Enabling CO2 Gas Exchange.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import traceback

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from hydrosheaf import Config, fit_network
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge
except ImportError:
    # If installed in editable mode or not found, try adding 'hydrosheaf' subdir
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import Config, fit_network
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge

def main():
    print("="*60)
    print("HYDROSHEAF MODEL IMPROVEMENT: ADDING CO2(g)")
    print("="*60)

    # 1. Load Data
    data_dir = Path("data/synthetic")
    chem_file = data_dir / "Hydrochem_CBE_Routine.csv"
    chem_df = pd.read_csv(chem_file)
    
    ion_map = {
        'Ca_mgL': 'Ca', 'Mg_mgL': 'Mg', 'Na_mgL': 'Na', 'K_mgL': 'K',
        'HCO3_mgL': 'HCO3', 'Cl_mgL': 'Cl', 'SO4_mgL': 'SO4', 'NO3_mgL': 'NO3'
    }
    
    # Extract L1 and BH2 for E1-DRY
    l1_row = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "L1")].iloc[0]
    bh2_row = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "BH2")].iloc[0]
    
    samples = []
    # Build sample list
    for row in [l1_row, bh2_row]:
        sample = {'site_id': row['Station'], 'sample_id': row['PairID'], 'pH': row.get('pH', 7.0), 'temp_C': 25.0}
        for csv, ion in ion_map.items():
            val_mg = row.get(csv, 0)
            if pd.isna(val_mg): val_mg = 0.0
            try:
                sample[ion] = mgL_to_mmolL(val_mg, ion)
            except:
                sample[ion] = 0.0
        # Fill
        for ion in ['F', 'Fe', 'PO4']: sample[ion] = 0.0
        samples.append(sample)

    # 2. Configure Model WITH CO2(g)
    print("\n[CONFIG] Enabling CO2(g) gas phase...")
    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
        weights=[1.0]*11,
        phreeqc_enabled=False,
        transport_models_enabled=['evap'],
        active_minerals=[
            "calcite", "dolomite", "gypsum", "pyrite_oxidation_aerobic",
            "CO2(g)" # <--- THE FIX: Allows degassing
        ],
        exchange_enabled=True,
        gibbs_enabled=True
    )
    
    edge = Edge(u="L1", v="BH2", edge_id="L1->BH2")
    
    print("\n[MODEL] Running Fit for L1 -> BH2...")
    try:
        results = fit_network(samples, [edge], config)
    except Exception as e:
        print(f"Fit failed: {e}")
        traceback.print_exc()
        return

    res = results[0]
    
    if not res:
        print("Model failed.")
        return

    # 3. Analyze Results
    print(f"\n1. New Objective Score: {res.objective_score:.4f} (Was ~6.7)")
    
    # Reconstruct Modeled V
    obs_v = []
    for ion in config.ion_order:
        obs_v.append(samples[1].get(ion, 0.0))
        
    residuals = res.residual_vector
    modeled_v = [o - r for o, r in zip(obs_v, residuals)]
    
    print(f"\n{'Ion':<10} {'Obs BH2':>10} {'Modeled BH2':>12} {'Residual':>10} {'% Error':>10}")
    print("-" * 60)
    
    total_sse = 0
    total_sst = 0
    mean_obs = sum(obs_v) / len(obs_v)
    
    for i, ion in enumerate(config.ion_order):
        obs_v_val = obs_v[i]
        mod_v_val = modeled_v[i]
        resid = residuals[i]
        
        pct_err = (resid / obs_v_val * 100) if obs_v_val > 1e-6 else 0.0
        print(f"{ion:<10} {obs_v_val:>10.3f} {mod_v_val:>12.3f} {resid:>10.3f} {pct_err:>9.1f}%")
        
        total_sse += resid**2
        total_sst += (obs_v_val - mean_obs)**2

    r2 = 1 - (total_sse / total_sst) if total_sst > 0 else 0
    print("-" * 60)
    print(f"New R-squared: {r2:.4f} (Was ~0.35)")
    
    # Check CO2 reaction
    print("\n2. Reaction Extents:")
    if res.z_labels and res.z_extents:
        for lbl, ext in zip(res.z_labels, res.z_extents):
            if abs(ext) > 0.01:
                print(f"  {lbl}: {ext:+.4f} mmol/L")
                if "CO2" in lbl and ext < 0:
                     print("  -> DEGASSING DETECTED! (CO2 Loss)")

if __name__ == "__main__":
    main()
