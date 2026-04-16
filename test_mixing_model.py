"""
Test Mixing Hypothesis: BH2 = f * L1 + (1-f) * BH1
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
    print("HYDROSHEAF CONCEPTUAL MODEL CHECK: MIXING")
    print("="*60)

    # 1. Load Data
    data_dir = Path("data/synthetic")
    chem_file = data_dir / "Hydrochem_CBE_Routine.csv"
    chem_df = pd.read_csv(chem_file)
    
    ion_map = {
        'Ca_mgL': 'Ca', 'Mg_mgL': 'Mg', 'Na_mgL': 'Na', 'K_mgL': 'K',
        'HCO3_mgL': 'HCO3', 'Cl_mgL': 'Cl', 'SO4_mgL': 'SO4', 'NO3_mgL': 'NO3'
    }
    
    # Get L1, BH1 (Upgradient), BH2 (Target)
    l1_row = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "L1")].iloc[0]
    bh1_row = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "BH1")].iloc[0]
    bh2_row = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "BH2")].iloc[0]
    
    samples = []
    
    # Use standard ion order for all operations
    ion_order = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4']
    
    # Build Samples for fit_network (U=L1, V=BH2)
    # The mixing model will use Config.mixing_endmembers for the second component.
    for row in [l1_row, bh2_row]: 
        sample = {'site_id': row['Station'], 'sample_id': row['PairID'], 'pH': row.get('pH', 7.0), 'temp_C': 25.0}
        for csv, ion in ion_map.items():
            val = mgL_to_mmolL(row.get(csv, 0), ion)
            sample[ion] = val
        for ion in ['F', 'Fe', 'PO4']: sample[ion] = 0.0
        samples.append(sample)
        
    # Build BH1 Vector (Endmember)
    bh1_vec = []
    # We need to manually convert BH1 values to mmol/L in the correct order
    bh1_sample = {} # Temp dict
    for csv, ion in ion_map.items():
        bh1_sample[ion] = mgL_to_mmolL(bh1_row.get(csv, 0), ion)
    
    for ion in ion_order:
        bh1_vec.append(bh1_sample.get(ion, 0.0))

    # print("Endmember BH1 (Upgradient):", bh1_vec)

    # 2. Configure Model for MIXING
    print("\n[CONFIG] Enabling Mixing with BH1 Endmember...")
    config = Config(
        ion_order=ion_order,
        weights=[1.0]*11,
        phreeqc_enabled=False,
        transport_models_enabled=['mix'], # Force mixing
        mixing_endmembers={
            "Upgradient_GW": bh1_vec
        },
        active_minerals=["calcite", "dolomite", "gypsum", "pyrite_oxidation_aerobic"],
        exchange_enabled=True,
        gibbs_enabled=True
    )
    
    edge = Edge(u="L1", v="BH2", edge_id="L1->BH2")
    
    print("\n[MODEL] Running Mixing Fit...")
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
    print(f"\n1. Objective Score: {res.objective_score:.4f}")
    
    # f is fraction of U (L1)
    if res.f is not None:
        print(f"\n2. Mixing Fraction:")
        print(f"   f (Fraction of L1 - Recharge): {res.f:.3f}")
        print(f"   1-f (Fraction of BH1 - Regional): {1.0 - res.f:.3f}")
    else:
        print("\n2. Mixing Fraction: None (Failed?)")
    
    # Residuals
    residuals = res.residual_vector
    # Reconstruct V (Observed)
    obs_v = []
    for ion in ion_order:
        obs_v.append(samples[1].get(ion, 0.0))
    
    obs_v_mean = np.mean(obs_v)
    sse = sum(r**2 for r in residuals)
    sst = sum((o - obs_v_mean)**2 for o in obs_v)
    r2 = 1 - (sse/sst) if sst > 0 else 0
    
    print(f"\n3. Goodness of Fit:")
    print(f"   SSE: {sse:.4f}")
    print(f"   R-squared: {r2:.4f} (Prev Evap R2: ~0.35)")
    
    print("\n4. Reactions in Mixed Water:")
    if res.z_labels and res.z_extents:
        for lbl, ext in zip(res.z_labels, res.z_extents):
            if abs(ext) > 0.01:
                print(f"  {lbl}: {ext:+.4f} mmol/L")

if __name__ == "__main__":
    main()
