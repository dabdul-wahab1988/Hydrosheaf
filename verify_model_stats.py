"""
Statistical Verification of Hydrosheaf Model Reliability.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import json

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
    print("HYDROSHEAF MODEL RELIABILITY VERIFICATION")
    print("="*60)

    # 1. Load Data (Focus on L1 -> BH2)
    data_dir = Path("data/synthetic")
    chem_file = data_dir / "Hydrochem_CBE_Routine.csv"
    chem_df = pd.read_csv(chem_file)
    
    ion_map = {
        'Ca_mgL': 'Ca', 'Mg_mgL': 'Mg', 'Na_mgL': 'Na', 'K_mgL': 'K',
        'HCO3_mgL': 'HCO3', 'Cl_mgL': 'Cl', 'SO4_mgL': 'SO4', 'NO3_mgL': 'NO3'
    }
    
    # Extract L1 and BH2 for E1-DRY
    l1_data = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "L1")]
    if l1_data.empty:
        print("L1 data not found.")
        return
    l1_row = l1_data.iloc[0]
    
    bh2_data = chem_df[(chem_df['EventCode'] == "E1-DRY") & (chem_df['Station'] == "BH2")]
    if bh2_data.empty:
        print("BH2 data not found.")
        return
    bh2_row = bh2_data.iloc[0]
    
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

    # 2. Configure with BOOTSTRAP UNCERTAINTY
    print("\n[CONFIG] Enabling Bootstrap Uncertainty Quantification (n=100)...")
    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
        weights=[1.0]*11,
        phreeqc_enabled=False, # Focus on stats/inverse engine
        transport_models_enabled=['evap'],
        active_minerals=["calcite", "dolomite", "gypsum", "pyrite_oxidation_aerobic"],
        exchange_enabled=True,
        gibbs_enabled=True,
        
        # ENABLE STATISTICS
        uncertainty_method="bootstrap",
        bootstrap_n_resamples=100,
        input_uncertainty_pct=5.0 # Assume 5% measurement error
    )
    
    edge = Edge(u="L1", v="BH2", edge_id="L1->BH2")
    
    print("\n[MODEL] Fitting Network...")
    try:
        results = fit_network(samples, [edge], config)
    except Exception as e:
        print(f"Fitting failed: {e}")
        return
        
    if not results:
        print("Model failed to fit (no results).")
        return
        
    res = results[0]
    
    # 3. Report Statistics
    print("\n" + "-"*60)
    print("STATISTICAL RELIABILITY REPORT")
    print("-" * 60)
    
    # A. Goodness of Fit (RMSE / R2)
    print("1. Goodness of Fit:")
    print(f"   Objective Score: {res.objective_score:.4f} (Likelihood Proxy)")
    
    # Manually calculate RMSE if needed, but assume residual_vector is present
    if res.residual_vector:
        residuals = np.array(res.residual_vector)
        rmse = np.sqrt(np.mean(residuals**2))
        print(f"   RMSE (Residual): {rmse:.4f} mmol/L")
    else:
        print("   RMSE: N/A")

    if hasattr(res, 'chemistry_r2') and res.chemistry_r2 is not None:
         print(f"   R-squared:       {res.chemistry_r2:.4f} (Ideal: 1.0)")
    
    # B. Parameter Uncertainty (Bootstrap)
    print("\n2. Parameter Confidence (Bootstrap n=100):")
    print("   Evaporation Factor (Gamma):")
    
    # Check if uncertainty fields are populated
    # Some versions put it in a separate object or directly on res
    gamma_std = getattr(res, 'gamma_std', None)
    gamma_ci_low = getattr(res, 'gamma_ci_low', None)
    gamma_ci_high = getattr(res, 'gamma_ci_high', None)
    
    if gamma_std is not None:
        print(f"     Mean: {res.gamma:.3f}")
        print(f"     Std Dev: {gamma_std:.3f}")
        print(f"     95% CI: [{gamma_ci_low:.3f}, {gamma_ci_high:.3f}]")
        
        if gamma_ci_low > 1.0:
            print("     -> Evaporation is Statistically Significant (CI > 1.0)")
        elif gamma_ci_high < 1.0:
            print("     -> Dilution is Statistically Significant (CI < 1.0)")
        else:
            print("     -> Evaporation is NOT Significant (CI overlaps 1.0)")
            
    print("\n   Reaction Extents (mmol/L):")
    
    extents_std = getattr(res, 'extents_std', None)
    extents_ci_low = getattr(res, 'extents_ci_low', None)
    extents_ci_high = getattr(res, 'extents_ci_high', None)
    
    if res.z_labels and extents_std:
        print(f"   {'Reaction':<25} {'Extent':>8} {'StdErr':>8} {'95% CI':>20}")
        print("   " + "-"*65)
        for i, lbl in enumerate(res.z_labels):
            val = res.z_extents[i]
            std = extents_std[i]
            low = extents_ci_low[i]
            high = extents_ci_high[i]
            
            # Filter for non-zero mean or significant CI
            if abs(val) > 0.01 or (low > 0 or high < 0):
                sig = "*" if (low > 0 or high < 0) else " " 
                print(f" {sig} {lbl:<25}: {val:>8.3f} +/- {std:.3f}   [{low:>6.3f}, {high:>6.3f}]")
    else:
        print("   No bootstrap statistics available for reactions.")
    
    print("-" * 60)
    print("INTERPRETATION:")
    print(" * = Statistically Significant (95% CI excludes zero)")
    print(" R-squared > 0.9 and RMSE < 0.5 mmol/L indicate excellent model fit.")

if __name__ == "__main__":
    main()
