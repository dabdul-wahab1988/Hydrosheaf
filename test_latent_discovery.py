"""
Test Automated Discovery of Latent Endmembers (Direct Call).
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from hydrosheaf import Config
    from hydrosheaf.models.latent import identify_latent_endmembers
    from hydrosheaf.data.units import mgL_to_mmolL
except ImportError:
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import Config
    from hydrosheaf.models.latent import identify_latent_endmembers
    from hydrosheaf.data.units import mgL_to_mmolL

def main():
    print("="*60)
    print("HYDROSHEAF INTELLIGENT DISCOVERY: LATENT ENDMEMBERS")
    print("="*60)

    # 1. Load Data
    data_dir = Path("data/synthetic")
    chem_file = data_dir / "Hydrochem_CBE_Routine.csv"
    chem_df = pd.read_csv(chem_file)
    
    ion_map = {
        'Ca_mgL': 'Ca', 'Mg_mgL': 'Mg', 'Na_mgL': 'Na', 'K_mgL': 'K',
        'HCO3_mgL': 'HCO3', 'Cl_mgL': 'Cl', 'SO4_mgL': 'SO4', 'NO3_mgL': 'NO3'
    }
    
    samples = []
    for _, row in chem_df.iterrows():
        sample = {'site_id': row['Station'], 'sample_id': row['PairID']}
        for csv, ion in ion_map.items():
            try:
                sample[ion] = mgL_to_mmolL(row.get(csv, 0), ion)
            except:
                sample[ion] = 0.0
        samples.append(sample)

    print(f"Loaded {len(samples)} samples.")

    # 2. Run Intelligent Discovery
    ion_order = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3']
    print("\n[AI] Identifying Latent Endmembers via PCA/CLR...")
    
    virtual_nodes = identify_latent_endmembers(samples, ion_order, n_endmembers=2)
    
    # 3. Inspect Results
    print(f"\n[RESULTS] Discovered {len(virtual_nodes)} Virtual Endmembers:")
    
    for vn in virtual_nodes:
        name = vn['site_id']
        print(f"\n  Source: {name}")
        
        no3 = vn.get('NO3', 0.0)
        cl = vn.get('Cl', 0.0)
        na = vn.get('Na', 0.0)
        hco3 = vn.get('HCO3', 0.0)
        ca = vn.get('Ca', 0.0)
        mg = vn.get('Mg', 0.0)
        
        print(f"    NO3 : {no3:.2f} mmol/L")
        print(f"    Cl  : {cl:.2f} mmol/L")
        print(f"    Na  : {na:.2f} mmol/L")
        print(f"    Ca  : {ca:.2f} mmol/L")
        print(f"    HCO3: {hco3:.2f} mmol/L")
        
        # Automated Interpretation Logic
        interpretation = []
        if no3 > 2.0 and cl > 2.0:
            interpretation.append("Manure/Pollution (High NO3, High Cl)")
        if na > 2.0 and hco3 > 4.0:
            interpretation.append("Sodic Recharge/Alkaline Source")
        if cl < 0.5 and no3 < 0.5 and hco3 < 1.0:
            interpretation.append("Rainwater/Dilute Background")
        if ca > 2.0 and mg > 1.0 and hco3 > 2.0:
            interpretation.append("Regional Groundwater (Carbonate/Hard)")
            
        if interpretation:
            print(f"    -> AI Interpretation: {', '.join(interpretation)}")
        else:
            print("    -> AI Interpretation: Unclassified Composition")

if __name__ == "__main__":
    main()
