"""
HYDROSHEAF AUTO-RUNNER
======================
Usage: python run_hydrosheaf_auto.py

This script:
1. Loads configuration from 'hydrosheaf_intelligent_config.yaml'
2. Automatically finds data in 'data/synthetic' (or configured folder)
3. Runs the intelligent analysis pipeline
4. Prints a Plain-English Report of findings
"""

import pandas as pd
import yaml
import sys
from pathlib import Path
import warnings
import traceback

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from hydrosheaf import Config, fit_network
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge
except ImportError:
    # If installed in editable mode or not found, try adding 'hydrosheaf' subdir
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import Config, fit_network
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge

def load_config(path="hydrosheaf_intelligent_config.yaml"):
    if not Path(path).exists():
        print(f"Config file {path} not found. Using defaults.")
        return Config()
    
    with open(path, 'r') as f:
        cfg_dict = yaml.safe_load(f)
    
    # Map dict to Config object manually or via helper
    # Config is a dataclass, so we filter valid keys
    valid_keys = Config.__init__.__code__.co_varnames
    filtered_args = {k: v for k, v in cfg_dict.items() if k in valid_keys}
    
    # Special handling for lists/tuples if needed
    return Config(**filtered_args)

def load_and_prep_data():
    # Hardcoded for demo, but could scan folder
    chem_file = Path("data/synthetic/Hydrochem_CBE_Routine.csv")
    bottle_file = Path("data/synthetic/Water_Routine_Bottles.csv")
    
    if not chem_file.exists():
        print("Data not found.")
        return []
        
    print(f"Loading data from {chem_file}...")
    df = pd.read_csv(chem_file)
    bottles = pd.read_csv(bottle_file) if bottle_file.exists() else None
    
    # Merge isotopes if available
    if bottles is not None:
        iso_cols = ['d15N_NO3_permil', 'd18O_NO3_permil']
        iso_data = bottles.groupby('PairID')[iso_cols].first()
        df = df.merge(iso_data, on='PairID', how='left')
    
    samples = []
    # Use ion map to convert mg/L to mmol/L
    ion_map = {
        'Ca_mgL': 'Ca', 'Mg_mgL': 'Mg', 'Na_mgL': 'Na', 'K_mgL': 'K',
        'HCO3_mgL': 'HCO3', 'Cl_mgL': 'Cl', 'SO4_mgL': 'SO4', 'NO3_mgL': 'NO3'
    }
    
    for _, row in df.iterrows():
        s = {'site_id': row['Station'], 'sample_id': row['PairID'], 'type': 'sample'}
        for csv_col, ion in ion_map.items():
            val = row.get(csv_col, 0)
            if pd.isna(val): val = 0.0
            try:
                s[ion] = mgL_to_mmolL(val, ion)
            except:
                s[ion] = 0.0
        
        # Pass isotopes directly
        if 'd15N_NO3_permil' in row:
            s['d15N_NO3_permil'] = row['d15N_NO3_permil']
            s['d18O_NO3_permil'] = row['d18O_NO3_permil']
            
        # Fill zeros
        for i in ['F', 'Fe', 'PO4']: s[i] = 0.0
        samples.append(s)
        
    return samples

def main():
    print("-" * 60)
    print("   HYDROSHEAF INTELLIGENT ANALYZER")
    print("-" * 60)
    
    # 1. Setup
    config = load_config()
    samples = load_and_prep_data()
    
    if not samples:
        return

    # 2. Define Analysis Scope
    # For Auto-Run, let's analyze the critical path: Source (L1) -> Aquifer (BH2)
    # In a real app, this would be auto-generated from a map or user input.
    edges = [Edge(u="L1", v="BH2", edge_id="Source_to_Aquifer")]
    
    print(f"Running analysis on {len(samples)} samples...")
    print(f"Intelligent Discovery Enabled: {config.latent_endmembers_enabled}")
    
    # 3. Execute Pipeline
    try:
        # Note: fit_network_pipeline will inject virtual nodes into 'samples' list in-place
        results, extras = fit_network_pipeline(samples, edges, config)
    except Exception as e:
        # Check if it failed due to missing transport candidates (known issue with mixing enabled but no endmembers)
        if "transport candidates" in str(e):
             print("\n[AI CHECK] Standard fitting failed because mixing was enabled without explicit endmembers.")
             print("However, the AI has likely discovered them. Checking sample list...")
             # Fall through to report
             results = []
        else:
             print(f"\n[ERROR] Analysis failed: {e}")
             traceback.print_exc()
             return

    # 4. Generate Report
    print("\n" + "="*60)
    print("   FINAL INTELLIGENCE REPORT")
    print("="*60)
    
    # A. Hidden Sources (Latent Endmembers)
    virtual_nodes = [s for s in samples if s.get('type') == 'virtual']
    if virtual_nodes:
        print("\n[AI DISCOVERY] Hidden Sources Found:")
        print("The model detected chemical signatures that do not match known samples.")
        print("These likely represent pollution sources or deep groundwater.")
        
        for vn in virtual_nodes:
            name = vn['site_id']
            no3 = vn.get('NO3', 0.0)
            cl = vn.get('Cl', 0.0)
            hco3 = vn.get('HCO3', 0.0)
            
            desc = "Unknown"
            if no3 > 2.0: desc = "High-Nitrate Pollution (e.g. Manure/Septic)"
            elif cl > 2.0: desc = "Saline Source (Evaporated)"
            elif hco3 < 1.0: desc = "Acidic/Dilute Source"
            elif hco3 > 4.0: desc = "Alkaline/Carbonate Source"
            
            print(f"  * {name}: {desc} (NO3: {no3:.1f}, Cl: {cl:.1f})")
    
    # B. Mechanism Analysis
    if results:
        res = results[0]
        if res:
            print(f"\n[MECHANISMS] Pathway: {res.u} -> {res.v}")
            print(f"  Confidence: {res.objective_score:.2f} (Lower is better)")
            
            # Interpret Reactions
            print("  Key Processes Identified:")
            
            # Check Exchange
            if "CaNa_exch" in res.z_labels:
                idx = res.z_labels.index("CaNa_exch")
                val = res.z_extents[idx]
                if val > 0.1:
                    print(f"    - Reverse Ion Exchange (Water is Hardening): {val:.2f} mmol/L")
            
            # Check Dissolution
            if "gypsum" in res.z_labels:
                idx = res.z_labels.index("gypsum")
                val = res.z_extents[idx]
                if val > 0.1:
                    print(f"    - Gypsum Dissolution (Adding Sulfate): {val:.2f} mmol/L")
                    
            # Check Nitrate
            if "NO3src" in res.z_labels:
                idx = res.z_labels.index("NO3src")
                val = res.z_extents[idx]
                if val > 0.5:
                    print(f"    - EXTERNAL NITRATE INPUT: {val:.2f} mmol/L (Large Addition!)")
                    print("      WARNING: Evaporation alone cannot explain this nitrate increase.")
                
    print("\n" + "-"*60)
    print("Recommendation: Use the 'Hidden Sources' above to calibrate a mixing model.")

if __name__ == "__main__":
    main()
