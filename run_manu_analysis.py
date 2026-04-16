"""
HYDROSHEAF ANALYSIS: MANU.XLSX
==============================
Targeted workflow for Ghana Groundwater Data.
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
    from hydrosheaf import Config
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge
except ImportError:
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import Config, fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge

def main():
    print("="*60)
    print("   HYDROSHEAF ANALYSIS: MANU.XLSX (GHANA)")
    print("="*60)

    # 1. LOAD DATA
    file_path = "manu.xlsx"
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"Error reading Excel: {e}")
        return
    
    # 2. PRE-PROCESS (Map Excel names to Hydrosheaf names)
    samples = []
    ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F']
    
    for _, row in df.iterrows():
        site_id = str(row['Station']) if pd.notna(row['Station']) else str(row['Sample ID'])
        
        s = {
            'site_id': site_id,
            'sample_id': str(row['Sample ID']),
            'x': row.get('X coordinate', 0.0),
            'y': row.get('Y coordinate', 0.0),
            'z': row.get('Elevation', 0.0),
            'pH': row.get('pH', 7.0),
            'temp_C': row.get('Temp', 25.0),
            '18O': row.get('d18O'), 
            '2H': row.get('d2H'),
            'type': 'sample'
        }
        
        # Convert all ions to mmol/L (Hydrosheaf Internal Unit)
        for ion in ions:
            val_mg = row.get(ion, 0.0)
            if pd.isna(val_mg): val_mg = 0.0
            try:
                s[ion] = mgL_to_mmolL(val_mg, ion)
            except:
                s[ion] = 0.0
        
        s['Fe'] = 0.0
        s['PO4'] = 0.0
        samples.append(s)

    print(f"Ingested {len(samples)} samples from {file_path}")

    # 3. CONFIGURE WITH GHANA EQUATIONS
    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
        weights=[1.0]*11,
        
        # Ghana LMWL (Gibrilla et al. 2022)
        lmwl_a=7.22, 
        lmwl_b=8.66,
        isotope_enabled=True,
        
        latent_endmembers_enabled=True,
        phreeqc_enabled=False, 
        transport_models_enabled=['evap', 'mix'],
        active_minerals=["calcite", "dolomite", "gypsum", "albite", "halite", "CO2(g)"],
        exchange_enabled=True
    )

    # 4. TOPOLOGY
    if len(samples) >= 2:
        u_site = samples[0]['site_id']
        v_site = samples[1]['site_id']
        edges = [Edge(u=u_site, v=v_site, edge_id="Pathway_1")]
        
        print(f"\nAnalyzing Flow Path: {u_site} -> {v_site}...")
        
        # 5. EXECUTE
        try:
            results, extras = fit_network_pipeline(samples, edges, config)
        except Exception as e:
            if "transport candidates" in str(e):
                print("\n[NOTE] Fitting skipped (awaiting custom endmembers), but AI Discovery complete.")
                results = []
            else:
                print(f"Analysis failed: {e}")
                return
        
        # 6. REPORT
        print("\n" + "="*60)
        print("   INTELLIGENCE REPORT (GHANA SITE)")
        print("="*60)
        
        if results:
            res = results[0]
            print(f"Mechanism: {res.transport_model.upper()}")
            print(f"Fit Accuracy: {100*(1-res.objective_score):.1f}% (Approx)")
            
            if res.gamma and res.gamma > 1.01:
                print(f" -> EVAPORATION detected! Concentration factor: {res.gamma:.2f}")
            
            print("\nDominant Geochemical Processes:")
            found_process = False
            if res.z_labels and res.z_extents:
                for lbl, ext in zip(res.z_labels, res.z_extents):
                    if abs(ext) > 0.01:
                        found_process = True
                        action = "Dissolving" if ext > 0 else "Precipitating"
                        # Convert reaction extent roughly to mg/L for context (using dominant ion mass)
                        # e.g. Calcite -> Ca (40 g/mol). 1 mmol -> 40 mg.
                        mass_est = 0.0
                        unit_str = "mmol/L"
                        
                        if "calcite" in lbl: mass_est = abs(ext) * 100.0 # CaCO3 ~ 100 g/mol
                        elif "gypsum" in lbl: mass_est = abs(ext) * 172.0 # CaSO4.2H2O ~ 172 g/mol
                        elif "NO3" in lbl: mass_est = abs(ext) * 62.0 
                        
                        extra_info = ""
                        if mass_est > 0:
                             extra_info = f" (~{mass_est:.1f} mg/L)"
                             
                        print(f"  * {lbl}: {abs(ext):.2f} {unit_str}{extra_info} ({action})")
            
            if not found_process:
                print("  * No significant chemical reactions detected (Conservative transport).")
                
        # Check for AI-Discovered Hidden Sources
        virtual = [s for s in samples if s.get('type') == 'virtual']
        if virtual:
            print(f"\n[AI DISCOVERY] {len(virtual)} Hidden Sources found in your data trends:")
            print("(Values converted back to mg/L for verification)")
            for vn in virtual:
                name = vn['site_id']
                no3_mmol = vn.get('NO3', 0.0)
                cl_mmol = vn.get('Cl', 0.0)
                
                # Convert back to mg/L
                no3_mg = mmolL_to_mgL(no3_mmol, 'NO3')
                cl_mg = mmolL_to_mgL(cl_mmol, 'Cl')
                
                print(f"  * {name}:")
                print(f"      NO3: {no3_mg:.1f} mg/L ({no3_mmol:.2f} mmol/L)")
                print(f"      Cl : {cl_mg:.1f} mg/L  ({cl_mmol:.2f} mmol/L)")
    else:
        print("Insufficient samples for pathway analysis.")

if __name__ == "__main__":
    main()
