"""
HYDROSHEAF ANALYSIS: SILICATE WEATHERING MODEL
==============================================
Optimized for Crystalline Basement / Granitic Aquifers.
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
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge
except ImportError:
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import Config, fit_network
    from hydrosheaf.api import fit_network_pipeline
    from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL
    from hydrosheaf.graph.types import Edge

def main():
    print("="*60)
    print("   HYDROSHEAF SILICATE AQUIFER DIAGNOSTIC")
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

    # 2. CONFIGURE FOR SILICATE TERRAIN
    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
        weights=[1.0]*11,
        
        lmwl_a=7.22, lmwl_b=8.66, isotope_enabled=True,
        latent_endmembers_enabled=True,
        phreeqc_enabled=False, 
        transport_models_enabled=['evap', 'mix'],
        
        # --- THE SILICATE SUITE ---
        active_minerals=[
            "albite", "anorthite", "k_feldspar", "biotite", # Silicates
            "calcite", "dolomite", "gypsum", "halite",      # Others
            "CO2(g)", "CaNa_exch", "NO3src"
        ],
        exchange_enabled=True
    )

    # 3. TOPOLOGY & EXECUTION
    edges = []
    import math
    for i, u in enumerate(samples):
        for j, v in enumerate(samples):
            if i == j: continue
            dx = u['x'] - v['x']
            dy = u['y'] - v['y']
            dist = math.sqrt(dx*dx + dy*dy)
            dz = u['z'] - v['z']
            if dz > -5.0:
                 edges.append({'u': u['site_id'], 'v': v['site_id'], 'dist': dist, 'dz': dz})
    
    final_edges = []
    u_groups = {}
    for e in edges: u_groups.setdefault(e['u'], []).append(e)
    for u, candidates in u_groups.items():
        candidates.sort(key=lambda x: x['dist'])
        for cand in candidates[:2]:
            final_edges.append(Edge(u=cand['u'], v=cand['v'], edge_id=f"{cand['u']}->{cand['v']}"))

    print(f"\nRunning Silicate Weathering Model on {len(final_edges)} paths...")
    
    try:
        results, extras = fit_network_pipeline(samples, final_edges, config)
    except Exception as e:
        print(f"Analysis failed: {e}")
        return

    # 4. REPORT: SILICATE VS CARBONATE
    print("\n" + "="*60)
    print("   SILICATE VS CARBONATE WEATHERING REPORT")
    print("="*60)
    
    silicate_flux = 0.0
    carbonate_flux = 0.0
    exchange_flux = 0.0
    biotite_flux = 0.0
    
    silicate_minerals = ["albite", "anorthite", "k_feldspar", "biotite"]
    carbonate_minerals = ["calcite", "dolomite"]
    
    for res in results:
        if res.z_labels and res.z_extents:
            for lbl, ext in zip(res.z_labels, res.z_extents):
                if lbl in silicate_minerals:
                    silicate_flux += abs(ext)
                    if lbl == "biotite":
                        biotite_flux += abs(ext)
                elif lbl in carbonate_minerals:
                    carbonate_flux += abs(ext)
                elif lbl == "CaNa_exch":
                    exchange_flux += abs(ext)

    print(f"Total Weathering Flux (mmol/L summed across system):")
    print(f"  Silicate Weathering:  {silicate_flux:.2f}")
    print(f"  Carbonate Dissolution: {carbonate_flux:.2f}")
    print(f"  Ion Exchange Activity: {exchange_flux:.2f}")
    
    ratio = silicate_flux / (carbonate_flux + 1e-6)
    print(f"\nSilicate/Carbonate Ratio: {ratio:.2f}")
    
    if ratio > 1.0:
        print("-> The aquifer chemistry is DOMINATED by Silicate Hydrolysis (Granitic Signature).")
    else:
        print("-> Despite being a silicate aquifer, Carbonate dissolution (fracture fillings?) dominates the chemistry.")

    if biotite_flux > 0.5:
        print(f"\nSignificant Biotite Weathering detected (Total Flux: {biotite_flux:.2f}).")
        print("Expect elevated Fluoride (F) and Iron (Fe) in these waters.")

if __name__ == "__main__":
    main()
