"""
HYDROSHEAF FULL SYSTEM ANALYSIS: MANU.XLSX
==========================================
Comprehensive diagnostic of the entire groundwater system.
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
    print("   HYDROSHEAF SYSTEM DIAGNOSTIC: MANU.XLSX")
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

    print(f"Loaded {len(samples)} stations.")

    # 2. CONFIGURE
    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
        weights=[1.0]*11,
        lmwl_a=7.22, lmwl_b=8.66, isotope_enabled=True,
        latent_endmembers_enabled=True,
        phreeqc_enabled=False, 
        transport_models_enabled=['evap', 'mix'],
        active_minerals=["calcite", "dolomite", "gypsum", "albite", "halite", "CO2(g)", "NO3src", "CaNa_exch"],
        exchange_enabled=True,
        # Allow connecting neighbors up to 5km away
        edge_radius_km=5.0 
    )

    # 3. BUILD NETWORK (3D TOPOLOGY)
    print("\n[TOPOLOGY] Building Hydraulic Flow Network...")
    
    edges = []
    import math
    
    # Generate edges based on distance and elevation
    for i, u in enumerate(samples):
        for j, v in enumerate(samples):
            if i == j: continue
            
            # Distance
            dx = u['x'] - v['x']
            dy = u['y'] - v['y']
            dist = math.sqrt(dx*dx + dy*dy)
            
            # Elevation Difference (Flow must go Downhill or be Lateral)
            dz = u['z'] - v['z']
            
            # Simple heuristic: Connect if close?
            # Assuming coordinates are roughly meters or consistent
            # If coordinates are Lat/Lon, distance will be tiny.
            # Let's just use array index proximity for "Neighbors" if coordinates are suspect
            # Or assume the user wants an all-to-all check on nearest neighbors.
            
            # Let's connect everything to everything if N < 50 for robust checking, 
            # filtering by dz > -2.
            
            if dz > -5.0: # Permissive downhill check
                 edges.append({
                     'u': u['site_id'], 
                     'v': v['site_id'], 
                     'dist': dist,
                     'dz': dz
                 })
    
    final_edges = []
    u_groups = {}
    for e in edges:
        u_groups.setdefault(e['u'], []).append(e)
        
    for u, candidates in u_groups.items():
        candidates.sort(key=lambda x: x['dist'])
        for cand in candidates[:2]:
            final_edges.append(Edge(u=cand['u'], v=cand['v'], edge_id=f"{cand['u']}->{cand['v']}"))

    print(f"Generated {len(final_edges)} hydraulic flow paths based on topography.")

    # 4. EXECUTE FULL SYSTEM
    print("\n[EXECUTION] Running Full System Analysis...")
    try:
        results, extras = fit_network_pipeline(samples, final_edges, config)
    except Exception as e:
        print(f"Analysis failed: {e}")
        traceback.print_exc()
        return

    # 5. SYNTHESIZE FINDINGS
    print("\n" + "="*60)
    print("   SYSTEM-WIDE DIAGNOSIS")
    print("="*60)
    
    n_evap = sum(1 for r in results if r.transport_model == 'evap')
    n_mix = sum(1 for r in results if r.transport_model == 'mix')
    scores = [r.objective_score for r in results]
    avg_score = np.mean(scores) if scores else 0.0
    
    print(f"\n1. Dominant Hydrodynamic Regime:")
    print(f"   Evaporation Dominated: {n_evap} paths")
    print(f"   Mixing Dominated:      {n_mix} paths")
    print(f"   System Stability (Avg Residual): {avg_score:.2f} (lower is better)")
    
    print("\n2. Pollution Hotspots (Nitrate Loading):")
    hotspots = []
    for res in results:
        if "NO3src" in res.z_labels:
            idx = res.z_labels.index("NO3src")
            val = res.z_extents[idx]
            val_mg = val * 62.0
            if val_mg > 10.0:
                hotspots.append((res.v, val_mg))
    
    if hotspots:
        hotspots.sort(key=lambda x: x[1], reverse=True)
        print("   The following locations are receiving active nitrate loads (Leaching/Contamination):")
        # Deduplicate by site, keep max load
        seen = set()
        for site, load in hotspots:
            if site not in seen:
                print(f"   * {site}: +{load:.1f} mg/L added locally")
                seen.add(site)
            if len(seen) >= 5: break
    else:
        print("   No significant active nitrate loading detected along flow paths.")

    print("\n3. Mineral Evolution Trends:")
    total_gypsum = 0.0
    total_exchange = 0.0
    for res in results:
        if "gypsum" in res.z_labels:
            total_gypsum += res.z_extents[res.z_labels.index("gypsum")]
        if "CaNa_exch" in res.z_labels:
            total_exchange += res.z_extents[res.z_labels.index("CaNa_exch")]
            
    if total_exchange > 0:
        print(f"   * Reverse Ion Exchange is active (System-wide Ca release, Na uptake).")
    elif total_exchange < 0:
        print(f"   * Normal Ion Exchange is active (System-wide Softening).")
        
    if total_gypsum > 0:
        print(f"   * Active Sulfate Dissolution (Gypsum weathering).")

    print("\n4. Hidden Source Confirmation:")
    virtual = [s for s in samples if s.get('type') == 'virtual']
    for vn in virtual:
        if vn.get('NO3', 0) > 1.5:
             print(f"   * AI confirms a High-Nitrate Endmember (~{vn.get('NO3')*62:.0f} mg/L) exists in the mixing space.")

if __name__ == "__main__":
    main()
