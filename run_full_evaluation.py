"""
Full Evaluation of Hydrosheaf Package using Synthetic Data.
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
    from hydrosheaf import (
        Config,
        fit_network,
        summarize_network,
    )
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge
except ImportError:
    # If installed in editable mode or not found, try adding 'hydrosheaf' subdir
    sys.path.insert(0, str(Path.cwd() / "hydrosheaf"))
    from hydrosheaf import (
        Config,
        fit_network,
        summarize_network,
    )
    from hydrosheaf.data.units import mgL_to_mmolL
    from hydrosheaf.graph.types import Edge

def main():
    print("Starting Hydrosheaf Evaluation...")
    
    # 1. Load Data
    data_dir = Path("data/synthetic")
    chem_file = data_dir / "Hydrochem_CBE_Routine.csv"
    bottles_file = data_dir / "Water_Routine_Bottles.csv"
    
    if not chem_file.exists() or not bottles_file.exists():
        print("Error: Data files not found.")
        return

    print("Loading data...")
    chem_df = pd.read_csv(chem_file)
    bottles_df = pd.read_csv(bottles_file)
    
    # 2. Merge Isotope Data
    # Group bottles by PairID to get isotope values (taking mean or first non-null)
    # Using 'first' because usually one bottle has the isotope data
    iso_cols = ['d15N_NO3_permil', 'd18O_NO3_permil', 'd2H_H2O_permil', 'd18O_H2O_permil']
    iso_data = bottles_df.groupby('PairID')[iso_cols].first()
    
    # Merge into chem_df
    # Note: PairID in bottles matches PairID in chem
    full_df = chem_df.merge(iso_data, on='PairID', how='left')
    
    print(f"Merged data shape: {full_df.shape}")
    
    # 3. Define Network Structure
    # Sources: L1-L4 (Lysimeters/Soil), BH1 (Upgradient)
    # Targets: BH2-BH4 (Downgradient)
    # We will run this for each EVENT (EventCode) separately.
    
    sources = ['L1', 'L2', 'L3', 'L4', 'BH1']
    targets = ['BH2', 'BH3', 'BH4']
    
    # Generate edges
    edges = []
    for s in sources:
        for t in targets:
            # Simple edge ID: s->t
            edges.append(Edge(u=s, v=t, edge_id=f"{s}->{t}"))
            
    print(f"Defined {len(edges)} potential edges.")
    
    # 4. Prepare Samples for Hydrosheaf
    # Need to convert mg/L to mmol/L and rename columns
    ion_map = {
        'Ca_mgL': 'Ca',
        'Mg_mgL': 'Mg',
        'Na_mgL': 'Na',
        'K_mgL': 'K',
        'HCO3_mgL': 'HCO3',
        'Cl_mgL': 'Cl',
        'SO4_mgL': 'SO4',
        'NO3_mgL': 'NO3'
    }
    
    events = full_df['EventCode'].unique()
    results_summary = []
    
    for event in events:
        print(f"\nProcessing Event: {event}")
        event_df = full_df[full_df['EventCode'] == event]
        
        samples = []
        for _, row in event_df.iterrows():
            sample = {
                'site_id': row['Station'],
                'sample_id': row['PairID'],
                'pH': row.get('pH', 7.0),
                'temp_C': 25.0, # Default
                # Isotopes
                'd2H': row.get('d2H_H2O_permil'),
                'd18O': row.get('d18O_H2O_permil'),
                'd15N_NO3': row.get('d15N_NO3_permil'),
                'd18O_NO3': row.get('d18O_NO3_permil'),
            }
            
            # Ions
            for csv_col, ion_name in ion_map.items():
                val_mg = row.get(csv_col)
                if pd.notna(val_mg):
                    # Convert to mmol/L
                    try:
                        val_mmol = mgL_to_mmolL(val_mg, ion_name)
                        sample[ion_name] = val_mmol
                    except Exception as e:
                        # print(f"  Warning: conversion failed for {ion_name}: {e}")
                        sample[ion_name] = 0.0
                else:
                    sample[ion_name] = 0.0
            
            # Add required but missing ions (F, Fe, PO4)
            for ion in ['F', 'Fe', 'PO4']:
                sample[ion] = 0.0
                
            samples.append(sample)
            
        # 5. Run Fit Network
        config = Config(
            # Standard config
            ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe', 'PO4'],
            weights=[1.0]*11,
            phreeqc_enabled=False, # Disable for speed/simplicity first
            transport_models_enabled=['mix', 'evap'], # Enable mixing and evaporation
            active_minerals=[], # Disable minerals effectively
            gibbs_enabled=True,
            exchange_enabled=True
        )
        
        try:
            # Filter edges to only those where nodes exist in samples
            available_sites = set(s['site_id'] for s in samples)
            valid_edges = [e for e in edges if e.u in available_sites and e.v in available_sites]
            
            print(f"  Sites: {len(available_sites)}, Valid Edges: {len(valid_edges)}")
            
            if not valid_edges:
                print("  No valid edges for this event.")
                continue
                
            fitted_edges = fit_network(samples, valid_edges, config)
            
            print(f"  Fitted {len(fitted_edges)} edges.")
            
            # Summarize
            for fe in fitted_edges:
                if fe:
                    results_summary.append({
                        'Event': event,
                        'Edge': f"{fe.u}->{fe.v}",
                        'Model': fe.transport_model,
                        'Objective': fe.objective_score,
                        'Gamma (Evap)': fe.gamma,
                        'Mixing (f)': fe.f
                    })
                
        except Exception as e:
            print(f"  Error fitting network: {e}")
            traceback.print_exc()

    # 6. Report
    print("\n" + "="*40)
    print("SUMMARY RESULTS")
    print("="*40)
    if results_summary:
        df_res = pd.DataFrame(results_summary)
        # Fix: observed that groupby might fail if columns are missing or empty
        try:
            print(df_res.groupby(['Event', 'Model']).size().unstack(fill_value=0))
        except Exception as e:
            print("Could not group by Event/Model:", e)
            
        print("\nTop 5 Best Fits (Lowest Objective):")
        print(df_res.sort_values('Objective').head(5))
        
        # Diagnostics: Check if any mixing happened
        mixing_results = df_res[df_res['Model'] == 'mix']
        if not mixing_results.empty:
            print(f"\nMixing edges found: {len(mixing_results)}")
            print(mixing_results[['Edge', 'Mixing (f)', 'Objective']].head())
    else:
        print("No results generated.")

if __name__ == "__main__":
    main()
