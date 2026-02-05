"""
Add L3 and L4 lysimeter data to VeaCatchment CSV files.
L3: Agricultural plot east (10.8850, -0.8680) - low_input zone
L4: Agricultural plot west (10.8910, -0.8900) - high_input zone (manure applied)
"""

import argparse
import warnings
import pandas as pd
import numpy as np

# By default do NOT use random augmentation; pass --allow-random to enable
parser = argparse.ArgumentParser(description='Add L3/L4 lysimeter data (optional random augmentation)')
parser.add_argument('--allow-random', action='store_true', help='Allow randomized augmentation of synthetic samples')
args = parser.parse_args()
ALLOW_RANDOM = bool(args.allow_random)

if not ALLOW_RANDOM:
    # Replace common numpy random functions with deterministic fallbacks that warn
    def _det_uniform(low, high, *args, **kwargs):
        warnings.warn(f"Deterministic fallback used for uniform({low}, {high}). Pass --allow-random to enable randomness.")
        return (low + high) / 2.0

    def _det_choice(seq, *args, **kwargs):
        warnings.warn("Deterministic fallback used for choice(). Pass --allow-random to enable randomness.")
        # return first element deterministically
        return seq[0]

    np.random.uniform = _det_uniform
    np.random.choice = _det_choice
else:
    np.random.seed(42)

base_dir = 'VeaCatchment_SyntheticDataset_CBEpass_v1_csv'

# Lysimeter coordinates and characteristics
lysimeter_info = {
    'L3': {'lat': 10.8850, 'lon': -0.8680, 'input_zone': 'low_input'},
    'L4': {'lat': 10.8910, 'lon': -0.8900, 'input_zone': 'high_input_manure'},
}

# ============================================================
# 1. Update Hydrochem_CBE_Routine.csv
# ============================================================
print("Updating Hydrochem_CBE_Routine.csv...")
hydrochem = pd.read_csv(f'{base_dir}/Hydrochem_CBE_Routine.csv')

# Get unique events
events = hydrochem['EventCode'].unique()
print(f"  Events: {events}")

# Generate L3 and L4 data based on L1 and L2 patterns
new_rows = []

for event in events:
    event_data = hydrochem[hydrochem['EventCode'] == event]
    
    # Get event metadata from existing L1 row
    l1_30 = event_data[(event_data['Station'] == 'L1') & (event_data['Depth_cm'] == 30.0)]
    if len(l1_30) == 0:
        continue
    
    l1_30 = l1_30.iloc[0]
    month = l1_30['Month']
    season = l1_30['Season']
    
    for lys, info in lysimeter_info.items():
        for depth in [30.0, 60.0]:
            # Base values - vary from L1/L2 patterns
            if lys == 'L3':
                # L3: Low input - lower NO3, lower EC, more natural signature
                base_l2 = event_data[(event_data['Station'] == 'L2') & (event_data['Depth_cm'] == depth)]
                if len(base_l2) == 0:
                    continue
                base = base_l2.iloc[0].copy()
                
                # Modify for low input zone
                Ca = base['Ca_mgL'] * np.random.uniform(0.7, 0.9)
                Mg = base['Mg_mgL'] * np.random.uniform(0.8, 1.0)
                Na = base['Na_mgL'] * np.random.uniform(0.6, 0.8)
                K = base['K_mgL'] * np.random.uniform(0.5, 0.7)
                HCO3 = base['HCO3_mgL'] * np.random.uniform(0.8, 1.1)
                Cl = base['Cl_mgL'] * np.random.uniform(0.5, 0.7)
                SO4 = base['SO4_mgL'] * np.random.uniform(0.6, 0.8)
                NO3 = base['NO3_mgL'] * np.random.uniform(0.3, 0.5)  # Much lower NO3
                pH = np.random.uniform(5.8, 6.8)
                EC = base['EC_uScm'] * np.random.uniform(0.6, 0.8)
                
            else:  # L4
                # L4: High input with manure - elevated NO3, higher K
                base_l1 = event_data[(event_data['Station'] == 'L1') & (event_data['Depth_cm'] == depth)]
                if len(base_l1) == 0:
                    continue
                base = base_l1.iloc[0].copy()
                
                # Modify for high manure input
                Ca = base['Ca_mgL'] * np.random.uniform(1.0, 1.3)
                Mg = base['Mg_mgL'] * np.random.uniform(1.1, 1.4)
                Na = base['Na_mgL'] * np.random.uniform(0.9, 1.1)
                K = base['K_mgL'] * np.random.uniform(1.5, 2.5)  # Higher K from manure
                HCO3 = base['HCO3_mgL'] * np.random.uniform(1.0, 1.3)
                Cl = base['Cl_mgL'] * np.random.uniform(1.1, 1.4)
                SO4 = base['SO4_mgL'] * np.random.uniform(0.8, 1.1)
                NO3 = base['NO3_mgL'] * np.random.uniform(1.3, 1.8)  # Higher NO3 from manure
                pH = np.random.uniform(5.5, 6.5)
                EC = base['EC_uScm'] * np.random.uniform(1.1, 1.4)
            
            # Calculate charge balance
            cations = Ca/20.04 + Mg/12.15 + Na/22.99 + K/39.10
            anions = HCO3/61.02 + Cl/35.45 + SO4/48.03 + NO3/62.00
            cbe = ((cations - anions) / (cations + anions)) * 100
            
            # Create PairID
            pair_id = f"{event}-{lys}-{int(depth)}"
            
            new_row = {
                'PairID': pair_id,
                'EventCode': event,
                'Month': month,
                'Season': season,
                'Matrix': 'soil_water',
                'Station': lys,
                'Latitude': info['lat'],
                'Longitude': info['lon'],
                'Depth_cm': depth,
                'InputZone': info['input_zone'],
                'pH': pH,
                'EC_uScm': EC,
                'Ca_mgL': Ca,
                'Mg_mgL': Mg,
                'Na_mgL': Na,
                'K_mgL': K,
                'HCO3_mgL': HCO3,
                'Cl_mgL': Cl,
                'SO4_mgL': SO4,
                'NO3_mgL': NO3,
                'Cations_meqL': cations,
                'Anions_meqL': anions,
                'CBE_pct': cbe,
                'CBE_pass': abs(cbe) < 5,
            }
            new_rows.append(new_row)

# Add new rows
new_df = pd.DataFrame(new_rows)
hydrochem_updated = pd.concat([hydrochem, new_df], ignore_index=True)

# Sort by EventCode, Station, Depth
hydrochem_updated = hydrochem_updated.sort_values(['EventCode', 'Station', 'Depth_cm']).reset_index(drop=True)

hydrochem_updated.to_csv(f'{base_dir}/Hydrochem_CBE_Routine.csv', index=False)
print(f"  Added {len(new_rows)} new rows. Total: {len(hydrochem_updated)}")

# ============================================================
# 2. Update Water_Routine_Bottles.csv
# ============================================================
print("\nUpdating Water_Routine_Bottles.csv...")
water_bottles = pd.read_csv(f'{base_dir}/Water_Routine_Bottles.csv')

new_water_rows = []

for event in events:
    event_data = water_bottles[water_bottles['EventCode'] == event]
    
    for lys, info in lysimeter_info.items():
        for depth in [30.0, 60.0]:
            # Get corresponding hydrochem row for values
            hc_row = hydrochem_updated[(hydrochem_updated['EventCode'] == event) & 
                                       (hydrochem_updated['Station'] == lys) & 
                                       (hydrochem_updated['Depth_cm'] == depth)]
            if len(hc_row) == 0:
                continue
            hc = hc_row.iloc[0]
            
            pair_id = f"{event}-{lys}-{int(depth)}"
            
            # Generate isotope values based on input zone
            if lys == 'L3':
                # Low input - more natural soil N signature
                d15N = np.random.uniform(2, 6)  # Soil organic N
                d18O_NO3 = np.random.uniform(0, 8)
                d2H = np.random.uniform(-25, -15)
                d18O_H2O = np.random.uniform(-4.5, -3.0)
            else:  # L4
                # High manure input - enriched δ15N
                d15N = np.random.uniform(8, 15)  # Manure signature
                d18O_NO3 = np.random.uniform(2, 10)
                d2H = np.random.uniform(-20, -10)
                d18O_H2O = np.random.uniform(-4.0, -2.5)
            
            # Get event metadata
            ref_row = event_data.iloc[0]
            month = ref_row['Month']
            season = ref_row['Season']
            
            # Create 4 bottle types: A, U, I, W
            bottle_data = [
                ('A', {'Ca_mgL': hc['Ca_mgL'], 'Mg_mgL': hc['Mg_mgL'], 'Na_mgL': hc['Na_mgL'], 'K_mgL': hc['K_mgL']}),
                ('U', {'HCO3_mgL': hc['HCO3_mgL'], 'Cl_mgL': hc['Cl_mgL'], 'SO4_mgL': hc['SO4_mgL'], 'NO3_mgL': hc['NO3_mgL']}),
                ('I', {'NO3_mgL': hc['NO3_mgL'] * np.random.uniform(0.98, 1.02), 'd15N_NO3_permil': d15N, 'd18O_NO3_permil': d18O_NO3}),
                ('W', {'d2H_H2O_permil': d2H, 'd18O_H2O_permil': d18O_H2O}),
            ]
            
            for bottle_type, specific_data in bottle_data:
                sample_id = f"{event}-{lys}-{int(depth)}cm-{bottle_type}"
                
                row = {
                    'SampleID': sample_id,
                    'PairID': pair_id,
                    'Category': 'Routine',
                    'EventCode': event,
                    'Month': month,
                    'Season': season,
                    'Matrix': 'soil_water',
                    'Station': lys,
                    'Latitude': info['lat'],
                    'Longitude': info['lon'],
                    'Depth_cm': depth,
                    'BottleType': bottle_type,
                    'InputZone': info['input_zone'],
                    'pH': hc['pH'],
                    'EC_uScm': hc['EC_uScm'],
                    'Ca_mgL': specific_data.get('Ca_mgL', np.nan),
                    'Mg_mgL': specific_data.get('Mg_mgL', np.nan),
                    'Na_mgL': specific_data.get('Na_mgL', np.nan),
                    'K_mgL': specific_data.get('K_mgL', np.nan),
                    'HCO3_mgL': specific_data.get('HCO3_mgL', np.nan),
                    'Cl_mgL': specific_data.get('Cl_mgL', np.nan),
                    'SO4_mgL': specific_data.get('SO4_mgL', np.nan),
                    'NO3_mgL': specific_data.get('NO3_mgL', np.nan),
                    'd15N_NO3_permil': specific_data.get('d15N_NO3_permil', np.nan),
                    'd18O_NO3_permil': specific_data.get('d18O_NO3_permil', np.nan),
                    'd2H_H2O_permil': specific_data.get('d2H_H2O_permil', np.nan),
                    'd18O_H2O_permil': specific_data.get('d18O_H2O_permil', np.nan),
                }
                new_water_rows.append(row)

new_water_df = pd.DataFrame(new_water_rows)
water_updated = pd.concat([water_bottles, new_water_df], ignore_index=True)
water_updated = water_updated.sort_values(['EventCode', 'Station', 'Depth_cm', 'BottleType']).reset_index(drop=True)

water_updated.to_csv(f'{base_dir}/Water_Routine_Bottles.csv', index=False)
print(f"  Added {len(new_water_rows)} new rows. Total: {len(water_updated)}")

# ============================================================
# 3. Update Hydrochem_CBE_QC.csv (add QC duplicates for L3, L4)
# ============================================================
print("\nUpdating Hydrochem_CBE_QC.csv...")
try:
    hydrochem_qc = pd.read_csv(f'{base_dir}/Hydrochem_CBE_QC.csv')
    
    # Add QC samples for L3 and L4 (randomly select some events)
    qc_events = np.random.choice(events, size=min(3, len(events)), replace=False)
    new_qc_rows = []
    
    for event in qc_events:
        for lys in ['L3', 'L4']:
            depth = np.random.choice([30.0, 60.0])
            hc_row = hydrochem_updated[(hydrochem_updated['EventCode'] == event) & 
                                       (hydrochem_updated['Station'] == lys) & 
                                       (hydrochem_updated['Depth_cm'] == depth)]
            if len(hc_row) == 0:
                continue
            
            hc = hc_row.iloc[0].copy()
            hc['PairID'] = f"{event}-{lys}-{int(depth)}-QC"
            # Add slight variation for QC duplicate
            for col in ['Ca_mgL', 'Mg_mgL', 'Na_mgL', 'K_mgL', 'HCO3_mgL', 'Cl_mgL', 'SO4_mgL', 'NO3_mgL']:
                if col in hc and pd.notna(hc[col]):
                    hc[col] = hc[col] * np.random.uniform(0.97, 1.03)
            
            new_qc_rows.append(hc.to_dict())
    
    if new_qc_rows:
        new_qc_df = pd.DataFrame(new_qc_rows)
        hydrochem_qc_updated = pd.concat([hydrochem_qc, new_qc_df], ignore_index=True)
        hydrochem_qc_updated.to_csv(f'{base_dir}/Hydrochem_CBE_QC.csv', index=False)
        print(f"  Added {len(new_qc_rows)} QC rows")
    else:
        print("  No QC rows added")
except Exception as e:
    print(f"  Error updating QC file: {e}")

# ============================================================
# 4. Update Water_QC_Bottles.csv
# ============================================================
print("\nUpdating Water_QC_Bottles.csv...")
try:
    water_qc = pd.read_csv(f'{base_dir}/Water_QC_Bottles.csv')
    
    new_water_qc_rows = []
    for event in qc_events:
        for lys in ['L3', 'L4']:
            info = lysimeter_info[lys]
            depth = np.random.choice([30.0, 60.0])
            pair_id = f"{event}-{lys}-{int(depth)}-QC"
            
            hc_row = hydrochem_updated[(hydrochem_updated['EventCode'] == event) & 
                                       (hydrochem_updated['Station'] == lys) & 
                                       (hydrochem_updated['Depth_cm'] == depth)]
            if len(hc_row) == 0:
                continue
            hc = hc_row.iloc[0]
            
            ref = water_qc.iloc[0]
            month = hc['Month']
            season = hc['Season']
            
            for bottle_type in ['A', 'U']:
                sample_id = f"{pair_id}-{bottle_type}"
                row = {
                    'SampleID': sample_id,
                    'PairID': pair_id,
                    'Category': 'QC',
                    'EventCode': event,
                    'Month': month,
                    'Season': season,
                    'Matrix': 'soil_water',
                    'Station': lys,
                    'Latitude': info['lat'],
                    'Longitude': info['lon'],
                    'Depth_cm': depth,
                    'BottleType': bottle_type,
                    'InputZone': info['input_zone'],
                    'pH': hc['pH'],
                    'EC_uScm': hc['EC_uScm'],
                    'Ca_mgL': hc['Ca_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'A' else np.nan,
                    'Mg_mgL': hc['Mg_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'A' else np.nan,
                    'Na_mgL': hc['Na_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'A' else np.nan,
                    'K_mgL': hc['K_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'A' else np.nan,
                    'HCO3_mgL': hc['HCO3_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'U' else np.nan,
                    'Cl_mgL': hc['Cl_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'U' else np.nan,
                    'SO4_mgL': hc['SO4_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'U' else np.nan,
                    'NO3_mgL': hc['NO3_mgL'] * np.random.uniform(0.97, 1.03) if bottle_type == 'U' else np.nan,
                    'd15N_NO3_permil': np.nan,
                    'd18O_NO3_permil': np.nan,
                    'd2H_H2O_permil': np.nan,
                    'd18O_H2O_permil': np.nan,
                }
                new_water_qc_rows.append(row)
    
    if new_water_qc_rows:
        new_water_qc_df = pd.DataFrame(new_water_qc_rows)
        water_qc_updated = pd.concat([water_qc, new_water_qc_df], ignore_index=True)
        water_qc_updated.to_csv(f'{base_dir}/Water_QC_Bottles.csv', index=False)
        print(f"  Added {len(new_water_qc_rows)} QC water bottle rows")
except Exception as e:
    print(f"  Error: {e}")

print("\n✅ Done! VeaCatchment CSV files now include L1, L2, L3, and L4 lysimeters.")
print("\nLysimeter characteristics:")
print("  L1: High input (fertilizer) - north recharge area")
print("  L2: Moderate input - south recharge area")
print("  L3: Low input (natural) - agricultural plot east")
print("  L4: High input (manure) - agricultural plot west")
