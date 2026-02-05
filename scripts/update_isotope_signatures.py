"""
Update nitrate isotope signatures in VeaCatchment CSV files to reflect realistic field conditions:

LYSIMETERS (soil water from farmland):
- Farmers apply fertilizers → δ¹⁵N should be in fertilizer range
- NH₄ fertilizer: δ¹⁵N = -6 to +6‰, δ¹⁸O = -5 to +5‰ (after nitrification)
- NO₃ fertilizer: δ¹⁵N = -4 to +4‰, δ¹⁸O = +18 to +25‰
- Mixed fertilizer use: δ¹⁵N = -3 to +5‰, δ¹⁸O = +5 to +15‰

BOREHOLES (groundwater):
- Mixed sources: fertilizer leaching + soil organic N + some denitrification
- δ¹⁵N = +5 to +15‰ (enriched due to mixing and denitrification)
- δ¹⁸O = +5 to +15‰ (mixed signature)
- Denitrification trend: both δ¹⁵N and δ¹⁸O increase along flowpath
"""

import argparse
import warnings
import pandas as pd
import numpy as np

# By default do NOT use random augmentation; pass --allow-random to enable
parser = argparse.ArgumentParser(description='Update isotope signatures (optional random augmentation)')
parser.add_argument('--allow-random', action='store_true', help='Allow randomized augmentation of synthetic isotope values')
args = parser.parse_args()
ALLOW_RANDOM = bool(args.allow_random)

if not ALLOW_RANDOM:
    def _det_uniform(low, high, *a, **k):
        warnings.warn(f"Deterministic fallback for uniform({low}, {high}) used. Pass --allow-random to enable randomness.")
        return (low + high) / 2.0
    np.random.uniform = _det_uniform
else:
    np.random.seed(123)

base_dir = 'VeaCatchment_SyntheticDataset_CBEpass_v1_csv'

# ============================================================
# Update Water_Routine_Bottles.csv
# ============================================================
print("Updating isotope signatures in Water_Routine_Bottles.csv...")
water_bottles = pd.read_csv(f'{base_dir}/Water_Routine_Bottles.csv')

print(f"  Total rows: {len(water_bottles)}")
print(f"  Stations: {water_bottles['Station'].unique()}")

# Update isotope values based on Matrix type
for idx, row in water_bottles.iterrows():
    station = row['Station']
    matrix = row['Matrix']
    bottle_type = row['BottleType']
    depth = row.get('Depth_cm', 0)
    
    # Only update I (isotope) and W (water isotope) bottles
    if bottle_type == 'I':  # Nitrate isotopes
        if matrix == 'soil_water':  # Lysimeters - FERTILIZER signature
            # Fertilizer application area - NH4 fertilizer that gets nitrified
            # Some NO3 fertilizer also applied
            # δ¹⁵N: -3 to +5‰ (fertilizer range)
            # δ¹⁸O: +5 to +15‰ (nitrified NH4 + some NO3 fertilizer)
            
            if station in ['L1', 'L4']:  # High input areas
                d15N = np.random.uniform(-2, 4)
                d18O = np.random.uniform(8, 18)
            elif station == 'L3':  # Lower fertilizer input
                d15N = np.random.uniform(0, 5)
                d18O = np.random.uniform(5, 12)
            else:  # L2 - moderate
                d15N = np.random.uniform(-1, 4)
                d18O = np.random.uniform(6, 15)
            
            # Deeper samples may show some transformation
            if depth and depth > 30:
                d15N += np.random.uniform(0, 2)  # Slight enrichment with depth
                d18O += np.random.uniform(-2, 2)
            
            water_bottles.at[idx, 'd15N_NO3_permil'] = d15N
            water_bottles.at[idx, 'd18O_NO3_permil'] = d18O
            
        elif matrix == 'groundwater':  # Boreholes - MIXED sources
            # Mixed signature: fertilizer leaching + soil organic N + denitrification
            # δ¹⁵N: +5 to +15‰ (enriched from original fertilizer)
            # δ¹⁸O: +5 to +15‰ (mixed)
            
            # Upgradient boreholes - more fertilizer influence
            if station == 'BH1':
                d15N = np.random.uniform(5, 10)
                d18O = np.random.uniform(6, 12)
            # Downgradient boreholes - more denitrification, more mixing
            elif station in ['BH2', 'BH3']:
                d15N = np.random.uniform(8, 14)
                d18O = np.random.uniform(8, 14)
            else:  # BH4 - furthest downgradient
                d15N = np.random.uniform(10, 16)
                d18O = np.random.uniform(10, 16)
            
            water_bottles.at[idx, 'd15N_NO3_permil'] = d15N
            water_bottles.at[idx, 'd18O_NO3_permil'] = d18O
    
    elif bottle_type == 'W':  # Water isotopes (δ²H, δ¹⁸O of H2O)
        # These should remain similar - reflect local precipitation and evaporation
        # Ghana LMWL: δ²H = 7.87 * δ¹⁸O + 13.61
        
        if matrix == 'soil_water':
            # Soil water - some evaporative enrichment
            d18O_H2O = np.random.uniform(-3.5, -1.5)
            d2H_H2O = 7.87 * d18O_H2O + 13.61 + np.random.uniform(-5, 0)  # Below LMWL due to evaporation
        else:  # groundwater
            # Groundwater - less evaporated, closer to LMWL
            d18O_H2O = np.random.uniform(-4.5, -2.5)
            d2H_H2O = 7.87 * d18O_H2O + 13.61 + np.random.uniform(-3, 2)
        
        water_bottles.at[idx, 'd18O_H2O_permil'] = d18O_H2O
        water_bottles.at[idx, 'd2H_H2O_permil'] = d2H_H2O

# Save updated file
water_bottles.to_csv(f'{base_dir}/Water_Routine_Bottles.csv', index=False)
print("  ✓ Updated Water_Routine_Bottles.csv")

# ============================================================
# Verify the updates
# ============================================================
print("\n" + "="*60)
print("VERIFICATION: Nitrate Isotope Signatures")
print("="*60)

# Reload and summarize
wb = pd.read_csv(f'{base_dir}/Water_Routine_Bottles.csv')
isotope_bottles = wb[wb['BottleType'] == 'I']

print("\nLYSIMETERS (soil water) - Expected: Fertilizer signature")
print("-" * 50)
lys_data = isotope_bottles[isotope_bottles['Matrix'] == 'soil_water']
for station in ['L1', 'L2', 'L3', 'L4']:
    st_data = lys_data[lys_data['Station'] == station]
    if len(st_data) > 0:
        d15N_mean = st_data['d15N_NO3_permil'].mean()
        d18O_mean = st_data['d18O_NO3_permil'].mean()
        print(f"  {station}: δ¹⁵N = {d15N_mean:.1f}‰, δ¹⁸O = {d18O_mean:.1f}‰")

print("\nBOREHOLES (groundwater) - Expected: Mixed sources")
print("-" * 50)
bh_data = isotope_bottles[isotope_bottles['Matrix'] == 'groundwater']
for station in ['BH1', 'BH2', 'BH3', 'BH4']:
    st_data = bh_data[bh_data['Station'] == station]
    if len(st_data) > 0:
        d15N_mean = st_data['d15N_NO3_permil'].mean()
        d18O_mean = st_data['d18O_NO3_permil'].mean()
        print(f"  {station}: δ¹⁵N = {d15N_mean:.1f}‰, δ¹⁸O = {d18O_mean:.1f}‰")

print("\n" + "="*60)
print("Expected pattern in nitrate source plot:")
print("="*60)
print("• Lysimeters: Cluster in FERTILIZER zones (NH₄/NO₃)")
print("  - δ¹⁵N: -3 to +5‰")
print("  - δ¹⁸O: +5 to +18‰")
print("• Boreholes: Spread across MIXED/DENITRIFICATION zones")
print("  - δ¹⁵N: +5 to +16‰ (enriched)")
print("  - δ¹⁸O: +6 to +16‰")
print("• Denitrification arrow shows enrichment trend from")
print("  lysimeters → boreholes")
