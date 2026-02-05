"""
Inspect EdgeResult structure to understand what data is available
from hydrosheaf fit_network analysis.
"""

from hydrosheaf import Config, fit_network
from hydrosheaf.data.units import mgL_to_mmolL
import pandas as pd

# Load sample data
data_dir = 'hydrosheaf_synthetic_csv'
water_chem = pd.read_csv(f'{data_dir}/water_chem_full.csv')
edges_df = pd.read_csv(f'{data_dir}/network_edges.csv')

# Prepare one sample for analysis
event = 'E1-DRY'
event_data = water_chem[water_chem['event_code'] == event]

samples = []
for _, row in event_data.iterrows():
    sample = {
        'site_id': row['station_code'],
        'sample_id': f"{row['event_code']}_{row['station_code']}",
        'Ca': mgL_to_mmolL(row['Ca_mg_L'], 'Ca'),
        'Mg': mgL_to_mmolL(row['Mg_mg_L'], 'Mg'),
        'Na': mgL_to_mmolL(row['Na_mg_L'], 'Na'),
        'K': mgL_to_mmolL(row['K_mg_L'], 'K'),
        'HCO3': mgL_to_mmolL(row['HCO3_mg_L'], 'HCO3'),
        'Cl': mgL_to_mmolL(row['Cl_mg_L'], 'Cl'),
        'SO4': mgL_to_mmolL(row['SO4_mg_L'], 'SO4'),
        'NO3': mgL_to_mmolL(row['NO3_mg_L'], 'NO3'),
        'F': 0, 'Fe': 0, 'PO4': 0,
        'pH': row.get('pH', 7.0),
    }
    samples.append(sample)

edges = [(row['from_station'], row['to_station']) for _, row in edges_df.iterrows()]
config = Config()

# Run analysis
results = fit_network(samples, edges, config)

print('='*70)
print('EDGE RESULT STRUCTURE ANALYSIS')
print('='*70)

if results:
    res = results[0]
    print(f'\nEdgeResult type: {type(res).__name__}')
    print(f'\nAvailable attributes and values:')
    print('-'*70)
    
    for attr in sorted(dir(res)):
        if not attr.startswith('_'):
            try:
                val = getattr(res, attr, None)
                if not callable(val):
                    val_str = str(val)[:60] if val is not None else 'None'
                    print(f'  {attr:25s}: {type(val).__name__:15s} = {val_str}')
            except:
                pass
    
    print('\n' + '='*70)
    print('DETAILED VALUES FOR KEY ATTRIBUTES')
    print('='*70)
    
    # Check for gamma (concentration factor)
    if hasattr(res, 'gamma'):
        print(f'\ngamma (concentration factor): {res.gamma}')
    
    # Check for objective/score
    for attr in ['objective', 'score', 'obj', 'phi', 'residual', 'fit_score', 'rmse']:
        if hasattr(res, attr):
            print(f'{attr}: {getattr(res, attr)}')
    
    # Check for reaction data
    for attr in ['reactions', 'mineral_reactions', 'reaction_extents', 'delta_minerals']:
        if hasattr(res, attr):
            print(f'{attr}: {getattr(res, attr)}')
    
    # Check for transport model info
    for attr in ['transport_model', 'best_model', 'model_type', 'model']:
        if hasattr(res, attr):
            print(f'{attr}: {getattr(res, attr)}')
    
    print('\n' + '='*70)
    print('ALL RESULTS SUMMARY')
    print('='*70)
    print(f'Total edges analyzed: {len(results)}')
    
    # Collect gamma values
    gammas = [r.gamma for r in results if hasattr(r, 'gamma') and r.gamma is not None]
    print(f'Edges with gamma values: {len(gammas)}')
    if gammas:
        print(f'  Gamma range: {min(gammas):.3f} to {max(gammas):.3f}')
        print(f'  Gamma mean: {sum(gammas)/len(gammas):.3f}')
    
    # Collect objective values
    objectives = []
    for r in results:
        for attr in ['objective', 'obj', 'score']:
            if hasattr(r, attr):
                val = getattr(r, attr)
                if val is not None:
                    objectives.append(val)
                    break
    print(f'Edges with objective values: {len(objectives)}')
    if objectives:
        print(f'  Objective range: {min(objectives):.6f} to {max(objectives):.6f}')

else:
    print('No results returned from fit_network')
