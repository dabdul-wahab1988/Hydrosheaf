import pandas as pd
import numpy as np

df = pd.read_csv(r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\results\edge_fit_results.csv')
base = df[(df['scenario']=='complete') & (df['topology_variant']=='full')]

print('=== PANEL B: Reaction Recovery ===')
for m in ['calcite', 'gypsum', 'halite']:
    col = f'reaction_{m}'
    vals = pd.to_numeric(base[col], errors='coerce').dropna()
    print(f'{m}: n={len(vals)}, mean={vals.mean():.3f}, std={vals.std():.3f}, max={vals.max():.3f}')
    print(f'  >0.01 count: {(vals > 0.01).sum()}, near-zero count: {(vals.abs() < 0.01).sum()}')

print()
print('=== PANEL D: Isotope Shift ===')
iso = pd.read_csv(r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\results\res_isotopes.csv')
print('n=', len(iso))
print('unique true=', sorted(iso['true_shift'].unique().tolist()))
print('inf range=', iso['inf_shift'].min(), 'to', iso['inf_shift'].max())

mask = iso['true_shift'] > 0.01
r2_all = np.corrcoef(iso['true_shift'], iso['inf_shift'])[0,1]**2
r2_pos = np.corrcoef(iso[mask]['true_shift'], iso[mask]['inf_shift'])[0,1]**2
print('R2 all=', round(r2_all, 3))
print('R2 pos=', round(r2_pos, 3))
print('true>0 count=', mask.sum())

mae = np.mean(np.abs(iso['true_shift'] - iso['inf_shift']))
mae_pos = np.mean(np.abs(iso[mask]['true_shift'] - iso[mask]['inf_shift']))
print('MAE all=', round(mae, 3))
print('MAE pos=', round(mae_pos, 3))

false_pos = iso[(iso['true_shift'].abs() < 0.01) & (iso['inf_shift'].abs() > 0.3)]
print('False positives (true=0, |inf|>0.3):', len(false_pos))
