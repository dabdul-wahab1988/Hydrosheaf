import pandas as pd

csv_path = r'c:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\hydrosheafupgraded\Hydrosheaf\hydrosheaf_synthetic_csv\water_chem_full.csv'
df = pd.read_csv(csv_path)

summary = df.groupby(['event_label', 'station_type']).size().unstack(fill_value=0)
summary['Total'] = summary.sum(axis=1)
print(summary)
print("\nColumn Totals:")
print(summary.sum())
