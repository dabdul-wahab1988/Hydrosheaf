import pandas as pd

csv_path = r'c:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\hydrosheafupgraded\Hydrosheaf\hydrosheaf_synthetic_csv\water_chem_full.csv'
df = pd.read_csv(csv_path)

print("Samples per Station Type:")
print(df['station_type'].value_counts())

print("\nSamples per Event:")
print(df['event_label'].value_counts().sort_index())

print("\nTotal Samples:", len(df))
