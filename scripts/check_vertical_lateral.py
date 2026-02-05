"""Check vertical vs lateral edge statistics."""
import pandas as pd

df = pd.read_csv('hydrosheaf_manual/figures/vea_catchment/data/edge_results.csv')

print("=== Vertical Edges (Lysimeter -> Borehole) ===")
vert = df[df['u'].str.startswith('L')]
print(f"Count: {len(vert)}")
print(f"Mean gamma: {vert['gamma'].mean():.3f}")
print(f"Mean isotope_penalty: {vert['isotope_penalty'].mean():.3f}")
print(f"Mean objective_score: {vert['objective_score'].mean():.3f}")
print()

print("=== Lateral Edges (Borehole/AW -> AW) ===")
lat = df[~df['u'].str.startswith('L')]
print(f"Count: {len(lat)}")
print(f"Mean gamma: {lat['gamma'].mean():.3f}")  
print(f"Mean isotope_penalty: {lat['isotope_penalty'].mean():.3f}")
print(f"Mean objective_score: {lat['objective_score'].mean():.3f}")
print()

print("=== Sample Edges ===")
print("Vertical edges (first 5):")
print(vert[['u', 'v', 'gamma', 'isotope_penalty']].head())
print()
print("Lateral edges (first 5):")
print(lat[['u', 'v', 'gamma', 'isotope_penalty']].head())
