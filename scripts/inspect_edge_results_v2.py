"""
Inspect EdgeResult structure in detail to understand what data is available
for replacing random placeholders in figures.
"""

from hydrosheaf import Config, fit_network
from hydrosheaf.data.units import mgL_to_mmolL
import pandas as pd
import numpy as np

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

print('=' * 70)
print('EXTRACTING DATA FOR FIGURE GENERATION')
print('=' * 70)

if not results:
    print("No results returned!")
    exit()

# 1. GAMMA VALUES
print('\n1. GAMMA (CONCENTRATION FACTOR) VALUES')
print('-' * 70)
gamma_data = {}
for r in results:
    edge_id = f"{r.u}->{r.v}"
    gamma_data[edge_id] = {
        'gamma': r.gamma,
        'transport_model': r.transport_model,
        'objective_score': r.objective_score if hasattr(r, 'objective_score') else None
    }
    
df_gamma = pd.DataFrame(gamma_data).T
print(df_gamma.to_string())

# 2. REACTION EXTENTS
print('\n2. REACTION EXTENTS (z_extents and z_labels)')
print('-' * 70)

# Get labels from first result
first = results[0]
labels = first.z_labels if first.z_labels else []
print(f"Reaction labels: {labels}")

# Collect extents for all edges
reaction_data = {}
for r in results:
    edge_id = f"{r.u}->{r.v}"
    if r.z_extents:
        reaction_data[edge_id] = dict(zip(labels[:len(r.z_extents)], r.z_extents))

if reaction_data:
    df_reactions = pd.DataFrame(reaction_data).T
    print(df_reactions.to_string())
    
    # Summary statistics
    print('\nReaction Summary Statistics:')
    for col in df_reactions.columns:
        vals = df_reactions[col].dropna()
        if len(vals) > 0:
            print(f"  {col}: mean={vals.mean():.4f}, range=[{vals.min():.4f}, {vals.max():.4f}]")
else:
    print("No reaction extent data available!")

# 3. OBJECTIVE SCORES (for candidate_scores)
print('\n3. CANDIDATE SCORES (Transport Model Selection)')
print('-' * 70)

for r in results[:5]:  # First 5 edges
    edge_id = f"{r.u}->{r.v}"
    print(f"\n{edge_id}:")
    if r.candidate_scores:
        for cs in r.candidate_scores[:3]:
            model = cs.get('transport_model', 'unknown')
            obj = cs.get('objective', 0)
            print(f"  {model}: objective={obj:.6f}")
    else:
        print("  No candidate scores")

# 4. REACTION FIT OBJECT
print('\n4. REACTION FIT OBJECT DETAILS')
print('-' * 70)

rf = first.reaction_fit
print(f"Type: {type(rf)}")
print(f"Attributes: {[a for a in dir(rf) if not a.startswith('_')]}")

# Access ReactionFit fields
if hasattr(rf, 'extents'):
    print(f"rf.extents: {rf.extents}")
if hasattr(rf, 'labels'):
    print(f"rf.labels: {rf.labels}")
if hasattr(rf, 'residuals'):
    print(f"rf.residuals: {rf.residuals}")
if hasattr(rf, 'converged'):
    print(f"rf.converged: {rf.converged}")

# 5. SUMMARY FOR FIGURES
print('\n' + '=' * 70)
print('SUMMARY: DATA AVAILABLE FOR FIGURES')
print('=' * 70)

print("""
FOR plot_reaction_summary():
  - Use EdgeResult.z_extents (list of floats for each mineral)
  - Use EdgeResult.z_labels (list of mineral names)
  - Available for ALL edges in the network

FOR plot_transport_parameters():
  - Gamma: EdgeResult.gamma (concentration factor, always available)
  - Objective: EdgeResult.objective_score (always available)
  - Transport model: EdgeResult.transport_model (always available)
  - Candidate scores: EdgeResult.candidate_scores (list of dicts)

NO NEED FOR RANDOM VALUES - all data is computed!
""")

# 6. CREATE SUMMARY TABLES FOR FIGURE GENERATOR
print('\n' + '=' * 70)
print('DATA EXTRACTION CODE FOR FIGURE GENERATOR')
print('=' * 70)

print("""
# Extract gamma values
gammas = {f"{e.u}->{e.v}": e.gamma for e in edge_results}

# Extract reaction extents as matrix
labels = edge_results[0].z_labels
reaction_matrix = np.array([e.z_extents for e in edge_results])
edge_ids = [f"{e.u}->{e.v}" for e in edge_results]

# Extract objective scores
objectives = {f"{e.u}->{e.v}": e.objective_score for e in edge_results}

# Get transport models
transport_models = {f"{e.u}->{e.v}": e.transport_model for e in edge_results}
""")
