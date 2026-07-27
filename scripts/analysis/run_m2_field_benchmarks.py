"""
M2 Field Benchmark Runner.
Processes LowerAnayari (Manu) and Talensi datasets using the Hydrosheaf pipeline.
Saves results to M2/m2_benchmark/results/field_discovery_results.csv
"""
import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf import Config
from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.graph.types import Edge


def configure_local_mixing_endmember(config, samples, label):
    """Use the lowest ionic-strength sample as a transparent local dilution endmember."""
    candidates = []
    for sample in samples:
        vector = []
        valid = True
        for ion in config.ion_order:
            value = sample.get(ion)
            if value is None or not np.isfinite(float(value)):
                valid = False
                break
            vector.append(float(value))
        if valid:
            ionic_strength_proxy = sum(max(value, 0.0) for value in vector)
            candidates.append((ionic_strength_proxy, str(sample["site_id"]), vector, sample))
    if not candidates:
        return

    _, site_id, vector, sample = min(candidates, key=lambda item: item[0])
    endmember_name = f"{label}_local_dilute_{site_id}"
    config.mixing_endmembers = {endmember_name: vector}
    if "18O" in sample and "2H" in sample:
        config.mixing_endmembers_isotopes = {
            endmember_name: (float(sample["18O"]), float(sample["2H"]))
        }


def process_manu():
    print("Processing Manu (LowerAnayari)...")
    df = pd.read_csv(ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv")
    samples = []
    ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F']

    for _, row in df.iterrows():
        s = {
            'site_id': str(row['Sample ID']),
            'x': row['X coordinate'],
            'y': row['Y coordinate'],
            'z': row['Elevation'],
            'pH': row['pH'],
            'temp_C': row['Temp'],
            '18O': row['d18O'],
            '2H': row['d2H'],
            'type': 'sample'
        }
        for ion in ions:
            val = row[ion]
            if isinstance(val, str) and '<' in val: val = 0.0005 # Detection limit half
            s[ion] = mgL_to_mmolL(float(val), ion)
        s['Fe'] = mgL_to_mmolL(float(row['Fe']), 'Fe')
        samples.append(s)

    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe'],
        weights=[1.0]*10,
        conservative_weights=[0.01]*10,
        lmwl_a=7.22, lmwl_b=8.66,
        isotope_enabled=True,
        active_minerals=["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"],
        exchange_enabled=True,
        honest_modeling=True,
        geologic_bias="crystalline"
    )
    configure_local_mixing_endmember(config, samples, "manu")

    # Use a more realistic topology: Connect to 2 nearest spatial neighbors
    from scipy.spatial import cKDTree
    coords = np.array([[s['x'], s['y']] for s in samples])
    tree = cKDTree(coords)
    edges = []
    for i, s in enumerate(samples):
        # Find 3 nearest (including self)
        dist, idx = tree.query([s['x'], s['y']], k=3)
        for j in idx:
            if i != j:
                edges.append(Edge(u=samples[i]['site_id'], v=samples[j]['site_id'], edge_id=f"Manu_{i}_{j}"))

    results, _ = fit_network_pipeline(samples, edges, config)
    return results

def process_talensi():
    print("Processing Talensi...")
    df = pd.read_csv(ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv")
    samples = []
    ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3']

    for _, row in df.iterrows():
        s = {
            'site_id': str(row['Code']),
            'x': row['Longitude'],
            'y': row['Latitude'],
            'z': row['Elevation'],
            'pH': row['pH'],
            'temp_C': row['Temp'],
            '18O': row['d18O'],
            '2H': row['d2H'],
            'type': 'sample'
        }
        for ion in ions:
            s[ion] = mgL_to_mmolL(float(row[ion]), ion)
        s['Fe'] = mgL_to_mmolL(float(row['Fe']), 'Fe')
        s['F'] = 0.0
        samples.append(s)

    config = Config(
        ion_order=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe'],
        weights=[1.0]*10,
        conservative_weights=[0.01]*10,
        lmwl_a=7.22, lmwl_b=8.66,
        isotope_enabled=True,
        active_minerals=["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"],
        exchange_enabled=True,
        honest_modeling=True,
        geologic_bias="crystalline"
    )
    configure_local_mixing_endmember(config, samples, "talensi")

    # Use a more realistic topology: Connect to 2 nearest spatial neighbors
    from scipy.spatial import cKDTree
    coords = np.array([[s['x'], s['y']] for s in samples])
    tree = cKDTree(coords)
    edges = []
    for i, s in enumerate(samples):
        dist, idx = tree.query([s['x'], s['y']], k=3)
        for j in idx:
            if i != j:
                edges.append(Edge(u=samples[i]['site_id'], v=samples[j]['site_id'], edge_id=f"Talensi_{i}_{j}"))

    results, _ = fit_network_pipeline(samples, edges, config)
    return results

def main():
    m_res = process_manu()
    t_res = process_talensi()

    out_rows = []
    for r in m_res + t_res:
        row = {
            "edge_id": r.edge_id,
            "u": r.u,
            "v": r.v,
            "chemistry_r2": r.chemistry_r2,
            "objective_score": r.objective_score,
            "transport_model": r.transport_model,
            "gamma": r.gamma,
            "f": r.f,
            "endmember_id": r.endmember_id,
        }
        # Flatten extents
        if r.z_labels and r.z_extents:
            for lbl, ext in zip(r.z_labels, r.z_extents):
                row[f"extent_{lbl}"] = ext
        out_rows.append(row)

    out_df = pd.DataFrame(out_rows)
    save_path = ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(save_path, index=False)
    print(f"Saved results to {save_path}")

if __name__ == "__main__":
    main()
