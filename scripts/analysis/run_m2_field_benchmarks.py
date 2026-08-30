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
from hydrosheaf.inference.network_fit import infer_edges


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


def _build_and_fit(samples, config, site_label):
    """Build the candidate graph and fit it using the documented pipeline.

    Candidate edges come from the probabilistic builder of Section 2.3 --
    directional confidence p_ij = Phi((h_i - h_j) / sigma_dh) with elevation as
    the head proxy, capped at edge_max_neighbors per node, filtered on
    edge_p_min / edge_gradient_min / edge_radius_km -- and the retained set is
    then chosen by the sheaf refinement, with the cohomology diagnostics
    (H0, H1, obstruction energy, per-edge leverage) attached.

    Earlier revisions of this script built the field graph with a plain 2-D
    Euclidean KD-tree over (lon, lat) and ran no refinement at all, so no
    directional test was applied to the Ghana edges and nothing was ever
    rejected. That construction is retained for comparison in
    M2/m2_benchmark/scripts/run_m2_field_documented_pipeline.py.
    """
    # No nuclear tracers exist at either field site. Leaving the age stage on
    # makes the pipeline try to treat Cl as a residence-time tracer and fail.
    config.sheaf_age_enabled = False
    config.sheaf_cohomology_enabled = True

    edges = infer_edges(samples, method="probabilistic", config=config)
    # The builder names edges "u->v"; downstream table and figure scripts key
    # the site off an edge_id prefix, so stamp it on here.
    for e in edges:
        e.edge_id = f"{site_label}_{e.edge_id}"
    print(f"  {site_label}: {len(edges)} candidate edges from the probabilistic builder")
    results, report = fit_network_pipeline(
        samples, edges, config, sheaf_refinement_enabled=True
    )
    stage = (report.get("stage_status", {}) or {}).get("sheaf_refinement", {})
    status = stage.get("status") if isinstance(stage, dict) else stage
    print(f"  {site_label}: sheaf refinement {status}; {len(results)} edges retained "
          f"({len(edges) - len(results)} rejected)")
    if status != "completed":
        raise RuntimeError(
            f"{site_label}: sheaf refinement did not complete (status={status!r}); "
            "the reported field topology depends on it"
        )
    # p_ij is written onto the candidate Edge attrs by the builder, not onto the
    # EdgeResult, so carry it across for reporting.
    conf = {e.edge_id: (e.attrs or {}).get("edge_confidence",
                                           (e.attrs or {}).get("p_uv"))
            for e in edges}
    for r in results:
        if r.edge_confidence is None:
            r.edge_confidence = conf.get(r.edge_id)
    return results


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
            # keys the documented graph builder reads (Section 2.3): coordinates
            # for the search radius, elevation as the head proxy inside p_ij
            'lon': float(row['X coordinate'] if 'X coordinate' in row else row['Longitude']),
            'lat': float(row['Y coordinate'] if 'Y coordinate' in row else row['Latitude']),
            'elevation': float(row['Elevation']),
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
    return _build_and_fit(samples, config, "Manu")

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
            # keys the documented graph builder reads (Section 2.3): coordinates
            # for the search radius, elevation as the head proxy inside p_ij
            'lon': float(row['X coordinate'] if 'X coordinate' in row else row['Longitude']),
            'lat': float(row['Y coordinate'] if 'Y coordinate' in row else row['Latitude']),
            'elevation': float(row['Elevation']),
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
    return _build_and_fit(samples, config, "Talensi")

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
            # directional confidence from the probabilistic builder, and the
            # sheaf cohomology diagnostics attached by the refinement
            "edge_confidence": r.edge_confidence,
            "sheaf_h0_dim": r.sheaf_h0_dim,
            "sheaf_h1_dim": r.sheaf_h1_dim,
            "sheaf_obstruction_energy": r.sheaf_obstruction_energy,
            "sheaf_obstruction_leverage": r.sheaf_obstruction_leverage,
            "sheaf_cycle_obstruction_max": r.sheaf_cycle_obstruction_max,
            "sheaf_cycle_count": r.sheaf_cycle_count,
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
