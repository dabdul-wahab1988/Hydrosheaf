"""Execute field network integration pipeline on the completed Upper East Region (UER) dataset.

Demonstrates:
1. Loading the harmonized, elevation-complete field dataset.
2. Building the directed hydraulic flow graph using completed DEM elevations.
3. Evaluating geochemical ratios and mass-transfer diagnostics along inferred flowpaths.
4. Exporting network edge tables, summary statistics, and integration report.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf.graph.build import infer_edges_from_coordinates, infer_edges_probabilistic
from hydrosheaf.models.ratios import compare_ratio_diagnostics, compute_geochemical_ratios

ROOT = Path(__file__).resolve().parents[2]
DATA_CSV = ROOT / "data" / "FieldData" / "derived" / "uer_field_integration_dataset.csv"
OUT_DIR = ROOT / "outputs" / "uer_field_integration"


def main() -> None:
    raise RuntimeError(
        "Retired: this legacy UER network promotes DEM elevation to hydraulic "
        "head and labels topographic candidates as groundwater flowpaths. The "
        "completed UER workbook contains no measured head series or independent "
        "flow-path truth, so this output is not regenerated."
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Loading field integration dataset from {DATA_CSV.name}...")
    df = pd.read_csv(DATA_CSV)
    print(f"Loaded {len(df)} samples across {df['community'].nunique()} communities.")

    # Prepare node records
    samples = []
    sample_map = {}
    for _, row in df.iterrows():
        s_id = str(row["node_id"])
        record = {
            "node_id": s_id,
            "site_id": s_id,
            "community": row["community"],
            "sample_type": row["sample_type"],
            "lat": float(row["latitude_dd"]),
            "lon": float(row["longitude_dd"]),
            "elevation": float(row["elevation_m"]),
            "hydraulic_head": float(row["elevation_m"]),
            "pH": float(row["pH"]) if pd.notna(row["pH"]) else None,
            "EC": float(row["ec_uS_cm"]) if pd.notna(row["ec_uS_cm"]) else None,
            "TDS": float(row["tds_mg_L"]) if pd.notna(row["tds_mg_L"]) else None,
            "Ca": float(row["ca_mmol_L"]) if pd.notna(row["ca_mmol_L"]) else None,
            "Mg": float(row["mg_mmol_L"]) if pd.notna(row["mg_mmol_L"]) else None,
            "Na": float(row["na_mmol_L"]) if pd.notna(row["na_mmol_L"]) else None,
            "K": float(row["k_mmol_L"]) if pd.notna(row["k_mmol_L"]) else None,
            "HCO3": float(row["hco3_mmol_L"]) if pd.notna(row["hco3_mmol_L"]) else None,
            "Cl": float(row["cl_mmol_L"]) if pd.notna(row["cl_mmol_L"]) else None,
            "SO4": float(row["so4_mmol_L"]) if pd.notna(row["so4_mmol_L"]) else None,
            "NO3": float(row["no3_mmol_L"]) if pd.notna(row["no3_mmol_L"]) else None,
            "F": float(row["f_mmol_L"]) if pd.notna(row["f_mmol_L"]) else None,
            "18O": float(row["d18O_permil"]) if pd.notna(row["d18O_permil"]) else None,
            "2H": float(row["d2H_permil"]) if pd.notna(row["d2H_permil"]) else None,
            "tritium": float(row["tritium_TU"]) if pd.notna(row["tritium_TU"]) else None,
            "geology": str(row["geology_symbol"]) if pd.notna(row["geology_symbol"]) else "Unknown",
        }
        samples.append(record)
        sample_map[s_id] = record

    # 1. Infer spatial-hydraulic edges (probabilistic Tier C head from DEM)
    print("\nInferring probabilistic flow edges (search radius = 25 km, p_min = 0.6)...")
    edges_prob = infer_edges_probabilistic(
        samples,
        radius_km=25.0,
        max_neighbors=3,
        p_min=0.6,
        sigma_topo=10.0,
        gradient_min=1e-4,
    )
    print(f"  Inferred {len(edges_prob)} directed flow edges.")

    # 2. Evaluate geochemical evolution along edges
    print("\nEvaluating geochemical ratio evolution along flowpaths...")
    edge_rows = []
    ratio_similarities = []

    for edge in edges_prob:
        u_rec = sample_map[edge.u]
        v_rec = sample_map[edge.v]

        u_ratios = compute_geochemical_ratios(u_rec)
        v_ratios = compute_geochemical_ratios(v_rec)
        comp = compare_ratio_diagnostics(u_rec, v_rec)

        sim = comp.get("similarity")
        if sim is not None:
            ratio_similarities.append(sim)

        # Delta TDS and Delta Cl
        delta_tds = (v_rec["TDS"] - u_rec["TDS"]) if (v_rec["TDS"] and u_rec["TDS"]) else None
        delta_cl = (v_rec["Cl"] - u_rec["Cl"]) if (v_rec["Cl"] and u_rec["Cl"]) else None
        delta_h = edge.attrs.get("delta_h", u_rec["elevation"] - v_rec["elevation"])
        dist_km = edge.attrs.get("distance_km", 0.0)
        grad = delta_h / (dist_km * 1000.0) if dist_km > 0 else 0.0

        edge_rows.append(
            {
                "edge_id": edge.edge_id,
                "u": edge.u,
                "v": edge.v,
                "u_community": u_rec["community"],
                "v_community": v_rec["community"],
                "u_type": u_rec["sample_type"],
                "v_type": v_rec["sample_type"],
                "u_elevation_m": u_rec["elevation"],
                "v_elevation_m": v_rec["elevation"],
                "delta_h_m": delta_h,
                "distance_km": dist_km,
                "hydraulic_gradient": grad,
                "p_uv": edge.attrs.get("p_uv"),
                "edge_confidence": edge.attrs.get("edge_confidence"),
                "ratio_similarity": sim,
                "delta_tds_mg_L": delta_tds,
                "delta_cl_mmol_L": delta_cl,
                "u_geology": u_rec["geology"],
                "v_geology": v_rec["geology"],
                "same_geology": u_rec["geology"] == v_rec["geology"],
            }
        )

    df_edges = pd.DataFrame(edge_rows)
    edges_csv = OUT_DIR / "uer_network_edges.csv"
    df_edges.to_csv(edges_csv, index=False)
    print(f"Wrote {len(df_edges)} flow edges to {edges_csv.name}")

    # Summary statistics
    summary = {
        "dataset": "northern_ghana_new",
        "total_nodes": len(samples),
        "elevation_completeness_percent": 100.0,
        "elevation_range_m": {
            "min": float(df["elevation_m"].min()),
            "mean": float(df["elevation_m"].mean()),
            "max": float(df["elevation_m"].max()),
        },
        "total_inferred_edges": len(df_edges),
        "mean_edge_distance_km": float(df_edges["distance_km"].mean()) if not df_edges.empty else 0.0,
        "mean_delta_h_m": float(df_edges["delta_h_m"].mean()) if not df_edges.empty else 0.0,
        "mean_hydraulic_gradient": float(df_edges["hydraulic_gradient"].mean()) if not df_edges.empty else 0.0,
        "mean_ratio_similarity": float(np.mean(ratio_similarities)) if ratio_similarities else 0.0,
        "same_geology_edge_fraction": float(df_edges["same_geology"].mean()) if not df_edges.empty else 0.0,
    }

    summary_json = OUT_DIR / "uer_integration_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote summary metrics to {summary_json.name}")

    # Markdown Report
    report_md = OUT_DIR / "uer_field_integration_report.md"
    report_content = f"""# Upper East Region (Northern Ghana New) Field Integration Report

## 1. Dataset Overview & Elevation Completion
- **Total sample locations:** {summary['total_nodes']}
- **Elevation completeness:** 100.0% (237/237)
- **Elevation distribution:** Min = {summary['elevation_range_m']['min']:.1f} m, Mean = {summary['elevation_range_m']['mean']:.1f} m, Max = {summary['elevation_range_m']['max']:.1f} m a.s.l.
- **Elevation sources:** 2 surveyed ground controls reconciled with NASA SRTM 30m DEM (RMSE = 2.24 m); 235 extracted from SRTM 30m DEM with ASTER GDEM v3 verification.

## 2. Directed Hydraulic Graph Construction
- **Total inferred directed edges:** {summary['total_inferred_edges']}
- **Mean flowpath length:** {summary['mean_edge_distance_km']:.2f} km
- **Mean elevation drop (\\Delta h):** {summary['mean_delta_h_m']:.2f} m
- **Mean hydraulic gradient (i):** {summary['mean_hydraulic_gradient']:.4f} m/m
- **Same-geology flowpath share:** {summary['same_geology_edge_fraction']:.1%}

## 3. Geochemical Transport & Ratio Diagnostics
- **Mean geochemical ratio similarity along flowpaths:** {summary['mean_ratio_similarity']:.3f}
- **Primary chemistry completeness:** 100% (237/237) for Ca, Mg, Na, K, HCO3, Cl, SO4, NO3.
- **Charge balance acceptability:** 97.5% of samples within acceptable/caution limits (|CBE| <= 10%).

## 4. Integration Artifacts
- **Completed Excel Workbook:** `data/FieldData/NorthernGhanaNew/compiled UER data_new_completed.xlsx`
- **Canonical Field CSV:** `data/FieldData/derived/uer_field_integration_dataset.csv`
- **DEM Elevation Sidecar:** `data/FieldData/derived/uer_elevations_dem.csv`
- **Network Edge Table:** `outputs/uer_field_integration/uer_network_edges.csv`
"""
    with open(report_md, "w", encoding="utf-8") as f:
        f.write(report_content)
    print(f"Wrote field integration report to {report_md.name}")
    print("\nField integration pipeline completed successfully!")


if __name__ == "__main__":
    main()
