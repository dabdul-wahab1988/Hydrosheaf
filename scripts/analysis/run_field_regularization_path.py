"""Field Regularization Path Analysis for Lower Anayari and Talensi.

Evaluates coordinate-descent sparse reaction inversion across regularisation
strengths (lambda in [1e-4, 1.0]) using exclusively the real field water
chemistry datasets and site-specific geology-aware mineral dictionaries.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PROJECT_ROOT = REPO_ROOT
DATA_DIR = PROJECT_ROOT / "data" / "FieldData"
BENCHMARK_RESULTS_DIR = PROJECT_ROOT / "M2" / "m2_benchmark" / "results"

from hydrosheaf.config import Config
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.models.reactions import build_reaction_dictionary, fit_reactions

ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]

def load_all_sites_data():
    raise RuntimeError(
        "Retired: this field regularization analysis depends on excluded legacy "
        "cohorts and zero-imputes missing chemistry. No approved field edge/reaction "
        "truth exists for a replacement run."
    )
    manu_path = PROJECT_ROOT / "manu.csv" if (PROJECT_ROOT / "manu.csv").exists() else DATA_DIR / "LowerAnayari" / "manu.csv"
    talensi_path = PROJECT_ROOT / "talensi.csv" if (PROJECT_ROOT / "talensi.csv").exists() else DATA_DIR / "Talensi_MiningArea" / "talensi.csv"
    manu_df = pd.read_csv(manu_path)
    talensi_df = pd.read_csv(talensi_path)
    uer_df = pd.read_csv(DATA_DIR / "derived" / "uer_field_integration_dataset.csv")
    cr_df = pd.read_csv(DATA_DIR / "derived" / "cr_field_integration_dataset.csv")

    manu_samples = {}
    for _, row in manu_df.iterrows():
        sid = str(row["Sample ID"])
        vec = []
        for ion in ION_ORDER:
            if ion in row and pd.notna(row[ion]):
                val = row[ion]
                if isinstance(val, str) and "<" in val:
                    val = 0.0005
                vec.append(mgL_to_mmolL(float(val), ion))
            else:
                vec.append(0.0)
        manu_samples[sid] = vec
        manu_samples[f"Manu_{sid}"] = vec

    talensi_samples = {}
    for _, row in talensi_df.iterrows():
        code = str(row["Code"])
        vec = []
        for ion in ION_ORDER:
            if ion in row and pd.notna(row[ion]):
                vec.append(mgL_to_mmolL(float(row[ion]), ion))
            else:
                vec.append(0.0)
        talensi_samples[code] = vec
        talensi_samples[f"Talensi_{code}"] = vec

    uer_samples = {}
    for _, row in uer_df.iterrows():
        nid = str(row["node_id"])
        uer_samples[nid] = [
            float(row["ca_mmol_L"]) if pd.notna(row["ca_mmol_L"]) else 0.0,
            float(row["mg_mmol_L"]) if pd.notna(row["mg_mmol_L"]) else 0.0,
            float(row["na_mmol_L"]) if pd.notna(row["na_mmol_L"]) else 0.0,
            float(row["k_mmol_L"]) if pd.notna(row["k_mmol_L"]) else 0.0,
            float(row["hco3_mmol_L"]) if pd.notna(row["hco3_mmol_L"]) else 0.0,
            float(row["cl_mmol_L"]) if pd.notna(row["cl_mmol_L"]) else 0.0,
            float(row["so4_mmol_L"]) if pd.notna(row["so4_mmol_L"]) else 0.0,
            float(row["no3_mmol_L"]) if pd.notna(row["no3_mmol_L"]) else 0.0,
            float(row["f_mmol_L"]) if pd.notna(row["f_mmol_L"]) else 0.0,
            0.0,
        ]

    cr_samples = {}
    cr_nodes = []
    for _, row in cr_df.iterrows():
        nid = str(row["node_id"])
        fe_mmol = (float(row["fe_mg_L"]) / 55.845) if pd.notna(row.get("fe_mg_L")) else 0.0
        cr_samples[nid] = [
            float(row["ca_mmol_L"]) if pd.notna(row["ca_mmol_L"]) else 0.0,
            float(row["mg_mmol_L"]) if pd.notna(row["mg_mmol_L"]) else 0.0,
            float(row["na_mmol_L"]) if pd.notna(row["na_mmol_L"]) else 0.0,
            float(row["k_mmol_L"]) if pd.notna(row["k_mmol_L"]) else 0.0,
            float(row["hco3_mmol_L"]) if pd.notna(row["hco3_mmol_L"]) else 0.0,
            float(row["cl_mmol_L"]) if pd.notna(row["cl_mmol_L"]) else 0.0,
            float(row["so4_mmol_L"]) if pd.notna(row["so4_mmol_L"]) else 0.0,
            float(row["no3_mmol_L"]) if pd.notna(row["no3_mmol_L"]) else 0.0,
            0.0,
            fe_mmol,
        ]
        cr_nodes.append({
            "site_id": nid,
            "lat": float(row["latitude_dd"]),
            "lon": float(row["longitude_dd"]),
            "head_meas": float(row["hydraulic_head_m"]) if pd.notna(row.get("hydraulic_head_m")) else None,
            "elevation": float(row["elevation_m"]) if pd.notna(row.get("elevation_m")) else 50.0,
            "dtw": float(row["swl_m"]) if pd.notna(row.get("swl_m")) else None,
        })

    return manu_samples, talensi_samples, uer_samples, cr_samples, cr_nodes


def main():
    BENCHMARK_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    manu_samples, talensi_samples, uer_samples, cr_samples, cr_nodes = load_all_sites_data()
    edges_df = pd.read_csv(BENCHMARK_RESULTS_DIR / "field_discovery_results.csv")

    # Load or infer UER and CR edges
    from hydrosheaf.graph.build import infer_edges_probabilistic

    uer_edges_path = PROJECT_ROOT / "outputs" / "uer_field_integration" / "uer_network_edges.csv"
    if uer_edges_path.exists():
        uer_edges_df = pd.read_csv(uer_edges_path)
    else:
        uer_edges_df = pd.DataFrame()

    cr_edges = infer_edges_probabilistic(
        cr_nodes,
        radius_km=25.0,
        max_neighbors=3,
        p_min=0.6,
        sigma_meas=0.5,
        sigma_dtw=1.0,
        sigma_elev=1.0,
        sigma_topo=10.0,
        gradient_min=1e-4,
    )

    site_configs = {
        "Lower Anayari": {
            "samples": manu_samples,
            "edges": [(str(r["u"]), str(r["v"]), float(r["gamma"]) if pd.notna(r.get("gamma")) else 1.0)
                      for _, r in edges_df[edges_df["edge_id"].str.startswith("Manu")].iterrows()],
            "minerals": ["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"],
            "bias": "crystalline",
            "delta_mode": "median",
        },
        "Talensi": {
            "samples": talensi_samples,
            "edges": [(str(r["u"]), str(r["v"]), float(r["gamma"]) if pd.notna(r.get("gamma")) else 1.0)
                      for _, r in edges_df[edges_df["edge_id"].str.startswith("Talensi")].iterrows()],
            "minerals": ["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"],
            "bias": "crystalline",
            "delta_mode": "median",
        },
        "Upper East Region": {
            "samples": uer_samples,
            "edges": [(str(r["u"]), str(r["v"]), 1.0) for _, r in uer_edges_df.iterrows()],
            "minerals": ["calcite", "dolomite", "fluorite", "albite", "halite", "gypsum"],
            "bias": "crystalline",
            "delta_mode": "median",
        },
        "Central Region": {
            "samples": cr_samples,
            "edges": [(e.u, e.v, 1.0) for e in cr_edges],
            "minerals": ["calcite", "dolomite", "albite", "halite", "gypsum"],
            "bias": "crystalline",
            "delta_mode": "mean",
        },
    }

    lambdas = np.logspace(-4, 0, 25)
    all_path_rows = []

    for site_name, cfg_info in site_configs.items():
        samples = cfg_info["samples"]
        edge_list = cfg_info["edges"]
        minerals = cfg_info["minerals"]
        bias = cfg_info["bias"]
        delta_mode = cfg_info["delta_mode"]

        config = Config(
            ion_order=ION_ORDER,
            weights=[1.0] * 10,
            active_minerals=minerals,
            exchange_enabled=True,
            honest_modeling=True,
            geologic_bias=bias,
        )
        matrix, labels, _, penalty_scales = build_reaction_dictionary(config)

        edge_residuals = []
        for u, v, gamma in edge_list:
            if u in samples and v in samples:
                u_vec = samples[u]
                v_vec = samples[v]
                resid = [v_val - gamma * u_val for u_val, v_val in zip(u_vec, v_vec)]
                edge_residuals.append(resid)

        if not edge_residuals:
            print(f"Warning: No valid edges found for {site_name}")
            continue

        if delta_mode == "mean":
            empirical_delta = np.mean(edge_residuals, axis=0)
        else:
            empirical_delta = np.median(edge_residuals, axis=0)

        print(f"\nEvaluating Regularization Path for {site_name} ({len(edge_residuals)} field edges)...")

        # Sample up to 100 edges for edge-level statistics to ensure fast convergence
        np.random.seed(42)
        sub_indices = np.random.choice(len(edge_residuals), size=min(100, len(edge_residuals)), replace=False)
        sample_residuals = [edge_residuals[i] for i in sub_indices]

        for l in lambdas:
            fit = fit_reactions(
                empirical_delta,
                matrix,
                [1.0] * 10,
                lambda_l1=l,
                lambda_l2=config.lambda_l2,
                penalty_scales=penalty_scales,
            )

            edge_aicc_list = []
            edge_rss_list = []
            edge_k_list = []
            for resid in sample_residuals:
                efit = fit_reactions(
                    resid,
                    matrix,
                    [1.0] * 10,
                    lambda_l1=l,
                    lambda_l2=config.lambda_l2,
                    penalty_scales=penalty_scales,
                )
                if np.isfinite(efit.aicc):
                    edge_aicc_list.append(efit.aicc)
                edge_rss_list.append(efit.residual_norm)
                edge_k_list.append(len([e for e in efit.extents if abs(e) > 1e-6]))

            row = {
                "site": site_name,
                "lambda": l,
                "residual_norm": fit.residual_norm,
                "aicc": fit.aicc,
                "active_k": len([e for e in fit.extents if abs(e) > 1e-6]),
                "mean_edge_aicc": np.mean(edge_aicc_list) if edge_aicc_list else np.nan,
                "median_edge_aicc": np.median(edge_aicc_list) if edge_aicc_list else np.nan,
                "mean_edge_rss": np.mean(edge_rss_list),
                "median_edge_rss": np.median(edge_rss_list),
                "mean_edge_k": np.mean(edge_k_list),
            }
            for lbl, extent in zip(labels, fit.extents):
                row[lbl] = extent
            all_path_rows.append(row)

    out_df = pd.DataFrame(all_path_rows)
    out_path = BENCHMARK_RESULTS_DIR / "field_regularization_path.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved field regularization path to {out_path}")

    for site_name in site_configs.keys():
        sub = out_df[out_df["site"] == site_name]
        finite_sub = sub[np.isfinite(sub["aicc"])]
        if not finite_sub.empty:
            valid_k = finite_sub[finite_sub["active_k"] >= 1]
            best_row = valid_k.loc[valid_k["aicc"].idxmin()] if not valid_k.empty else finite_sub.loc[finite_sub["aicc"].idxmin()]
            print(f"{site_name} (Empirical Delta): Optimal lambda = {best_row['lambda']:.4f}, min AICc = {best_row['aicc']:.2f}, active_k = {best_row['active_k']}")


if __name__ == "__main__":
    main()
