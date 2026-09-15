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

def load_samples():
    manu_df = pd.read_csv(DATA_DIR / "LowerAnayari" / "manu.csv")
    talensi_df = pd.read_csv(DATA_DIR / "Talensi_MiningArea" / "talensi.csv")

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

    return manu_samples, talensi_samples

def main():
    BENCHMARK_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    manu_samples, talensi_samples = load_samples()
    edges_df = pd.read_csv(BENCHMARK_RESULTS_DIR / "field_discovery_results.csv")

    site_configs = {
        "Lower Anayari": {
            "prefix": "Manu",
            "samples": manu_samples,
            "minerals": ["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"],
            "bias": "crystalline",
        },
        "Talensi": {
            "prefix": "Talensi",
            "samples": talensi_samples,
            "minerals": ["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"],
            "bias": "crystalline",
        },
    }

    lambdas = np.logspace(-4, 0, 25)
    all_path_rows = []

    for site_name, cfg_info in site_configs.items():
        prefix = cfg_info["prefix"]
        samples = cfg_info["samples"]
        minerals = cfg_info["minerals"]
        bias = cfg_info["bias"]

        config = Config(
            ion_order=ION_ORDER,
            weights=[1.0] * 10,
            active_minerals=minerals,
            exchange_enabled=True,
            honest_modeling=True,
            geologic_bias=bias,
        )
        matrix, labels, _, penalty_scales = build_reaction_dictionary(config)
        site_edges = edges_df[edges_df["edge_id"].str.startswith(prefix)].copy()

        edge_residuals = []
        for _, row in site_edges.iterrows():
            u, v = str(row["u"]), str(row["v"])
            if u in samples and v in samples:
                u_vec = samples[u]
                v_vec = samples[v]
                gamma = float(row["gamma"]) if pd.notna(row["gamma"]) else 1.0
                resid = [v_val - gamma * u_val for u_val, v_val in zip(u_vec, v_vec)]
                edge_residuals.append(resid)

        empirical_delta = np.median(edge_residuals, axis=0)

        print(f"\nEvaluating Regularization Path for {site_name} ({len(site_edges)} field edges)...")

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
            for resid in edge_residuals:
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
            # We look for the non-trivial parsimonious minimum (k >= 1)
            valid_k = finite_sub[finite_sub["active_k"] >= 1]
            best_row = valid_k.loc[valid_k["aicc"].idxmin()] if not valid_k.empty else finite_sub.loc[finite_sub["aicc"].idxmin()]
            print(f"{site_name} (Empirical Delta): Optimal lambda = {best_row['lambda']:.4f}, min AICc = {best_row['aicc']:.2f}, active_k = {best_row['active_k']}")

if __name__ == "__main__":
    main()
