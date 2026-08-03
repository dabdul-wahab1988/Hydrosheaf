"""
M3 Cross-Validation & Graph Benchmark Script

This script performs true tracer-withholding cross-validation on the USGS datasets.
It drops a specific target tracer from the inversion, recalculates the single-node age,
applies physical (hydraulic proxy & depth-constrained) and randomized graph regularization,
and then forward-models the withheld tracer to compare predictive accuracy against ground truth observations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import concurrent.futures
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np

# Setup paths to ensure hydrosheaf is importable
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# Import Hydrosheaf modules
from hydrosheaf.nuclear.joint_lpm import predict_lpm_tracers
import run_m3_usgs_benchmark as usgs
import run_m3_real_usgs_graph_benchmark as graph_bench

BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"

TRACER_MAP = {
    "SF6": ("SF6", "sf6_pptv"),
    "3H": ("3H", "tritium_TU"),
    "14C": ("14C", "c14_pmc"),
    "CFC11": ("CFC11", "cfc11_pptv"),
    "CFC12": ("CFC12", "cfc12_pptv"),
}

CV_FACTORS = {
    "gas_correction_mode": "usgs_dgm",
    "tracer_set": "reported",
    "age_target_mode": "saturated_only",
    "uztt_mode": "ignore",
    "c14_correction_mode": "selected",
    "he4_mode": "calibrated",
}

def mask_tracer_in_row(row: dict[str, Any], tracer_to_withhold: str) -> dict[str, Any]:
    """Remove a target tracer from the reported list so the inversion ignores it."""
    new_row = dict(row)
    if "LPM_TracersMod" in new_row and isinstance(new_row["LPM_TracersMod"], str):
        parts = [p.strip() for p in new_row["LPM_TracersMod"].replace(",", "|").replace(";", "|").split("|")]
        new_parts = []
        for p in parts:
            if not p: continue
            up = p.upper()
            if tracer_to_withhold == "SF6" and up == "SF6": continue
            if tracer_to_withhold == "3H" and up in ["3H", "3HE(TRIT)", "3H/3HE"]: continue
            if tracer_to_withhold == "14C" and up == "14C": continue
            if tracer_to_withhold == "CFC11" and up in ["CFC11", "CFC-11"]: continue
            if tracer_to_withhold == "CFC12" and up in ["CFC12", "CFC-12"]: continue
            new_parts.append(p)
        new_row["LPM_TracersMod"] = "|".join(new_parts)
    return new_row

def _finite_float(val: Any) -> float | None:
    try:
        f = float(val)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None

def fit_row_task(args_tuple: tuple[int, dict[str, Any], str, int]) -> dict[str, Any] | None:
    idx, row, tracer_to_withhold, age_steps = args_tuple
    masked_row = mask_tracer_in_row(row, tracer_to_withhold)
    try:
        res = usgs.fit_benchmark_row(
            masked_row,
            age_steps=age_steps,
            factors=CV_FACTORS,
            model_strategy="selection",
        )
        res["site_id"] = row["site_id"]
        res["cv_original_row"] = dict(row)
        return res
    except Exception as e:
        print(f"Failed to fit {row['site_id']}: {e}")
        return None


def _withheld_observation(
    row: dict[str, Any] | pd.Series, tracer: str
) -> tuple[float | None, str]:
    """Return the observation on the same correction scale used by the fit."""
    _, canonical_key = TRACER_MAP[tracer]
    if tracer == "14C":
        corrected = _finite_float(row.get("corrected_c14_pmc"))
        if corrected is not None and corrected > 0:
            return corrected, "selected_corrected_c14_pmc"
    elif tracer in {"SF6", "CFC11", "CFC12"}:
        corrected_key = f"dgm_{canonical_key}"
        corrected = _finite_float(row.get(corrected_key))
        if corrected is not None:
            return corrected, corrected_key
    return _finite_float(row.get(canonical_key)), canonical_key


def _is_reportable_fit(row: dict[str, Any] | pd.Series) -> bool:
    """Apply the same reference-free scalar-age guard used by the main benchmark."""
    identifiable = str(row.get("fit_identifiable", False)).lower() == "true"
    tau = _finite_float(row.get("fit_mean_age_years"))
    return identifiable and tau is not None and tau > 0


def _filter_edges_to_nodes(
    edges: list[tuple[str, str]], valid_nodes: set[str]
) -> list[tuple[str, str]]:
    """Prevent target-tracer leakage through non-withheld neighbouring fits."""
    return [(u, v) for u, v in edges if u in valid_nodes and v in valid_nodes]

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="M3 Cross-Validation Benchmark")
    parser.add_argument("--withhold-tracer", required=True, choices=["SF6", "3H", "14C", "CFC11", "CFC12"], help="Tracer to withhold and predict")
    parser.add_argument("--study-unit", default=None, help="Filter to a specific StudyUnit for faster testing")
    parser.add_argument("--age-steps", type=int, default=usgs.M3_DEFAULT_AGE_STEPS)
    parser.add_argument("--max-workers", type=int, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    print("Loading USGS data...")
    df = usgs.load_benchmark_dataset()
    if args.study_unit:
        df = df[df["StudyUnit"] == args.study_unit].copy()
    
    print(f"Data loaded: {len(df)} rows. Withholding {args.withhold_tracer}...")

    # Step 1: Run single node inference without the tracer in parallel
    tasks = []
    for idx, row in df.iterrows():
        orig_tracers_mod = str(row.get("LPM_TracersMod", "")).upper()
        if args.withhold_tracer not in orig_tracers_mod.replace("-", "") and args.withhold_tracer != "3H":
             continue
        if args.withhold_tracer == "3H" and "3H" not in orig_tracers_mod:
             continue
        true_value, _ = _withheld_observation(row, args.withhold_tracer)
        if true_value is None:
            continue
        tasks.append((idx, dict(row), args.withhold_tracer, args.age_steps))

    if not tasks:
        print("No valid rows for CV found.")
        return 0

    print(f"Submitting {len(tasks)} fitting tasks to process pool...")
    cv_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        results = executor.map(fit_row_task, tasks)
        for res in results:
            if res is not None:
                cv_results.append(res)

    if not cv_results:
        print("No valid fits completed.")
        return 0

    all_fit_df = pd.DataFrame(cv_results)
    reportable_mask = all_fit_df.apply(_is_reportable_fit, axis=1)
    res_df = all_fit_df[reportable_mask].copy()
    res_df["fit_mean_age_years"] = pd.to_numeric(res_df["fit_mean_age_years"], errors="coerce")
    
    cv_ages_dict = {str(r["site_id"]): float(r["fit_mean_age_years"]) for _, r in res_df.iterrows()}
    print(
        f"Single-node inference complete: {len(res_df)} reportable withheld fits "
        f"out of {len(all_fit_df)} completed fits ({len(tasks)} eligible rows)."
    )

    # Step 2: Build Graphs using the entire study unit/dataset for full context
    print("Building physical and randomized graphs on the full dataset...")
    mock_df = pd.DataFrame()
    mock_df["site_id"] = df["site_id"]
    mock_df["StudyUnit"] = df["StudyUnit"]
    mock_df["lat"] = df.apply(lambda r: r.get("lat", r.get("LatDD83")), axis=1)
    mock_df["lon"] = df.apply(lambda r: r.get("lon", r.get("LongDD83")), axis=1)
    mock_df["depth_m"] = df.apply(lambda r: r.get("depth_m", r.get("Depth_m")), axis=1)
    
    # Only target-withheld, reportable fits may supply ages. Including full-fit
    # neighbours here would leak the withheld tracer back through the graph.
    mock_df["est_age_multi"] = mock_df["site_id"].astype(str).map(cv_ages_dict)
    mock_df["ref_age"] = mock_df["est_age_multi"] # placeholder for build_graph_families
    
    family_edges, edge_rows = graph_bench.build_graph_families(mock_df, min_unit_size=5)
    
    valid_nodes = set(cv_ages_dict)
    phys_hydraulic_edges = _filter_edges_to_nodes(
        family_edges.get("hydraulic_proxy_constrained", []), valid_nodes
    )
    phys_depth_edges = _filter_edges_to_nodes(
        family_edges.get("depth_constrained", []), valid_nodes
    )
    random_edges = _filter_edges_to_nodes(
        family_edges.get("randomized_negative_control", []), valid_nodes
    )

    print(f"Edges built. Hydraulic Proxy: {len(phys_hydraulic_edges)}, Depth Constrained: {len(phys_depth_edges)}, Random: {len(random_edges)}")

    # Step 3: Regularize Taus across the entire graph
    taus = {
        str(sid): float(age)
        for sid, age in zip(mock_df["site_id"], mock_df["est_age_multi"])
        if math.isfinite(age)
    }
    
    phys_hyd_taus = graph_bench.graph_regularize(taus, phys_hydraulic_edges, strength=0.85)
    phys_dep_taus = graph_bench.graph_regularize(taus, phys_depth_edges, strength=0.85)
    rand_taus = graph_bench.graph_regularize(taus, random_edges, strength=0.85)
    
    # Step 4: Forward Modeling & Performance Evaluation (only for the CV nodes)
    print("Predicting withheld tracers...")
    cv_records = []
    
    pred_key, _ = TRACER_MAP[args.withhold_tracer]

    for _, r in res_df.iterrows():
        sid = str(r["site_id"])
        orig = r["cv_original_row"]
        model = str(r.get("multi_model", "PFM"))
        c14_a0 = _finite_float(r.get("c14_initial_pmc")) or 100.0
        sec_name = r.get("fit_secondary_param_name")
        sec_val = _finite_float(r.get("fit_secondary_param_value"))

        def get_pred(tau_val):
            if not math.isfinite(tau_val): return np.nan
            try:
                sample_year = _finite_float(orig.get("sample_year")) or 2020.0
                histories, _ = usgs._get_site_histories(orig)
                max_age = max(200.0, 5.0 * tau_val)
                params = {"mean_age_years": tau_val}
                if sec_name and sec_val is not None:
                    params[sec_name] = sec_val
                preds = predict_lpm_tracers(
                    model=model,
                    parameters=params,
                    sample_year=sample_year,
                    tracers=[pred_key],
                    histories=histories,
                    initial_c14_pmc=c14_a0,
                    max_age_years=max_age
                )
                return preds.get(pred_key, np.nan)
            except Exception:
                return np.nan
                
        base_pred = get_pred(taus[sid])
        phys_hyd_pred = get_pred(phys_hyd_taus.get(sid, taus[sid]))
        phys_dep_pred = get_pred(phys_dep_taus.get(sid, taus[sid]))
        rand_pred = get_pred(rand_taus.get(sid, taus[sid]))
        
        true_val, target_source = _withheld_observation(orig, args.withhold_tracer)
        
        if true_val is not None and math.isfinite(true_val):
            cv_records.append({
                "site_id": sid,
                "StudyUnit": r["StudyUnit"],
                "tracer": args.withhold_tracer,
                "true_val": true_val,
                "target_source": target_source,
                "tau_single": taus[sid],
                "tau_phys_hyd": phys_hyd_taus.get(sid, taus[sid]),
                "tau_phys_dep": phys_dep_taus.get(sid, taus[sid]),
                "tau_rand_graph": rand_taus.get(sid, taus[sid]),
                "pred_single": base_pred,
                "pred_phys_hyd": phys_hyd_pred,
                "pred_phys_dep": phys_dep_pred,
                "pred_rand_graph": rand_pred,
                "err_single_abs": abs(true_val - base_pred) if math.isfinite(base_pred) else np.nan,
                "err_phys_hyd_abs": abs(true_val - phys_hyd_pred) if math.isfinite(phys_hyd_pred) else np.nan,
                "err_phys_dep_abs": abs(true_val - phys_dep_pred) if math.isfinite(phys_dep_pred) else np.nan,
                "err_rand_abs": abs(true_val - rand_pred) if math.isfinite(rand_pred) else np.nan,
                "fit_rmse_standardized": r.get("fit_rmse_standardized"),
                "fit_effective_n_params": r.get("fit_effective_n_params"),
                "fit_age_profile_log10_span": r.get("fit_age_profile_log10_span"),
                "fit_identifiability_reason": r.get("fit_identifiability_reason", ""),
            })

    if not cv_records:
        print("No valid predictions could be made.")
        return 0
        
    out_df = pd.DataFrame(cv_records)
    
    def rmse(errs):
        errs = errs.dropna()
        if len(errs) == 0: return np.nan
        return float(np.sqrt(np.mean(errs**2)))
        
    summary = {
        "withheld_tracer": args.withhold_tracer,
        "n_samples": len(out_df),
        "rmse_single": rmse(out_df["err_single_abs"]),
        "rmse_phys_hyd": rmse(out_df["err_phys_hyd_abs"]),
        "rmse_phys_dep": rmse(out_df["err_phys_dep_abs"]),
        "rmse_rand_graph": rmse(out_df["err_rand_abs"]),
    }
    
    print("\n--- Cross-Validation Summary ---")
    print(json.dumps(summary, indent=2))
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = args.output or (RESULTS_DIR / f"m3_cv_benchmark_{args.withhold_tracer}.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    manifest = {
        "schema_version": "1.0",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "withheld_tracer": args.withhold_tracer,
        "age_steps": args.age_steps,
        "model_strategy": "selection",
        "factors": CV_FACTORS,
        "reportability_guard": "reference-free M3 scalar-age guard",
        "leakage_control": "all graph ages were fitted with the target tracer withheld",
        "eligible_rows": len(tasks),
        "completed_fits": len(all_fit_df),
        "reportable_fits": len(res_df),
        "prediction_rows": len(out_df),
        "graph_edges": {
            "hydraulic_proxy_constrained": len(phys_hydraulic_edges),
            "depth_constrained": len(phys_depth_edges),
            "randomized_negative_control": len(random_edges),
        },
        "summary": summary,
        "output": str(out_path),
        "output_sha256": hashlib.sha256(out_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
    }
    manifest_path = out_path.with_name(f"{out_path.stem}_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nWrote full CV results to {out_path}")
    
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
