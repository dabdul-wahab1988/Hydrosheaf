
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error
import os
import sys
from typing import Callable, Dict, cast
from pathlib import Path

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).parents[1]))

from hydrosheaf.calibration.glm import run_pestglm_multistart
from hydrosheaf.calibration.pestpp.runner import run_pestpp
from hydrosheaf.transport.validation import TransportValidationProblem

def evaluate_fit(problem, param_values):
    obs_targets = problem.get_observation_targets()
    obs_names = list(obs_targets.keys())
    obs_vals = np.array([obs_targets.get(name, -9999) for name in obs_names])
    sim_vals = np.array([problem.run_model_raw(param_values).get(name, -9999) for name in obs_names])

    valid = sim_vals != -9999
    if np.any(valid):
        obs_vals = obs_vals[valid]
        sim_vals = sim_vals[valid]
        r2 = r2_score(obs_vals, sim_vals)
        rmse = np.sqrt(mean_squared_error(obs_vals, sim_vals))
        nse = 1 - (np.sum((obs_vals - sim_vals)**2) / np.sum((obs_vals - np.mean(obs_vals))**2))
    else:
        r2, rmse, nse = 0, 0, 0
    return r2, rmse, nse


def run_validation(
    use_pestpp=True,
    use_second_stage=True,
    use_seasonal=True,
    use_hybrid=True,
    fit_lag=True,
    max_nfev=80,
    n_starts=6,
    pestpp_engine="pestpp-ies"
):
    engine_label = "PESTGLM"
    if use_pestpp:
        engine_label = f"PEST++ ({pestpp_engine})"
    print(f"Running Calibrated Validation (1D ADE) using {engine_label}...")
    print("Running Calibrated Validation (1D ADE)...")
    if use_second_stage or use_seasonal or use_hybrid:
        flags = []
        if use_second_stage:
            flags.append("two-stage")
        if use_seasonal:
            flags.append("seasonal")
        if use_hybrid:
            flags.append("hybrid")
        print(f"Model options: {', '.join(flags)}")

    # 1. Load Data
    chem_df = pd.read_csv("hydrosheaf_synthetic_csv/water_chem_full.csv")
    chem_df["collection_date"] = pd.to_datetime(chem_df["collection_date"])
    
    # 2. Define Pairs (Source -> Receptor)
    # Based on Vea Dataset layout: L1 -> BH1, L2 -> BH2
    # Distance x = 500m (L=0, BH=500)
    pairs = [
        ("L1", "BH1"),
        ("L2", "BH2"),
    ]
    
    DISTANCE_X = 500.0 # meters

    validation_results = []
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, (source, target) in enumerate(pairs):
        ax = axes[idx]
        print(f"\nCalibrating Pair: {source} -> {target}")

        # Get Data
        # Group by date to handle multiple depths (30cm/60cm)
        src_df_raw = chem_df[chem_df["station_code"] == source]
        src_df = src_df_raw.groupby("collection_date")["NO3_mg_L"].mean().reset_index().sort_values("collection_date")
        
        tgt_df_raw = chem_df[chem_df["station_code"] == target]
        tgt_df = tgt_df_raw.groupby("collection_date")["NO3_mg_L"].mean().reset_index().sort_values("collection_date")
        
        # Prepare Time Grid (Daily)
        start_date = min(src_df["collection_date"].min(), tgt_df["collection_date"].min())
        end_date = max(src_df["collection_date"].max(), tgt_df["collection_date"].max())
        
        days_total = (end_date - start_date).days + 100 # buffer
        t_grid = np.arange(days_total)
        
        # Interpolate Source on Grid
        src_df["days"] = (src_df["collection_date"] - start_date).dt.days
        # Use linear interpolation for input signal
        # Fill NaNs with mean or 0? Linear interp handles gaps.
        # Create full index
        full_idx = pd.Index(t_grid, name="days")
        
        # Resample source to daily
        # First make sure we have valid floats
        src_valid = src_df.dropna(subset=["NO3_mg_L"])
        
        # We need a continuous input signal C_in(t)
        # Use numpy interp
        C_in = np.interp(t_grid, src_valid["days"], src_valid["NO3_mg_L"])
        
        # Target Observations
        tgt_valid = tgt_df.dropna(subset=["NO3_mg_L"]).copy()
        tgt_valid["days"] = (tgt_valid["collection_date"] - start_date).dt.days
        
        # --- Run Calibration ---
        
        # Create Problem Instance
        problem = TransportValidationProblem(
            t_grid,
            C_in,
            tgt_valid,
            DISTANCE_X,
            use_log=not use_hybrid,
            fit_lag=fit_lag,
            use_second_stage=use_second_stage,
            use_seasonal=use_seasonal,
            use_hybrid=use_hybrid
        )
        
        opt_p = {}
        
        if use_pestpp:
            # Run PEST++
            # We need to ensure the class is pickleable and available in worker
            # The runner handles generating worker.py that imports hydrosheaf
            
            pestpp_options = {}
            if pestpp_engine == "pestpp-ies":
                pestpp_options["ies_num_reals"] = 4

            pest_res = run_pestpp(
                problem,
                engine=pestpp_engine,
                work_dir=f"pest_workspace_{source}_{target}",
                case_name="transport",
                max_nfev=max_nfev,
                pestpp_options=pestpp_options if pestpp_options else None,
            )
            
            if pest_res["success"]:
                opt_p = pest_res.get("optimal_parameters", {})
                print(f"  Optimization Success: {pest_res['success']}")
                # If optimal_parameters is empty (e.g. failed to read par file), fallback?
                if not opt_p:
                     print("  Warning: No optimal parameters returned from PEST++.")
                     # Fallback to initial values from problem to avoid crash
                     opt_p = {p.name: p.value for p in problem.get_parameters()}
            else:
                print("  PEST++ Optimization Failed. Using initial parameters.")
                opt_p = {p.name: p.value for p in problem.get_parameters()}
                
        else:
            # Use Internal PESTGLM with multi-start
            def score_fn(params: Dict[str, float]) -> float:
                return float(evaluate_fit(problem, params)[2])
            score_fn_typed = cast(Callable[[Dict[str, float]], float], score_fn)
            multistart = run_pestglm_multistart(
                problem,
                n_starts=n_starts,
                max_nfev=max_nfev,
                seed=42,
                n_workers=1,
                worker_type="thread",
                score_fn=score_fn_typed,
            )

            for start in multistart["start_results"]:
                start_id = start.get("start", 0)
                if not start.get("success"):
                    error = start.get("error", "failed")
                    print(f"  Start {start_id}/{n_starts}: failed ({error})")
                    continue

                score = start.get("score", float("nan"))
                if score is None or not np.isfinite(score):
                    print(f"  Start {start_id}/{n_starts}: NSE=nan")
                    continue
                print(f"  Start {start_id}/{n_starts}: success=True, NSE={score:.4f}")

            opt_p = multistart["best_parameters"]
            best_result = multistart.get("best_result")
            print(f"  Optimization Success: {bool(best_result and best_result.get('success'))}")

        # Check if we got parameters
        if not opt_p:
             opt_p = {p.name: p.value for p in problem.get_parameters()}

        v = opt_p.get("v", problem.default_params["v"])
        alpha = opt_p.get("alpha", problem.default_params["alpha"])
        k = opt_p.get("k", problem.default_params["k"])
        lag = opt_p.get("lag", problem.default_params["lag"])
        tau = opt_p.get("tau", problem.default_params["tau"])
        dilution = opt_p.get("dilution", problem.default_params["dilution"])

        params_text = (
            f"v={v:.2f} m/d, alpha={alpha:.1f} m, k={k:.4f}, "
            f"lag={lag:.1f} d, tau={tau:.1f} d, dil={dilution:.2f}"
        )
        if use_second_stage:
            params_text += (
                f", v2={opt_p.get('v2',0):.2f} m/d, alpha2={opt_p.get('alpha2',0):.1f} m, "
                f"k2={opt_p.get('k2',0):.4f}"
            )
        if use_seasonal:
            params_text += (
                f", dil_amp={opt_p.get('dilution_amp',0):.2f}, "
                f"dil_phase={opt_p.get('dilution_phase',0):.1f} d, "
                f"base_amp={opt_p.get('baseflow_amp',0):.2f}, "
                f"base_phase={opt_p.get('baseflow_phase',0):.1f} d"
            )
        print(f"  Optimal Params: {params_text}")
        
        # --- Final Prediction ---
        # Rerun model with optimal params for plotting
        final_preds_raw = problem.run_model_raw(opt_p)
        
        # Get full time series for plotting (Manual reconstruction using optimal params)
        C_out_full = problem.simulate_series(opt_p)
        
        # Dates for plotting
        plot_dates = [start_date + pd.Timedelta(days=d) for d in t_grid]
        
        # Stats
        r2, rmse, nse = evaluate_fit(problem, opt_p)
        
        print(f"  NSE: {nse:.4f}, RMSE: {rmse:.2f}")
        
        # --- Plotting ---
        ax.plot(src_valid["collection_date"], src_valid["NO3_mg_L"], 'k:', alpha=0.3, label=f"Input ({source})")
        ax.plot(tgt_valid["collection_date"], tgt_valid["NO3_mg_L"], 'ko', label=f"Observed ({target})")
        ax.plot(plot_dates, C_out_full, 'r-', linewidth=2, label="Calibrated ADE")
        
        stats_text = (
            f"NSE: {nse:.2f}\nRMSE: {rmse:.1f}\n"
            f"v: {v:.1f} m/d\n$\\alpha$: {alpha:.0f} m\n"
            f"lag: {lag:.0f} d\n$\\tau$: {tau:.0f} d"
        )
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_title(f"Calibrated Transport: {source} -> {target}")
        ax.set_ylabel("NO3 (mg/L)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.setp(ax.get_xticklabels(), rotation=30)
        
        validation_results.append({
            "pair": f"{source}->{target}",
            "rmse": rmse,
            "r2": r2,
            "nse": nse,
            "v": v,
            "alpha": alpha
        })

    # Save
    os.makedirs("analysis_results_complete/validation", exist_ok=True)
    plt.tight_layout()
    plt.savefig("analysis_results_complete/validation/split_sample_validation.png")
    
    stats_df = pd.DataFrame(validation_results)
    stats_df.to_csv("analysis_results_complete/validation/stats.csv", index=False)
    print("\nFinal Validation Statistics:")
    print(stats_df)

if __name__ == "__main__":
    run_validation()
