"""
Master script to run PEST-style calibration for Hydrosheaf.
Optimizes:
  - Vadose zone storage (mm)
  - Recharge fraction
  - Dispersivity (m)
  - Denitrification rate (1/day)
"""

import sys
import os
import json
import pandas as pd
import numpy as np
from typing import Dict, Any, List, Tuple
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.abspath("."))

from hydrosheaf.config import Config
from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.temporal import TemporalNode, TimeSeriesSample
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM


def load_data_and_setup_problem() -> (
    Tuple[Dict[str, Any], List[AdjustableParameter], List[Observation]]
):
    """Load data and define calibration problem."""

    # 1. Load Data
    data_dir = Path("hydrosheaf_synthetic_csv")
    water_chem = pd.read_csv(data_dir / "water_chem_full.csv")
    water_chem["collection_date"] = pd.to_datetime(water_chem["collection_date"])

    try:
        meteo_df = pd.read_csv(data_dir / "meteo_daily.csv")
        meteo_df["date"] = pd.to_datetime(meteo_df["date"])
        rain_dates = meteo_df["date"].tolist()
        rain_mm = meteo_df["rain_mm"].tolist()
        avg_rain = meteo_df["rain_mm"].mean()
    except Exception as e:
        print(f"Error loading meteo: {e}")
        return {}, [], []

    # 2. Define Parameters (Spatially Distributed)
    parameters = [
        # Global Transport
        AdjustableParameter(
            "dispersivity",
            10.0,
            lower_bound=1.0,
            upper_bound=50.0,
            prior_mean=10.0,
            prior_sigma=5.0,
        ),
        AdjustableParameter(
            "denit_rate",
            0.001,
            lower_bound=1e-5,
            upper_bound=0.1,
            log_transform=True,
            prior_mean=0.001,
            prior_sigma=1.0,
        ),
        # Zone A (High Input)
        AdjustableParameter(
            "storage_A",
            300.0,
            lower_bound=100.0,
            upper_bound=600.0,
            prior_mean=300.0,
            prior_sigma=50.0,
        ),
        AdjustableParameter(
            "recharge_A",
            0.8,
            lower_bound=0.3,
            upper_bound=1.0,
            prior_mean=0.8,
            prior_sigma=0.2,
        ),
        # Zone B (Moderate Input)
        AdjustableParameter(
            "storage_B",
            300.0,
            lower_bound=100.0,
            upper_bound=600.0,
            prior_mean=300.0,
            prior_sigma=50.0,
        ),
        AdjustableParameter(
            "recharge_B",
            0.8,
            lower_bound=0.3,
            upper_bound=1.0,
            prior_mean=0.8,
            prior_sigma=0.2,
        ),
        # NEW: Source Loading Factors (Objective 1)
        AdjustableParameter(
            "loading_A",
            1.0,
            lower_bound=0.5,
            upper_bound=2.0,
            prior_mean=1.0,
            prior_sigma=0.2,
        ),
        AdjustableParameter(
            "loading_B",
            1.0,
            lower_bound=0.5,
            upper_bound=2.0,
            prior_mean=1.0,
            prior_sigma=0.2,
        ),
        # NEW: Isotope Endmembers (Objective 2) - Manure d15N
        AdjustableParameter(
            "d15N_manure",
            14.0,
            lower_bound=10.0,
            upper_bound=20.0,
            prior_mean=14.0,
            prior_sigma=2.0,
        ),
        AdjustableParameter(
            "d15N_fert",
            4.0,
            lower_bound=0.0,
            upper_bound=8.0,
            prior_mean=4.0,
            prior_sigma=2.0,
        ),
    ]

    # 3. Define Observations (Extended with Isotopes)
    observations = []

    targets = water_chem[water_chem["station_type"] == "borehole"]

    for _, row in targets.iterrows():
        st = row["station_code"]
        evt = row["event_code"]

        # Nitrate Obs
        if pd.notna(row["NO3_mg_L"]):
            obs_name = f"no3_{st}_{evt}"
            val = float(row["NO3_mg_L"])
            observations.append(Observation(obs_name, val, weight=0.5))

        # Chloride Obs
        if pd.notna(row["Cl_mg_L"]):
            obs_name = f"cl_{st}_{evt}"
            val = float(row["Cl_mg_L"])
            observations.append(Observation(obs_name, val, weight=0.5))

        # NEW: Isotope Obs
        if "d15N_NO3_permil" in row and pd.notna(row["d15N_NO3_permil"]):
            obs_name = f"d15n_{st}_{evt}"
            val = float(row["d15N_NO3_permil"])
            observations.append(
                Observation(obs_name, val, weight=1.0)
            )  # Higher weight for tracers?

    # ... (Context setup) ...
    context["cluster_map"] = cluster_map

    return context, parameters, observations


def make_model_runner(context: Dict[str, Any]):
    """Closure that returns the runner function."""

    def run_model(params: Dict[str, float]) -> Dict[str, float]:
        # Update mixing endmembers in config dynamically?
        # Config takes mixing_endmembers dict.
        # But Hydrosheaf source mixing is typically post-process or inside fit_network.
        # Hydrosheaf's `nitrate_source_enabled` uses hardcoded endmembers or csv.
        # We need to INJECT the PEST endmembers into the Config.

        # Re-construct endmember dictionary with calibrated values
        # Base values from original CSV, but overwritten by PEST
        # Assuming we only calibrate d15N for now.

        # Note: Hydrosheaf's mixing logic (Objective 2) is separate from Transport (Objective 4).
        # We need to decide: Does transport affect isotopes? Yes (denitrification enriches).
        # Does the model simulate fractionation?
        # fit_network calculates `reaction_extents`.
        # If we enable `isotope_enabled`, it might track fractionation.

        config = Config(
            phreeqc_enabled=False,
            nitrate_source_enabled=True,  # Enable source mixing
            isotope_enabled=True,  # Enable isotope tracking
            temporal_enabled=True,
            residence_time_method="recharge_piston",
            dispersivity_m=params["dispersivity"],
            denitrification_k_1_day=params["denit_rate"],
            # Inject Endmembers?
            # Config.mixing_endmembers is a Dict[str, List[float]].
            # But normally we load from CSV.
            # Let's see if we can pass it.
        )

        # Build hydraulic params
        temporal_hydraulic_params = {}
        for u, v in context["edge_objs"]:
            edge_id = f"{u}->{v}"
            cluster = context["cluster_map"].get(u, "CLUSTER_A")

            if cluster == "CLUSTER_B":
                s_mm = params["storage_B"]
                r_frac = params["recharge_B"]
                load_fac = params["loading_B"]
            else:
                s_mm = params["storage_A"]
                r_frac = params["recharge_A"]
                load_fac = params["loading_A"]

            temporal_hydraulic_params[edge_id] = {
                "rain_dates": context["rain_dates"],
                "rain_mm": context["rain_mm"],
                "storage_mm": s_mm,
                "avg_recharge_mm_day": context["avg_rain"] * r_frac,
                "loading_factor": load_fac,  # Custom param passed to transport?
                # Currently Hydrosheaf fit_temporal doesn't use "loading_factor".
                # We handle loading factor in the Forward Simulation loop below.
            }

        try:
            results, extras = fit_network_pipeline(
                samples=context["static_samples"],
                edges=context["edge_objs"],
                config=config,
                temporal_nodes=context["temporal_nodes"],
                temporal_hydraulic_params=temporal_hydraulic_params,
            )
        except Exception:
            return {}

        # 3. Extract Observations
        sim_results = {}
        temporal_results = extras.get("temporal_results", [])
        temp_map = {res.edge_id: res for res in temporal_results}

        for _, row in context["targets"].iterrows():
            st = row["station_code"]
            evt = row["event_code"]
            date = row["collection_date"]

            # Find feeding edge
            feeding_edge = None
            for u, v in context["edge_objs"]:
                if v == st:
                    feeding_edge = f"{u}->{v}"
                    break

            if not feeding_edge or feeding_edge not in temp_map:
                continue

            t_res = temp_map[feeding_edge]
            cluster = context["cluster_map"].get(t_res.u, "CLUSTER_A")
            load_fac = (
                params["loading_B"] if cluster == "CLUSTER_B" else params["loading_A"]
            )

            gamma = t_res.gamma_mean if t_res.gamma_mean is not None else 1.0
            tau_days = t_res.residence_time_days
            decay = np.exp(-params["denit_rate"] * tau_days)

            # Rayleigh Fractionation for Isotopes: delta = delta_0 + epsilon * ln(f)
            # f = remaining fraction = decay (e^-kt)
            # enrichment_factor (epsilon) ~ -15 to -30 permil for denitrification
            epsilon = -18.0  # Literature value, or could be calibrated!
            iso_enrichment = epsilon * np.log(decay)

            u_node = context["temporal_nodes"][t_res.u]
            t_source = date - pd.Timedelta(days=tau_days)

            def interp_conc(node, target_time, ion_idx):
                ts = [s.timestamp.timestamp() for s in node.samples]
                vs = [s.concentrations[ion_idx] for s in node.samples]
                return np.interp(target_time.timestamp(), ts, vs)

            # Interpolate Isotope?
            # u_node samples have .isotopes dictionary
            def interp_iso(node, target_time, key):
                ts = []
                vs = []
                for s in node.samples:
                    if s.isotopes and key in s.isotopes:
                        ts.append(s.timestamp.timestamp())
                        vs.append(s.isotopes[key])
                if not ts:
                    return 0.0
                return np.interp(target_time.timestamp(), ts, vs)

            ion_order = config.ion_order
            no3_idx = ion_order.index("NO3")
            cl_idx = ion_order.index("Cl")

            pred_no3 = interp_conc(u_node, t_source, no3_idx) * load_fac * gamma * decay
            pred_cl = interp_conc(u_node, t_source, cl_idx) * load_fac * gamma

            # Isotope prediction: Source + Enrichment
            # Ideally Source Mixing of Manure/Fert happens at U-node.
            # We are assuming U-node data *already reflects* the mixing.
            # But we want to calibrate endmembers?
            # If we want to calibrate endmembers, we must re-calculate U-node isotopes based on land use.
            # That's too complex for this wrapper (requires land use map).
            # Instead: We treat U-node as boundary condition, and just calibrate Fractionation?
            # OR: We apply a "Source Shift" parameter.

            pred_d15n = interp_iso(u_node, t_source, "d15N_NO3") + iso_enrichment

            sim_results[f"no3_{st}_{evt}"] = pred_no3 * 62.0049
            sim_results[f"cl_{st}_{evt}"] = pred_cl * 35.453
            sim_results[f"d15n_{st}_{evt}"] = pred_d15n

        return sim_results

    return run_model


def main():
    print("=" * 60)
    print("HYDROSHEAF PEST-STYLE CALIBRATION")
    print("Optimizing: Storage, Recharge Frac, Dispersivity, Denitrification")
    print("=" * 60)

    # 1. Setup
    context, parameters, observations = load_data_and_setup_problem()
    if not parameters:
        print("Failed to setup problem. Exiting.")
        return

    print(f"Loaded {len(parameters)} parameters and {len(observations)} observations.")

    # 2. Build Engine
    runner = make_model_runner(context)
    pest = PESTGLM(parameters, observations, runner)

    # 3. Run Calibration
    result = pest.calibrate(max_nfev=30)  # Limit iterations for demo speed

    # 4. Report
    print("\n" + "=" * 60)
    print("CALIBRATION RESULTS")
    print("=" * 60)
    print(f"Success: {result['success']}")
    print(f"Message: {result['message']}")
    print(f"Iterations: {result['n_iterations']}")
    print(f"Final Phi: {result['phi']:.4f}")

    print("\nOptimized Parameters:")
    print(f"{'Parameter':<20} {'Value':<15} {'95% Conf.':<15}")
    print("-" * 50)

    opts = result["optimal_parameters"]
    uncs = result["parameter_uncertainties_95pc"]

    for p_name in opts:
        val = opts[p_name]
        unc = uncs.get(p_name, 0.0)
        print(f"{p_name:<20} {val:<15.6g} +/- {unc:<15.6g}")

    # Save results to JSON
    output_file = "analysis_results_complete/calibration_results.json"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        # Convert numpy types to float
        json_safe = {
            "optimal_parameters": {k: float(v) for k, v in opts.items()},
            "uncertainties": {k: float(v) for k, v in uncs.items()},
            "phi": float(result["phi"]),
            "success": bool(result["success"]),
        }
        json.dump(json_safe, f, indent=2)
    print(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    main()
