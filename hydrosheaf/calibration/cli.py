"""
Calibration CLI Module.
"""

import argparse
import os
import json
from pathlib import Path
from typing import Any, Optional, Dict

from .config import load_calibration_config, load_observations_from_csv
from .glm import PESTGLM
from .pestpp.runner import run_pestpp
from .adapters import (
    TransportCalibrationAdapter,
    TransportExperiment,
    VadoseCalibrationAdapter,
    CompositeCalibrationAdapter,
    KineticCalibrationAdapter,
    KineticExperiment
)
from .adapters_iso import WaterIsotopeMixingAdapter, WaterEndmember
from ..log import setup_logging, get_logger


def setup_transport_adapter(config, settings=None):
    """Helper to setup transport adapter from config settings."""
    s = settings if settings else config.adapter_settings

    # Check if we have observations file
    obs_file = s.get("observations_file", config.observations_file)
    if not obs_file:
        # Maybe embedded observations?
        pass

    # ... (existing transport logic adapted) ...
    # Reuse the logic from run_calibration_cli but make it reusable

    dist = float(s.get("distance_m", 10.0))

    # Observations
    times = []
    vals = []
    weights = []

    if obs_file:
        import pandas as pd
        df = pd.read_csv(obs_file)
        times = df["time"].tolist()
        vals = df["value"].tolist()
        weights = df["weight"].tolist() if "weight" in df.columns else [1.0] * len(df)

    exp = TransportExperiment(
        id=str(s.get("id", "exp1")),
        times=times,
        observed_concentrations=vals,
        distance_m=dist,
        weights=weights,
        velocity_m_day=float(s.get("velocity", 0.1)),
        source_concentration=float(s.get("source_conc", 1.0)),
    )

    # Params to fit - if in composite, they are defined in global params or sub-model specific?
    # Usually in composite, we define parameters globally in config.parameters,
    # and the adapter maps them.
    # TransportCalibrationAdapter expects 'params_to_fit' list of names.

    # We can infer from config.parameters names?
    # e.g. if 'dispersivity' is in config.parameters, we fit it.

    params_to_fit = []
    known_transport_params = ["dispersivity", "decay", "velocity"]

    # If settings has explicit list
    if "fit_parameters" in s:
        params_to_fit = s["fit_parameters"]
    else:
        # Check global params
        for p in config.parameters:
            if p.name in known_transport_params:
                params_to_fit.append(p.name)

    return TransportCalibrationAdapter(
        experiments=[exp],
        params_to_fit=params_to_fit,
        base_dispersivity=float(s.get("dispersivity", 1.0)),
        base_decay=float(s.get("decay", 0.0)),
    )


def setup_vadose_adapter(config, settings=None):
    """Helper to setup vadose adapter."""
    raise NotImplementedError(
        "Vadose calibration is not implemented in hydrosheaf.calibration.cli yet. "
        "Use the main CLI vadose workflow in hydrosheaf/cli.py (--vadose-calibrate), "
        "or implement setup_vadose_adapter to build VadoseCalibrationAdapter from vadose/io loaders."
    )


def setup_kinetic_adapter(config, settings=None):
    """Helper to setup kinetic PHREEQC adapter."""
    s = settings if settings else config.adapter_settings

    from ..config import Config as HydrosheafConfig
    from ..reactive_transport import KineticParameters

    # 1. Base Parameters (Minerals to simulate)
    # usually defined in config.parameters? Or specific kinetic section?
    # If using advanced config, we might have a 'minerals' section in settings.

    minerals = s.get("minerals", ["calcite"])  # Default list
    base_params = {}

    # Construct base parameters with reasonable defaults
    for m in minerals:
        # Check if we have specific overrides in config
        # This is a simplification. Real implementation needs robust dict parsing.
        # We use default kinetic params if possible, or fallback
        # Ideally we'd use get_default_kinetic_params from rate_laws but we need to import it
        from ..reactive_transport.rate_laws import get_default_kinetic_params, DEFAULT_KINETIC_PARAMS

        if m in DEFAULT_KINETIC_PARAMS:
            kp = get_default_kinetic_params(m)
        else:
            kp = KineticParameters(reaction_name=m, rate_constant=1e-9, surface_area=1.0)

        base_params[m] = kp

    # 2. Experiments
    # We need to load observations and map them to experiments.
    # Logic similar to transport: Load CSV, create KineticExperiment objects.

    obs_file = s.get("observations_file", config.observations_file)
    experiments = []

    default_cfg = HydrosheafConfig()
    hydro_config = HydrosheafConfig(
        ion_order=s.get("ion_order", default_cfg.ion_order),
        phreeqc_mode=s.get("phreeqc_mode", default_cfg.phreeqc_mode),
        phreeqc_database=s.get("phreeqc_database", default_cfg.phreeqc_database),
        phreeqc_executable=s.get("phreeqc_executable", default_cfg.phreeqc_executable),
        temp_default_c=float(s.get("temp_default_c", default_cfg.temp_default_c)),
    )

    if obs_file:
        import pandas as pd
        df = pd.read_csv(obs_file)

        # Group by experiment ID?
        # If CSV has 'experiment_id' column
        if "experiment_id" in df.columns:
            groups = df.groupby("experiment_id")
        else:
            # Treat all as one experiment? Or row-based?
            # Kinetic batch: each row could be a time point or a final state of a different reactor.
            # Let's assume grouping by 'experiment_id' if present, else one exp.
            groups = [("exp1", df)]

        for exp_id, group in groups:
            # Gather observations for this experiment
            obs_dict = {}
            for _, row in group.iterrows():
                # ion name in 'id' or 'ion' column?
                ion = row.get("ion", row.get("id"))  # fallback
                val = float(row.get("value"))
                obs_dict[ion] = val

            # Experiment conditions (residence time, etc)
            # Should be in settings or in a separate experiment config file.
            # Simplified: Use global settings for all experiments
            res_time = float(s.get("residence_time_days", 10.0))

            exp = KineticExperiment(
                id=str(exp_id),
                initial_solution=s.get("initial_solution", {"pH": 7.0}),
                residence_time_days=res_time,
                reaction_extents=[0.0] * len(minerals),  # Placeholder
                reaction_labels=minerals,
                observations=obs_dict,
                units=str(s.get("units", "mg/L")),
                edge_id=s.get("edge_id"),
                default_extent_mmol_L=float(s.get("default_extent_mmol_L", 10.0)),
            )
            experiments.append(exp)

    # 3. Params to fit
    params_to_fit = []
    # If settings has explicit list
    if "fit_parameters" in s:
        params_to_fit = s["fit_parameters"]
    else:
        # Default: fit k for all minerals
        for m in minerals:
            params_to_fit.append(f"{m}:k")

    return KineticCalibrationAdapter(
        base_params=base_params,
        experiments=experiments,
        config=hydro_config,
        params_to_fit=params_to_fit
    )


def setup_hydroiso_adapter(config, settings=None):
    """Helper to setup Water Isotope Mixing adapter."""
    s = settings if settings else config.adapter_settings

    # 1. Endmembers
    # List of {name: str, d18O: float, d2H: float, ...}
    endmembers = []
    for em_def in s.get("endmembers", []):
        endmembers.append(
            WaterEndmember(
                name=str(em_def["name"]),
                d18O=float(em_def["d18O"]),
                d2H=float(em_def["d2H"]),
                d18O_sigma=float(em_def.get("d18O_sigma", 0.5)),
                d2H_sigma=float(em_def.get("d2H_sigma", 2.0))
            )
        )

    # 2. Observations and Groups
    # Load from obs file
    obs_file = s.get("observations_file", config.observations_file)
    observations = {}
    group_map = {}

    if obs_file:
        import pandas as pd
        df = pd.read_csv(obs_file)
        # Expected cols: sample_id, d18O, d2H, group_id

        for _, row in df.iterrows():
            sid = str(row["sample_id"])
            group = str(row.get("group_id", "default_group"))

            obs_dict = {}
            if "d18O" in row and not pd.isna(row["d18O"]):
                obs_dict["18O"] = float(row["d18O"])
            if "d2H" in row and not pd.isna(row["d2H"]):
                obs_dict["2H"] = float(row["d2H"])

            observations[sid] = obs_dict
            group_map[sid] = group

    # 3. Weights
    # Defaults or from config
    weights = s.get("weights", {"18O": 5.0, "2H": 1.0})

    return WaterIsotopeMixingAdapter(
        endmembers=endmembers,
        observations=observations,
        group_map=group_map,
        weights=weights
    )


def run_calibration_cli(args):
    """
    Main entry point for calibration commands.
    """
    # 1. Setup Logging
    log_file = args.log if args.log else "calibration.log"
    setup_logging(verbose=args.verbose, log_file=log_file)
    logger = get_logger("calibration.cli")

    logger.info("Hydrosheaf PEST Calibration Started")
    logger.info(f"Loading config: {args.config}")

    # 2. Load Config
    try:
        config = load_calibration_config(args.config)
    except Exception as e:
        logger.error(f"Error loading config: {e}")
        return

    # 3. Setup Adapter
    problem = None

    if config.problem_type == "composite":
        logger.info("Setting up Composite Calibration...")
        sub_problems = []

        for model_def in config.sub_models:
            m_type = model_def.get("type")
            logger.info(f"Adding sub-model: {m_type}")

            if m_type == "transport":
                prob = setup_transport_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)
            elif m_type == "vadose":
                prob = setup_vadose_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)
            elif m_type == "kinetic":
                prob = setup_kinetic_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)
            elif m_type == "hydroiso":
                prob = setup_hydroiso_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)

        if sub_problems:
            problem = CompositeCalibrationAdapter(sub_problems)
        else:
            logger.error("No valid sub-models found for composite problem.")
            return

    elif config.problem_type == "transport":
        logger.info("Setting up Transport Calibration...")
        problem = setup_transport_adapter(config)

    elif config.problem_type == "vadose":
        logger.info("Setting up Vadose Calibration...")
        try:
            problem = setup_vadose_adapter(config)
        except NotImplementedError as e:
            logger.error(str(e))
            return

    elif config.problem_type == "kinetic":
        logger.info("Setting up Kinetic (PHREEQC) Calibration...")
        problem = setup_kinetic_adapter(config)

    elif config.problem_type == "hydroiso":
        logger.info("Setting up Water Isotope Mixing Calibration...")
        problem = setup_hydroiso_adapter(config)

    else:
        logger.error(f"Unknown problem type: {config.problem_type}")
        return

    if problem is None:
        logger.error("Failed to initialize problem adapter.")
        return

    # 4. Create PEST or Run PEST++
    if config.engine.startswith("pestpp"):
        # Use external PEST++ binary
        # Updated to pass IES settings
        result = run_pestpp(
            problem=problem,
            engine=config.engine,
            work_dir=config.work_dir or "pest_workspace",
            case_name="calibration",
            max_nfev=config.max_nfev,
            n_workers=config.n_workers,
            pestpp_options=config.ies_settings
        )
    else:
        # Internal PESTGLM (Python)
        # Override parameters with config definitions (bounds, log, etc)
        pest_params = config.parameters

        pest = PESTGLM(
            parameters=pest_params,
            observations=problem.get_observations(),
            model_runner=problem.run_model,
            n_workers=config.n_workers,
            worker_type="thread",  # Safe default
        )

        # 5. Run
        result = pest.calibrate(max_nfev=config.max_nfev)

    # 6. Save Results
    os.makedirs(config.output_dir, exist_ok=True)
    out_path = Path(config.output_dir) / "results.json"

    import numpy as np  # Import locally if not at top

    # Convert numpy types
    def convert(o):
        if isinstance(o, np.int64):
            return int(o)
        if isinstance(o, np.float64):
            return float(o)
        raise TypeError

    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=convert)

    logger.info(f"Calibration Complete. Results saved to {out_path}")
    logger.info(f"Success: {result['success']}")
    logger.info(f"Phi: {result['phi']:.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Path to calibration.yaml")
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose console logging"
    )
    parser.add_argument("--log", type=str, help="Path to write debug log file")
    args = parser.parse_args()
    run_calibration_cli(args)
