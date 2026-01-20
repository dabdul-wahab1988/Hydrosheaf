"""
Calibration CLI Module.
"""

import argparse
import os
import json
from pathlib import Path

from .config import load_calibration_config, load_observations_from_csv
from .glm import PESTGLM
from .adapters import TransportCalibrationAdapter, TransportExperiment
from ..log import setup_logging, get_logger


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

    if config.problem_type == "transport":
        logger.info("Setting up Transport Calibration...")
        load_observations_from_csv(config.observations_file)

        # In the simple transport adapter, we need experiments.
        # This mapping from generic obs to TransportExperiment needs logic.
        # For simplicity, we assume 1 experiment "exp1" and match obs ids.

        # Group obs by experiment prefix if present?
        # Or just assume all obs belong to one experiment for now.

        # Get experiment settings from config
        dist = config.adapter_settings.get("distance_m", 10.0)
        times = []
        vals = []
        weights = []

        # Parse times from observation names? Or need time column in CSV.
        # Let's assume the OBS CSV has a 'time' column which our loader ignored,
        # or we need to update loader.
        # Let's read CSV directly here to get times.
        import pandas as pd

        df = pd.read_csv(config.observations_file)
        times = df["time"].tolist()
        vals = df["value"].tolist()
        weights = df.get("weight", [1.0] * len(df)).tolist()

        exp = TransportExperiment(
            id="exp1",
            times=times,
            observed_concentrations=vals,
            distance_m=dist,
            weights=weights,
            velocity_m_day=config.adapter_settings.get("velocity", 0.1),
            source_concentration=config.adapter_settings.get("source_conc", 1.0),
        )

        # Params to fit
        p_names = [p.name for p in config.parameters]
        # Map generic parameter names to adapter keys (dispersivity, decay, velocity)
        # We assume user named them correctly in YAML.

        problem = TransportCalibrationAdapter(
            experiments=[exp],
            params_to_fit=p_names,  # e.g. ["dispersivity"]
            base_dispersivity=config.adapter_settings.get("dispersivity", 1.0),
            base_decay=config.adapter_settings.get("decay", 0.0),
        )

    elif config.problem_type == "vadose":
        logger.info("Setting up Vadose Calibration...")
        pass

    else:
        logger.error(f"Unknown problem type: {config.problem_type}")
        return

    if problem is None:
        logger.error("Failed to initialize problem adapter.")
        return

    # 4. Create PEST

    # Override parameters with config definitions (bounds, log, etc)
    # The adapter creates default parameters, but we want the ones from YAML.
    # We can inject them or wrap the adapter?
    # PESTGLM.from_problem calls problem.get_parameters().
    # We can patch the problem object to return OUR parameters.

    # Better: Adapter should initialize parameters based on config, but generic config loader
    # created AdjustableParameter objects.
    # Let's just pass these parameters to PESTGLM directly if possible.
    # PESTGLM constructor takes parameters list.

    # We need to ensure parameter names match what the model runner expects.
    pest_params = config.parameters

    pest = PESTGLM(
        parameters=pest_params,
        observations=problem.get_observations(),  # Use adapter's obs generation or config?
        # Adapter's obs might differ from config's generic list if adapter does processing.
        # But for TransportAdapter, we passed the data.
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
