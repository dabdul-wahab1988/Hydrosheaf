"""
Calibration CLI Module.
"""

import argparse
import os
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, cast

from .config import load_calibration_config
from .glm import PESTGLM
from .pestpp.runner import run_pestpp
from .definitions import AdjustableParameter
from .factory import (
    build_calibration_problem,
    setup_transport_adapter,
    setup_vadose_adapter,
    setup_kinetic_adapter,
    setup_hydroiso_adapter,
    setup_nitrate_adapter,
    setup_age_temporal_adapter,
    setup_topology_adapter,
)
from ..log import setup_logging, get_logger


def _resolve_internal_parameters(
    problem: Any,
    config_parameters: List[AdjustableParameter],
    logger: Any,
) -> List[AdjustableParameter]:
    """
    Resolve parameters for internal GLM runs.

    Default source is adapter/problem parameters. Config parameters with matching
    names override bounds/transforms/priors/initial values. Unknown config
    parameter names are ignored with a warning.
    """
    problem_parameters = [replace(param) for param in problem.get_parameters()]
    if not config_parameters:
        return problem_parameters

    if not problem_parameters:
        logger.warning(
            "Adapter returned no parameters; using config.parameters directly for internal calibration."
        )
        return [replace(param) for param in config_parameters]

    overrides_by_name: Dict[str, AdjustableParameter] = {}
    duplicate_names = set()
    for param in config_parameters:
        if param.name in overrides_by_name:
            duplicate_names.add(param.name)
        overrides_by_name[param.name] = param

    if duplicate_names:
        logger.warning(
            "Duplicate parameter definitions found in config.parameters; using the last value for: %s",
            sorted(duplicate_names),
        )

    problem_names = {param.name for param in problem_parameters}
    unknown_names = sorted(name for name in overrides_by_name if name not in problem_names)
    if unknown_names:
        logger.warning(
            "Ignoring config.parameters not present in adapter parameter set: %s",
            unknown_names,
        )

    resolved_parameters: List[AdjustableParameter] = []
    for base in problem_parameters:
        override = overrides_by_name.get(base.name)
        if override is None:
            resolved_parameters.append(base)
            continue
        resolved_parameters.append(
            replace(
                base,
                value=override.value,
                lower_bound=override.lower_bound,
                upper_bound=override.upper_bound,
                log_transform=override.log_transform,
                group=override.group,
                prior_mean=override.prior_mean,
                prior_sigma=override.prior_sigma,
                description=override.description,
                fixed=override.fixed,
                tied_to=override.tied_to,
            )
        )
    return resolved_parameters


def run_calibration_cli(args):
    """
    Main entry point for calibration commands.
    """
    import math
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
        raise ValueError(f"Failed to load config: {e}")

    # 3. Setup Adapter
    problem = build_calibration_problem(config, logger)

    # Expose dry-run mode
    if getattr(args, "dry_run", False):
        logger.info("Dry-run validation complete.")
        pest_params = _resolve_internal_parameters(problem, config.parameters, logger)
        obs = problem.get_observations()
        print("\n" + "=" * 60)
        print("DRY RUN: RESOLVED CALIBRATION TARGETS")
        print("=" * 60)
        print(f"Problem Type: {config.problem_type}")
        print(f"Engine:       {config.engine}")
        print(f"Output Dir:   {config.output_dir}")
        print("-" * 60)
        print(f"Adjustable Parameters ({len(pest_params)}):")
        for p in pest_params:
            log_str = " (log-transformed)" if p.log_transform else ""
            print(f"  - {p.name:<30} Initial: {p.value:<10} Bounds: [{p.lower_bound}, {p.upper_bound}]{log_str}")
        print("-" * 60)
        print(f"Observations ({len(obs)}):")
        for o in obs:
            print(f"  - {o.name:<30} Observed: {o.value:<10} Weight: {o.weight}")
        print("=" * 60)
        return

    # 4. Create PEST or Run PEST++ or run Outer-Loop Topology Search
    if config.problem_type == "topology" and config.adapter_settings.get("outer_loop", False):
        logger.info("Executing Outer-Loop Topology Search...")
        from .adapters import TopologyOuterLoopCalibrator, TopologyCalibrationAdapter
        
        samples_file = (
            config.adapter_settings.get("samples_file")
            or config.adapter_settings.get("nodes_file")
            or config.adapter_settings.get("observations_file")
            or config.observations_file
        )
        if not samples_file:
            raise ValueError("Outer-loop topology search requires model.samples_file.")
            
        import pandas as pd
        samples_df = pd.read_csv(samples_file)
        
        topo_adapter = cast(TopologyCalibrationAdapter, problem)
        calibrator = TopologyOuterLoopCalibrator(
            samples_df=samples_df,
            candidate_edges=topo_adapter.candidate_edges,
            observations=topo_adapter.observations,
            config=config.adapter_settings.get("model_config"),
            max_iterations=int(config.adapter_settings.get("max_iterations", 5)),
            max_neighbors=int(config.adapter_settings.get("max_neighbors", 2)),
            head_key=config.adapter_settings.get("head_key", "hydraulic_head"),
            elevation_key=config.adapter_settings.get("elevation_key", "elevation"),
            layer_key=config.adapter_settings.get("layer_key", "aquifer_layer"),
        )
        result = calibrator.search()
    elif config.engine.startswith("pestpp"):
        # Use external PEST++ binary
        pestpp_opts = config.pestpp_options.copy()
        if config.engine == "pestpp-ies":
            pestpp_opts.update(config.ies_settings)
        elif config.engine == "pestpp-sen":
            pestpp_opts.update(config.sen_settings)
        elif config.engine == "pestpp-swp":
            pestpp_opts.update(config.swp_settings)
        elif config.engine in ("pestpp-mou", "pestpp-opt"):
            pestpp_opts.update(config.opt_settings)
        elif config.engine == "pestpp-da":
            pestpp_opts.update(config.da_settings)

        result = run_pestpp(
            problem=problem,
            engine=config.engine,
            work_dir=config.work_dir or "pest_workspace",
            case_name="calibration",
            max_nfev=config.max_nfev,
            n_workers=config.n_workers,
            pestpp_options=pestpp_opts,
            pestpp_version=config.pestpp_version
        )
    else:
        # Internal PESTGLM (Python)
        # Start from adapter parameters, then overlay matching config definitions.
        pest_params = _resolve_internal_parameters(problem, config.parameters, logger)
        if not pest_params:
            logger.error(
                "No calibration parameters were resolved. Configure fit_parameters or calibration.parameters."
            )
            return

        pest = PESTGLM(
            parameters=pest_params,
            observations=problem.get_observations(),
            model_runner=problem.run_model,
            n_workers=config.n_workers,
            worker_type="thread",  # Safe default
            loss=config.loss,
        )

        # 5. Run
        result = pest.calibrate(max_nfev=config.max_nfev)

    # 7. Post-calibration diagnostics & topology details
    # Expose edge probabilities & selected edges for topology models
    for sub_prob in getattr(problem, "sub_problems", [problem]):
        if hasattr(sub_prob, "edge_probabilities") and hasattr(sub_prob, "selected_edges"):
            try:
                edge_probs = sub_prob.edge_probabilities(result.get('optimal_parameters', {}))
                sel_edges = sub_prob.selected_edges(result.get('optimal_parameters', {}))
                result['edge_probabilities'] = edge_probs
                result['selected_edges'] = [
                    {"edge_id": edge.edge_id, "u": edge.u, "v": edge.v}
                    for edge in sel_edges
                ]
            except Exception as ex:
                logger.warning(f"Failed to extract topology post-calibration details: {ex}")

    # Calculate AIC/BIC
    try:
        nobs = len(problem.get_observations())
        npar = len(result.get('optimal_parameters', {}))
        phi = result.get('phi', 0.0)
        if nobs > 0 and phi > 0:
            mse = phi / nobs
            log_lik = -0.5 * nobs * (math.log(2.0 * math.pi * mse) + 1.0)
            result["aic"] = 2.0 * npar - 2.0 * log_lik
            result["bic"] = math.log(nobs) * npar - 2.0 * log_lik
    except Exception as ex:
        logger.warning(f"Failed to calculate AIC/BIC diagnostics: {ex}")

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
    logger.info(f"Success: {result.get('success', False)}")
    logger.info(f"Phi: {result.get('phi', 0.0):.4f}")

    # ── 8. Post-calibration validation workflow (topology only) ───────
    val_file = config.adapter_settings.get(
        "validation_observations_file",
        config.validation_observations_file,
    )
    if config.problem_type == "topology" and val_file:
        logger.info(
            "Validation observations file detected — running independent "
            "validation workflow."
        )
        from .validation_workflow import run_assumption_calibration_validation

        # Pass the calibration result so validation reuses the same
        # optimal parameters instead of re-running calibration.
        val_report = run_assumption_calibration_validation(
            config, calibration_result=result,
        )
        result["validation_metrics"] = val_report["validation_metrics"]
        result["independent_validation"] = val_report["independent_validation"]
        logger.info(
            "Validation metrics appended to results.json — "
            "precision=%.4f recall=%.4f f1=%.4f",
            val_report["validation_metrics"]["precision"],
            val_report["validation_metrics"]["recall"],
            val_report["validation_metrics"]["f1"],
        )
        # Re-write results.json with validation fields included
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2, default=convert)

    # ── 9. Post-calibration benchmark workflow (topology only) ──────────
    run_benchmark = config.adapter_settings.get(
        "run_assumption_benchmark", False
    )
    if config.problem_type == "topology" and run_benchmark:
        logger.info(
            "run_assumption_benchmark=true — running benchmark workflow."
        )
        from .benchmark import run_assumption_benchmark

        benchmark_report = run_assumption_benchmark(
            config, calibration_result=result,
        )
        result["benchmark"] = {
            "improvement_summary": benchmark_report["improvement_summary"],
            "variants": benchmark_report["variants"],
            "independent_validation": benchmark_report["independent_validation"],
            "manuscript_claim_allowed": benchmark_report["manuscript_claim_allowed"],
            "uncertainty": benchmark_report.get("uncertainty"),
        }
        logger.info(
            "Benchmark complete — delta_f1=%.4f",
            benchmark_report["improvement_summary"]["delta_f1"],
        )
        # Re-write results.json with benchmark fields included
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2, default=convert)

        # ── 10. Active-learning measurement recommendations (topology only) ──
        run_active_learning = config.adapter_settings.get(
            "active_learning", False
        )
        if run_active_learning:
            logger.info(
                "active_learning=true — ranking next measurement priorities."
            )
            from .active_learning import rank_next_measurements

            # Read full benchmark + validation reports from disk
            bench_json_path = Path(config.output_dir) / "assumption_benchmark_results.json"
            val_json_path = Path(config.output_dir) / "assumption_validation_results.json"

            full_benchmark_report = {}
            if bench_json_path.exists():
                with open(bench_json_path) as fh:
                    full_benchmark_report = json.load(fh)

            full_validation_report = None
            if val_json_path.exists():
                with open(val_json_path) as fh:
                    full_validation_report = json.load(fh)

            # Extract samples / candidate_edges / config from the adapter
            topo_samples = None
            topo_edges = None
            topo_config = None
            for sub_prob in getattr(problem, "sub_problems", [problem]):
                if hasattr(sub_prob, "samples"):
                    topo_samples = sub_prob.samples
                if hasattr(sub_prob, "candidate_edges"):
                    topo_edges = sub_prob.candidate_edges
                if hasattr(sub_prob, "config"):
                    topo_config = sub_prob.config

            # Normalize samples to list-of-dicts if DataFrame
            if topo_samples is not None:
                try:
                    import pandas as pd
                    if isinstance(topo_samples, pd.DataFrame):
                        topo_samples = topo_samples.to_dict("records")
                except Exception:
                    pass

            top_k = int(config.adapter_settings.get("active_learning_top_k", 20))

            al_result = rank_next_measurements(
                benchmark_report=full_benchmark_report,
                validation_report=full_validation_report,
                samples=topo_samples,
                candidate_edges=topo_edges,
                config=topo_config,
                top_k=top_k,
                output_dir=config.output_dir,
            )

            result["active_learning"] = {
                "summary": al_result["summary"],
                "n_recommendations": al_result["summary"]["n_recommendations"],
                "top_priority_score": al_result["summary"]["top_priority_score"],
            }
            logger.info(
                "Active learning complete — %d recommendations, top score=%.4f",
                al_result["summary"]["n_recommendations"],
                al_result["summary"]["top_priority_score"],
            )
            # Re-write results.json with active_learning fields included
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2, default=convert)

    # ── Guard: active_learning requires benchmark ──────────────────────
    run_active_learning = config.adapter_settings.get("active_learning", False)
    if run_active_learning and "benchmark" not in result:
        logger.error(
            "active_learning=true requires run_assumption_benchmark=true. "
            "Set both in adapter_settings to enable measurement recommendations."
        )
        raise ValueError(
            "active_learning requires run_assumption_benchmark. "
            "No benchmark report available for ranking."
        )

    # Print human-readable summary
    print("\n" + "=" * 60)
    print("HYDROSHEAF CALIBRATION SUMMARY")
    print("=" * 60)
    print(f"Success:      {result.get('success', False)}")
    print(f"Final Phi:    {result.get('phi', 0.0):.6f}")
    if "aic" in result:
        print(f"AIC:          {result['aic']:.4f}")
    if "bic" in result:
        print(f"BIC:          {result['bic']:.4f}")
    print("-" * 60)
    print(f"{'Parameter':<30} | {'Initial':<12} | {'Optimal':<12}")
    print("-" * 60)
    
    init_vals = {p.name: p.value for p in problem.get_parameters()}
    opt_vals = result.get('optimal_parameters', {})
    for p_name, opt_val in opt_vals.items():
        init_val = init_vals.get(p_name, "N/A")
        if isinstance(init_val, (int, float)):
            init_str = f"{init_val:.6f}"
        else:
            init_str = str(init_val)
        print(f"{p_name:<30} | {init_str:<12} | {opt_val:.6f}")
    print("=" * 60)
    
    if "selected_edges" in result:
        print("\nCalibrated Selected Edges:")
        for edge in result["selected_edges"]:
            p = result["edge_probabilities"].get(edge["edge_id"], 1.0)
            print(f"  - {edge['edge_id']} ({edge['u']} -> {edge['v']}) [p={p:.4f}]")
        print("=" * 60)

    if "active_learning" in result:
        al = result["active_learning"]
        print("\n" + "=" * 60)
        print("ACTIVE LEARNING RECOMMENDATIONS")
        print("=" * 60)
        print(f"Recommendations: {al['n_recommendations']}")
        print(f"Top Priority:    {al['top_priority_score']:.4f}")
        print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Hydrosheaf Calibration CLI")
    parser.add_argument("config", nargs="?", help="Path to calibration.yaml")
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose console logging"
    )
    parser.add_argument("--log", type=str, help="Path to write debug log file")
    parser.add_argument("--dry-run", action="store_true", help="Validate config and show parameters/observations")
    parser.add_argument("--write-template", type=str, help="Write starter YAML config of given type")
    args = parser.parse_args()

    if args.write_template:
        from .templates import TEMPLATES
        t_type = args.write_template.lower()
        if t_type not in TEMPLATES:
            print(f"Error: Unknown template type '{args.write_template}'.")
            print("Available types: " + ", ".join(TEMPLATES.keys()))
            return 1
        print(TEMPLATES[t_type])
        return 0

    if not args.config:
        parser.print_help()
        return 1

    try:
        run_calibration_cli(args)
    except Exception as e:
        print(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
