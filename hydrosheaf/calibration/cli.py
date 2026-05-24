"""
Calibration CLI Module.
"""

import argparse
import os
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List

from .config import load_calibration_config
from .glm import PESTGLM
from .pestpp.runner import run_pestpp
from .definitions import AdjustableParameter
from .adapters import (
    TransportCalibrationAdapter,
    TransportExperiment,
    CompositeCalibrationAdapter,
    KineticCalibrationAdapter,
    KineticExperiment,
    VadoseCalibrationAdapter,
    NitrateSourceCalibrationAdapter,
    NitrateSourceCalibrationObservation,
    AgeTemporalCalibrationAdapter,
    AgeTemporalExperiment,
    TopologyCalibrationAdapter,
    TopologyCalibrationObservation,
)
from .adapters_iso import WaterIsotopeMixingAdapter, WaterEndmember
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
            )
        )
    return resolved_parameters


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
    s = settings if settings else config.adapter_settings

    profile_file = (
        s.get("profile_file")
        or s.get("profile")
        or s.get("config_file")
        or config.model_config_file
    )
    forcing_file = s.get("forcing_file") or s.get("forcing")
    obs_file = s.get("observations_file", config.observations_file)

    if not profile_file:
        raise ValueError("Vadose calibration requires model.profile_file or model.config_file.")
    if not forcing_file:
        raise ValueError("Vadose calibration requires model.forcing_file.")
    if not obs_file:
        raise ValueError("Vadose calibration requires observations.file.")

    from ..vadose.io import load_forcing_csv, load_profile
    from ..vadose.contracts import VadoseRunConfig
    import pandas as pd

    profile = load_profile(str(profile_file))
    forcing = load_forcing_csv(
        str(forcing_file),
        time_column=str(s.get("time_column", "timestamp")),
        time_format=str(s.get("time_format", "%Y-%m-%d")),
    )

    run_config = VadoseRunConfig(
        dz_m=float(s.get("dz_m", s.get("dz", 0.1))),
        kc=float(s.get("kc", 1.0)),
        evaporation_fraction=float(s.get("evaporation_fraction", s.get("evap_frac", 0.3))),
        lower_boundary=s.get("lower_boundary", "free_drainage"),
    )

    obs_df = pd.read_csv(obs_file)
    forcing_times = [sample.timestamp for sample in forcing]

    def nearest_time_index(raw_value):
        if raw_value is None or raw_value == "":
            return 0
        try:
            return int(raw_value)
        except (TypeError, ValueError):
            pass
        ts = pd.to_datetime(raw_value).to_pydatetime()
        if not forcing_times:
            return 0
        diffs = [abs((t - ts).total_seconds()) for t in forcing_times]
        return int(min(range(len(diffs)), key=lambda idx: diffs[idx]))

    observations = []
    for _, row in obs_df.iterrows():
        obs_type = str(row.get("type", s.get("observation_type", "theta")))
        if "time_idx" in row and not pd.isna(row["time_idx"]):
            time_idx = nearest_time_index(row["time_idx"])
        else:
            time_idx = nearest_time_index(row.get("timestamp", row.get("time")))
        value = row.get("value", row.get(obs_type))
        if pd.isna(value):
            continue
        observations.append(
            {
                "time_idx": int(time_idx),
                "depth_m": float(row.get("depth_m", s.get("depth_m", 0.0))),
                "type": obs_type,
                "value": float(value),
            }
        )

    layers_to_fit = s.get("layers_to_fit", [0])
    if isinstance(layers_to_fit, int):
        layers_to_fit = [layers_to_fit]

    return VadoseCalibrationAdapter(
        profile=profile,
        forcing=forcing,
        observations=observations,
        config=run_config,
        layers_to_fit=[int(layer) for layer in layers_to_fit],
        params_to_fit=s.get("fit_parameters"),
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
            first = group.iloc[0] if len(group) else {}

            def _setting_or_column(key, default=None):
                if key in s:
                    return s.get(key)
                try:
                    value = first.get(key, default)
                except AttributeError:
                    return default
                if value is None:
                    return default
                try:
                    import pandas as pd
                    if pd.isna(value):
                        return default
                except Exception:
                    pass
                return value

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
            res_time = float(_setting_or_column("residence_time_days", 10.0))
            default_extent = _setting_or_column("default_extent_mmol_L", 10.0)

            exp = KineticExperiment(
                id=str(exp_id),
                initial_solution=s.get("initial_solution", {"pH": 7.0}),
                residence_time_days=res_time,
                reaction_extents=[0.0] * len(minerals),  # Placeholder
                reaction_labels=minerals,
                observations=obs_dict,
                units=str(s.get("units", "mg/L")),
                edge_id=_setting_or_column("edge_id"),
                default_extent_mmol_L=(
                    float(default_extent) if default_extent is not None else None
                ),
                geological_layer=_setting_or_column("geological_layer"),
                aquifer_unit=_setting_or_column("aquifer_unit"),
                site_group=_setting_or_column("site_group"),
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


def setup_nitrate_adapter(config, settings=None):
    """Helper to setup nitrate-source calibration adapter."""
    s = settings if settings else config.adapter_settings
    import pandas as pd

    samples_file = (
        s.get("samples_file")
        or s.get("nodes_file")
        or s.get("observations_file")
        or config.observations_file
    )
    obs_file = s.get("target_file") or s.get("targets_file") or config.observations_file
    if not samples_file:
        raise ValueError("Nitrate calibration requires model.samples_file or observations.file.")

    nodes_df = pd.read_csv(samples_file)
    obs_df = pd.read_csv(obs_file) if obs_file else nodes_df
    node_col = s.get("node_id_column")
    if node_col is None:
        for candidate in ("site_id", "node_id", "sample_id", "id"):
            if candidate in obs_df.columns:
                node_col = candidate
                break
    if node_col is None:
        raise ValueError("Nitrate calibration targets need a site_id/node_id/sample_id column.")

    target_col = s.get("target_column")
    if target_col is None:
        for candidate in ("observed_p_manure", "target_p_manure", "p_manure", "value"):
            if candidate in obs_df.columns:
                target_col = candidate
                break
    if target_col is None:
        raise ValueError("Nitrate calibration targets need observed_p_manure or value column.")

    observations = []
    for _, row in obs_df.iterrows():
        if pd.isna(row.get(node_col)) or pd.isna(row.get(target_col)):
            continue
        observations.append(
            NitrateSourceCalibrationObservation(
                node_id=str(row[node_col]),
                value=float(row[target_col]),
                target=str(row.get("target", s.get("target", "p_manure"))),
                weight=float(row.get("weight", 1.0)),
            )
        )
    if not observations:
        raise ValueError("No nitrate calibration observations were loaded.")

    return NitrateSourceCalibrationAdapter(
        nodes_df=nodes_df,
        observations=observations,
        edge_results=[],
        config=None,
        params_to_fit=s.get("fit_parameters"),
        base_overrides=s.get("overrides", {}),
    )


def setup_age_temporal_adapter(config, settings=None):
    """Helper to setup groundwater age / temporal calibration adapter."""
    s = settings if settings else config.adapter_settings
    obs_file = s.get("observations_file", config.observations_file)
    if not obs_file:
        raise ValueError("Age/temporal calibration requires observations.file.")

    import pandas as pd

    df = pd.read_csv(obs_file)
    experiments = []
    
    tracer = str(s.get("tracer", "Cl"))

    for idx, row in df.iterrows():
        exp_id = str(row.get("id", row.get("sample_id", row.get("node_id", f"age_{idx}"))))
        
        u_times = None
        u_vals = None
        v_times = None
        v_vals = None
        
        if "u_file" in row and not pd.isna(row["u_file"]):
            u_df = pd.read_csv(row["u_file"])
            t_col = "time" if "time" in u_df.columns else ("timestamp" if "timestamp" in u_df.columns else u_df.columns[0])
            v_col = tracer if tracer in u_df.columns else ("value" if "value" in u_df.columns else u_df.columns[1])
            try:
                ts = pd.to_datetime(u_df[t_col])
                t0 = ts.min()
                u_times = [(t - t0).total_seconds() / 86400.0 for t in ts]
            except Exception:
                u_times = u_df[t_col].astype(float).tolist()
            u_vals = u_df[v_col].astype(float).tolist()

        if "v_file" in row and not pd.isna(row["v_file"]):
            v_df = pd.read_csv(row["v_file"])
            t_col = "time" if "time" in v_df.columns else ("timestamp" if "timestamp" in v_df.columns else v_df.columns[0])
            v_col = tracer if tracer in v_df.columns else ("value" if "value" in v_df.columns else v_df.columns[1])
            try:
                ts = pd.to_datetime(v_df[t_col])
                t0 = ts.min()
                v_times = [(t - t0).total_seconds() / 86400.0 for t in ts]
            except Exception:
                v_times = v_df[t_col].astype(float).tolist()
            v_vals = v_df[v_col].astype(float).tolist()

        age_value = row.get("observed_age_days", row.get("age_days", row.get("value")))
        if pd.isna(age_value) and v_vals is None:
            continue
            
        experiments.append(
            AgeTemporalExperiment(
                id=exp_id,
                observed_age_days=float(age_value) if not pd.isna(age_value) else None,
                weight=float(row.get("weight", 1.0)),
                default_tau_days=float(row.get("default_tau_days", s.get("default_tau_days", age_value if not pd.isna(age_value) else 30.0))),
                distance_m=(
                    float(row["distance_m"]) if "distance_m" in row and not pd.isna(row["distance_m"]) else None
                ),
                tracer_initial=(
                    float(row["tracer_initial"]) if "tracer_initial" in row and not pd.isna(row["tracer_initial"]) else None
                ),
                tracer_observed=(
                    float(row["tracer_observed"]) if "tracer_observed" in row and not pd.isna(row["tracer_observed"]) else None
                ),
                node_u_times=u_times,
                node_u_values=u_vals,
                node_v_times=v_times,
                node_v_values=v_vals,
                tracer_ion=tracer,
            )
        )
    if not experiments:
        raise ValueError("No age/temporal calibration observations were loaded.")

    return AgeTemporalCalibrationAdapter(
        experiments=experiments,
        params_to_fit=s.get("fit_parameters"),
        base_decay_rate_1_day=float(s.get("decay_rate_1_day", s.get("decay", 1e-3))),
        base_velocity_m_day=float(s.get("velocity_m_day", s.get("velocity", 0.1))),
        base_porosity=float(s.get("porosity", 0.2)),
        base_ttd_cv=float(s.get("ttd_cv", 0.7)),
    )


def setup_topology_adapter(config, settings=None):
    """Helper to setup flow-topology calibration adapter."""
    s = settings if settings else config.adapter_settings
    import pandas as pd
    from ..graph.types import Edge

    candidates_file = s.get("candidates_file") or s.get("edges_file")
    obs_file = s.get("observations_file", config.observations_file)
    if not obs_file:
        raise ValueError("Topology calibration requires observations.file.")

    obs_df = pd.read_csv(obs_file)
    if candidates_file:
        edge_df = pd.read_csv(candidates_file)
    else:
        edge_df = obs_df

    candidate_edges = []
    for idx, row in edge_df.iterrows():
        edge_id = str(row.get("edge_id", f"{row.get('u', 'u')}_{row.get('v', 'v')}_{idx}"))
        attrs = {}
        for key in ("edge_confidence", "p_uv"):
            if key in row and not pd.isna(row[key]):
                attrs[key] = float(row[key])
        candidate_edges.append(
            Edge(
                edge_id=edge_id,
                u=str(row.get("u", "")),
                v=str(row.get("v", "")),
                attrs=attrs,
            )
        )

    observations = []
    for _, row in obs_df.iterrows():
        edge_id = str(row.get("edge_id"))
        if not edge_id or edge_id == "nan":
            continue
        value = row.get("observed_present", row.get("present", row.get("value")))
        if pd.isna(value):
            continue
        observations.append(
            TopologyCalibrationObservation(
                edge_id=edge_id,
                observed_present=float(value),
                weight=float(row.get("weight", 1.0)),
            )
        )
    if not candidate_edges or not observations:
        raise ValueError("Topology calibration requires candidate edges and observations.")

    return TopologyCalibrationAdapter(
        candidate_edges=candidate_edges,
        observations=observations,
        prior_sigma=float(s.get("prior_sigma", 2.0)),
        normalize_by_upstream=bool(s.get("normalize_by_upstream", False)),
    )


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
            elif m_type == "nitrate":
                prob = setup_nitrate_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)
            elif m_type in {"age", "temporal", "age_temporal"}:
                prob = setup_age_temporal_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)
            elif m_type == "topology":
                prob = setup_topology_adapter(config, settings=model_def)
                if prob:
                    sub_problems.append(prob)

        if sub_problems:
            problem = CompositeCalibrationAdapter(sub_problems)
        else:
            logger.error("No valid sub-models found for composite problem.")
            raise ValueError("No valid sub-models found for composite problem.")

    elif config.problem_type == "transport":
        logger.info("Setting up Transport Calibration...")
        problem = setup_transport_adapter(config)

    elif config.problem_type == "vadose":
        logger.info("Setting up Vadose Calibration...")
        problem = setup_vadose_adapter(config)

    elif config.problem_type == "kinetic":
        logger.info("Setting up Kinetic (PHREEQC) Calibration...")
        problem = setup_kinetic_adapter(config)

    elif config.problem_type == "hydroiso":
        logger.info("Setting up Water Isotope Mixing Calibration...")
        problem = setup_hydroiso_adapter(config)

    elif config.problem_type == "nitrate":
        logger.info("Setting up Nitrate Source Calibration...")
        problem = setup_nitrate_adapter(config)

    elif config.problem_type in {"age", "temporal", "age_temporal"}:
        logger.info("Setting up Age/Temporal Calibration...")
        problem = setup_age_temporal_adapter(config)

    elif config.problem_type == "topology":
        logger.info("Setting up Flow Topology Calibration...")
        problem = setup_topology_adapter(config)

    else:
        logger.error(f"Unknown problem type: {config.problem_type}")
        raise ValueError(f"Unknown problem type: {config.problem_type}")

    if problem is None:
        logger.error("Failed to initialize problem adapter.")
        raise ValueError("Failed to initialize problem adapter.")

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
