"""Orchestrate the evidence-bounded graph-TTD virtual benchmark components.

This runner persists truth-blind inference artifacts before it writes their
sealed scoring truth.  A successful run is an execution/reproducibility result;
the common claim gate deliberately remains ``ABSTAIN`` until all preregistered
generator families and comparator arms are available.
"""

from __future__ import annotations

from dataclasses import fields, is_dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .independent_particle_ttd import (
    PARTICLE_TTD_GENERATOR_FAMILY,
    ParticleTTDConfig,
    ParticleTTDRecovery,
    ParticleTTDScore,
    assert_particle_observation_view_is_truth_blind,
    generate_independent_particle_ttd_case,
    particle_ttd_public_manifest,
    recover_particle_ttd_baseline,
    score_particle_ttd_recovery,
)
from .ttd_graph_dynamic import (
    DynamicTTDGraphConfig,
    evaluate_dynamic_ttd_recovery,
    generate_dynamic_ttd_case,
    recover_dynamic_ttd_baseline,
)
from .ttd_graph_evidence import (
    BenchmarkRecord,
    assess_claim_readiness,
    assert_observation_view_is_truth_blind,
    write_virtual_benchmark_artifacts,
)
from .ttd_graph_static import (
    GENERATOR_FAMILY as STATIC_GENERATOR_FAMILY,
    default_static_ttd_graph_protocol,
    generate_static_ttd_graph_case,
    run_all_pre_registered_controls,
    score_static_ttd_graph_submission,
    write_static_ttd_graph_case,
)


CONFIG_SCHEMA = "hydrosheaf-ttd-graph-virtual-benchmark-config-v1"
STATIC_METHODS = {
    "local": "local_lumped_ttd_baseline",
    "correct": "source_mixing_graph_correct_topology_control",
    "reversed": "reversed_graph_control",
    "random": "random_graph_control",
    "edge_removed": "edge_removed_graph_control",
}
DYNAMIC_GENERATOR_FAMILY = "independent_dynamic_forward_network"
DYNAMIC_METHOD = "seasonal_phase_lag_baseline"
PARTICLE_GENERATOR_FAMILY = PARTICLE_TTD_GENERATOR_FAMILY
PARTICLE_METHOD = "particle_lag_mixing_baseline"


def _json_safe(value: Any) -> Any:
    """Convert NumPy/dataclass values to strict JSON without writing NaN."""

    if is_dataclass(value):
        return {item.name: _json_safe(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_safe(payload), indent=2, ensure_ascii=False, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_virtual_benchmark_config(path: Path | str) -> dict[str, Any]:
    """Load and minimally validate the frozen graph-TTD benchmark config."""

    config_path = Path(path)
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Benchmark config is not valid JSON: {config_path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Benchmark config root must be a JSON object.")
    if payload.get("schema") != CONFIG_SCHEMA:
        raise ValueError(
            f"Expected config schema {CONFIG_SCHEMA!r}, got {payload.get('schema')!r}."
        )
    for field in ("static", "dynamic", "required_comparators", "claim_gate"):
        if field not in payload:
            raise ValueError(f"Benchmark config is missing {field!r}.")
    if not isinstance(payload["static"], Mapping) or not isinstance(payload["dynamic"], Mapping):
        raise ValueError("Benchmark static and dynamic configuration must be objects.")
    return payload


def _static_metric_record(score: Mapping[str, Any]) -> dict[str, float | int | None]:
    interval = score["interval_recovery"]
    prediction = score["held_out_prediction"]
    selected = score["topology"]["selected_edges"]
    return {
        "interval_coverage": interval["conditional_coverage"],
        "mean_interval_width": interval["mean_interval_width"],
        "abstention_rate": interval["abstention_rate"],
        "held_out_prediction_mae": prediction["mae"],
        "held_out_rmse": prediction["rmse"],
        "prediction_rate": prediction["prediction_rate"],
        "topology_precision": selected["precision"],
        "topology_recall": selected["recall"],
        "false_edge_acceptance_rate": (
            None
            if selected["n_predicted"] == 0
            else float(selected["false_positive"]) / float(selected["n_predicted"])
        ),
    }


def _run_static_component(
    output_dir: Path, config: Mapping[str, Any]
) -> tuple[list[BenchmarkRecord], list[dict[str, object]], dict[str, object]]:
    component_dir = output_dir / "static"
    records: list[BenchmarkRecord] = []
    public_index: list[dict[str, object]] = []
    scoring_index: dict[str, object] = {}
    seeds = tuple(int(seed) for seed in config.get("seeds", ()))
    if not seeds:
        raise ValueError("Static benchmark configuration must declare at least one seed.")
    n_steps = int(config.get("n_time_steps", 300))
    protocol = default_static_ttd_graph_protocol(n_steps=n_steps)
    protocol_overrides: dict[str, object] = {}
    field_mapping = {
        "young_water_threshold_days": "young_water_cutoff_days",
        "input_sampling_fraction": "input_sampling_fraction",
        "output_sampling_fraction": "output_sampling_fraction",
        "input_missing_fraction": "input_missing_fraction",
        "output_missing_fraction": "output_missing_fraction",
        "held_out_fraction": "held_out_fraction",
    }
    for config_key, protocol_key in field_mapping.items():
        if config_key in config:
            protocol_overrides[protocol_key] = config[config_key]
    if protocol_overrides:
        protocol = replace(protocol, **protocol_overrides)

    for seed in seeds:
        case = generate_static_ttd_graph_case(seed=seed, protocol=protocol)
        if not case.verify_manifest():
            raise ValueError(f"Static case {case.observations.benchmark_id} failed its manifest check.")
        case_dir = component_dir / case.observations.benchmark_id
        paths = write_static_ttd_graph_case(case, case_dir)
        public_payload = case.observations.to_dict()
        assert_observation_view_is_truth_blind([public_payload])
        submissions = run_all_pre_registered_controls(case.observations, case.controls)
        scores = {
            control_id: score_static_ttd_graph_submission(
                case.truth, case.observations, submission
            )
            for control_id, submission in submissions.items()
        }
        _write_json(
            case_dir / "truth_blind_submissions.json",
            {control_id: submission.to_dict() for control_id, submission in submissions.items()},
        )
        _write_json(
            case_dir / "sealed_scores.json",
            {control_id: score.to_dict() for control_id, score in scores.items()},
        )
        public_index.append(
            {
                "component": "static",
                "case_id": case.observations.benchmark_id,
                "observation_artifact": str(paths["observations"].relative_to(output_dir)),
                "control_plan_artifact": str(paths["control_plan"].relative_to(output_dir)),
                "manifest_artifact": str(paths["manifest"].relative_to(output_dir)),
                "observation_sha256": _sha256(paths["observations"]),
            }
        )
        scoring_index[case.observations.benchmark_id] = {
            "sealed_score_artifact": str((case_dir / "sealed_scores.json").relative_to(output_dir)),
            "sealed_truth_artifact": str(paths["sealed_truth"].relative_to(output_dir)),
            "scores": {control_id: score.to_dict() for control_id, score in scores.items()},
        }
        for control_id, submission in submissions.items():
            if control_id not in STATIC_METHODS:
                raise ValueError(f"Static control {control_id!r} has no declared method mapping.")
            score = scores[control_id].metrics
            records.append(
                BenchmarkRecord(
                    case_id=case.observations.benchmark_id,
                    generator_family=STATIC_GENERATOR_FAMILY,
                    regime="static",
                    scenario="nominal",
                    method=STATIC_METHODS[control_id],
                    held_out=True,
                    truth_blind=bool(submission.diagnostics.get("truth_used") is False),
                    metrics=_static_metric_record(score),
                    abstention_reason=(
                        ";".join(
                            sorted(
                                {
                                    str(interval.reason)
                                    for interval in submission.node_intervals.values()
                                    if interval.status == "ABSTAIN" and interval.reason
                                }
                            )
                        )
                        or None
                    ),
                    notes=(
                        ("correct-topology arm is a harness-only conditional oracle control",)
                        if control_id == "correct"
                        else ()
                    ),
                )
            )
    return records, public_index, scoring_index


def _dynamic_public_payload(observations: object) -> dict[str, object]:
    payload = _json_safe(observations)
    if not isinstance(payload, dict):  # pragma: no cover - defensive dataclass gate
        raise TypeError("Dynamic observation payload was not serializable as an object.")
    return payload


def _dynamic_metrics(
    recoveries: Sequence[object], evaluation: object
) -> dict[str, float | int | None]:
    summary = getattr(evaluation, "summary")
    forecast_rmse = [
        float(getattr(recovery, "forecast_rmse"))
        for recovery in recoveries
        if getattr(recovery, "forecast_rmse") is not None
    ]
    n_recoveries = int(summary["n_recoveries"])
    n_abstained = sum(getattr(recovery, "status") == "ABSTAIN" for recovery in recoveries)
    return {
        "justified_recoveries": int(summary["justified_recoveries"]),
        "identifiable_recovery_misses": int(summary["identifiable_recovery_misses"]),
        "correct_abstentions": int(summary["correct_abstentions"]),
        "false_abstentions": int(summary["false_abstentions"]),
        "unsupported_point_estimates": int(summary["unsupported_point_estimates"]),
        "abstention_rate": n_abstained / n_recoveries if n_recoveries else None,
        "appropriate_abstention_rate": (
            int(summary["correct_abstentions"]) / n_abstained if n_abstained else None
        ),
        "held_out_rmse": (sum(forecast_rmse) / len(forecast_rmse) if forecast_rmse else None),
        "held_out_forecast_r2": summary["mean_heldout_forecast_r2"],
        "held_out_forecast_n": int(summary["heldout_forecast_n"]),
    }


def _run_dynamic_component(
    output_dir: Path, config: Mapping[str, Any]
) -> tuple[list[BenchmarkRecord], list[dict[str, object]], dict[str, object]]:
    component_dir = output_dir / "dynamic"
    records: list[BenchmarkRecord] = []
    public_index: list[dict[str, object]] = []
    scoring_index: dict[str, object] = {}
    seeds = tuple(int(seed) for seed in config.get("seeds", ()))
    if not seeds:
        raise ValueError("Dynamic benchmark configuration must declare at least one seed.")
    scenarios = tuple(str(item) for item in config.get("scenarios", ()))
    if not scenarios:
        raise ValueError("Dynamic benchmark configuration must declare scenarios.")

    for seed in seeds:
        for scenario in scenarios:
            dynamic_config = DynamicTTDGraphConfig(
                seed=seed,
                scenario=scenario,
                n_steps=int(config.get("n_time_steps", 208)),
                step_days=float(config.get("time_step_days", 7.0)),
                season_period_steps=int(config.get("season_period_days", 52)),
                observation_probability=float(config.get("observation_keep_probability", 0.82)),
                observation_noise_std=float(config.get("observation_noise_sd", 0.035)),
                training_fraction=1.0 - float(config.get("holdout_fraction", 0.30)),
            )
            truth, observations, manifest = generate_dynamic_ttd_case(dynamic_config)
            public_payload = _dynamic_public_payload(observations)
            assert_observation_view_is_truth_blind([public_payload])
            recoveries = recover_dynamic_ttd_baseline(observations)
            evaluation = evaluate_dynamic_ttd_recovery(
                truth,
                observations,
                recoveries,
                season_period_steps=dynamic_config.season_period_steps,
                n_phase_bins=dynamic_config.n_phase_bins,
            )
            case_dir = component_dir / observations.case_id
            _write_json(case_dir / "observations.json", public_payload)
            _write_json(case_dir / "truth_scoring_only.json", truth)
            _write_json(case_dir / "generator_manifest.json", manifest)
            _write_json(case_dir / "truth_blind_recoveries.json", recoveries)
            _write_json(case_dir / "sealed_evaluation.json", evaluation)
            public_index.append(
                {
                    "component": "dynamic",
                    "case_id": observations.case_id,
                    "observation_artifact": str((case_dir / "observations.json").relative_to(output_dir)),
                    "manifest_artifact": str((case_dir / "generator_manifest.json").relative_to(output_dir)),
                    "observation_sha256": _sha256(case_dir / "observations.json"),
                }
            )
            scoring_index[observations.case_id] = {
                "sealed_truth_artifact": str((case_dir / "truth_scoring_only.json").relative_to(output_dir)),
                "sealed_evaluation_artifact": str((case_dir / "sealed_evaluation.json").relative_to(output_dir)),
                "evaluation": _json_safe(evaluation),
            }
            records.append(
                BenchmarkRecord(
                    case_id=observations.case_id,
                    generator_family=DYNAMIC_GENERATOR_FAMILY,
                    regime="dynamic",
                    scenario=scenario,
                    method=DYNAMIC_METHOD,
                    held_out=True,
                    truth_blind=True,
                    metrics=_dynamic_metrics(recoveries, evaluation),
                    abstention_reason=(
                        ";".join(
                            sorted(
                                {
                                    str(getattr(recovery, "reason"))
                                    for recovery in recoveries
                                    if getattr(recovery, "status") == "ABSTAIN"
                                }
                            )
                        )
                        or None
                    ),
                    notes=(
                        "candidate topology is a declared conditional stress scenario; "
                        "it is not independent field-edge truth",
                    ),
                )
            )
    return records, public_index, scoring_index


def _particle_metric_record(
    recovery: ParticleTTDRecovery, score: ParticleTTDScore
) -> dict[str, float | int | None]:
    return {
        "young_water_absolute_error": score.young_water_absolute_error,
        "mean_lag_absolute_error_days": score.mean_lag_absolute_error_days,
        "held_out_rmse": score.held_out_rmse,
        "n_held_out": score.n_held_out,
        "abstained": 1 if recovery.status == "ABSTAIN" else 0,
        "appropriate_abstention": 1 if score.appropriate_abstention else 0,
        "calibration_rmse": recovery.calibration_rmse,
        "n_calibration_samples": recovery.n_calibration_samples,
    }


def _run_particle_component(
    output_dir: Path, config: Mapping[str, Any]
) -> tuple[list[BenchmarkRecord], list[dict[str, object]], dict[str, object]]:
    component_dir = output_dir / "particle"
    records: list[BenchmarkRecord] = []
    public_index: list[dict[str, object]] = []
    scoring_index: dict[str, object] = {}
    seeds = tuple(int(seed) for seed in config.get("seeds", (1729,)))
    scenarios = tuple(
        str(item)
        for item in config.get(
            "scenarios",
            (
                "nominal",
                "sparse",
                "wrong_forcing",
                "local_recharge_misspecification",
                "wrong_topology",
            ),
        )
    )
    n_particles = int(config.get("n_particles", 600))
    n_steps = int(config.get("n_steps", 160))
    max_lag_steps = int(config.get("max_lag_steps", 36))

    for seed in seeds:
        for scenario in scenarios:
            particle_config = ParticleTTDConfig(
                seed=seed,
                scenario=scenario,
                n_particles=n_particles,
                n_steps=n_steps,
                max_lag_steps=max_lag_steps,
            )
            case = generate_independent_particle_ttd_case(particle_config)
            assert_particle_observation_view_is_truth_blind(case.observations)
            recovery = recover_particle_ttd_baseline(case.observations)
            score = score_particle_ttd_recovery(case, recovery)

            case_id = f"particle_{scenario}_seed_{seed}"
            case_dir = component_dir / case_id
            public_payload = _json_safe(case.observations)
            truth_payload = _json_safe(case.truth)
            recovery_payload = _json_safe(recovery)
            score_payload = _json_safe(score)

            public_manifest = particle_ttd_public_manifest(case.observations)
            _write_json(case_dir / "observations.json", public_payload)
            _write_json(case_dir / "public_manifest.json", public_manifest)
            _write_json(case_dir / "truth_scoring_only.json", truth_payload)
            _write_json(case_dir / "truth_blind_recovery.json", recovery_payload)
            _write_json(case_dir / "sealed_score.json", score_payload)

            public_index.append(
                {
                    "component": "particle",
                    "case_id": case_id,
                    "observation_artifact": str((case_dir / "observations.json").relative_to(output_dir)),
                    "manifest_artifact": str((case_dir / "public_manifest.json").relative_to(output_dir)),
                    "observation_sha256": _sha256(case_dir / "observations.json"),
                }
            )
            scoring_index[case_id] = {
                "sealed_truth_artifact": str((case_dir / "truth_scoring_only.json").relative_to(output_dir)),
                "sealed_score_artifact": str((case_dir / "sealed_score.json").relative_to(output_dir)),
                "score": score_payload,
            }
            records.append(
                BenchmarkRecord(
                    case_id=case_id,
                    generator_family=PARTICLE_TTD_GENERATOR_FAMILY,
                    regime="particle",
                    scenario=scenario,
                    method=PARTICLE_METHOD,
                    held_out=True,
                    truth_blind=True,
                    metrics=_particle_metric_record(recovery, score),
                    abstention_reason=recovery.reason,
                    notes=(
                        "independent particle Monte-Carlo transport network generator family",
                    ),
                )
            )
    return records, public_index, scoring_index


def run_ttd_graph_virtual_benchmark(
    output_dir: Path | str,
    *,
    config_path: Path | str = "configs/ttd_graph_virtual_benchmark_v1.json",
    overwrite: bool = False,
) -> dict[str, object]:
    """Execute the static/dynamic virtual benchmark and write auditable artifacts.

    The runner intentionally does not download WATRES or use Ghana field data.
    WATRES remains an opt-in external component panel, while Ghana remains a
    field-readiness/campaign-design panel under the protocol.
    """

    config = load_virtual_benchmark_config(config_path)
    target = Path(output_dir)
    if target.exists() and any(target.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite populated benchmark directory: {target}"
        )
    target.mkdir(parents=True, exist_ok=True)

    static_records, static_public, static_scores = _run_static_component(
        target, config["static"]
    )
    dynamic_records, dynamic_public, dynamic_scores = _run_dynamic_component(
        target, config["dynamic"]
    )
    particle_config = config.get(
        "particle",
        {
            "seeds": [1729],
            "n_particles": 600,
            "n_steps": 160,
            "max_lag_steps": 36,
            "scenarios": [
                "nominal",
                "sparse",
                "wrong_forcing",
                "local_recharge_misspecification",
                "wrong_topology",
            ],
        },
    )
    particle_records, particle_public, particle_scores = _run_particle_component(
        target, particle_config
    )

    all_records = static_records + dynamic_records + particle_records
    static_readiness = assess_claim_readiness(
        static_records,
        required_generator_families=[STATIC_GENERATOR_FAMILY],
        required_comparators=list(STATIC_METHODS.values()),
        required_scenarios=["nominal"],
    )
    dynamic_readiness = assess_claim_readiness(
        dynamic_records,
        required_generator_families=[DYNAMIC_GENERATOR_FAMILY],
        required_comparators=[DYNAMIC_METHOD],
        required_scenarios=list(config["dynamic"]["scenarios"]),
    )
    particle_readiness = assess_claim_readiness(
        particle_records,
        required_generator_families=[PARTICLE_TTD_GENERATOR_FAMILY],
        required_comparators=[PARTICLE_METHOD],
        required_scenarios=list(particle_config.get("scenarios", ["nominal"])),
    )
    programme_readiness = assess_claim_readiness(
        all_records,
        required_generator_families=[
            STATIC_GENERATOR_FAMILY,
            PARTICLE_TTD_GENERATOR_FAMILY,
        ],
        required_comparators=list(STATIC_METHODS.values()) + [DYNAMIC_METHOD, PARTICLE_METHOD],
        required_scenarios=["nominal", "sparse_sampling", "wrong_forcing", "wrong_topology"],
    )
    readiness = {
        "static_component": static_readiness,
        "dynamic_component": dynamic_readiness,
        "particle_component": particle_readiness,
        "programme": programme_readiness,
        "external_watres": {
            "status": "NOT_RUN",
            "reason": "WATRES is opt-in and is an external catchment-level temporal-TTD component panel only.",
        },
        "field_transfer": {
            "status": "DEFERRED",
            "reason": "Current Ghana data are a readiness/campaign-design panel, not dynamic graph-TTD field validation.",
        },
    }
    public_index = static_public + dynamic_public + particle_public
    scoring_index = {
        "static": static_scores,
        "dynamic": dynamic_scores,
        "particle": particle_scores,
    }
    generator_provenance = {
        "independent_from_hydrosheaf_inference": True,
        "static_generator_family": STATIC_GENERATOR_FAMILY,
        "dynamic_generator_family": DYNAMIC_GENERATOR_FAMILY,
        "particle_generator_family": PARTICLE_TTD_GENERATOR_FAMILY,
        "source_files": {
            "static": str(Path(__file__).with_name("ttd_graph_static.py")),
            "dynamic": str(Path(__file__).with_name("ttd_graph_dynamic.py")),
            "particle": str(Path(__file__).with_name("independent_particle_ttd.py")),
        },
        "second_independent_generator_family": PARTICLE_TTD_GENERATOR_FAMILY,
    }
    artifacts = write_virtual_benchmark_artifacts(
        target / "programme",
        run_id=f"TTD-GRAPH-VIRTUAL-{len(all_records):03d}",
        protocol_path=config["protocol"],
        config=config,
        observations=public_index,
        records=all_records,
        truth_for_scoring=scoring_index,
        generator_provenance=generator_provenance,
        readiness=readiness,
        overwrite=overwrite,
    )
    return {
        "artifacts": artifacts,
        "static_record_count": len(static_records),
        "dynamic_record_count": len(dynamic_records),
        "particle_record_count": len(particle_records),
        "readiness": readiness,
    }


__all__ = [
    "CONFIG_SCHEMA",
    "DYNAMIC_GENERATOR_FAMILY",
    "DYNAMIC_METHOD",
    "PARTICLE_GENERATOR_FAMILY",
    "PARTICLE_METHOD",
    "STATIC_METHODS",
    "load_virtual_benchmark_config",
    "run_ttd_graph_virtual_benchmark",
]
