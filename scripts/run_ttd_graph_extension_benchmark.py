"""Master execution runner for HydroSheaf M9.1–M9.3 extension benchmark suite.

Usage:
    python scripts/run_ttd_graph_extension_benchmark.py --config configs/ttd_graph_extension_v1.json --output .codex_work/m9_extension_verification_v1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from hydrosheaf.benchmarks.ttd_graph_dynamic_kernel import (
    DYNAMIC_SCENARIOS,
    DynamicKernelBenchmarkConfig,
    evaluate_dynamic_kernel_benchmark,
    generate_dynamic_kernel_case,
)
from hydrosheaf.benchmarks.ttd_graph_multitracer import (
    TRACER_PANELS,
    MultiTracerBenchmarkConfig,
    evaluate_multitracer_recovery,
    generate_multitracer_case,
)
from hydrosheaf.models.ttd_losses import LossConfig, evaluate_composite_ttd_loss
from hydrosheaf.nuclear.dynamic_kernel_inversion import (
    DynamicTTDInversionConfig,
    DynamicTTDRecovery,
    solve_dynamic_node_inversion,
)
from hydrosheaf.nuclear.multi_tracer_graph_inversion import (
    MultiTracerGraphConfig,
    solve_joint_multitracer_node_inversion,
)
from hydrosheaf.nuclear.graph_tracer_forward import validate_candidate_graph


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def run_benchmark(config_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(config_path, "r", encoding="utf-8") as f:
        cfg_dict = json.load(f)

    start_time = datetime.now(timezone.utc).isoformat()
    config_hash = hashlib.sha256(json.dumps(cfg_dict, sort_keys=True).encode()).hexdigest()

    seeds = cfg_dict["global_parameters"]["seeds"]
    scenarios = cfg_dict.get("scenarios", DYNAMIC_SCENARIOS)
    panels = cfg_dict.get("multitracer_panels", TRACER_PANELS)

    # Inversion config for dynamic kernel benchmark
    inv_params = cfg_dict["regularization"]
    max_lag_steps = int(cfg_dict["global_parameters"].get("max_lag_steps", 18))
    step_days = float(cfg_dict["global_parameters"]["step_days"])
    dyn_lag_grid = tuple(float(x) for x in np.arange(max_lag_steps + 1) * step_days)

    dyn_inv_cfg = DynamicTTDInversionConfig(
        lag_grid_days=dyn_lag_grid,
        young_water_cutoff_days=float(
            cfg_dict["global_parameters"].get("young_water_cutoff_days", 4.0 * step_days)
        ),
        step_days=step_days,
        kernel_mode="phase",
        season_period_steps=cfg_dict["global_parameters"]["season_period_steps"],
        n_phase_bins=cfg_dict["global_parameters"]["n_phase_bins"],
        lambda_lag=inv_params["lambda_lag"],
        lambda_time=inv_params["lambda_time"],
        lambda_edge_sparsity=inv_params["lambda_edge_sparsity"],
        max_condition_number=inv_params["max_condition_number"],
        residual_cone_tolerance=inv_params["residual_cone_tolerance"],
        min_observed_samples=inv_params["min_observed_samples"],
        min_pairs_per_phase=inv_params["min_pairs_per_phase"],
        min_identified_phases=inv_params["min_identified_phases"],
        min_r2=inv_params["min_r2"],
        min_forecast_r2=inv_params.get("min_forecast_r2", 0.0),
        require_local_input=bool(inv_params.get("require_local_input", False)),
    )

    loss_cfg = LossConfig(
        time_loss=cfg_dict["loss_weights"]["time_loss"],
        huber_delta=cfg_dict["loss_weights"]["huber_delta"],
        wasserstein_ttd_weight=cfg_dict["loss_weights"]["wasserstein_ttd_weight"],
        wasserstein_series_weight=cfg_dict["loss_weights"]["wasserstein_series_weight"],
        spectral_weight=cfg_dict["loss_weights"]["spectral_weight"],
        wavelet_weight=cfg_dict["loss_weights"]["wavelet_weight"],
        coherence_weight=cfg_dict["loss_weights"].get("coherence_weight", 0.0),
    )

    inference_records: list[dict[str, Any]] = []
    truth_scoring_records: list[dict[str, Any]] = []
    loss_decomposition_records: list[dict[str, Any]] = []
    dynamic_kernel_summaries: list[dict[str, Any]] = []
    abstention_diagnostics: list[dict[str, Any]] = []
    held_out_metrics: list[dict[str, Any]] = []

    print(f"Running M9.1 Dynamic Kernel Benchmarks ({len(scenarios)} scenarios x {len(seeds)} seeds)...")

    # 1. Run Dynamic Kernel Benchmarks
    for scen in scenarios:
        for seed in seeds:
            b_cfg = DynamicKernelBenchmarkConfig(
                seed=seed,
                scenario=scen,
                n_steps=cfg_dict["global_parameters"]["n_steps"],
                step_days=cfg_dict["global_parameters"]["step_days"],
                season_period_steps=cfg_dict["global_parameters"]["season_period_steps"],
                observation_probability=cfg_dict["global_parameters"]["observation_probability"],
                observation_noise_std=cfg_dict["global_parameters"]["observation_noise_std"],
                training_fraction=cfg_dict["global_parameters"]["training_fraction"],
            )
            truth, obs = generate_dynamic_kernel_case(b_cfg)

            # Solve for each target node
            all_recovs = {}
            graph_ok, graph_reason = validate_candidate_graph(
                ("R", "A", "B", "C"), obs.candidate_edges, source_node="R"
            )
            for target_node in ("A", "B", "C"):
                target_obs = obs.node_observations[target_node]
                cal_mask = np.zeros(len(obs.time_steps), dtype=bool)
                cal_mask[: obs.training_end_step] = obs.observation_masks[target_node][: obs.training_end_step]
                cal_times = np.where(cal_mask)[0]
                cal_vals = target_obs[cal_times]

                ho_mask = np.zeros(len(obs.time_steps), dtype=bool)
                ho_mask[obs.training_end_step :] = obs.observation_masks[target_node][obs.training_end_step :]
                ho_times = np.where(ho_mask)[0]
                ho_vals = target_obs[ho_times]

                # Identify parents from candidate_edges
                parents = {}
                for u, v in obs.candidate_edges:
                    if v == target_node:
                        if u == "R":
                            parents[u] = obs.regional_forcing
                        elif u in obs.node_observations:
                            parents[u] = obs.node_observations[u]

                # Local recharge is supplied only when it is part of the
                # visible observation contract.  No scenario label is passed
                # to the estimator and no input is altered based on the label.
                loc_in = obs.local_inputs.get(target_node)

                if not graph_ok:
                    node_recovs = {
                        f"{parent}->{target_node}": DynamicTTDRecovery(
                            method_id="hydrosheaf_dynamic_ttd_inversion_v2",
                            status="ABSTAIN",
                            reason_codes=(f"candidate_graph_validation_failed:{graph_reason}",),
                            edge_id=f"{parent}->{target_node}",
                            target_node=target_node,
                            configuration_hash="",
                        )
                        for parent in parents
                    }
                else:
                    node_recovs = solve_dynamic_node_inversion(
                        target_node=target_node,
                        target_times=cal_times,
                        target_values=cal_vals,
                        candidate_parents=parents,
                        local_input=loc_in,
                        config=dyn_inv_cfg,
                        holdout_times=ho_times,
                        holdout_values=ho_vals,
                        loss_config=loss_cfg,
                    )
                all_recovs.update(node_recovs)

            # Serialize inference records (truth-blind)
            for eid, r in all_recovs.items():
                r_dict = r.to_dict()
                r_dict["case_id"] = obs.case_id
                r_dict["scenario"] = obs.scenario
                r_dict["seed"] = seed
                inference_records.append(r_dict)
                if r.status == "ABSTAIN":
                    abstention_diagnostics.append(
                        {
                            "case_id": obs.case_id,
                            "scenario": obs.scenario,
                            "seed": seed,
                            "edge_id": eid,
                            "reason_codes": list(r.reason_codes),
                            "condition_number": r.condition_number,
                        }
                    )
                else:
                    if r.forecast_r2 is not None:
                        held_out_metrics.append(
                            {
                                "case_id": obs.case_id,
                                "scenario": obs.scenario,
                                "seed": seed,
                                "edge_id": eid,
                                "forecast_r2": r.forecast_r2,
                                "forecast_rmse": r.forecast_rmse,
                            }
                        )

            # Evaluate against sealed truth
            score = evaluate_dynamic_kernel_benchmark(truth, obs, all_recovs)
            truth_scoring_records.append(score)
            dynamic_kernel_summaries.append(
                {
                    "case_id": obs.case_id,
                    "scenario": obs.scenario,
                    "seed": seed,
                    "justified_recoveries": score["justified_recoveries"],
                    "correct_abstentions": score["correct_abstentions"],
                    "false_abstentions": score["false_abstentions"],
                    "unsupported_estimates": score["unsupported_estimates"],
                    "mean_w1_distance": score["mean_w1_distance"],
                    "empirical_coverage": score["empirical_coverage"],
                }
            )

            # Use the fitted node-B prediction, never the sealed truth, for
            # time/frequency loss decomposition.  TTD W1 is scored separately
            # below against revealed truth and is explicitly labelled scorer-only.
            node_b_records = [
                rec for rec in all_recovs.values() if rec.target_node == "B" and rec.status == "RECOVERED"
            ]
            if node_b_records:
                decomp = dict(node_b_records[0].loss_decomposition)
                scorer_w1 = []
                for rec in node_b_records:
                    if rec.edge_id in truth.true_edge_kernels and rec.estimated_kernel is not None:
                        true_kernel = truth.true_edge_kernels[rec.edge_id]
                        lags = np.asarray(obs.metadata["lag_grid_days"], dtype=float)
                        scorer_w1.append(
                            float(
                                np.mean(
                                    [
                                        evaluate_composite_ttd_loss(
                                            predicted_signal=np.array([0.0, 0.0]),
                                            observed_signal=np.array([0.0, 0.0]),
                                            time_grid=np.array([0.0, 1.0]),
                                            predicted_ttd=rec.estimated_kernel[t],
                                            observed_ttd=true_kernel[t],
                                            lag_grid=lags,
                                            config=LossConfig(wasserstein_ttd_weight=1.0),
                                        ).decomposition["wasserstein_ttd"]
                                        for t in range(min(len(true_kernel), len(rec.estimated_kernel)))
                                    ]
                                )
                            )
                        )
                decomp["scorer_only_wasserstein_ttd_mean"] = float(np.mean(scorer_w1)) if scorer_w1 else None
                decomp["scorer_only_truth_revealed"] = True
                loss_decomposition_records.append(
                    {
                        "case_id": obs.case_id,
                        "scenario": obs.scenario,
                        "seed": seed,
                        "status": "FITTED_PREDICTION",
                        "decomposition": decomp,
                    }
                )
            else:
                loss_decomposition_records.append(
                    {
                        "case_id": obs.case_id,
                        "scenario": obs.scenario,
                        "seed": seed,
                        "status": "NO_RECOVERED_NODE_B_EDGE",
                        "decomposition": {"status": "not_available", "reason": "all_node_B_edges_abstained"},
                    }
                )

    print(f"Running M9.3 Multi-Tracer Benchmarks ({len(panels)} panels x {len(seeds)} seeds)...")

    # 2. Run Multi-Tracer Benchmarks
    per_tracer_records: list[dict[str, Any]] = []

    for panel in panels:
        for seed in seeds:
            mt_cfg = MultiTracerBenchmarkConfig(
                seed=seed,
                panel=panel,
                n_steps=156,
                step_days=14.0,
                season_period_steps=26,
            )
            truth_mt, obs_mt = generate_multitracer_case(mt_cfg)

            cal_times = obs_mt.time_steps[: obs_mt.training_end_step]
            ho_times = obs_mt.time_steps[obs_mt.training_end_step :]

            target_cal_obs = {
                t: obs_mt.node_tracer_observations["B"][t][: obs_mt.training_end_step]
                for t in obs_mt.available_tracers
            }
            target_ho_obs = {
                t: obs_mt.node_tracer_observations["B"][t][obs_mt.training_end_step :]
                for t in obs_mt.available_tracers
            }

            parents_mt = {"A": obs_mt.node_tracer_observations["A"]}
            local_mt = obs_mt.local_tracer_inputs.get("B", {})

            mt_inv_cfg = DynamicTTDInversionConfig(
                lag_grid_days=tuple(float(v) for v in obs_mt.metadata["lag_grid_days"]),
                season_period_steps=26,
                step_days=14.0,
                n_phase_bins=4,
                min_r2=0.10,
            )

            res_mt = solve_joint_multitracer_node_inversion(
                target_node="B",
                target_times=cal_times,
                tracer_observations=target_cal_obs,
                candidate_parents=parents_mt,
                local_tracer_inputs=local_mt,
                config=MultiTracerGraphConfig(inversion_config=mt_inv_cfg, enforce_holdout_gate=True),
                holdout_times=ho_times,
                holdout_tracer_observations=target_ho_obs,
            )

            score_mt = evaluate_multitracer_recovery(truth_mt, obs_mt, res_mt)
            per_tracer_records.append(
                {
                    "case_id": obs_mt.case_id,
                    "panel": panel,
                    "seed": seed,
                    "status": res_mt.status,
                    "decision": score_mt.get("decision", "UNKNOWN"),
                    "conflict_detected": res_mt.conflict_detected,
                    "per_tracer_rmse": dict(res_mt.per_tracer_rmse),
                    "held_out_per_tracer_rmse": dict(res_mt.held_out_per_tracer_rmse),
                    "held_out_per_tracer_r2": dict(res_mt.held_out_per_tracer_r2),
                    "held_out_per_tracer_normalized_rmse": dict(res_mt.held_out_per_tracer_normalized_rmse),
                    "loto_sensitivity": dict(res_mt.loto_sensitivity),
                    "effective_rank": res_mt.effective_rank,
                    "condition_number": res_mt.condition_number,
                    "diagnostics": dict(res_mt.diagnostics),
                    "truth_commitment": score_mt.get("truth_commitment"),
                    "observation_commitment": score_mt.get("observation_commitment"),
                }
            )

    # 3. Write All Artifact Files
    files_to_write = {
        "inference_records.json": inference_records,
        "truth_scoring_only.json": truth_scoring_records,
        "loss_decomposition.json": loss_decomposition_records,
        "dynamic_kernel_summaries.json": dynamic_kernel_summaries,
        "per_tracer_results.json": per_tracer_records,
        "held_out_metrics.json": held_out_metrics,
        "abstention_diagnostics.json": abstention_diagnostics,
    }

    for fname, data in files_to_write.items():
        with open(output_dir / fname, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, default=str)

    # Provenance and claim readiness
    total_false_abstentions = sum(s["false_abstentions"] for s in dynamic_kernel_summaries)
    total_dynamic_unsupported = sum(s["unsupported_estimates"] for s in dynamic_kernel_summaries)
    non_conflict_mt = [r for r in per_tracer_records if r["panel"] != "conflicting_tracers"]
    mt_gate_passed = bool(non_conflict_mt) and all(
        r["decision"] == "JUSTIFIED_RECOVERY"
        and bool(r.get("diagnostics", {}).get("holdout_gate_passed", True))
        for r in non_conflict_mt
    )
    loss_gate_passed = len(loss_decomposition_records) == len(scenarios) * len(seeds) and all(
        isinstance(record.get("decomposition"), dict) and bool(record["decomposition"]) for record in loss_decomposition_records
    )
    truth_blind_gate_passed = all(
        record.get("truth_blindness_declared") is True
        for record in inference_records
        if record.get("status") == "RECOVERED"
    )
    # This is deliberately a strict protocol gate.  A completed runner is not
    # called claim-ready when it has false abstentions, unsupported estimates,
    # missing loss records, or failed multi-tracer holdout gates.
    controlled_ready = bool(
        loss_gate_passed
        and truth_blind_gate_passed
        and total_dynamic_unsupported == 0
        and total_false_abstentions == 0
        and mt_gate_passed
    )
    claim_readiness = {
        "protocol_id": "HS-M9-EXTENSION-V1",
        "claim_boundary": cfg_dict["claim_boundary"],
        "controlled_synthetic_ready": controlled_ready,
        "protocol_completed": True,
        "truth_blindness_enforced": truth_blind_gate_passed,
        "loss_decomposition_complete": loss_gate_passed,
        "multitracer_holdout_gate_passed": mt_gate_passed,
        "field_validation_claimed": False,
        "universal_superiority_claimed": False,
        "total_dynamic_cases": len(scenarios) * len(seeds),
        "total_multitracer_cases": len(panels) * len(seeds),
        "total_justified_recoveries": sum(s["justified_recoveries"] for s in dynamic_kernel_summaries),
        "total_correct_abstentions": sum(s["correct_abstentions"] for s in dynamic_kernel_summaries),
        "total_false_abstentions": total_false_abstentions,
        "total_unsupported_estimates": total_dynamic_unsupported,
        "multitracer_conflicts_detected": sum(1 for r in per_tracer_records if r["conflict_detected"]),
    }
    with open(output_dir / "claim_readiness.json", "w", encoding="utf-8") as f:
        json.dump(claim_readiness, f, indent=2)

    repo_root = Path(__file__).resolve().parents[1]
    source_files = [
        Path(__file__).resolve(),
        repo_root / "hydrosheaf/benchmarks/ttd_graph_dynamic_kernel.py",
        repo_root / "hydrosheaf/benchmarks/ttd_graph_multitracer.py",
        repo_root / "hydrosheaf/nuclear/dynamic_kernel_inversion.py",
        repo_root / "hydrosheaf/nuclear/multi_tracer_graph_inversion.py",
        repo_root / "hydrosheaf/models/ttd_losses.py",
    ]
    source_hashes = {str(path.relative_to(repo_root)): sha256_of_file(path) for path in source_files if path.exists()}
    gen_provenance = {
        "generator_family": "independent_dynamic_and_multitracer_synthetic_v1",
        "scenarios": scenarios,
        "multitracer_panels": panels,
        "config_hash": config_hash,
        "source_hashes": source_hashes,
        "truth_commitments_are_recorded_in_scoring_artifact": True,
        "scenario_labels_passed_to_estimator": False,
        "timestamp": start_time,
    }
    with open(output_dir / "generator_provenance.json", "w", encoding="utf-8") as f:
        json.dump(gen_provenance, f, indent=2)

    # Manifest and artifact hashes
    hashes = {}
    for p in output_dir.glob("*.json"):
        hashes[p.name] = sha256_of_file(p)

    with open(output_dir / "artifact_hashes.json", "w", encoding="utf-8") as f:
        json.dump(hashes, f, indent=2)

    hashes["artifact_hashes.json"] = sha256_of_file(output_dir / "artifact_hashes.json")

    run_manifest = {
        "run_id": f"RUN-M9-EXTENSION-{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}",
        "config_path": str(config_path),
        "config_hash": config_hash,
        "start_time": start_time,
        "end_time": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "status": "COMPLETED",
        "claim_boundary": cfg_dict["claim_boundary"],
        "artifact_hashes": hashes,
    }
    with open(output_dir / "run_manifest.json", "w", encoding="utf-8") as f:
        json.dump(run_manifest, f, indent=2)

    print(f"[OK] Master M9 extension benchmark completed successfully. Artifacts in {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run HydroSheaf M9 Extension Benchmark Suite")
    parser.add_argument("--config", required=True, type=Path, help="Path to ttd_graph_extension_v1.json")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    args = parser.parse_args()

    run_benchmark(args.config, args.output)


if __name__ == "__main__":
    main()
