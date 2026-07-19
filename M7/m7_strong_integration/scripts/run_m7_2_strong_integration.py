"""Run the locked M7.2 strong-integration validation programme."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import sys
from typing import Dict, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from field_prequential import run_field_prequential  # noqa: E402
from independent_modflow_generator import (  # noqa: E402
    IndependentAquifer,
    generate_independent_aquifer,
)
from strong_inference import (  # noqa: E402
    FUSION_FEATURES,
    StrongInferenceResult,
    run_strong_inference,
)


DEV_SEEDS = (2101, 2102, 2103, 2104, 2105, 2106)
TEST_SEEDS = (
    3101,
    3102,
    3103,
    3104,
    3105,
    3106,
    3107,
    3108,
    3109,
    3110,
    3111,
    3112,
)
CONFIRMATORY_TEST_SEEDS = tuple(range(4101, 4113))
DIRECTION_GATE_AGE_COST_MAX = 0.05
BASELINE_FEATURES = (
    "hydraulic_logit",
    "negative_chemistry_log_objective",
)


def _json_default(value: object) -> object:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    raise TypeError(f"Cannot JSON encode {type(value)!r}")


def _truth_family(process: str) -> str:
    if process in {"carbonate_weathering", "carbonate_precipitation"}:
        return "carbonate"
    if process == "silicate_weathering":
        return "silicate_exchange"
    if process == "sulfate_reduction":
        return "sulfate_reduction"
    if process == "denitrification":
        return "denitrification"
    if process == "iron_reduction":
        return "iron_reduction"
    raise ValueError(f"Unknown truth process {process}")


def _label_rows(
    case: IndependentAquifer,
    result: StrongInferenceResult,
    split: str,
) -> list[Dict[str, object]]:
    truth = {f"{u}->{v}" for u, v in case.true_edges}
    rows: list[Dict[str, object]] = []
    for feature in result.edge_features:
        row = dict(feature)
        edge_id = str(row["edge_id"])
        row.update(
            {
                "seed": case.seed,
                "split": split,
                "is_true_edge": int(edge_id in truth),
                "true_process": case.true_processes.get(edge_id),
                "true_process_family": (
                    _truth_family(case.true_processes[edge_id])
                    if edge_id in case.true_processes
                    else None
                ),
            }
        )
        rows.append(row)
    return rows


def _fit_fusion_model(
    frame: pd.DataFrame,
    feature_names: Sequence[str],
) -> Dict[str, object]:
    x = frame[list(feature_names)].to_numpy(float)
    y = frame["is_true_edge"].to_numpy(int)
    means = x.mean(axis=0)
    scales = x.std(axis=0)
    scales[scales < 1e-9] = 1.0
    standardized = (x - means) / scales
    model = LogisticRegression(
        C=1.0,
        class_weight="balanced",
        solver="lbfgs",
        max_iter=5000,
        random_state=1947,
    )
    model.fit(standardized, y)
    return {
        "feature_names": list(feature_names),
        "means": means.tolist(),
        "scales": scales.tolist(),
        "coefficients": model.coef_[0].tolist(),
        "intercept": float(model.intercept_[0]),
        "training_split": "development_external_simulator_cases_only",
        "class_weight": "balanced",
        "regularization_C": 1.0,
    }


def _predict_fusion(
    frame: pd.DataFrame,
    model: Mapping[str, object],
) -> np.ndarray:
    names = list(model["feature_names"])
    x = frame[names].to_numpy(float)
    means = np.asarray(model["means"], float)
    scales = np.asarray(model["scales"], float)
    coefficients = np.asarray(model["coefficients"], float)
    linear = float(model["intercept"]) + ((x - means) / scales) @ coefficients
    probability = 1.0 / (
        1.0 + np.exp(-np.clip(linear, -40.0, 40.0))
    )
    if model.get("kind") == "age_compatibility_gate":
        incompatible = (
            frame["age_cost"].to_numpy(float)
            > float(model["age_cost_max"])
        )
        probability[incompatible] = np.minimum(
            probability[incompatible],
            float(model.get("incompatible_probability", 1.0e-6)),
        )
    return probability


def _age_permuted_frame(
    frame: pd.DataFrame,
    rng: np.random.Generator,
    age_evidence_mode: str,
) -> pd.DataFrame:
    permuted = frame.copy()
    if age_evidence_mode == "direction_gate":
        permuted["age_cost"] = frame.groupby("seed")["age_cost"].transform(
            lambda values: rng.permutation(values.to_numpy())
        )
    else:
        permuted["negative_age_cost"] = frame.groupby("seed")[
            "negative_age_cost"
        ].transform(lambda values: rng.permutation(values.to_numpy()))
    return permuted


def _select_threshold(
    truth: np.ndarray,
    probability: np.ndarray,
) -> float:
    candidates = np.unique(
        np.concatenate(
            [
                np.linspace(0.05, 0.95, 181),
                np.asarray(probability, dtype=float),
            ]
        )
    )
    best = (-np.inf, 0.5)
    for threshold in candidates:
        prediction = probability >= threshold
        mcc = matthews_corrcoef(truth, prediction)
        f1 = f1_score(truth, prediction, zero_division=0)
        objective = mcc + 0.1 * f1
        if objective > best[0]:
            best = (objective, float(threshold))
    return best[1]


def _metrics(
    truth: np.ndarray,
    probability: np.ndarray,
    threshold: float,
) -> Dict[str, float]:
    prediction = probability >= float(threshold)
    return {
        "pr_auc": float(average_precision_score(truth, probability)),
        "roc_auc": float(roc_auc_score(truth, probability)),
        "brier": float(brier_score_loss(truth, probability)),
        "precision": float(precision_score(truth, prediction, zero_division=0)),
        "recall": float(recall_score(truth, prediction, zero_division=0)),
        "f1": float(f1_score(truth, prediction, zero_division=0)),
        "mcc": float(matthews_corrcoef(truth, prediction)),
        "threshold_locked_on_development": float(threshold),
    }


def _paired_case_bootstrap(
    test_frame: pd.DataFrame,
    methods: Mapping[str, tuple[str, float]],
    *,
    n_bootstrap: int = 10_000,
) -> pd.DataFrame:
    seeds = np.asarray(sorted(test_frame["seed"].unique()), dtype=int)
    rng = np.random.default_rng(7781)
    rows = []
    full_column, full_threshold = methods["hydraulic_chemistry_age"]
    for baseline_name in ("hydraulic_chemistry", "age_permuted_control"):
        baseline_column, baseline_threshold = methods[baseline_name]
        differences = {"f1": [], "pr_auc": []}
        for _ in range(int(n_bootstrap)):
            sampled = rng.choice(seeds, size=len(seeds), replace=True)
            pieces = [
                test_frame[test_frame["seed"] == seed] for seed in sampled
            ]
            boot = pd.concat(pieces, ignore_index=True)
            truth = boot["is_true_edge"].to_numpy(int)
            full = _metrics(
                truth,
                boot[full_column].to_numpy(float),
                full_threshold,
            )
            baseline = _metrics(
                truth,
                boot[baseline_column].to_numpy(float),
                baseline_threshold,
            )
            for metric in differences:
                differences[metric].append(full[metric] - baseline[metric])
        for metric, values in differences.items():
            rows.append(
                {
                    "contrast": f"hydraulic_chemistry_age_minus_{baseline_name}",
                    "metric": metric,
                    "mean_difference": float(np.mean(values)),
                    "ci95_low": float(np.quantile(values, 0.025)),
                    "ci95_high": float(np.quantile(values, 0.975)),
                    "resampling_unit": "independent_MODFLOW_case",
                    "n_cases": len(seeds),
                    "n_bootstrap": int(n_bootstrap),
                }
            )
    return pd.DataFrame(rows)


def _topology_summary(
    case: IndependentAquifer,
    result: StrongInferenceResult,
) -> Dict[str, object]:
    truth = {f"{u}->{v}" for u, v in case.true_edges}
    inferred = set(result.posterior_diagnostics.get("map_edges", []))
    tp = len(truth & inferred)
    precision = tp / len(inferred) if inferred else 0.0
    recall = tp / len(truth) if truth else 0.0
    edge_rhat = [
        float(value)
        for value in result.posterior_diagnostics.get("edge_r_hat", {}).values()
        if value is not None
    ]
    edge_ess = [
        float(value)
        for value in result.posterior_diagnostics.get("edge_ess", {}).values()
    ]
    return {
        "seed": case.seed,
        "map_n_edges": len(inferred),
        "map_precision": precision,
        "map_recall": recall,
        "map_f1": (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        ),
        "n_edges_r_hat": result.posterior_diagnostics.get("n_edges_r_hat"),
        "n_edges_ess_bulk": result.posterior_diagnostics.get(
            "n_edges_ess_bulk"
        ),
        "n_edges_ess_tail": result.posterior_diagnostics.get(
            "n_edges_ess_tail"
        ),
        "edge_r_hat_max": max(edge_rhat) if edge_rhat else None,
        "edge_ess_min": min(edge_ess) if edge_ess else None,
        "acceptance_rate": result.posterior_diagnostics.get("acceptance_rate"),
        "converged": bool(
            result.posterior_diagnostics.get("n_edges_r_hat") is not None
            and float(result.posterior_diagnostics["n_edges_r_hat"]) <= 1.01
            and float(
                result.posterior_diagnostics.get("n_edges_ess_bulk", 0.0)
            )
            >= 400.0
            and float(
                result.posterior_diagnostics.get("n_edges_ess_tail", 0.0)
            )
            >= 400.0
            and (max(edge_rhat) if edge_rhat else float("inf")) <= 1.01
            and (min(edge_ess) if edge_ess else 0.0) >= 400.0
        ),
    }


def _age_summary(
    case: IndependentAquifer,
    result: StrongInferenceResult,
) -> Dict[str, object]:
    errors = []
    covered = []
    for node, truth in case.true_ages_years.items():
        estimate = result.age_results[node]
        errors.append(float(estimate["mean_age_years"]) - float(truth))
        covered.append(
            float(estimate["age_95_low"])
            <= float(truth)
            <= float(estimate["age_95_high"])
        )
    diagnostic = result.age_diagnostics
    return {
        "seed": case.seed,
        "mae_years": float(np.mean(np.abs(errors))),
        "bias_years": float(np.mean(errors)),
        "interval95_coverage": float(np.mean(covered)),
        "r_hat_max": diagnostic.get("age_r_hat_max"),
        "ess_bulk_min": diagnostic.get("age_ess_bulk_min"),
        "ess_tail_min": diagnostic.get("age_ess_tail_min"),
        "divergences": diagnostic.get("divergences"),
        "n_identifiable_nodes": diagnostic.get("n_nodes_tracer_identifiable"),
        "converged": diagnostic.get("converged"),
        "sampler": diagnostic.get("sampler"),
    }


def _write_case_artifacts(
    output: Path,
    case: IndependentAquifer,
    result: StrongInferenceResult,
    split: str,
) -> None:
    case_dir = output / "cases" / f"{split}_{case.seed}"
    case_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(case.observations).to_csv(
        case_dir / "blind_observations.csv", index=False
    )
    pd.DataFrame(case.pathline_rows).to_csv(
        case_dir / "modpath_pathline_truth.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "edge_id": f"{u}->{v}",
                "u": u,
                "v": v,
                "true_process": case.true_processes[f"{u}->{v}"],
            }
            for u, v in case.true_edges
        ]
    ).to_csv(case_dir / "heldout_truth.csv", index=False)
    (case_dir / "provenance.json").write_text(
        json.dumps(case.provenance, indent=2, default=_json_default),
        encoding="utf-8",
    )
    (case_dir / "diagnostics.json").write_text(
        json.dumps(
            {
                "age": result.age_diagnostics,
                "phreeqc": result.phreeqc_diagnostics,
                "topology": result.posterior_diagnostics,
                "module_calls": result.module_calls,
            },
            indent=2,
            default=_json_default,
        ),
        encoding="utf-8",
    )


def run_benchmark(
    *,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    dev_seeds: Sequence[int] = DEV_SEEDS,
    test_seeds: Sequence[int] = TEST_SEEDS,
    age_draws: int = 500,
    age_evidence_mode: str = "logistic",
    topology_updates_per_sample: int = 16,
    protocol_stage: str = "initial_locked",
) -> Dict[str, object]:
    if age_evidence_mode not in {"logistic", "direction_gate"}:
        raise ValueError(
            "age_evidence_mode must be 'logistic' or 'direction_gate'."
        )
    age_travel_cost_weight = (
        0.0 if age_evidence_mode == "direction_gate" else 0.1
    )
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    simulator_workspace.mkdir(parents=True, exist_ok=True)

    dev_cases: Dict[int, IndependentAquifer] = {}
    dev_results: Dict[int, StrongInferenceResult] = {}
    dev_rows = []
    for seed in dev_seeds:
        case = generate_independent_aquifer(
            seed,
            simulator_workspace / f"dev_{seed}",
            mf6_executable,
            mp7_executable,
        )
        result = run_strong_inference(
            case.observations,
            seed,
            age_draws=age_draws,
            age_chains=4,
            age_travel_cost_weight=age_travel_cost_weight,
        )
        dev_cases[seed] = case
        dev_results[seed] = result
        dev_rows.extend(_label_rows(case, result, "development"))
        _write_case_artifacts(output, case, result, "development")
    development = pd.DataFrame(dev_rows)

    baseline_model = _fit_fusion_model(development, BASELINE_FEATURES)
    if age_evidence_mode == "direction_gate":
        full_model = {
            **baseline_model,
            "kind": "age_compatibility_gate",
            "age_cost_max": DIRECTION_GATE_AGE_COST_MAX,
            "incompatible_probability": 1.0e-6,
            "age_evidence": (
                "one_sided_downstream_older_probability_with_uncertainty"
            ),
            "development_true_edge_retention": float(
                np.mean(
                    development.loc[
                        development["is_true_edge"] == 1, "age_cost"
                    ]
                    <= DIRECTION_GATE_AGE_COST_MAX
                )
            ),
        }
    else:
        full_model = _fit_fusion_model(development, FUSION_FEATURES)
    development["probability_hydraulic_chemistry"] = _predict_fusion(
        development, baseline_model
    )
    development["probability_hydraulic_chemistry_age"] = _predict_fusion(
        development, full_model
    )
    rng = np.random.default_rng(9991)
    permuted_dev = _age_permuted_frame(
        development, rng, age_evidence_mode
    )
    development["probability_age_permuted_control"] = _predict_fusion(
        permuted_dev, full_model
    )
    truth_dev = development["is_true_edge"].to_numpy(int)
    thresholds = {
        "hydraulic_chemistry": _select_threshold(
            truth_dev,
            development["probability_hydraulic_chemistry"].to_numpy(float),
        ),
        "hydraulic_chemistry_age": _select_threshold(
            truth_dev,
            development[
                "probability_hydraulic_chemistry_age"
            ].to_numpy(float),
        ),
        "age_permuted_control": _select_threshold(
            truth_dev,
            development["probability_age_permuted_control"].to_numpy(float),
        ),
    }

    test_cases: Dict[int, IndependentAquifer] = {}
    test_results: Dict[int, StrongInferenceResult] = {}
    test_rows = []
    topology_rows = []
    age_rows = []
    phreeqc_rows = []
    for seed in test_seeds:
        case = generate_independent_aquifer(
            seed,
            simulator_workspace / f"test_{seed}",
            mf6_executable,
            mp7_executable,
        )
        result = run_strong_inference(
            case.observations,
            seed,
            fusion_model=full_model,
            run_posterior=True,
            age_draws=age_draws,
            age_chains=4,
            topology_updates_per_sample=topology_updates_per_sample,
            age_travel_cost_weight=age_travel_cost_weight,
        )
        test_cases[seed] = case
        test_results[seed] = result
        test_rows.extend(_label_rows(case, result, "locked_test"))
        topology_rows.append(_topology_summary(case, result))
        age_rows.append(_age_summary(case, result))
        phreeqc_rows.append({"seed": seed, **result.phreeqc_diagnostics})
        _write_case_artifacts(output, case, result, "locked_test")
    test = pd.DataFrame(test_rows)
    test["probability_hydraulic_chemistry"] = _predict_fusion(
        test, baseline_model
    )
    test["probability_hydraulic_chemistry_age"] = _predict_fusion(
        test, full_model
    )
    rng = np.random.default_rng(9992)
    permuted_test = _age_permuted_frame(test, rng, age_evidence_mode)
    test["probability_age_permuted_control"] = _predict_fusion(
        permuted_test, full_model
    )

    truth_test = test["is_true_edge"].to_numpy(int)
    methods = {
        "hydraulic_chemistry": (
            "probability_hydraulic_chemistry",
            thresholds["hydraulic_chemistry"],
        ),
        "hydraulic_chemistry_age": (
            "probability_hydraulic_chemistry_age",
            thresholds["hydraulic_chemistry_age"],
        ),
        "age_permuted_control": (
            "probability_age_permuted_control",
            thresholds["age_permuted_control"],
        ),
    }
    method_rows = []
    for name, (column, threshold) in methods.items():
        method_rows.append(
            {
                "method": name,
                **_metrics(
                    truth_test,
                    test[column].to_numpy(float),
                    threshold,
                ),
            }
        )
    method_summary = pd.DataFrame(method_rows)
    bootstrap = _paired_case_bootstrap(test, methods)

    reaction_true = test[test["is_true_edge"] == 1].copy()
    reaction_true["constrained_family_correct"] = (
        reaction_true["dominant_reaction_family"]
        == reaction_true["true_process_family"]
    )
    reaction_true["unconstrained_family_correct"] = (
        reaction_true["dominant_reaction_family_unconstrained"]
        == reaction_true["true_process_family"]
    )
    reaction_summary = (
        reaction_true.groupby("true_process", as_index=False)
        .agg(
            n=("edge_id", "size"),
            constrained_family_accuracy=(
                "constrained_family_correct",
                "mean",
            ),
            unconstrained_family_accuracy=(
                "unconstrained_family_correct",
                "mean",
            ),
            phreeqc_material_change_fraction=(
                "phreeqc_objective_change",
                lambda values: float(np.mean(np.abs(values) > 1e-6)),
            ),
            thermodynamic_bound_hit_fraction=(
                "thermodynamic_bound_hit_count",
                lambda values: float(np.mean(np.asarray(values) > 0)),
            ),
        )
    )
    reaction_overall = pd.DataFrame(
        [
            {
                "true_process": "ALL",
                "n": len(reaction_true),
                "constrained_family_accuracy": float(
                    reaction_true["constrained_family_correct"].mean()
                ),
                "unconstrained_family_accuracy": float(
                    reaction_true["unconstrained_family_correct"].mean()
                ),
                "phreeqc_material_change_fraction": float(
                    np.mean(
                        np.abs(reaction_true["phreeqc_objective_change"]) > 1e-6
                    )
                ),
                "thermodynamic_bound_hit_fraction": float(
                    np.mean(
                        reaction_true["thermodynamic_bound_hit_count"] > 0
                    )
                ),
            }
        ]
    )
    reaction_summary = pd.concat(
        [reaction_overall, reaction_summary], ignore_index=True
    )

    field = run_field_prequential(
        REPO_ROOT
        / "data"
        / "NorthenGhana"
        / "Aquifers_Dataset_Mendeley.xlsx"
    )
    field.predictions.to_csv(
        output / "field_prequential_predictions.csv", index=False
    )
    field.summary.to_csv(output / "field_prequential_summary.csv", index=False)
    (output / "field_prequential_audit.json").write_text(
        json.dumps(field.audit, indent=2, default=_json_default),
        encoding="utf-8",
    )

    development.to_csv(output / "development_edge_features.csv", index=False)
    test.to_csv(output / "locked_test_edge_results.csv", index=False)
    method_summary.to_csv(output / "method_summary.csv", index=False)
    bootstrap.to_csv(output / "age_incremental_bootstrap.csv", index=False)
    pd.DataFrame(topology_rows).to_csv(
        output / "topology_posterior_diagnostics.csv", index=False
    )
    pd.DataFrame(age_rows).to_csv(
        output / "bayesian_age_diagnostics.csv", index=False
    )
    pd.DataFrame(phreeqc_rows).to_csv(
        output / "phreeqc_constraint_audit.csv", index=False
    )
    reaction_true.to_csv(output / "reaction_edge_recovery.csv", index=False)
    reaction_summary.to_csv(output / "reaction_summary.csv", index=False)
    (output / "frozen_fusion_models.json").write_text(
        json.dumps(
            {
                "baseline": baseline_model,
                "full_with_age": full_model,
                "thresholds": thresholds,
            },
            indent=2,
            default=_json_default,
        ),
        encoding="utf-8",
    )

    full_metrics = method_summary.set_index("method").loc[
        "hydraulic_chemistry_age"
    ]
    baseline_metrics = method_summary.set_index("method").loc[
        "hydraulic_chemistry"
    ]
    topology_frame = pd.DataFrame(topology_rows)
    age_frame = pd.DataFrame(age_rows)
    manifest: Dict[str, object] = {
        "benchmark": "M7.2 strong integration",
        "protocol_stage": protocol_stage,
        "age_evidence_mode": age_evidence_mode,
        "age_travel_cost_weight": age_travel_cost_weight,
        "topology_updates_per_sample": topology_updates_per_sample,
        "development_seeds": list(dev_seeds),
        "locked_test_seeds": list(test_seeds),
        "n_development_cases": len(dev_seeds),
        "n_locked_test_cases": len(test_seeds),
        "external_generator": (
            "official MODFLOW 6 + MODPATH 7 + independent nonlinear chemistry"
        ),
        "candidate_recall_test": float(
            np.mean(
                [
                    len(
                        {f"{u}->{v}" for u, v in test_cases[seed].true_edges}
                        & set(test_results[seed].candidate_edges)
                    )
                    / len(test_cases[seed].true_edges)
                    for seed in test_seeds
                ]
            )
        ),
        "age_incremental_f1": float(
            full_metrics["f1"] - baseline_metrics["f1"]
        ),
        "age_incremental_pr_auc": float(
            full_metrics["pr_auc"] - baseline_metrics["pr_auc"]
        ),
        "all_age_cases_converged": bool(age_frame["converged"].all()),
        "all_topology_cases_converged": bool(
            topology_frame["converged"].all()
        ),
        "phreeqc_all_samples_successful": bool(
            all(row["success_fraction"] == 1.0 for row in phreeqc_rows)
        ),
        "phreeqc_materially_active_all_cases": bool(
            all(
                row["n_edges_with_material_objective_change"] > 0
                for row in phreeqc_rows
            )
        ),
        "reaction_truth_families_covered": [
            "carbonate",
            "silicate_exchange",
            "sulfate_reduction",
            "denitrification",
            "iron_reduction",
        ],
        "field_prequential_claim": field.audit["field_claim"],
        "field_n_pairs": field.audit["n_complete_quantitative_pairs"],
        "claim_guardrail": (
            "External simulators provide model-conditioned synthetic truth. "
            "Northern Ghana supplies within-campaign chemistry hold-forward, "
            "not field topology, age, or reaction truth."
        ),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, default=_json_default),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--simulator-workspace",
        type=Path,
        default=REPO_ROOT / ".codex_work" / "m7_2_simulators",
    )
    parser.add_argument(
        "--bin-dir",
        type=Path,
        default=REPO_ROOT / ".codex_work" / "modflow-bin",
    )
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--confirmatory", action="store_true")
    parser.add_argument("--age-draws", type=int, default=500)
    args = parser.parse_args()
    dev = DEV_SEEDS[:2] if args.quick else DEV_SEEDS
    available_test = (
        CONFIRMATORY_TEST_SEEDS if args.confirmatory else TEST_SEEDS
    )
    test = available_test[:3] if args.quick else available_test
    output = args.output or (
        REPO_ROOT
        / "M7"
        / "m7_strong_integration"
        / "results"
        / (
            "m7_2_strong_confirmatory"
            if args.confirmatory
            else "m7_2_strong"
        )
    )
    manifest = run_benchmark(
        output=output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        dev_seeds=dev,
        test_seeds=test,
        age_draws=args.age_draws,
        age_evidence_mode=(
            "direction_gate" if args.confirmatory else "logistic"
        ),
        topology_updates_per_sample=(
            32 if args.confirmatory else 16
        ),
        protocol_stage=(
            "fresh_seed_confirmatory_after_initial_age_contrast_failure"
            if args.confirmatory
            else "initial_locked"
        ),
    )
    print(json.dumps(manifest, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
