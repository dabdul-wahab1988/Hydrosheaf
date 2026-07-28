"""Prospective closed-loop benchmark for Bayesian active learning."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss
from sklearn.linear_model import LogisticRegression

PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
M7_SCRIPTS = REPO / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
for path in (REPO, M7_SCRIPTS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hydrosheaf.calibration.bayesian_active_learning import (  # noqa: E402
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    rank_measurement_options,
    shannon_entropy,
    update_hypothesis_posterior,
)
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402
from independent_modflow_generator import (  # noqa: E402
    ION_ORDER,
    generate_independent_aquifer,
)
from strong_inference import strong_config  # noqa: E402


RUN_ID = "RUN-M8-FRONTIER-AL-20260728-01"
DEVELOPMENT_RUN_ID = "DEV-M8-FRONTIER-AL-20260728-01"
CALIBRATION_SEEDS = tuple(range(7401, 7409))
TUNING_SEEDS = tuple(range(7451, 7459))
LOCKED_TEST_SEEDS = tuple(range(7601, 7625))
MEASUREMENT_TYPES = ("chemistry_panel", "age_tracer", "connectivity_tracer")
COSTS = {"chemistry_panel": 2.0, "age_tracer": 5.0, "connectivity_tracer": 9.0}
STRATEGIES = (
    "robust_information_decision_per_cost",
    "mean_information_decision_per_cost",
    "legacy_uncertainty_chemistry",
    "random_feasible",
    "realised_oracle",
)
N_PARTICLES = 256
BUDGET = 10.0
MAX_ACTIONS = 5
BOOTSTRAP_SAMPLES = 5000
ROBUSTNESS_WEIGHT = 0.75
DECISION_WEIGHT = 0.95
MINIMUM_EIG = 0.002
LIKELIHOOD_COEFFICIENT_SHRINKAGE = 0.25
LEGACY_BRIER_NONINFERIORITY_MARGIN = 0.01
LEGACY_INFORMATION_NONINFERIORITY_MARGIN = -0.01
FEATURE_NAMES = ("intercept", "distance_km", "abs_head_drop_m", "hydraulic_prior")
STRESS = {
    "separation_multiplier": 0.65,
    "separation_sd_multiplier": 1.25,
    "noise_sd_multiplier": 1.60,
}


@dataclass(frozen=True)
class CaseRecord:
    seed: int
    edge_ids: tuple[str, ...]
    sources: tuple[str, ...]
    targets: tuple[str, ...]
    truth: np.ndarray
    hydraulic_prior: np.ndarray
    features: np.ndarray
    outcomes: Mapping[str, np.ndarray]
    candidate_recall: float
    provenance: Mapping[str, Any]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def _safe_output(output: Path, overwrite: bool) -> None:
    resolved = output.resolve()
    if not resolved.is_relative_to(PROJECT.resolve()):
        raise ValueError("Output must stay within the M8 benchmark directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {resolved}")
        if resolved == PROJECT.resolve():
            raise ValueError("Refusing to remove the M8 benchmark root.")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)


def _hydraulic_config():
    config = strong_config(phreeqc_enabled=False)
    config.edge_max_neighbors = 3
    config.sheaf_age_enabled = False
    config.validate()
    return config


def _distance_km(left: Mapping[str, object], right: Mapping[str, object]) -> float:
    dx = float(right["x_m"]) - float(left["x_m"])
    dy = float(right["y_m"]) - float(left["y_m"])
    return math.hypot(dx, dy) / 1000.0


def _hidden_outcomes(
    seed: int,
    edge_ids: Sequence[str],
    sources: Sequence[str],
    targets: Sequence[str],
    truth: np.ndarray,
    sample_map: Mapping[str, Mapping[str, object]],
    true_ages: Mapping[str, float],
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(int(seed) + 9_100_003)
    chemistry = []
    ages = []
    connectivity = []
    for edge_index, (edge_id, source, target) in enumerate(
        zip(edge_ids, sources, targets)
    ):
        left = sample_map[source]
        right = sample_map[target]
        left_chemistry = np.asarray([float(left[ion]) for ion in ION_ORDER])
        right_chemistry = np.asarray([float(right[ion]) for ion in ION_ORDER])
        chemistry.append(
            float(
                np.linalg.norm(np.log1p(right_chemistry) - np.log1p(left_chemistry))
                + rng.normal(0.0, 0.012)
            )
        )
        distance = _distance_km(left, right)
        age_noise = 0.25 + 0.12 * distance
        ages.append(
            float(
                true_ages[target]
                - true_ages[source]
                + rng.normal(0.0, age_noise)
            )
        )
        present = bool(truth[edge_index])
        base = 0.88 if present else 0.12
        # Case- and edge-specific discrepancy prevents a perfect direct oracle.
        bias = 0.035 * math.sin(0.17 * seed + 0.31 * edge_index)
        connectivity.append(float(rng.normal(base + bias, 0.11 + 0.02 * distance)))
        if not edge_id:
            raise ValueError("Empty edge identifier in hidden outcome generator.")
    return {
        "chemistry_panel": np.asarray(chemistry, dtype=float),
        "age_tracer": np.asarray(ages, dtype=float),
        "connectivity_tracer": np.asarray(connectivity, dtype=float),
    }


def generate_case(
    seed: int,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
) -> CaseRecord:
    case = generate_independent_aquifer(
        int(seed),
        simulator_workspace / f"case_{int(seed)}",
        mf6_executable,
        mp7_executable,
    )
    rows = [dict(row) for row in case.observations]
    forbidden = {"true_age", "true_edges", "true_process"}
    leaked = forbidden & set().union(*(set(row) for row in rows))
    if leaked:
        raise ValueError(f"Truth leaked into acquisition inputs: {sorted(leaked)}")
    candidates = infer_edges(rows, method="probabilistic", config=_hydraulic_config())
    edge_ids = tuple(edge.edge_id for edge in candidates)
    sources = tuple(edge.u for edge in candidates)
    targets = tuple(edge.v for edge in candidates)
    truth_ids = {f"{source}->{target}" for source, target in case.true_edges}
    truth = np.asarray([edge_id in truth_ids for edge_id in edge_ids], dtype=int)
    candidate_set = set(edge_ids)
    recall = float(len(candidate_set & truth_ids) / max(1, len(truth_ids)))
    sample_map = {str(row["site_id"]): row for row in rows}

    hydraulic_prior = []
    feature_rows = []
    for edge in candidates:
        attrs = edge.attrs or {}
        probability = float(
            attrs.get(
                "prior_edge_probability",
                attrs.get("p_uv", attrs.get("edge_confidence", 0.5)),
            )
        )
        probability = float(np.clip(probability, 0.02, 0.98))
        left = sample_map[edge.u]
        right = sample_map[edge.v]
        distance = _distance_km(left, right)
        head_drop = abs(float(left["hydraulic_head"]) - float(right["hydraulic_head"]))
        hydraulic_prior.append(probability)
        feature_rows.append((1.0, distance, head_drop, probability))

    outcomes = _hidden_outcomes(
        int(seed),
        edge_ids,
        sources,
        targets,
        truth,
        sample_map,
        case.true_ages_years,
    )
    return CaseRecord(
        seed=int(seed),
        edge_ids=edge_ids,
        sources=sources,
        targets=targets,
        truth=truth,
        hydraulic_prior=np.asarray(hydraulic_prior, dtype=float),
        features=np.asarray(feature_rows, dtype=float),
        outcomes=outcomes,
        candidate_recall=recall,
        provenance=dict(case.provenance),
    )


def fit_measurement_model(cases: Sequence[CaseRecord]) -> dict[str, Any]:
    if not cases:
        raise ValueError("At least one development case is required.")
    features = np.vstack([case.features for case in cases])
    truth = np.concatenate([case.truth for case in cases]).astype(int)
    model: dict[str, Any] = {
        "schema_version": "1.0",
        "development_run_id": DEVELOPMENT_RUN_ID,
        "development_seeds": [case.seed for case in cases],
        "feature_names": list(FEATURE_NAMES),
        "ridge_penalty": 1.0e-3,
        "stress": dict(STRESS),
        "measurement_types": {},
    }
    prior_model = LogisticRegression(
        C=1.0,
        class_weight="balanced",
        max_iter=2000,
        random_state=20260728,
    )
    prior_model.fit(features[:, 1:], truth)
    model["topology_prior_calibration"] = {
        "feature_names": list(FEATURE_NAMES[1:]),
        "intercept": float(prior_model.intercept_[0]),
        "coefficients": prior_model.coef_[0].tolist(),
        "class_weight": "balanced",
        "regularisation_C": 1.0,
    }
    penalty = np.eye(features.shape[1]) * 1.0e-3
    penalty[0, 0] = 0.0
    for measurement_type in MEASUREMENT_TYPES:
        outcomes = np.concatenate(
            [case.outcomes[measurement_type] for case in cases]
        )
        states: dict[str, Any] = {}
        for state in (0, 1):
            mask = truth == state
            if int(mask.sum()) < features.shape[1] + 2:
                raise ValueError(
                    f"Insufficient development examples for {measurement_type} "
                    f"state {state}."
                )
            design = features[mask]
            response = outcomes[mask]
            raw_coefficients = np.linalg.solve(
                design.T @ design + penalty,
                design.T @ response,
            )
            pooled_coefficients = np.zeros_like(raw_coefficients)
            pooled_coefficients[0] = float(np.mean(response))
            coefficients = pooled_coefficients + LIKELIHOOD_COEFFICIENT_SHRINKAGE * (
                raw_coefficients - pooled_coefficients
            )
            residuals = response - design @ coefficients
            floor = max(1.0e-3, 0.05 * float(np.std(response, ddof=1)))
            residual_sd = max(floor, float(np.sqrt(np.mean(residuals**2))))
            states[str(state)] = {
                "coefficients": coefficients.tolist(),
                "unshrunk_coefficients": raw_coefficients.tolist(),
                "coefficient_shrinkage": LIKELIHOOD_COEFFICIENT_SHRINKAGE,
                "residual_sd": residual_sd,
                "n": int(mask.sum()),
                "response_mean": float(np.mean(response)),
                "response_sd": float(np.std(response, ddof=1)),
            }
        model["measurement_types"][measurement_type] = states
    return model


def _calibrated_edge_prior(
    case: CaseRecord,
    measurement_model: Mapping[str, Any] | None,
) -> np.ndarray:
    if measurement_model is None or "topology_prior_calibration" not in measurement_model:
        return np.asarray(case.hydraulic_prior, dtype=float)
    calibration = measurement_model["topology_prior_calibration"]
    coefficients = np.asarray(calibration["coefficients"], dtype=float)
    linear = float(calibration["intercept"]) + case.features[:, 1:] @ coefficients
    return np.clip(1.0 / (1.0 + np.exp(-linear)), 0.01, 0.99)


def _topology_particles(
    case: CaseRecord,
    measurement_model: Mapping[str, Any] | None = None,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    rng = np.random.default_rng(case.seed + 4_200_013)
    particles = np.zeros((N_PARTICLES, len(case.edge_ids)), dtype=int)
    by_source: dict[str, list[int]] = {}
    for edge_index, source in enumerate(case.sources):
        by_source.setdefault(source, []).append(edge_index)
    calibrated_prior = _calibrated_edge_prior(case, measurement_model)
    for indices in by_source.values():
        probabilities_for_edges = np.asarray(
            [calibrated_prior[index] for index in indices], dtype=float
        )
        # Convert binary calibrated probabilities to odds and include a unit
        # no-edge reference state. This produces a coherent source-level
        # categorical prior rather than renormalising raw probabilities.
        edge_odds = probabilities_for_edges / (1.0 - probabilities_for_edges)
        probabilities = np.append(edge_odds, 1.0)
        probabilities /= probabilities.sum()
        draws = rng.choice(len(indices) + 1, size=N_PARTICLES, p=probabilities)
        for particle_index, draw in enumerate(draws):
            if draw < len(indices):
                particles[particle_index, indices[int(draw)]] = 1
    hypothesis_ids = [f"G{index:04d}" for index in range(N_PARTICLES)]
    weights = np.full(N_PARTICLES, 1.0 / N_PARTICLES)
    return hypothesis_ids, particles, weights


def _class_predictions(
    case: CaseRecord,
    measurement_model: Mapping[str, Any],
    measurement_type: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    states = measurement_model["measurement_types"][measurement_type]
    coefficients_0 = np.asarray(states["0"]["coefficients"], dtype=float)
    coefficients_1 = np.asarray(states["1"]["coefficients"], dtype=float)
    mean_0 = case.features @ coefficients_0
    mean_1 = case.features @ coefficients_1
    sd_0 = np.full(len(case.edge_ids), float(states["0"]["residual_sd"]))
    sd_1 = np.full(len(case.edge_ids), float(states["1"]["residual_sd"]))
    return mean_0, mean_1, sd_0, sd_1


def build_measurement_options(
    case: CaseRecord,
    particles: np.ndarray,
    measurement_model: Mapping[str, Any],
) -> tuple[list[MeasurementOption], dict[str, tuple[str, int]]]:
    options = []
    lookup: dict[str, tuple[str, int]] = {}
    for measurement_type in MEASUREMENT_TYPES:
        mean_0, mean_1, sd_0, sd_1 = _class_predictions(
            case, measurement_model, measurement_type
        )
        for edge_index, edge_id in enumerate(case.edge_ids):
            states = particles[:, edge_index].astype(bool)
            nominal_means = np.where(states, mean_1[edge_index], mean_0[edge_index])
            nominal_sds = np.where(states, sd_1[edge_index], sd_0[edge_index])
            midpoint = 0.5 * (mean_0[edge_index] + mean_1[edge_index])
            stressed_0 = midpoint + STRESS["separation_multiplier"] * (
                mean_0[edge_index] - midpoint
            )
            stressed_1 = midpoint + STRESS["separation_multiplier"] * (
                mean_1[edge_index] - midpoint
            )
            separation_means = np.where(states, stressed_1, stressed_0)
            separation_sds = nominal_sds * STRESS["separation_sd_multiplier"]
            noise_sds = nominal_sds * STRESS["noise_sd_multiplier"]
            option_id = f"{measurement_type}:{edge_id}"
            options.append(
                MeasurementOption(
                    option_id=option_id,
                    measurement_type=measurement_type,
                    target_id=edge_id,
                    cost=COSTS[measurement_type],
                    scenarios=(
                        PredictiveScenario(
                            "nominal", nominal_means, nominal_sds, weight=0.50
                        ),
                        PredictiveScenario(
                            "separation_stress",
                            separation_means,
                            separation_sds,
                            weight=0.25,
                        ),
                        PredictiveScenario(
                            "noise_stress", nominal_means, noise_sds, weight=0.25
                        ),
                    ),
                    metadata={"edge_index": edge_index},
                )
            )
            lookup[option_id] = (measurement_type, edge_index)
    return options, lookup


def _edge_probabilities(particles: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.asarray(weights @ particles, dtype=float)


def _marginal_edge_entropy(probabilities: np.ndarray) -> float:
    clipped = np.clip(probabilities, 1.0e-12, 1.0 - 1.0e-12)
    return float(np.sum(-clipped * np.log(clipped) - (1.0 - clipped) * np.log(1.0 - clipped)))


def _choose_legacy(
    case: CaseRecord,
    particles: np.ndarray,
    weights: np.ndarray,
    available: Sequence[MeasurementOption],
) -> MeasurementOption | None:
    chemistry = {
        option.target_id: option
        for option in available
        if option.measurement_type == "chemistry_panel"
    }
    probabilities = _edge_probabilities(particles, weights)
    rows = []
    for edge_index, edge_id in enumerate(case.edge_ids):
        if edge_id not in chemistry:
            continue
        probability = float(np.clip(probabilities[edge_index], 1.0e-12, 1 - 1.0e-12))
        entropy = -probability * math.log(probability) - (1 - probability) * math.log(
            1 - probability
        )
        rows.append((-entropy, edge_id, chemistry[edge_id]))
    return min(rows, key=lambda row: (row[0], row[1]))[2] if rows else None


def _choose_oracle(
    case: CaseRecord,
    hypothesis_ids: Sequence[str],
    particles: np.ndarray,
    weights: np.ndarray,
    available: Sequence[MeasurementOption],
    lookup: Mapping[str, tuple[str, int]],
) -> MeasurementOption | None:
    current = brier_score_loss(case.truth, _edge_probabilities(particles, weights))
    candidates = []
    for option in available:
        measurement_type, edge_index = lookup[option.option_id]
        update = update_hypothesis_posterior(
            hypothesis_ids,
            weights,
            option,
            float(case.outcomes[measurement_type][edge_index]),
        )
        posterior = np.asarray(update["posterior_probabilities"], dtype=float)
        brier = brier_score_loss(case.truth, _edge_probabilities(particles, posterior))
        improvement_per_cost = (current - brier) / float(option.cost)
        candidates.append((-improvement_per_cost, option.option_id, option))
    if not candidates:
        return None
    best = min(candidates, key=lambda row: (row[0], row[1]))
    return best[2] if -best[0] > 0.0 else None


def run_strategy(
    case: CaseRecord,
    measurement_model: Mapping[str, Any],
    strategy: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    hypothesis_ids, particles, weights = _topology_particles(case, measurement_model)
    options, lookup = build_measurement_options(case, particles, measurement_model)
    initial_joint_entropy = shannon_entropy(weights)
    initial_edge_entropy = _marginal_edge_entropy(
        _edge_probabilities(particles, weights)
    )
    available = list(options)
    cost_spent = 0.0
    action_rows = []
    rng = np.random.default_rng(case.seed + 8_800_101)

    for round_index in range(MAX_ACTIONS):
        affordable = [
            option
            for option in available
            if cost_spent + float(option.cost) <= BUDGET + 1.0e-12
        ]
        if not affordable:
            break
        option: MeasurementOption | None
        acquisition_status = "ACTIONABLE"
        if strategy in {
            "robust_information_decision_per_cost",
            "mean_information_decision_per_cost",
        }:
            robustness = ROBUSTNESS_WEIGHT if strategy.startswith("robust") else 0.0
            decision = rank_measurement_options(
                hypothesis_ids,
                weights,
                affordable,
                decision_values=particles,
                config=AcquisitionConfig(
                    quadrature_order=15,
                    robustness_weight=robustness,
                    decision_weight=DECISION_WEIGHT,
                    cost_exponent=1.0,
                    minimum_expected_information_gain=MINIMUM_EIG,
                ),
            )
            acquisition_status = str(decision["status"])
            selected_id = decision["selected_option_id"]
            option = (
                next(item for item in affordable if item.option_id == selected_id)
                if selected_id is not None
                else None
            )
        elif strategy == "legacy_uncertainty_chemistry":
            option = _choose_legacy(case, particles, weights, affordable)
        elif strategy == "random_feasible":
            option = affordable[int(rng.integers(0, len(affordable)))]
        elif strategy == "realised_oracle":
            option = _choose_oracle(
                case,
                hypothesis_ids,
                particles,
                weights,
                affordable,
                lookup,
            )
            if option is None:
                acquisition_status = "ABSTAIN"
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
        if option is None:
            break

        measurement_type, edge_index = lookup[option.option_id]
        observed = float(case.outcomes[measurement_type][edge_index])
        before_entropy = shannon_entropy(weights)
        update = update_hypothesis_posterior(
            hypothesis_ids,
            weights,
            option,
            observed,
        )
        weights = np.asarray(update["posterior_probabilities"], dtype=float)
        cost_spent += float(option.cost)
        available = [item for item in available if item.option_id != option.option_id]
        action_rows.append(
            {
                "seed": case.seed,
                "strategy": strategy,
                "round": round_index + 1,
                "option_id": option.option_id,
                "measurement_type": measurement_type,
                "target_edge_id": case.edge_ids[edge_index],
                "target_is_true": int(case.truth[edge_index]),
                "observed_value": observed,
                "cost": float(option.cost),
                "cumulative_cost": cost_spent,
                "entropy_before": before_entropy,
                "entropy_after": shannon_entropy(weights),
                "realised_information_gain": update["realised_information_gain"],
                "acquisition_status": acquisition_status,
            }
        )

    edge_probabilities = np.clip(
        _edge_probabilities(particles, weights), 1.0e-9, 1.0 - 1.0e-9
    )
    final_joint_entropy = shannon_entropy(weights)
    final_edge_entropy = _marginal_edge_entropy(edge_probabilities)
    joint_reduction = initial_joint_entropy - final_joint_entropy
    row = {
        "seed": case.seed,
        "strategy": strategy,
        "candidate_recall": case.candidate_recall,
        "n_candidates": len(case.edge_ids),
        "n_true_candidates": int(case.truth.sum()),
        "n_actions": len(action_rows),
        "actionable": int(bool(action_rows)),
        "cost_spent": cost_spent,
        "brier": float(brier_score_loss(case.truth, edge_probabilities)),
        "pr_auc": float(average_precision_score(case.truth, edge_probabilities)),
        "log_loss": float(log_loss(case.truth, edge_probabilities, labels=[0, 1])),
        "initial_joint_entropy": initial_joint_entropy,
        "final_joint_entropy": final_joint_entropy,
        "joint_entropy_reduction": joint_reduction,
        "joint_entropy_reduction_per_cost": (
            joint_reduction / cost_spent if cost_spent > 0.0 else 0.0
        ),
        "initial_marginal_edge_entropy": initial_edge_entropy,
        "final_marginal_edge_entropy": final_edge_entropy,
        "measurement_types": ";".join(
            action["measurement_type"] for action in action_rows
        ),
    }
    return row, action_rows


def _bootstrap_contrast(
    metrics: pd.DataFrame,
    comparator: str,
    metric: str,
    *,
    samples: int,
    seed: int,
) -> dict[str, Any]:
    pivot = metrics.pivot(index="seed", columns="strategy", values=metric).dropna()
    differences = (
        pivot["robust_information_decision_per_cost"] - pivot[comparator]
    ).to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(int(samples), dtype=float)
    for index in range(int(samples)):
        bootstrap[index] = float(
            np.median(rng.choice(differences, size=len(differences), replace=True))
        )
    return {
        "strategy": "robust_information_decision_per_cost",
        "comparator": comparator,
        "metric": metric,
        "n_cases": len(differences),
        "median_paired_difference": float(np.median(differences)),
        "ci95_low": float(np.quantile(bootstrap, 0.025)),
        "ci95_high": float(np.quantile(bootstrap, 0.975)),
    }


def _stable_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column in result.select_dtypes(include=["float"]).columns:
        result[column] = result[column].map(
            lambda value: float(f"{value:.12g}") if math.isfinite(value) else value
        )
    return result


def evaluate_cases(
    cases: Sequence[CaseRecord],
    measurement_model: Mapping[str, Any],
    *,
    bootstrap_samples: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    action_rows = []
    for case in cases:
        for strategy in STRATEGIES:
            metric, actions = run_strategy(case, measurement_model, strategy)
            metric_rows.append(metric)
            action_rows.extend(actions)
    metrics = _stable_frame(pd.DataFrame(metric_rows))
    actions = _stable_frame(pd.DataFrame(action_rows))
    summary = _stable_frame(
        metrics.groupby("strategy", sort=True)
        .agg(
            n_cases=("seed", "count"),
            actionability_rate=("actionable", "mean"),
            median_brier=("brier", "median"),
            median_pr_auc=("pr_auc", "median"),
            median_log_loss=("log_loss", "median"),
            median_entropy_reduction_per_cost=(
                "joint_entropy_reduction_per_cost",
                "median",
            ),
            median_cost=("cost_spent", "median"),
            median_actions=("n_actions", "median"),
        )
        .reset_index()
    )
    comparators = (
        "random_feasible",
        "legacy_uncertainty_chemistry",
        "mean_information_decision_per_cost",
        "realised_oracle",
    )
    contrast_metrics = (
        "brier",
        "pr_auc",
        "joint_entropy_reduction_per_cost",
    )
    contrast_rows = [
        _bootstrap_contrast(
            metrics,
            comparator,
            metric,
            samples=bootstrap_samples,
            seed=2026072900 + 100 * comparator_index + metric_index,
        )
        for comparator_index, comparator in enumerate(comparators)
        for metric_index, metric in enumerate(contrast_metrics)
    ]
    contrasts = _stable_frame(pd.DataFrame(contrast_rows))
    return metrics, actions, summary, contrasts


def _gate_decision(metrics: pd.DataFrame, contrasts: pd.DataFrame) -> dict[str, Any]:
    proposed = metrics[
        metrics["strategy"] == "robust_information_decision_per_cost"
    ]
    gates: dict[str, bool] = {
        "candidate_recall": float(proposed["candidate_recall"].mean()) >= 0.80,
        "actionability": float(proposed["actionable"].mean()) >= 0.90,
    }
    for comparator in ("random_feasible", "legacy_uncertainty_chemistry"):
        rows = contrasts[contrasts["comparator"] == comparator].set_index("metric")
        brier_limit = (
            LEGACY_BRIER_NONINFERIORITY_MARGIN
            if comparator == "legacy_uncertainty_chemistry"
            else 0.0
        )
        gate_name = (
            f"brier_noninferiority_vs_{comparator}"
            if comparator == "legacy_uncertainty_chemistry"
            else f"brier_superiority_vs_{comparator}"
        )
        gates[gate_name] = float(rows.loc["brier", "ci95_high"]) < brier_limit
        information_limit = (
            LEGACY_INFORMATION_NONINFERIORITY_MARGIN
            if comparator == "legacy_uncertainty_chemistry"
            else 0.0
        )
        information_gate_name = (
            f"information_efficiency_noninferiority_vs_{comparator}"
            if comparator == "legacy_uncertainty_chemistry"
            else f"information_efficiency_superiority_vs_{comparator}"
        )
        gates[information_gate_name] = (
            float(
                rows.loc[
                    "joint_entropy_reduction_per_cost",
                    "ci95_low",
                ]
            )
            > information_limit
        )
        gates[f"pr_auc_nonharm_vs_{comparator}"] = (
            float(rows.loc["pr_auc", "ci95_low"]) > -0.01
        )
    supported = bool(all(gates.values()))
    return {
        "run_id": RUN_ID,
        "frontier_active_learning_claim_supported": supported,
        "gates": gates,
        "allowed_claim": (
            "On untouched, code-independent controlled-synthetic aquifer cases, "
            "scenario-robust Bayesian acquisition improved Brier score and joint "
            "topology information per declared relative cost over random acquisition, "
            "while remaining noninferior to the strong legacy uncertainty policy."
            if supported
            else "The locked evidence did not satisfy every preregistered gate for "
            "a topology active-learning performance claim."
        ),
        "guardrail": (
            "Controlled-synthetic joint-topology evidence under a declared relative-"
            "cost and likelihood model; not prospective field validation."
        ),
    }


def _verify_lock(calibration_model: Path) -> dict[str, Any]:
    lock_path = PROJECT / "m8_frontier_active_learning_protocol.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("run_id") != RUN_ID or not lock.get("locked_before_test_execution"):
        raise RuntimeError("Frontier active-learning protocol is not validly locked.")
    if lock.get("locked_test_seeds") != list(LOCKED_TEST_SEEDS):
        raise RuntimeError("Locked-test seed family differs from the protocol lock.")
    expected_model = lock.get("calibration_model")
    if expected_model != calibration_model.relative_to(PROJECT).as_posix():
        raise RuntimeError("Calibration-model path differs from the protocol lock.")
    for relative, expected_hash in lock.get("sha256", {}).items():
        path = REPO / relative
        if not path.exists() or _sha256(path) != expected_hash:
            raise RuntimeError(f"Locked source mismatch: {relative}")
    return lock


def _write_outputs(
    output: Path,
    metrics: pd.DataFrame,
    actions: pd.DataFrame,
    summary: pd.DataFrame,
    contrasts: pd.DataFrame,
    provenance: Mapping[str, Any],
    claim_decision: Mapping[str, Any] | None,
    run_id: str,
    phase: str,
) -> dict[str, Any]:
    metrics.to_csv(output / "case_metrics.csv", index=False, lineterminator="\n")
    actions.to_csv(output / "action_log.csv", index=False, lineterminator="\n")
    summary.to_csv(output / "strategy_summary.csv", index=False, lineterminator="\n")
    contrasts.to_csv(output / "paired_bootstrap_contrasts.csv", index=False, lineterminator="\n")
    (output / "generator_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if claim_decision is not None:
        (output / "claim_decision.json").write_text(
            json.dumps(claim_decision, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    artifact_names = [
        "case_metrics.csv",
        "action_log.csv",
        "strategy_summary.csv",
        "paired_bootstrap_contrasts.csv",
        "generator_provenance.json",
    ]
    if claim_decision is not None:
        artifact_names.append("claim_decision.json")
    manifest = {
        "run_id": run_id,
        "phase": phase,
        "status": "PASS",
        "git_revision": _git_revision(),
        "n_cases": int(metrics["seed"].nunique()),
        "n_particles": N_PARTICLES,
        "budget": BUDGET,
        "max_actions": MAX_ACTIONS,
        "bootstrap_samples": BOOTSTRAP_SAMPLES,
        "artifacts": {name: _sha256(output / name) for name in artifact_names},
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def run_development(
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    *,
    overwrite: bool,
) -> dict[str, Any]:
    _safe_output(output, overwrite)
    calibration_cases = [
        generate_case(
            seed,
            simulator_workspace / "calibration",
            mf6_executable,
            mp7_executable,
        )
        for seed in CALIBRATION_SEEDS
    ]
    measurement_model = fit_measurement_model(calibration_cases)
    model_path = output / "development_calibration_model.json"
    model_path.write_text(
        json.dumps(measurement_model, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    tuning_cases = [
        generate_case(
            seed,
            simulator_workspace / "tuning",
            mf6_executable,
            mp7_executable,
        )
        for seed in TUNING_SEEDS
    ]
    metrics, actions, summary, contrasts = evaluate_cases(
        tuning_cases,
        measurement_model,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
    )
    provenance = {str(case.seed): dict(case.provenance) for case in tuning_cases}
    manifest = _write_outputs(
        output,
        metrics,
        actions,
        summary,
        contrasts,
        provenance,
        None,
        DEVELOPMENT_RUN_ID,
        "development",
    )
    manifest["calibration_model"] = {
        "path": model_path.name,
        "sha256": _sha256(model_path),
        "seeds": list(CALIBRATION_SEEDS),
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def run_locked_test(
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    calibration_model_path: Path,
    *,
    overwrite: bool,
) -> dict[str, Any]:
    calibration_model_path = calibration_model_path.resolve()
    if not calibration_model_path.is_relative_to(PROJECT.resolve()):
        raise ValueError("Calibration model must stay within the M8 benchmark.")
    _verify_lock(calibration_model_path)
    _safe_output(output, overwrite)
    measurement_model = json.loads(calibration_model_path.read_text(encoding="utf-8"))
    cases = [
        generate_case(
            seed,
            simulator_workspace / "locked_test",
            mf6_executable,
            mp7_executable,
        )
        for seed in LOCKED_TEST_SEEDS
    ]
    metrics, actions, summary, contrasts = evaluate_cases(
        cases,
        measurement_model,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
    )
    decision = _gate_decision(metrics, contrasts)
    provenance = {str(case.seed): dict(case.provenance) for case in cases}
    return _write_outputs(
        output,
        metrics,
        actions,
        summary,
        contrasts,
        provenance,
        decision,
        RUN_ID,
        "locked_test",
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("development", "locked-test"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--simulator-workspace", type=Path, required=True)
    parser.add_argument(
        "--bin-dir",
        type=Path,
        default=REPO / ".codex_work" / "modflow-bin",
    )
    parser.add_argument("--calibration-model", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.phase == "development":
        result = run_development(
            args.output,
            args.simulator_workspace,
            args.bin_dir / "mf6.exe",
            args.bin_dir / "mp7.exe",
            overwrite=args.overwrite,
        )
    else:
        if args.calibration_model is None:
            raise ValueError("--calibration-model is required for locked-test phase.")
        result = run_locked_test(
            args.output,
            args.simulator_workspace,
            args.bin_dir / "mf6.exe",
            args.bin_dir / "mp7.exe",
            args.calibration_model,
            overwrite=args.overwrite,
        )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
