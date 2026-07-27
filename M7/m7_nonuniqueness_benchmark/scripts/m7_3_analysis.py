"""Analysis primitives for the M7.3 conditional-integration benchmark."""

from __future__ import annotations

from collections import Counter
import hashlib
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.special import expit, gammaln, logsumexp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)

EVIDENCE_FEATURES: Dict[str, tuple[str, ...]] = {
    "H": ("hydraulic_logit",),
    "A": ("negative_age_cost",),
    "C": ("negative_chemistry_log_objective",),
    "HA": ("hydraulic_logit", "negative_age_cost"),
    "HC": ("hydraulic_logit", "negative_chemistry_log_objective"),
    "AC": ("negative_age_cost", "negative_chemistry_log_objective"),
    "HAC": (
        "hydraulic_logit",
        "negative_age_cost",
        "negative_chemistry_log_objective",
    ),
}

CONDITIONS = (
    "native",
    "age_permuted",
    "hydraulic_permuted",
    "joint_misspecified",
)

CORE_IONS = ("Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3")
ENHANCED_IONS = (
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
    "SiO2",
)


def _bounded_probability(values: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(values, dtype=float), 1.0e-8, 1.0 - 1.0e-8)


def fit_evidence_models(
    development: pd.DataFrame,
) -> Dict[str, Dict[str, object]]:
    """Fit unweighted, development-only logistic evidence models."""

    truth = development["is_true_edge"].to_numpy(int)
    if set(np.unique(truth)) != {0, 1}:
        raise ValueError("Development data must contain true and false edges.")
    models: Dict[str, Dict[str, object]] = {}
    for panel, feature_names in EVIDENCE_FEATURES.items():
        x = development[list(feature_names)].to_numpy(float)
        means = x.mean(axis=0)
        scales = x.std(axis=0)
        scales[scales < 1.0e-8] = 1.0
        standardized = (x - means) / scales
        estimator = LogisticRegression(
            C=1.0,
            class_weight=None,
            max_iter=2000,
            random_state=7300 + len(panel),
            solver="lbfgs",
        )
        estimator.fit(standardized, truth)
        models[panel] = {
            "panel": panel,
            "feature_names": list(feature_names),
            "means": means.tolist(),
            "scales": scales.tolist(),
            "coefficients": estimator.coef_[0].tolist(),
            "intercept": float(estimator.intercept_[0]),
            "training_split": "development_external_simulator_cases_only",
            "class_weight": None,
            "regularization_C": 1.0,
        }
    return models


def predict_evidence_model(
    frame: pd.DataFrame,
    model: Mapping[str, object],
) -> np.ndarray:
    names = list(model["feature_names"])
    means = np.asarray(model["means"], dtype=float)
    scales = np.asarray(model["scales"], dtype=float)
    coefficients = np.asarray(model["coefficients"], dtype=float)
    x = frame[names].to_numpy(float)
    linear = (x - means) / scales @ coefficients + float(model["intercept"])
    return _bounded_probability(expit(linear))


def _permute_within_cases(
    frame: pd.DataFrame,
    column: str,
    *,
    salt: int,
) -> pd.Series:
    output = frame[column].copy()
    for case_seed, indexes in frame.groupby("seed", sort=True).groups.items():
        indexes = np.asarray(list(indexes))
        rng = np.random.default_rng(int(case_seed) * 1009 + int(salt))
        values = frame.loc[indexes, column].to_numpy(copy=True)
        output.loc[indexes] = values[rng.permutation(len(values))]
    return output


def apply_evidence_condition(
    frame: pd.DataFrame,
    condition: str,
) -> pd.DataFrame:
    """Apply a locked within-case negative-control transformation."""

    if condition not in CONDITIONS:
        raise ValueError(f"Unknown condition: {condition}")
    transformed = frame.copy()
    if condition in {"age_permuted", "joint_misspecified"}:
        transformed["negative_age_cost"] = _permute_within_cases(
            transformed,
            "negative_age_cost",
            salt=17,
        )
    if condition in {"hydraulic_permuted", "joint_misspecified"}:
        transformed["hydraulic_logit"] = _permute_within_cases(
            transformed,
            "hydraulic_logit",
            salt=29,
        )
    return transformed


def normalized_binary_entropy(probability: np.ndarray) -> np.ndarray:
    p = _bounded_probability(probability)
    return -(p * np.log(p) + (1.0 - p) * np.log(1.0 - p)) / math.log(2.0)


def expected_calibration_error(
    truth: np.ndarray,
    probability: np.ndarray,
    *,
    n_bins: int = 10,
) -> float:
    truth = np.asarray(truth, dtype=int)
    probability = _bounded_probability(probability)
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    total = max(1, len(truth))
    error = 0.0
    for index in range(int(n_bins)):
        lower = edges[index]
        upper = edges[index + 1]
        selected = (probability >= lower) & (
            probability <= upper if index == n_bins - 1 else probability < upper
        )
        if not np.any(selected):
            continue
        error += (
            float(np.sum(selected))
            / total
            * abs(
                float(np.mean(probability[selected])) - float(np.mean(truth[selected]))
            )
        )
    return float(error)


def evidence_metric_row(
    frame: pd.DataFrame,
    probability: np.ndarray,
) -> Dict[str, float]:
    truth = frame["is_true_edge"].to_numpy(int)
    probability = _bounded_probability(probability)
    prediction = probability >= 0.5
    overconfident = (probability <= 0.1) | (probability >= 0.9)
    wrong = prediction != truth
    expected_by_case = (
        pd.Series(probability, index=frame.index).groupby(frame["seed"]).sum()
    )
    return {
        "pr_auc": float(average_precision_score(truth, probability)),
        "roc_auc": (
            float(roc_auc_score(truth, probability))
            if len(np.unique(truth)) == 2
            else float("nan")
        ),
        "brier": float(brier_score_loss(truth, probability)),
        "log_loss": float(log_loss(truth, probability, labels=[0, 1])),
        "mean_edge_entropy": float(np.mean(normalized_binary_entropy(probability))),
        "ambiguous_edge_fraction": float(
            np.mean((probability > 0.1) & (probability < 0.9))
        ),
        "mean_expected_edges_per_case": float(expected_by_case.mean()),
        "expected_calibration_error": expected_calibration_error(truth, probability),
        "overconfident_error_fraction": float(np.mean(overconfident & wrong)),
    }


def evaluate_evidence_conditions(
    test: pd.DataFrame,
    models: Mapping[str, Mapping[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluate all evidence panels under every locked condition."""

    aggregate_rows = []
    case_rows = []
    conflict_rows = []
    for condition in CONDITIONS:
        conditioned = apply_evidence_condition(test, condition)
        probabilities = {
            panel: predict_evidence_model(conditioned, model)
            for panel, model in models.items()
        }
        for panel, probability in probabilities.items():
            aggregate_rows.append(
                {
                    "condition": condition,
                    "panel": panel,
                    **evidence_metric_row(conditioned, probability),
                }
            )
            for case_seed, indexes in conditioned.groupby(
                "seed", sort=True
            ).groups.items():
                subset = conditioned.loc[indexes]
                subset_probability = pd.Series(
                    probability,
                    index=conditioned.index,
                ).loc[indexes]
                case_rows.append(
                    {
                        "seed": int(case_seed),
                        "condition": condition,
                        "panel": panel,
                        **evidence_metric_row(
                            subset,
                            subset_probability.to_numpy(float),
                        ),
                    }
                )

        univariate = np.column_stack(
            [probabilities["H"], probabilities["A"], probabilities["C"]]
        )
        span = univariate.max(axis=1) - univariate.min(axis=1)
        conflict = span >= 0.50
        integrated = probabilities["HAC"]
        wrong = (integrated >= 0.5) != conditioned["is_true_edge"].to_numpy(int)
        conflict_rows.append(
            {
                "condition": condition,
                "n_edges": len(conditioned),
                "conflict_fraction": float(np.mean(conflict)),
                "mean_univariate_probability_span": float(np.mean(span)),
                "integrated_error_rate_conflict": (
                    float(np.mean(wrong[conflict]))
                    if np.any(conflict)
                    else float("nan")
                ),
                "integrated_error_rate_concordant": (
                    float(np.mean(wrong[~conflict]))
                    if np.any(~conflict)
                    else float("nan")
                ),
                "integrated_overconfident_error_fraction": float(
                    np.mean(wrong & ((integrated <= 0.1) | (integrated >= 0.9)))
                ),
            }
        )
    return (
        pd.DataFrame(aggregate_rows),
        pd.DataFrame(case_rows),
        pd.DataFrame(conflict_rows),
    )


def bootstrap_evidence_contrasts(
    case_metrics: pd.DataFrame,
    *,
    n_bootstrap: int = 10_000,
    random_seed: int = 7331,
) -> pd.DataFrame:
    """Case-block bootstrap for the predeclared evidence contrasts."""

    contrasts = (
        ("native", "HAC", "HC", "native_incremental_age"),
        ("native", "HAC", "HA", "native_incremental_chemistry"),
        ("native", "HAC", "AC", "native_incremental_hydraulics"),
        ("age_permuted", "HAC", "HC", "permuted_age_increment"),
        (
            "hydraulic_permuted",
            "HAC",
            "AC",
            "permuted_hydraulic_increment",
        ),
        ("joint_misspecified", "HAC", "C", "joint_misspecification"),
    )
    metrics = (
        "pr_auc",
        "brier",
        "log_loss",
        "mean_edge_entropy",
        "overconfident_error_fraction",
    )
    rng = np.random.default_rng(int(random_seed))
    rows = []
    for condition, full_panel, baseline_panel, label in contrasts:
        selected = case_metrics[case_metrics["condition"] == condition]
        full = selected[selected["panel"] == full_panel].set_index("seed")
        baseline = selected[selected["panel"] == baseline_panel].set_index("seed")
        common = sorted(set(full.index) & set(baseline.index))
        if not common:
            continue
        for metric in metrics:
            paired = full.loc[common, metric].to_numpy(float) - baseline.loc[
                common, metric
            ].to_numpy(float)
            sampled = rng.choice(
                paired,
                size=(int(n_bootstrap), len(paired)),
                replace=True,
            ).mean(axis=1)
            rows.append(
                {
                    "contrast": label,
                    "condition": condition,
                    "full_panel": full_panel,
                    "baseline_panel": baseline_panel,
                    "metric": metric,
                    "mean_difference": float(np.mean(paired)),
                    "ci95_low": float(np.quantile(sampled, 0.025)),
                    "ci95_high": float(np.quantile(sampled, 0.975)),
                    "n_cases": len(paired),
                    "n_bootstrap": int(n_bootstrap),
                    "resampling_unit": "independent_MODFLOW_case",
                }
            )
    return pd.DataFrame(rows)


def _student_t_logpdf(
    residual: np.ndarray,
    sigma: float,
    degrees_of_freedom: float = 4.0,
) -> np.ndarray:
    nu = float(degrees_of_freedom)
    scale = max(1.0e-12, float(sigma))
    constant = (
        gammaln((nu + 1.0) / 2.0)
        - gammaln(nu / 2.0)
        - 0.5 * math.log(nu * math.pi)
        - math.log(scale)
    )
    return constant - 0.5 * (nu + 1.0) * np.log1p(
        (np.asarray(residual, dtype=float) / scale) ** 2 / nu
    )


def local_age_probability_grid(
    observations: Sequence[Mapping[str, object]],
    *,
    regime: str,
    age_grid: np.ndarray | None = None,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Return exact local age probabilities under the M7 tracer likelihood."""

    if regime not in {"informative", "tritium_only"}:
        raise ValueError("regime must be informative or tritium_only")
    grid = (
        np.asarray(age_grid, dtype=float)
        if age_grid is not None
        else np.linspace(0.25, 200.0, 800)
    )
    positive = np.maximum(grid, 1.0e-9)
    log_sigma = 1.2
    log_prior = (
        -np.log(positive * log_sigma * math.sqrt(2.0 * math.pi))
        - 0.5 * ((np.log(positive) - math.log(20.0)) / log_sigma) ** 2
    )
    tritium_expected = 6.2 * np.exp(-math.log(2.0) * grid / 12.32)
    argon_expected = 97.0 * np.exp(-math.log(2.0) * grid / 269.0)

    nodes = [str(row["site_id"]) for row in observations]
    probabilities = np.empty((len(nodes), len(grid)), dtype=float)
    for node_index, row in enumerate(observations):
        tritium = float(row["tritium_TU"])
        log_posterior = log_prior + _student_t_logpdf(
            tritium - tritium_expected,
            max(0.10, 0.12 * tritium),
        )
        if regime == "informative":
            argon = float(row["argon39_pmc"])
            log_posterior += _student_t_logpdf(
                argon - argon_expected,
                1.8,
            )
        relative = np.exp(log_posterior - np.max(log_posterior))
        probabilities[node_index] = relative / relative.sum()
    return nodes, grid, probabilities


def _weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: Sequence[float],
) -> np.ndarray:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    cumulative /= cumulative[-1]
    return np.interp(np.asarray(quantiles, dtype=float), cumulative, sorted_values)


def topology_age_sensitivity(
    *,
    observations: Sequence[Mapping[str, object]],
    true_ages: Mapping[str, float],
    true_edges: Sequence[tuple[str, str]],
    seed: int,
    regime: str,
    n_particles: int = 50_000,
    order_scale_years: float = 5.0,
) -> pd.DataFrame:
    """Condition local age posteriors on correct, partial and reversed graphs."""

    nodes, age_grid, probabilities = local_age_probability_grid(
        observations,
        regime=regime,
    )
    node_index = {node: index for index, node in enumerate(nodes)}
    rng = np.random.default_rng(
        int(seed) * 7919 + (0 if regime == "informative" else 1)
    )
    sampled_indexes = np.column_stack(
        [
            rng.choice(
                len(age_grid),
                size=int(n_particles),
                replace=True,
                p=probabilities[index],
            )
            for index in range(len(nodes))
        ]
    )
    draws = age_grid[sampled_indexes]
    partial_edges = tuple(
        edge for index, edge in enumerate(true_edges) if index % 2 == 0
    )
    graph_conditions: Dict[str, tuple[tuple[str, str], ...]] = {
        "none": (),
        "partial_true": partial_edges,
        "complete_true": tuple(true_edges),
        "reversed": tuple((v, u) for u, v in true_edges),
    }

    rows = []
    for graph_condition, assumed_edges in graph_conditions.items():
        log_weight = np.zeros(int(n_particles), dtype=float)
        for u, v in assumed_edges:
            delta = (draws[:, node_index[v]] - draws[:, node_index[u]]) / float(
                order_scale_years
            )
            log_weight += -np.logaddexp(0.0, -delta)
        log_mean_weight = float(logsumexp(log_weight) - math.log(len(log_weight)))
        weight = np.exp(log_weight - logsumexp(log_weight))
        ess = float(1.0 / np.sum(weight**2))

        posterior_means = {}
        lower = {}
        upper = {}
        widths = {}
        entropy = {}
        for index, node in enumerate(nodes):
            values = draws[:, index]
            posterior_means[node] = float(np.sum(weight * values))
            quantile = _weighted_quantile(values, weight, (0.025, 0.975))
            lower[node], upper[node] = map(float, quantile)
            widths[node] = float(quantile[1] - quantile[0])
            mass = np.bincount(
                sampled_indexes[:, index],
                weights=weight,
                minlength=len(age_grid),
            )
            positive_mass = mass[mass > 0.0]
            raw_entropy = -float(np.sum(positive_mass * np.log(positive_mass)))
            entropy[node] = raw_entropy / math.log(len(age_grid))

        errors = np.asarray(
            [posterior_means[node] - float(true_ages[node]) for node in nodes],
            dtype=float,
        )
        coverage = np.asarray(
            [lower[node] <= float(true_ages[node]) <= upper[node] for node in nodes],
            dtype=bool,
        )
        truth_order_violation = [
            float(np.sum(weight * (draws[:, node_index[v]] < draws[:, node_index[u]])))
            for u, v in true_edges
        ]
        assumed_order_violation = [
            float(np.sum(weight * (draws[:, node_index[v]] < draws[:, node_index[u]])))
            for u, v in assumed_edges
        ]
        rows.append(
            {
                "seed": int(seed),
                "tracer_regime": regime,
                "graph_condition": graph_condition,
                "n_nodes": len(nodes),
                "n_assumed_edges": len(assumed_edges),
                "n_particles": int(n_particles),
                "order_scale_years": float(order_scale_years),
                "age_mae_years": float(np.mean(np.abs(errors))),
                "age_bias_years": float(np.mean(errors)),
                "age_95_coverage": float(np.mean(coverage)),
                "mean_interval_width_years": float(np.mean(list(widths.values()))),
                "mean_normalized_age_entropy": float(np.mean(list(entropy.values()))),
                "true_edge_order_violation_probability": float(
                    np.mean(truth_order_violation)
                ),
                "assumed_edge_order_violation_probability": (
                    float(np.mean(assumed_order_violation))
                    if assumed_order_violation
                    else float("nan")
                ),
                "importance_ess": ess,
                "importance_ess_fraction": ess / int(n_particles),
                "log_mean_topology_weight": log_mean_weight,
                "importance_stable_ess_ge_400": bool(ess >= 400.0),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_topology_age_contrasts(
    rows: pd.DataFrame,
    *,
    n_bootstrap: int = 10_000,
    random_seed: int = 7341,
) -> pd.DataFrame:
    comparisons = (
        ("complete_true", "none", "correct_minus_no_topology"),
        ("partial_true", "none", "partial_minus_no_topology"),
        ("reversed", "complete_true", "reversed_minus_correct"),
    )
    metrics = (
        "age_mae_years",
        "age_95_coverage",
        "mean_interval_width_years",
        "mean_normalized_age_entropy",
        "true_edge_order_violation_probability",
        "log_mean_topology_weight",
    )
    rng = np.random.default_rng(int(random_seed))
    output = []
    for regime in sorted(rows["tracer_regime"].unique()):
        selected = rows[rows["tracer_regime"] == regime]
        for left_name, right_name, label in comparisons:
            left = selected[selected["graph_condition"] == left_name].set_index("seed")
            right = selected[selected["graph_condition"] == right_name].set_index(
                "seed"
            )
            common = sorted(set(left.index) & set(right.index))
            for metric in metrics:
                paired = left.loc[common, metric].to_numpy(float) - right.loc[
                    common, metric
                ].to_numpy(float)
                sampled = rng.choice(
                    paired,
                    size=(int(n_bootstrap), len(paired)),
                    replace=True,
                ).mean(axis=1)
                output.append(
                    {
                        "tracer_regime": regime,
                        "contrast": label,
                        "left": left_name,
                        "right": right_name,
                        "metric": metric,
                        "mean_difference": float(np.mean(paired)),
                        "ci95_low": float(np.quantile(sampled, 0.025)),
                        "ci95_high": float(np.quantile(sampled, 0.975)),
                        "n_cases": len(paired),
                        "n_bootstrap": int(n_bootstrap),
                    }
                )
    return pd.DataFrame(output)


def reaction_support_summary(
    bootstrap_predictions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize bootstrap family probabilities without hiding carbonate."""

    required = {"seed", "edge_id", "tier", "true_family", "predicted_family"}
    missing = required - set(bootstrap_predictions)
    if missing:
        raise ValueError(f"Missing reaction bootstrap columns: {sorted(missing)}")
    family_rows = []
    edge_rows = []
    group_columns = ["seed", "edge_id", "tier", "true_family", "true_process"]
    for keys, group in bootstrap_predictions.groupby(group_columns, sort=True):
        seed, edge_id, tier, true_family, true_process = keys
        counts = Counter(group["predicted_family"].astype(str))
        total = sum(counts.values())
        probabilities = {
            family: count / total for family, count in sorted(counts.items())
        }
        positive = np.asarray(list(probabilities.values()), dtype=float)
        entropy = -float(np.sum(positive * np.log(positive)))
        modal_family = max(probabilities, key=probabilities.get)
        edge_rows.append(
            {
                "seed": int(seed),
                "edge_id": edge_id,
                "tier": tier,
                "true_family": true_family,
                "true_process": true_process,
                "n_bootstrap": total,
                "modal_family": modal_family,
                "modal_family_correct": modal_family == true_family,
                "true_family_probability": float(
                    probabilities.get(str(true_family), 0.0)
                ),
                "family_support_entropy": entropy,
                "effective_supported_families": float(math.exp(entropy)),
            }
        )
        for family, probability in probabilities.items():
            family_rows.append(
                {
                    "seed": int(seed),
                    "edge_id": edge_id,
                    "tier": tier,
                    "true_family": true_family,
                    "true_process": true_process,
                    "predicted_family": family,
                    "support_probability": probability,
                    "n_bootstrap": total,
                }
            )
    edges = pd.DataFrame(edge_rows)
    summary = edges.groupby(["tier", "true_process"], as_index=False).agg(
        n_edges=("edge_id", "size"),
        modal_family_accuracy=("modal_family_correct", "mean"),
        mean_true_family_probability=("true_family_probability", "mean"),
        mean_family_support_entropy=("family_support_entropy", "mean"),
        mean_effective_supported_families=(
            "effective_supported_families",
            "mean",
        ),
    )
    overall = (
        edges.groupby("tier", as_index=False)
        .agg(
            n_edges=("edge_id", "size"),
            modal_family_accuracy=("modal_family_correct", "mean"),
            mean_true_family_probability=("true_family_probability", "mean"),
            mean_family_support_entropy=("family_support_entropy", "mean"),
            mean_effective_supported_families=(
                "effective_supported_families",
                "mean",
            ),
        )
        .assign(true_process="ALL")
    )
    summary = pd.concat([overall, summary], ignore_index=True)
    return pd.DataFrame(family_rows), edges, summary


def audit_ghana_workbook(workbook_path: Path) -> Dict[str, object]:
    """Create a machine-readable scope audit for the canonical Ghana workbook
    (data/FieldData/NorthenGhana/NorthernGhana.xlsx, Dry/Wet sheets).

    An earlier revision audited a different, antecedent study's own derived
    workbook (Wells_Nodes/Hydrochemistry_Seasonal/Graph_Edges/
    Coordinate_Masking_Note sheets, including a fabricated per-record
    Sampling_Date field with no equivalent in the real raw data); that
    workbook has been removed (DECISIONS.md). This audit reads the two
    canonical Dry/Wet sheets directly and reports what they do and do not
    contain, including that no graph edges or per-record sampling dates
    exist in the canonical source.
    """

    workbook_path = Path(workbook_path)
    dry = pd.read_excel(workbook_path, sheet_name="Dry").assign(Season="Dry")
    wet = pd.read_excel(workbook_path, sheet_name="Wet").assign(Season="Wet")
    chemistry = pd.concat([dry, wet], ignore_index=True)
    all_columns = {str(column).strip().lower() for column in chemistry.columns}
    age_tokens = ("tritium", "3h", "cfc", "sf6", "argon39", "39ar", "14c")
    age_columns = sorted(
        column
        for column in all_columns
        if any(token in column.replace("_", "") for token in age_tokens)
    )
    screen_columns = sorted(
        column
        for column in all_columns
        if "screen" in column and ("top" in column or "bottom" in column)
    )
    # Static_Water_Level_m is a per-well attribute duplicated across both
    # season rows, not two independent measurement occasions; count it once
    # per well from the Dry sheet alone.
    per_well_static_count = (
        dry.groupby("Well_ID")["Static_Water_Level_m"].count()
        if "Static_Water_Level_m" in dry
        else pd.Series(dtype=float)
    )
    sha256 = hashlib.sha256(workbook_path.read_bytes()).hexdigest()
    return {
        "workbook": str(workbook_path.resolve()),
        "workbook_sha256": sha256,
        "n_wells": int(chemistry["Well_ID"].nunique()),
        "n_hydrochemistry_rows": int(len(chemistry)),
        "n_seasons": int(chemistry["Season"].nunique()),
        "sampling_date_field_available": False,
        "sampling_granularity": (
            "one dry-season and one wet-season observation per well; no "
            "intra-season sampling-date field exists in the canonical source"
        ),
        "major_chemistry_available": all(
            f"{ion.lower()}_mg_l" in all_columns
            for ion in ("Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3")
        ),
        "stable_water_isotopes_available": {
            "d18O": "d18o_permil" in all_columns,
            "d2H": "d2h_permil" in all_columns,
            "interpretation": "recharge/source evidence, not residence-time truth",
        },
        "age_tracer_columns": age_columns,
        "environmental_age_tracer_panel_available": bool(age_columns),
        "screen_interval_columns": screen_columns,
        "screen_intervals_available": bool(screen_columns),
        "total_borehole_depth_available": "borehole_depth_m" in all_columns,
        "elevation_available": "elevation_m" in all_columns,
        "single_static_water_level_available": ("static_water_level_m" in all_columns),
        "static_measurements_per_well_max": (
            int(per_well_static_count.max()) if not per_well_static_count.empty else 0
        ),
        "time_varying_head_series_available": False,
        "single_occasion_head_proxy_possible": (
            "elevation_m" in all_columns and "static_water_level_m" in all_columns
        ),
        "coordinates_masked": True,
        "masking_distance_statement": "approximately 1-3 km (SI.pdf)",
        "processed_graph_edges_available": False,
        "independent_field_connectivity_truth_available": False,
        "independent_field_reaction_truth_available": False,
        "independent_aquifer_type_classification_available": False,
        "supportable_interpretations": [
            "data readiness and measurement-value audit",
            "within-campaign seasonal chemistry hold-forward under a "
            "disclosed arbitrary well-revelation order",
            "reaction-family plausibility and equivalence classes",
            "sensitivity to alternative assumed edge sets",
            "null explanations and non-identifiability classification",
        ],
        "non_identifiable_or_unsupported": [
            "field residence-time inversion",
            "exact directed field flow paths",
            "screen-resolved vertical connectivity",
            "unique field reaction mechanisms",
            "fully observed field digital twin",
            "true chronological sampling sequence within a season",
        ],
        "objective_6": (
            "Apply the framework and its component diagnostics to Ghanaian "
            "aquifer datasets to determine which integrated interpretations "
            "are supportable under available data and which remain "
            "non-identifiable."
        ),
    }
