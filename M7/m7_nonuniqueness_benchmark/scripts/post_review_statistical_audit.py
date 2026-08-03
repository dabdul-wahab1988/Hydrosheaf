"""Post-review statistical audits for the immutable benchmark result bundles.

This script never regenerates a locked test case. It verifies the stored inputs,
recomputes case-level summaries, applies a full-family max-z bootstrap correction
to every representation-benchmark and estimator-diagnostic contrast, and
produces transparent post-review precision
and public-pipeline selection audits.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from m7_3_analysis import evaluate_evidence_conditions  # noqa: E402
from run_m7_sheaf_vs_graph import (  # noqa: E402
    _case_metric_rows as m74_case_metric_rows,
    _fit_models as m74_fit_models,
    _predict_models as m74_predict_models,
)


OUTPUT = PROJECT / "results" / "post_review_audit_20260730"
TABLES = PROJECT / "tables" / "publication"
M73 = PROJECT / "results" / "m7_3_locked"
M74 = PROJECT / "results" / "RUN-M7-SHEAF-VS-GRAPH-20260729-01"
M75_DEV = (
    PROJECT
    / "results"
    / "RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01"
    / "development"
)
M75 = (
    PROJECT
    / "results"
    / "RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01"
    / "locked_test"
)
SYSTEM = PROJECT / "results" / "RUN-M7-SYSTEM-20260728-01"
BOOTSTRAP_SAMPLES = 10_000
POWER_SIMULATIONS = 20_000

HIGHER_IS_BETTER = {
    "pr_auc",
    "roc_auc",
    "selected_f1",
    "abstention_accuracy",
    "conflict_localisation_pr_auc",
    "age_95_coverage",
    "modal_family_accuracy",
}
LOWER_IS_BETTER = {
    "brier",
    "log_loss",
    "ece",
    "false_confidence_rate",
    "age_mae_years",
    "mean_interval_width_years",
}

# These margins were selected after review for replication planning. They are
# not field decision thresholds and are never described as preregistered.
PLANNING_MARGINS = {
    "pr_auc": 0.02,
    "brier": 0.01,
    "log_loss": 0.02,
    "age_mae_years": 0.25,
    "mean_interval_width_years": 0.50,
    "age_95_coverage": 0.05,
    "modal_family_accuracy": 0.10,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_locked_bundle(run: Path) -> dict[str, str]:
    manifest = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
    verified: dict[str, str] = {}
    for relative, expected in manifest["artifacts"].items():
        path = run / relative
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"Locked artifact hash mismatch: {path}")
        verified[path.relative_to(PROJECT).as_posix()] = actual
    return verified


def case_metrics_from_predictions(
    predictions: pd.DataFrame,
    *,
    group_columns: Sequence[str],
    model_column: str = "model",
) -> pd.DataFrame:
    rows: list[dict] = []
    grouping = [*group_columns, model_column]
    for keys, group in predictions.groupby(grouping, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        record = dict(zip(grouping, keys))
        truth = group["is_true_edge"].to_numpy(int)
        probability = group["probability"].to_numpy(float)
        record.update(
            {
                "pr_auc": float(average_precision_score(truth, probability)),
                "brier": float(brier_score_loss(truth, probability)),
                "log_loss": float(log_loss(truth, probability, labels=[0, 1])),
            }
        )
        rows.append(record)
    return pd.DataFrame(rows)


def development_metrics() -> dict[str, pd.DataFrame]:
    development73 = pd.read_csv(M73 / "development_edge_features.csv")
    models73 = json.loads((M73 / "frozen_evidence_models.json").read_text(encoding="utf-8"))
    _, metrics73, _ = evaluate_evidence_conditions(development73, models73)

    development74 = pd.read_csv(M74 / "development_edge_features.csv")
    fitted74 = m74_fit_models(development74)
    metrics74 = m74_case_metric_rows(m74_predict_models(development74, fitted74))

    oof75 = pd.read_csv(M75_DEV / "development_oof_predictions.csv")
    metrics75 = case_metrics_from_predictions(
        oof75[oof75["calibration_regime"] == "separate_crossfit"],
        group_columns=("seed", "scenario", "calibration_regime"),
    )
    return {"M7.3": metrics73, "M7.4": metrics74, "M7.5": metrics75}


def _difference_vector(
    metrics: pd.DataFrame,
    row: Mapping[str, object],
    *,
    model_column: str,
    regime_column: str | None,
) -> pd.Series:
    subset = metrics
    if regime_column is not None:
        subset = subset[
            subset[regime_column].astype(str) == str(row[regime_column])
        ]
    scenario = str(row["scenario"])
    if scenario != "all":
        subset = subset[subset["scenario"].astype(str) == scenario]
    pivot = subset.pivot(
        index="seed", columns=model_column, values=str(row["metric"])
    ).dropna()
    left, right = str(row["left"]), str(row["right"])
    if left not in pivot or right not in pivot:
        return pd.Series(dtype=float)
    return (pivot[left] - pivot[right]).sort_index()


def simultaneous_family(
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
    *,
    family_name: str,
    rng_seed: int,
    regime_column: str | None = None,
) -> tuple[pd.DataFrame, float]:
    rng = np.random.default_rng(rng_seed)
    scenario_by_seed = (
        metrics[["seed", "scenario"]]
        .drop_duplicates()
        .set_index("seed")["scenario"]
        .astype(str)
        .to_dict()
    )
    index_cache: dict[tuple[int, ...], np.ndarray] = {}

    def indices(seeds: tuple[int, ...]) -> np.ndarray:
        if seeds not in index_cache:
            index_cache[seeds] = rng.integers(
                0, len(seeds), size=(BOOTSTRAP_SAMPLES, len(seeds))
            )
        return index_cache[seeds]

    records: list[dict] = []
    max_z = np.zeros(BOOTSTRAP_SAMPLES, dtype=float)
    for source_row in contrasts.to_dict("records"):
        row = dict(source_row)
        differences = _difference_vector(
            metrics,
            row,
            model_column="model",
            regime_column=regime_column,
        )
        record = dict(row)
        if differences.empty:
            record.update(
                {
                    "standard_error": np.nan,
                    "simultaneous_ci95_low": np.nan,
                    "simultaneous_ci95_high": np.nan,
                    "simultaneous_excludes_zero": False,
                    "supports_benefit_fwer95": False,
                    "family": family_name,
                    "status": "NOT_ESTIMABLE",
                }
            )
            records.append(record)
            continue
        seed_values = tuple(map(int, differences.index))
        by_scenario: dict[str, list[int]] = {}
        for position, seed in enumerate(seed_values):
            by_scenario.setdefault(scenario_by_seed[seed], []).append(position)
        boot_total = np.zeros(BOOTSTRAP_SAMPLES, dtype=float)
        total_n = 0
        for positions in by_scenario.values():
            local = differences.iloc[positions].to_numpy(float)
            local_seeds = tuple(seed_values[position] for position in positions)
            boot_total += local[indices(local_seeds)].mean(axis=1) * len(local)
            total_n += len(local)
        boot_mean = boot_total / total_n
        observed = float(differences.mean())
        se = float(differences.std(ddof=1) / np.sqrt(len(differences)))
        if np.isfinite(se) and se > 0:
            max_z = np.maximum(max_z, np.abs((boot_mean - observed) / se))
        record.update(
            {
                "standard_error": se,
                "_observed": observed,
                "_boot_mean": boot_mean,
                "family": family_name,
                "status": "ESTIMABLE",
            }
        )
        records.append(record)

    critical = float(np.quantile(max_z, 0.95))
    output: list[dict] = []
    for record in records:
        boot = record.pop("_boot_mean", None)
        observed = record.pop("_observed", np.nan)
        se = float(record["standard_error"])
        if boot is None:
            output.append(record)
            continue
        if not np.isfinite(se) or se == 0:
            low = high = float(observed)
        else:
            low = float(observed - critical * se)
            high = float(observed + critical * se)
        metric = str(record["metric"])
        if metric in HIGHER_IS_BETTER:
            benefit = low > 0
        elif metric in LOWER_IS_BETTER:
            benefit = high < 0
        else:
            benefit = False
        record.update(
            {
                "simultaneous_ci95_low": low,
                "simultaneous_ci95_high": high,
                "simultaneous_excludes_zero": bool(low > 0 or high < 0),
                "supports_benefit_fwer95": bool(benefit),
            }
        )
        output.append(record)
    frame = pd.DataFrame(output)
    ordered = [column for column in contrasts.columns if column in frame.columns]
    ordered += [
        "standard_error",
        "simultaneous_ci95_low",
        "simultaneous_ci95_high",
        "simultaneous_excludes_zero",
        "supports_benefit_fwer95",
        "family",
        "status",
    ]
    return frame[ordered], critical


def _paired(
    metrics: pd.DataFrame,
    *,
    left: str,
    right: str,
    metric: str,
    condition_column: str,
) -> np.ndarray:
    pivot = metrics.pivot(index="seed", columns=condition_column, values=metric).dropna()
    return (pivot[left] - pivot[right]).to_numpy(float)


def _power_record(
    *,
    design: str,
    contrast: str,
    metric: str,
    differences: np.ndarray,
    planned_n: int,
    source: str,
    rng: np.random.Generator,
) -> dict:
    direction = 1.0 if metric in HIGHER_IS_BETTER else -1.0
    margin = float(PLANNING_MARGINS[metric])
    oriented = direction * np.asarray(differences, dtype=float)
    draws = rng.choice(
        oriented,
        size=(POWER_SIMULATIONS, int(planned_n)),
        replace=True,
    )
    means = draws.mean(axis=1)
    standard_errors = draws.std(axis=1, ddof=1) / np.sqrt(int(planned_n))
    lower = means - 1.96 * standard_errors
    return {
        "design": design,
        "contrast": contrast,
        "metric": metric,
        "development_or_planning_source": source,
        "source_n_cases": len(oriented),
        "planned_n_cases": int(planned_n),
        "planning_margin_benefit_scale": margin,
        "source_mean_benefit": float(oriented.mean()),
        "source_sd_benefit": float(oriented.std(ddof=1)),
        "probability_ci_excludes_zero_in_favourable_direction": float(
            np.mean(lower > 0)
        ),
        "probability_ci_clears_planning_margin": float(np.mean(lower > margin)),
        "simulations": POWER_SIMULATIONS,
        "interpretation_status": (
            "DEVELOPMENT_ONLY_EMPIRICAL_PLANNING"
            if source.startswith("development")
            else "POST_TEST_REPLICATION_PLANNING_ONLY"
        ),
    }


def precision_audit(development: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rng = np.random.default_rng(2026073003)
    rows: list[dict] = []
    m73 = development["M7.3"]
    for metric in ("pr_auc", "brier", "log_loss"):
        rows.append(
            _power_record(
                design="Process-based evidence-panel locked test",
                contrast="native HAC minus native HA",
                metric=metric,
                differences=_paired(
                    m73[m73["condition"] == "native"],
                    left="HAC",
                    right="HA",
                    metric=metric,
                    condition_column="panel",
                ),
                planned_n=12,
                source=(
                    "development cases evaluated with development-fitted models; "
                    "resubstitution may be optimistic"
                ),
                rng=rng,
            )
        )
    m74 = development["M7.4"]
    for metric in ("pr_auc", "brier", "log_loss"):
        rows.append(
            _power_record(
                design="Locked sheaf-versus-weighted-graph representation test",
                contrast="affine sheaf minus edge-local graph",
                metric=metric,
                differences=_paired(
                    m74,
                    left="affine_sheaf",
                    right="weighted_graph",
                    metric=metric,
                    condition_column="model",
                ),
                planned_n=64,
                source=(
                    "development cases evaluated with development-fitted models; "
                    "resubstitution may be optimistic"
                ),
                rng=rng,
            )
        )
    m75 = development["M7.5"]
    for metric in ("pr_auc", "brier", "log_loss"):
        rows.append(
            _power_record(
                design="Locked local-first/global-fallback estimator diagnostic",
                contrast="local-first/global-fallback minus edge-local",
                metric=metric,
                differences=_paired(
                    m75,
                    left="robust_hybrid",
                    right="edge_local",
                    metric=metric,
                    condition_column="model",
                ),
                planned_n=128,
                source="development eight-fold out-of-fold predictions",
                rng=rng,
            )
        )

    # M7.3 did not archive development topology/reaction metrics. The locked
    # case vectors are used only to plan a future independent replication.
    topology = pd.read_csv(M73 / "topology_age_sensitivity.csv")
    for regime in ("informative", "tritium_only"):
        subset = topology[topology["tracer_regime"] == regime]
        for metric in ("age_mae_years", "mean_interval_width_years"):
            rows.append(
                _power_record(
                    design=f"Topology-conditioned-age test ({regime})",
                    contrast="complete true topology minus no topology",
                    metric=metric,
                    differences=_paired(
                        subset,
                        left="complete_true",
                        right="none",
                        metric=metric,
                        condition_column="graph_condition",
                    ),
                    planned_n=12,
                    source=(
                        "locked-test case distribution because development topology "
                        "metrics were not archived; future replication planning only"
                    ),
                    rng=rng,
                )
            )
    reactions = pd.read_csv(M73 / "reaction_edge_nonuniqueness.csv")
    per_case = (
        reactions.groupby(["seed", "tier"], as_index=False)["modal_family_correct"]
        .mean()
        .rename(columns={"modal_family_correct": "modal_family_accuracy"})
    )
    rows.append(
        _power_record(
            design="Reaction-family recovery test",
            contrast="enhanced minus core indicator panel",
            metric="modal_family_accuracy",
            differences=_paired(
                per_case,
                left="enhanced",
                right="core",
                metric="modal_family_accuracy",
                condition_column="tier",
            ),
            planned_n=12,
            source=(
                "locked-test case distribution because development reaction-bootstrap "
                "metrics were not archived; future replication planning only"
            ),
            rng=rng,
        )
    )
    return pd.DataFrame(rows)


def practical_effects() -> pd.DataFrame:
    rows: list[dict] = []

    def add(
        design: str,
        left: str,
        right: str,
        metric: str,
        metrics: pd.DataFrame,
        filters: Mapping[str, str] | None = None,
    ) -> None:
        subset = metrics
        if filters:
            for column, value in filters.items():
                subset = subset[subset[column].astype(str) == value]
        pivot = subset.pivot(index="seed", columns="model", values=metric).dropna()
        difference = float((pivot[left] - pivot[right]).mean())
        baseline = float(pivot[right].mean())
        direction = 1.0 if metric in HIGHER_IS_BETTER else -1.0
        margin = PLANNING_MARGINS.get(metric, np.nan)
        rows.append(
            {
                "design": design,
                "contrast": f"{left} minus {right}",
                "metric": metric,
                "left_mean": float(pivot[left].mean()),
                "comparator_mean": baseline,
                "mean_difference": difference,
                "benefit_magnitude": direction * difference,
                "relative_effect_percent_of_comparator": (
                    100.0 * direction * difference / abs(baseline)
                    if baseline != 0
                    else np.nan
                ),
                "post_review_planning_margin": margin,
                "mean_clears_planning_margin": bool(
                    np.isfinite(margin) and direction * difference >= margin
                ),
                "margin_status": "POST_REVIEW_NOT_PRESPECIFIED",
            }
        )

    metrics74 = pd.read_csv(M74 / "case_metrics.csv")
    for metric in ("pr_auc", "brier", "log_loss"):
        add("M7.4", "affine_sheaf", "weighted_graph", metric, metrics74)
    metrics75 = pd.read_csv(M75 / "case_metrics.csv")
    for metric in ("pr_auc", "brier", "log_loss"):
        add(
            "M7.5",
            "robust_hybrid",
            "edge_local",
            metric,
            metrics75,
            {"calibration_regime": "separate_crossfit"},
        )
    topology = pd.read_csv(M73 / "topology_age_sensitivity.csv").rename(
        columns={"graph_condition": "model"}
    )
    for regime in ("informative", "tritium_only"):
        for metric in ("age_mae_years", "mean_interval_width_years"):
            add(
                f"Topology-conditioned age ({regime})",
                "complete_true",
                "none",
                metric,
                topology,
                {"tracer_regime": regime},
            )
    return pd.DataFrame(rows)


def public_pipeline_audit() -> pd.DataFrame:
    predictions = pd.read_csv(SYSTEM / "edge_predictions.csv")
    metrics = pd.read_csv(SYSTEM / "case_metrics.csv")
    rows: list[dict] = []
    total_truth = int(metrics.groupby("seed")["n_true_edges"].first().sum())
    for condition, group in predictions.groupby("condition", sort=True):
        truth = group["is_true_edge"].to_numpy(int)
        selected = group["selected"].to_numpy(int)
        true_positive = int(((truth == 1) & (selected == 1)).sum())
        false_positive = int(((truth == 0) & (selected == 1)).sum())
        conditional_false_negative = int(((truth == 1) & (selected == 0)).sum())
        missed_before_scoring = total_truth - int((truth == 1).sum())
        end_to_end_false_negative = conditional_false_negative + missed_before_scoring
        rows.append(
            {
                "condition": condition,
                "selection_rule": "pipeline-retained edge; no scalar probability threshold",
                "n_cases": int(group["seed"].nunique()),
                "candidate_edges": len(group),
                "selected_edges": int(selected.sum()),
                "minimum_selected_probability": float(
                    group.loc[group["selected"] == 1, "probability"].min()
                ),
                "maximum_unselected_probability": (
                    float(group.loc[group["selected"] == 0, "probability"].max())
                    if (group["selected"] == 0).any()
                    else np.nan
                ),
                "conditional_true_positive": true_positive,
                "conditional_false_positive": false_positive,
                "conditional_false_negative": conditional_false_negative,
                "truth_edges_missed_by_candidate_generation": missed_before_scoring,
                "end_to_end_false_negative": end_to_end_false_negative,
                "conditional_selected_f1": float(
                    2 * true_positive
                    / max(2 * true_positive + false_positive + conditional_false_negative, 1)
                ),
                "end_to_end_selected_f1": float(
                    2 * true_positive
                    / max(2 * true_positive + false_positive + end_to_end_false_negative, 1)
                ),
                "candidate_recall": float((truth == 1).sum() / total_truth),
            }
        )
    return pd.DataFrame(rows)


def write_markdown(
    m74: pd.DataFrame,
    critical74: float,
    m75: pd.DataFrame,
    critical75: float,
    precision: pd.DataFrame,
    effects: pd.DataFrame,
    public: pd.DataFrame,
) -> None:
    unadjusted74 = int(((m74["ci95_low"] > 0) | (m74["ci95_high"] < 0)).sum())
    unadjusted75 = int(((m75["ci95_low"] > 0) | (m75["ci95_high"] < 0)).sum())
    simultaneous74 = int(m74["simultaneous_excludes_zero"].fillna(False).sum())
    simultaneous75 = int(m75["simultaneous_excludes_zero"].fillna(False).sum())
    selected74 = m74[
        (m74["scenario"].isin(["incompatible_cycles", "noisy_missing", "all"]))
        & (m74["left"] == "affine_sheaf")
        & (m74["right"].isin(["weighted_graph", "graph_laplacian", "affine_sheaf_permuted"]))
        & (m74["metric"].isin(["pr_auc", "brier", "log_loss", "conflict_localisation_pr_auc"]))
    ]
    selected75 = m75[
        (m75["calibration_regime"] == "separate_crossfit")
        & (m75["left"] == "robust_hybrid")
        & (m75["right"].isin(["edge_local", "robust_hybrid_permuted"]))
        & (m75["metric"].isin(["pr_auc", "brier", "log_loss", "conflict_localisation_pr_auc"]))
    ]
    text = f"""# M7 post-review statistical audit

## Status and scope

This analysis was specified after the 29 July 2026 manuscript review. It reads
the immutable stored case tables and does not regenerate or tune any locked
test. The planning margins are transparent replication-planning values, not
field decision thresholds and not retrospectively labelled as prespecified.

## Full-family multiplicity correction

The family contains every row of the published representation-benchmark
contrast matrix (120 contrasts) or estimator-diagnostic matrix (560 contrasts).
Ten thousand stratified case-block
bootstrap resamples were shared across contrasts. A max-z critical value was
used to construct two-sided simultaneous 95% intervals, controlling the
family-wise error rate across the complete matrix.

- Representation-benchmark max-z critical value: {critical74:.4f}; unadjusted intervals excluding
  zero: {unadjusted74}; simultaneous intervals excluding zero: {simultaneous74}.
- Estimator-diagnostic max-z critical value: {critical75:.4f}; unadjusted intervals excluding
  zero: {unadjusted75}; simultaneous intervals excluding zero: {simultaneous75}.

The manuscript may treat a scenario statement as confirmatory only where
supports_benefit_fwer95 is true. Other scenario results are labelled
exploratory diagnostics and are excluded from the abstract and conclusion.

### Claim-facing representation-benchmark rows

{selected74.round(5).to_markdown(index=False)}

### Claim-facing estimator-diagnostic rows

{selected75.round(5).to_markdown(index=False)}

## Precision and practical magnitude

The original protocols did not contain prospective minimum important
differences. The following audit therefore reports post-review planning margins,
development-based empirical operating characteristics where development
metrics exist, and locked-test distributions only for planning a future
independent replication where development topology/reaction metrics were never
archived.

{precision.round(4).to_markdown(index=False)}

Relative effects and planning-margin checks are reported separately:

{effects.round(4).to_markdown(index=False)}

## Public-pipeline selection audit

All 33 candidates per case were retained in every arm. Hence selected F1 was
identical across arms even though probability ranking and calibration changed.
There was no scalar post-processing probability threshold; selection was the
pipeline's retained-edge set. One of 54 planted truth edges was absent before
scoring, so the scoring audit was conditional on candidate generation.

{public.round(6).to_markdown(index=False)}
"""
    (OUTPUT / "post_review_statistical_audit.md").write_text(text, encoding="utf-8")


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    input_hashes = {}
    input_hashes.update(verify_locked_bundle(M74))
    input_hashes.update(verify_locked_bundle(M75))
    input_hashes.update(verify_locked_bundle(SYSTEM))

    metrics74 = pd.read_csv(M74 / "case_metrics.csv")
    contrasts74 = pd.read_csv(M74 / "paired_bootstrap_contrasts.csv")
    adjusted74, critical74 = simultaneous_family(
        metrics74,
        contrasts74,
        family_name="Representation benchmark: all 120 published contrasts",
        rng_seed=2026073001,
    )
    metrics75 = pd.read_csv(M75 / "case_metrics.csv")
    contrasts75 = pd.read_csv(M75 / "paired_bootstrap_contrasts.csv")
    adjusted75, critical75 = simultaneous_family(
        metrics75,
        contrasts75,
        family_name="Estimator diagnostic: all 560 published contrasts",
        rng_seed=2026073002,
        regime_column="calibration_regime",
    )
    development = development_metrics()
    precision = precision_audit(development)
    effects = practical_effects()
    public = public_pipeline_audit()

    artifacts = {
        "m7_4_full_family_simultaneous.csv": adjusted74,
        "m7_5_full_family_simultaneous.csv": adjusted75,
        "precision_and_power_audit.csv": precision,
        "practical_effects.csv": effects,
        "public_pipeline_selection_audit.csv": public,
    }
    for name, frame in artifacts.items():
        frame.to_csv(OUTPUT / name, index=False, lineterminator="\n")
    adjusted74.to_csv(
        TABLES / "tableS10_m7_4_multiplicity_adjusted.csv",
        index=False,
        lineterminator="\n",
    )
    adjusted75.to_csv(
        TABLES / "tableS11_m7_5_multiplicity_adjusted.csv",
        index=False,
        lineterminator="\n",
    )
    precision.to_csv(
        TABLES / "tableS12_precision_and_power.csv", index=False, lineterminator="\n"
    )
    public.to_csv(
        TABLES / "tableS13_public_pipeline_selection.csv",
        index=False,
        lineterminator="\n",
    )
    write_markdown(
        adjusted74,
        critical74,
        adjusted75,
        critical75,
        precision,
        effects,
        public,
    )
    output_hashes = {
        path.name: sha256(path)
        for path in sorted(OUTPUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }
    manifest = {
        "schema_version": "1.0",
        "analysis": "M7 post-review multiplicity, precision and selection audit",
        "status": "PASS",
        "specified_after_review": True,
        "does_not_regenerate_locked_tests": True,
        "bootstrap_samples": BOOTSTRAP_SAMPLES,
        "power_simulations": POWER_SIMULATIONS,
        "m7_4_family_size": len(adjusted74),
        "m7_5_family_size": len(adjusted75),
        "m7_4_max_z_critical": critical74,
        "m7_5_max_z_critical": critical75,
        "planning_margins": PLANNING_MARGINS,
        "input_hashes": input_hashes,
        "output_hashes": output_hashes,
    }
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
