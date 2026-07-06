"""Generate all M5 main-text and supplementary tables from locked results."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from m5_common import ION_ORDER, REACTIONS, reaction_matrix  # noqa: E402


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
METHOD_LABELS = {
    "bounded_ls": "Bounded LS",
    "lasso": "Lasso",
    "elastic_net": "Elastic Net",
    "thermo_elastic_net": "Thermo Elastic Net",
    "hydrosheaf_guarded": "Hydrosheaf Guarded",
    "hydrosheaf_core": "Hydrosheaf-Core",
}


def mean_ci(series: pd.Series) -> str:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return ""
    mean = values.mean()
    if len(values) == 1:
        return f"{mean:.3f}"
    half = 1.96 * values.std(ddof=1) / np.sqrt(len(values))
    return f"{mean:.3f} [{mean - half:.3f}, {mean + half:.3f}]"


def write(frame: pd.DataFrame, name: str) -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    frame.to_csv(TABLES_DIR / name, index=False)


def main_table() -> None:
    fits = pd.read_csv(RESULTS_DIR / "benchmark_fits.csv")
    subset = fits[(fits["panel"] == "full_11") & (fits["noise_level"] == 0.03)]
    rows = []
    for method, group in subset.groupby("method", sort=False):
        rows.append(
            {
                "Method": METHOD_LABELS.get(method, method.replace("_", " ").title()),
                "Phase F1, mean [95% CI]": mean_ci(group["phase_f1"]),
                "Equivalence-class F1, mean [95% CI]": mean_ci(
                    group["class_f1"]
                ),
                "Extent RMSE (mmol/L), mean [95% CI]": mean_ci(
                    group["extent_rmse_mmolL"]
                ),
                "Reconstruction RMSE (mmol/L), mean [95% CI]": mean_ci(
                    group["reconstruction_rmse_mmolL"]
                ),
                "False-discovery rate, mean [95% CI]": mean_ci(
                    group["false_discovery_rate"]
                ),
                "Thermodynamic violations, mean": group[
                    "thermodynamic_violations"
                ].mean(),
                "Convergence (%)": 100.0 * group["converged"].mean(),
                "Runtime per fit (ms), median": group["runtime_ms"].median(),
                "n": len(group),
            }
        )
    baseline_path = RESULTS_DIR / "phreeqc_inverse_baseline.csv"
    if baseline_path.exists():
        baseline = pd.read_csv(baseline_path)
        rows.append(
            {
                "Method": "Conventional PHREEQC inverse (strict+fallback; first minimal)",
                "Phase F1, mean [95% CI]": mean_ci(
                    baseline["first_minimal_phase_f1"]
                ),
                "Equivalence-class F1, mean [95% CI]": mean_ci(
                    baseline["first_minimal_class_f1"]
                ),
                "Extent RMSE (mmol/L), mean [95% CI]": "",
                "Reconstruction RMSE (mmol/L), mean [95% CI]": "",
                "False-discovery rate, mean [95% CI]": mean_ci(
                    baseline["first_minimal_phase_false_discovery_rate"]
                ),
                "Thermodynamic violations, mean": "",
                "Convergence (%)": 100.0 * baseline[
                    "phreeqc_inverse_success"
                ].mean(),
                "Runtime per fit (ms), median": "",
                "n": len(baseline),
            }
        )
    write(pd.DataFrame(rows), "table1_comparative_inverse_performance.csv")


def supplementary_tables() -> None:
    matrix = reaction_matrix()
    stoich_rows = []
    for reaction, values in zip(REACTIONS, matrix):
        row = {
            "Reaction": reaction.label,
            "Family": reaction.family,
            "Signed extent": reaction.signed,
            "PHREEQC SI phase": reaction.si_phase or "",
        }
        row.update({ion: value for ion, value in zip(ION_ORDER, values)})
        stoich_rows.append(row)
    write(pd.DataFrame(stoich_rows), "tableS1_reaction_stoichiometry.csv")

    classes = pd.read_csv(RESULTS_DIR / "equivalence_classes.csv")
    write(classes, "tableS2_reaction_equivalence_classes.csv")

    truth = pd.read_csv(RESULTS_DIR / "phreeqc_ground_truth.csv")
    scenario_columns = [
        "scenario_id",
        "archetype",
        "replicate",
        "n_true_reactions",
        "true_support",
        "transport_error_level",
        "topology_confidence",
        "residence_time_confidence",
        "generation_rmse_mmolL",
        "upstream_pH",
        "downstream_pH",
    ]
    write(truth[scenario_columns], "tableS3_phreeqc_scenario_parameters.csv")

    hyper = pd.read_csv(RESULTS_DIR / "hyperparameter_selection.csv")
    write(hyper, "tableS4_hyperparameter_grid.csv")

    fits = pd.read_csv(RESULTS_DIR / "benchmark_fits.csv")
    grouped_rows = []
    for keys, group in fits.groupby(
        ["method", "noise_level", "panel", "archetype"], sort=False
    ):
        method, noise, panel, archetype = keys
        grouped_rows.append(
            {
                "Method": method,
                "Noise level": noise,
                "Panel": panel,
                "Archetype": archetype,
                "Phase F1": group["phase_f1"].mean(),
                "Class F1": group["class_f1"].mean(),
                "Extent RMSE (mmol/L)": group["extent_rmse_mmolL"].mean(),
                "Reconstruction RMSE (mmol/L)": group[
                    "reconstruction_rmse_mmolL"
                ].mean(),
                "Held-out RMSE (mmol/L)": group["heldout_rmse_mmolL"].mean(),
                "False-discovery rate": group["false_discovery_rate"].mean(),
                "Direction accuracy": group["direction_accuracy"].mean(),
                "Convergence fraction": group["converged"].mean(),
                "Median runtime (ms)": group["runtime_ms"].median(),
                "n": len(group),
            }
        )
    write(pd.DataFrame(grouped_rows), "tableS5_complete_model_metrics.csv")

    recovery = pd.read_csv(RESULTS_DIR / "reaction_recovery.csv")
    recovery = recovery[
        (recovery["panel"] == "full_11")
        & (recovery["noise_level"] == 0.03)
    ]
    confusion_rows = []
    for keys, group in recovery.groupby(["method", "reaction"]):
        method, reaction = keys
        truth_active = group["true_active"].astype(bool)
        predicted_active = group["recovered_active"].astype(bool)
        tp = int((truth_active & predicted_active).sum())
        fp = int((~truth_active & predicted_active).sum())
        fn = int((truth_active & ~predicted_active).sum())
        tn = int((~truth_active & ~predicted_active).sum())
        confusion_rows.append(
            {
                "Method": method,
                "Reaction": reaction,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "TN": tn,
                "Precision": tp / (tp + fp) if tp + fp else np.nan,
                "Recall": tp / (tp + fn) if tp + fn else np.nan,
            }
        )
    write(pd.DataFrame(confusion_rows), "tableS6_reaction_confusion_matrices.csv")

    next_best = pd.read_csv(RESULTS_DIR / "next_best_measurement.csv")
    ion_value = (
        next_best.groupby("candidate_measurement")
        .agg(
            mean_value_score=("measurement_value_score", "mean"),
            mean_class_f1_gain=("realised_class_f1_gain", "mean"),
            mean_support_change=("support_change", "mean"),
            mean_heldout_error_mmolL=("heldout_absolute_error_mmolL", "mean"),
            n=("scenario_id", "count"),
        )
        .reset_index()
        .sort_values("mean_value_score", ascending=False)
    )
    write(ion_value, "tableS7_missing_ion_and_measurement_value.csv")

    thermo = pd.read_csv(
        RESULTS_DIR / "thermodynamic_threshold_sensitivity.csv"
    )
    thermo_summary = (
        thermo.groupby(["si_threshold", "archetype"])
        .agg(
            phase_f1=("phase_f1", "mean"),
            class_f1=("class_f1", "mean"),
            extent_rmse_mmolL=("extent_rmse_mmolL", "mean"),
            thermodynamic_violations=("thermodynamic_violations", "mean"),
            n=("scenario_id", "count"),
        )
        .reset_index()
    )
    write(thermo_summary, "tableS8_thermodynamic_bounds.csv")

    summary = json.loads(
        (RESULTS_DIR / "analysis_summary.json").read_text(encoding="utf-8")
    )
    environment = [
        {"Item": "Python", "Value": summary["software"]["python"]},
        {"Item": "NumPy", "Value": summary["software"]["numpy"]},
        {"Item": "pandas", "Value": summary["software"]["pandas"]},
        {"Item": "PHREEQC", "Value": summary["phreeqc_version"]},
        {"Item": "PHREEQC executable", "Value": summary["phreeqc_executable"]},
        {"Item": "PHREEQC database", "Value": summary["phreeqc_database"]},
        {"Item": "Random seed", "Value": summary["random_seed"]},
        {
            "Item": "Selected lambda_l1",
            "Value": summary["selected_hyperparameters"]["lambda_l1"],
        },
        {
            "Item": "Selected lambda_l2",
            "Value": summary["selected_hyperparameters"]["lambda_l2"],
        },
        {
            "Item": "Selected support threshold (mmol/L)",
            "Value": summary["selected_hyperparameters"].get(
                "support_threshold_mmolL", ""
            ),
        },
    ]
    write(pd.DataFrame(environment), "tableS9_software_environment.csv")

    field = pd.read_csv(RESULTS_DIR / "ghana_field_pairs.csv")
    field_rows = []
    for aquifer, group in field.groupby("aquifer"):
        field_rows.append(
            {
                "Aquifer": aquifer,
                "Wet-dry pairs": len(group),
                "Lithologies": group["lithology"].nunique(),
                "Median MRS": group["mechanism_resolution_score"].median(),
                "Median held-out RMSE (mmol/L)": group[
                    "mean_heldout_rmse_mmolL"
                ].median(),
                "Median support stability": group["support_stability"].median(),
                "Identifiable (%)": 100.0
                * group["resolution_class"].eq("identifiable").mean(),
                "Equivalence-class (%)": 100.0
                * group["resolution_class"].eq("equivalence_class").mean(),
                "Partially identifiable (%)": 100.0
                * group["resolution_class"].eq(
                    "partially_identifiable"
                ).mean(),
                "Non-identifiable (%)": 100.0
                * group["resolution_class"].eq("non_identifiable").mean(),
            }
        )
    write(pd.DataFrame(field_rows), "tableS10_northern_ghana_summary.csv")

    baseline_path = RESULTS_DIR / "phreeqc_inverse_baseline.csv"
    if baseline_path.exists():
        baseline = pd.read_csv(baseline_path)
        baseline_summary = (
            baseline.groupby("archetype")
            .agg(
                success_fraction=("phreeqc_inverse_success", "mean"),
                strict_success_fraction=(
                    "phreeqc_inverse_setup",
                    lambda values: values.eq("strict_5pct").mean(),
                ),
                relaxed_fallback_fraction=(
                    "phreeqc_inverse_setup",
                    lambda values: values.eq("relaxed_20pct").mean(),
                ),
                mean_models_found=("models_found", "mean"),
                mean_minimal_models_found=("minimal_models_found", "mean"),
                first_minimal_phase_f1=("first_minimal_phase_f1", "mean"),
                first_minimal_class_f1=("first_minimal_class_f1", "mean"),
                minimal_union_class_f1=("minimal_union_class_f1", "mean"),
                oracle_minimal_class_f1=("oracle_minimal_class_f1", "mean"),
                n=("scenario_id", "count"),
            )
            .reset_index()
        )
        write(baseline_summary, "tableS11_phreeqc_inverse_baseline.csv")

    evidence_path = RESULTS_DIR / "hydrosheaf_core_evidence.csv"
    if evidence_path.exists():
        evidence = pd.read_csv(evidence_path)
        evidence_summary = (
            evidence[
                (evidence["panel"] == "full_11")
                & (evidence["noise_level"] == 0.03)
            ]
            .groupby(["reaction", "family", "equivalence_class"])
            .agg(
                mean_evidence_score=(
                    "hydrosheaf_core_evidence_score",
                    "mean",
                ),
                median_penalty_scale=("penalty_scale", "median"),
                mean_residual_alignment=("residual_alignment", "mean"),
                true_active_fraction=("true_active", "mean"),
                recovered_active_fraction=("recovered_active", "mean"),
                n=("scenario_id", "count"),
            )
            .reset_index()
            .sort_values("mean_evidence_score", ascending=False)
        )
        write(evidence_summary, "tableS12_hydrosheaf_core_evidence_gates.csv")

    field_evidence_path = RESULTS_DIR / "ghana_field_hydrosheaf_core_evidence.csv"
    if field_evidence_path.exists():
        field_evidence = pd.read_csv(field_evidence_path)
        field_evidence_summary = (
            field_evidence.groupby(["aquifer", "reaction", "family"])
            .agg(
                support_frequency=("selected", "mean"),
                mean_evidence_score=(
                    "hydrosheaf_core_evidence_score",
                    "mean",
                ),
                median_penalty_scale=("penalty_scale", "median"),
                n=("well_id", "count"),
            )
            .reset_index()
            .sort_values(["aquifer", "support_frequency"], ascending=[True, False])
        )
        write(
            field_evidence_summary,
            "tableS13_ghana_hydrosheaf_core_evidence.csv",
        )

    data_tier_path = RESULTS_DIR / "data_tier_experiment.csv"
    if data_tier_path.exists():
        data_tiers = pd.read_csv(data_tier_path)
        tier_summary = (
            data_tiers.groupby(["data_tier", "archetype"])
            .agg(
                phase_f1=("phase_f1", "mean"),
                class_f1=("class_f1", "mean"),
                false_discovery_rate=("false_discovery_rate", "mean"),
                extent_rmse_mmolL=("extent_rmse_mmolL", "mean"),
                reconstruction_rmse_mmolL=("reconstruction_rmse_mmolL", "mean"),
                n_optional_diagnostics=("n_optional_diagnostics", "max"),
                n=("scenario_id", "count"),
            )
            .reset_index()
        )
        write(tier_summary, "tableS14_data_tier_experiment.csv")

    data_tier_evidence_path = RESULTS_DIR / "data_tier_reaction_evidence.csv"
    if data_tier_evidence_path.exists():
        tier_evidence = pd.read_csv(data_tier_evidence_path)
        evidence_summary = (
            tier_evidence.groupby(["data_tier", "reaction", "family"])
            .agg(
                true_active_fraction=("true_active", "mean"),
                recovered_active_fraction=("recovered_active", "mean"),
                mean_core_evidence_score=("core_evidence_score", "mean"),
                mean_optional_evidence_score=("optional_evidence_score", "mean"),
                mean_combined_evidence_score=("combined_evidence_score", "mean"),
                median_combined_penalty_scale=("combined_penalty_scale", "median"),
                n=("scenario_id", "count"),
            )
            .reset_index()
            .sort_values(["data_tier", "mean_optional_evidence_score"], ascending=[True, False])
        )
        write(evidence_summary, "tableS15_data_tier_reaction_evidence.csv")

    data_tier_resolution_path = RESULTS_DIR / "data_tier_evidence_lifted_resolution.csv"
    if data_tier_resolution_path.exists():
        tier_resolution = pd.read_csv(data_tier_resolution_path)
        tier_resolution["class_active"] = (
            tier_resolution["n_true_active_members"] > 0
        )

        def active_top_fraction(values: pd.Series) -> float:
            active = tier_resolution.loc[values.index, "class_active"]
            subset = values[active]
            return float(subset.mean()) if len(subset) else np.nan

        resolution_summary = (
            tier_resolution.groupby(["data_tier", "class_id", "members"])
            .agg(
                mean_elri=("evidence_lifted_resolution_index", "mean"),
                median_elri=("evidence_lifted_resolution_index", "median"),
                mean_top_probability=("top_probability", "mean"),
                class_active_fraction=("class_active", "mean"),
                conditionally_preferred_or_resolved_fraction=(
                    "resolution_status",
                    lambda values: values.isin(
                        [
                            "conditionally_preferred",
                            "evidence_lifted_resolved",
                        ]
                    ).mean(),
                ),
                top_member_true_active_fraction_all=(
                    "top_member_true_active",
                    "mean",
                ),
                top_member_true_active_fraction_when_class_active=(
                    "top_member_true_active",
                    active_top_fraction,
                ),
                n=("scenario_id", "count"),
            )
            .reset_index()
            .sort_values(["data_tier", "class_id"])
        )
        write(
            resolution_summary,
            "tableS16_evidence_lifted_resolution.csv",
        )

    external_resolution_path = (
        RESULTS_DIR / "external_field_evidence_lifted_resolution.csv"
    )
    if external_resolution_path.exists():
        external_resolution = pd.read_csv(external_resolution_path)
        external_summary = (
            external_resolution.groupby(["dataset", "data_tier", "class_id", "members"])
            .agg(
                mean_elri=("evidence_lifted_resolution_index", "mean"),
                median_elri=("evidence_lifted_resolution_index", "median"),
                mean_top_probability=("top_probability", "mean"),
                conditionally_preferred_or_resolved_fraction=(
                    "resolution_status",
                    lambda values: values.isin(
                        [
                            "conditionally_preferred",
                            "evidence_lifted_resolved",
                        ]
                    ).mean(),
                ),
                n_edges=("edge_id", "nunique"),
                n_class_evaluations=("edge_id", "count"),
            )
            .reset_index()
            .sort_values(["dataset", "data_tier", "class_id"])
        )
        write(
            external_summary,
            "tableS17_external_field_evidence_lifted_resolution.csv",
        )


def main() -> None:
    main_table()
    supplementary_tables()
    print("Generated M5 main and supplementary tables.")


if __name__ == "__main__":
    main()
