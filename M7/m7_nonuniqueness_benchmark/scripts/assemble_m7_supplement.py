"""Assemble a readable M7 supplement from immutable result tables.

The complete machine-readable CSVs remain the authoritative supplementary
records.  The Word/PDF supplement contains compact, claim-bearing views so
that every printed table remains legible on portrait pages.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
SUPPLEMENT = PROJECT / "manuscript" / "supplementary"
TABLES = PROJECT / "tables" / "publication"
FIGURES_TABLES = SUPPLEMENT / "Supplementary-Figures-and-Tables.md"
OUTPUT = SUPPLEMENT / "Supplementary-Information.md"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name)


def _table(frame: pd.DataFrame, digits: int = 4) -> str:
    display = frame.copy()
    numeric = display.select_dtypes(include="number").columns
    display[numeric] = display[numeric].round(digits)
    display = display.fillna("NA")
    return display.to_markdown(index=False)


def _compact_figures_and_tables() -> str:
    s1 = _read("tableS1_all_evidence_conditions.csv")[
        ["condition", "panel", "pr_auc", "brier", "log_loss", "mean_edge_entropy"]
    ].rename(
        columns={
            "condition": "Condition",
            "panel": "Panel",
            "pr_auc": "PR-AUC",
            "brier": "Brier",
            "log_loss": "Log loss",
            "mean_edge_entropy": "Mean entropy",
        }
    )

    s2 = _read("tableS2_case_block_contrasts.csv")
    s2 = s2[s2["metric"].isin(["pr_auc", "brier", "log_loss"])][
        ["contrast", "metric", "mean_difference", "ci95_low", "ci95_high", "n_cases"]
    ].rename(
        columns={
            "contrast": "Contrast",
            "metric": "Metric",
            "mean_difference": "Difference",
            "ci95_low": "CI low",
            "ci95_high": "CI high",
            "n_cases": "Cases",
        }
    )
    s2["Contrast"] = s2["Contrast"].str.replace("_", " ")
    s2["Metric"] = s2["Metric"].replace(
        {"pr_auc": "PR-AUC", "brier": "Brier", "log_loss": "Log loss"}
    )

    s3_raw = _read("tableS3_topology_age_sensitivity.csv")
    stable = s3_raw["importance_stable_ess_ge_400"].astype(str).str.lower().eq("true")
    s3_raw = s3_raw.assign(stable=stable.astype(float))
    s3 = (
        s3_raw.groupby(["tracer_regime", "graph_condition"], sort=False)
        .agg(
            cases=("seed", "nunique"),
            age_mae=("age_mae_years", "mean"),
            coverage=("age_95_coverage", "mean"),
            interval_width=("mean_interval_width_years", "mean"),
            median_ess=("importance_ess", "median"),
            stable_fraction=("stable", "mean"),
        )
        .reset_index()
        .rename(
            columns={
                "tracer_regime": "Regime",
                "graph_condition": "Topology",
                "cases": "n",
                "age_mae": "MAE (yr)",
                "coverage": "Coverage",
                "interval_width": "Width (yr)",
                "median_ess": "ESS",
                "stable_fraction": "Stable",
            }
        )
    )
    s3["Regime"] = s3["Regime"].replace(
        {"informative": "Informative", "tritium_only": "3H only"}
    )
    s3["Topology"] = s3["Topology"].replace(
        {
            "none": "None",
            "partial_true": "Partial",
            "complete_true": "Complete",
            "reversed": "Reversed",
        }
    )

    s4_raw = _read("tableS4_reaction_edge_nonuniqueness.csv")
    correct = s4_raw["modal_family_correct"].astype(str).str.lower().eq("true")
    s4_raw = s4_raw.assign(modal_correct=correct.astype(float))
    s4 = (
        s4_raw.groupby(["tier", "true_process"], sort=False)
        .agg(
            edges=("edge_id", "count"),
            modal_accuracy=("modal_correct", "mean"),
            true_probability=("true_family_probability", "mean"),
            support_entropy=("family_support_entropy", "mean"),
            supported_families=("effective_supported_families", "mean"),
        )
        .reset_index()
        .rename(
            columns={
                "tier": "Tier",
                "true_process": "Process",
                "edges": "n",
                "modal_accuracy": "Accuracy",
                "true_probability": "True prob.",
                "support_entropy": "Entropy",
                "supported_families": "Eff. families",
            }
        )
    )
    s4["Tier"] = s4["Tier"].replace({"enhanced": "E", "core": "C"})
    s4["Process"] = s4["Process"].replace(
        {
            "carbonate_weathering": "Carbonate weather.",
            "carbonate_precipitation": "Carbonate precip.",
            "silicate_weathering": "Silicate weather.",
            "denitrification": "Denitrification",
            "sulfate_reduction": "Sulfate reduction",
            "iron_reduction": "Iron reduction",
        }
    )

    s5 = _read("tableS5_conflict_diagnostics.csv").rename(
        columns={
            "condition": "Condition",
            "n_edges": "N",
            "conflict_fraction": "Flagged",
            "mean_univariate_probability_span": "Span",
            "integrated_error_rate_concordant": "Error",
            "integrated_overconfident_error_fraction": "Overconf.",
        }
    )[
        ["Condition", "N", "Flagged", "Span", "Error", "Overconf."]
    ]
    s5["Condition"] = s5["Condition"].replace(
        {
            "native": "Native",
            "age_permuted": "Age perm.",
            "hydraulic_permuted": "Hyd. perm.",
            "joint_misspecified": "Joint misspec.",
        }
    )

    s6 = _read("tableS6_multiplicity_correction.csv").rename(
        columns={
            "Cases": "n",
            "Mean difference": "Difference",
            "Exact permutation p (two-sided)": "Exact p",
            "Benjamini-Hochberg adjusted p": "BH p",
        }
    ).drop(columns=["Survives q=0.05"])
    s6["Contrast"] = s6["Contrast"].replace(
        {
            "Native age added": "Native +age",
            "Native chemistry added": "Native +chem.",
            "Native hydraulics added": "Native +hyd.",
            "Permuted age added": "Permuted +age",
            "Permuted hydraulics added": "Permuted +hyd.",
            "Age + hydraulics misspecified": "Joint misspec.",
        }
    )

    s7_raw = _read("tableS7_sheaf_vs_graph_contrasts.csv")
    primary = ["pr_auc", "brier", "log_loss", "selected_f1"]
    keep = (
        (s7_raw["scenario"] == "all")
        & (s7_raw["left"] == "affine_sheaf")
        & s7_raw["metric"].isin(primary + ["conflict_localisation_pr_auc"])
    ) | (
        (s7_raw["scenario"] == "identity_limit")
        & (s7_raw["left"] == "affine_sheaf")
        & (s7_raw["right"] == "graph_laplacian")
        & s7_raw["metric"].isin(primary)
    ) | (
        (s7_raw["scenario"] == "incompatible_cycles")
        & (s7_raw["left"] == "affine_sheaf")
        & (s7_raw["right"] == "weighted_graph")
        & s7_raw["metric"].isin(["pr_auc", "conflict_localisation_pr_auc"])
    )
    s7 = s7_raw.loc[
        keep,
        ["scenario", "right", "metric", "n_cases", "mean_difference", "ci95_low", "ci95_high"],
    ].rename(
        columns={
            "scenario": "Scenario",
            "right": "Comparator",
            "metric": "Metric",
            "n_cases": "Cases",
            "mean_difference": "Difference",
            "ci95_low": "CI low",
            "ci95_high": "CI high",
        }
    )
    s7["Scenario"] = s7["Scenario"].replace(
        {
            "all": "All",
            "identity_limit": "Identity limit",
            "incompatible_cycles": "Incompatible",
        }
    )
    s7["Comparator"] = s7["Comparator"].replace(
        {
            "weighted_graph": "Edge-local",
            "graph_laplacian": "Identity",
            "affine_sheaf_permuted": "Permuted-map",
        }
    )
    s7["Metric"] = s7["Metric"].replace(
        {
            "pr_auc": "PR-AUC",
            "brier": "Brier",
            "log_loss": "Log loss",
            "selected_f1": "Selected F1",
            "conflict_localisation_pr_auc": "Conflict PR-AUC",
        }
    )

    s8_raw = _read("tableS8_public_pipeline_acceptance.csv")
    s8_summary = s8_raw[s8_raw["record_type"] == "condition_summary"][
        [
            "condition",
            "n_cases",
            "candidate_recall_mean",
            "pr_auc_mean",
            "brier_mean",
            "log_loss_mean",
            "selected_f1_mean",
        ]
    ].rename(
        columns={
            "condition": "Condition",
            "n_cases": "Cases",
            "candidate_recall_mean": "Recall",
            "pr_auc_mean": "PR-AUC",
            "brier_mean": "Brier",
            "log_loss_mean": "Log loss",
            "selected_f1_mean": "Selected F1",
        }
    )
    s8_summary["Condition"] = s8_summary["Condition"].str.replace("_", " ")
    s8_contrasts = s8_raw[s8_raw["record_type"] == "paired_contrast"][
        ["left", "right", "metric", "n_cases", "mean_difference", "ci95_low", "ci95_high"]
    ].rename(
        columns={
            "left": "Left",
            "right": "Comparator",
            "metric": "Metric",
            "n_cases": "Cases",
            "mean_difference": "Difference",
            "ci95_low": "CI low",
            "ci95_high": "CI high",
        }
    )
    for column in ["Left", "Comparator"]:
        s8_contrasts[column] = s8_contrasts[column].str.replace("_", " ")
    s8_contrasts["Metric"] = s8_contrasts["Metric"].replace(
        {"pr_auc": "PR-AUC", "brier": "Brier", "selected_f1": "Selected F1"}
    )

    s9_raw = _read("tableS9_robust_hybrid_contrasts.csv")
    primary_metrics = ["pr_auc", "brier", "log_loss"]
    mechanism_pairs = [
        ("original_hybrid", "original_affine_global"),
        ("robust_affine_global", "original_affine_global"),
        ("robust_hybrid", "original_hybrid"),
    ]
    pair_mask = pd.Series(False, index=s9_raw.index)
    for left, right in mechanism_pairs:
        pair_mask |= (s9_raw["left"] == left) & (s9_raw["right"] == right)
    keep_s9 = (
        (s9_raw["calibration_regime"] == "separate_crossfit")
        & s9_raw["metric"].isin(primary_metrics)
        & (
            (
                (s9_raw["scenario"] == "all")
                & (
                    (
                        (s9_raw["left"] == "robust_hybrid")
                        & s9_raw["right"].isin(
                            ["edge_local", "robust_hybrid_permuted"]
                        )
                    )
                    | pair_mask
                )
            )
            | (
                (s9_raw["scenario"] != "all")
                & (s9_raw["left"] == "robust_hybrid")
                & (s9_raw["right"] == "edge_local")
            )
        )
    )
    s9 = s9_raw.loc[
        keep_s9,
        [
            "scenario",
            "left",
            "right",
            "metric",
            "n_cases",
            "mean_difference",
            "ci95_low",
            "ci95_high",
        ],
    ].rename(
        columns={
            "scenario": "Scenario",
            "left": "Left",
            "right": "Comparator",
            "metric": "Metric",
            "n_cases": "n",
            "mean_difference": "Difference",
            "ci95_low": "CI low",
            "ci95_high": "CI high",
        }
    )
    model_labels = {
        "edge_local": "EL",
        "original_affine_global": "OG",
        "original_hybrid": "OH",
        "robust_affine_global": "RG",
        "robust_hybrid": "RH",
        "robust_hybrid_permuted": "PRH",
    }
    scenario_labels = {
        "all": "All",
        "identity_limit": "Identity limit",
        "heterogeneous_affine": "Heterog. affine",
        "incompatible_cycles": "Incompat. cycles",
        "noisy_missing": "Noisy/missing",
    }
    s9["Scenario"] = s9["Scenario"].map(scenario_labels)
    s9["Left"] = s9["Left"].map(model_labels)
    s9["Comparator"] = s9["Comparator"].map(model_labels)
    s9["Metric"] = s9["Metric"].replace(
        {"pr_auc": "PRAUC", "brier": "Brier", "log_loss": "LogLoss"}
    )
    s9["Comparison"] = s9["Left"] + " vs " + s9["Comparator"]
    s9["Difference [95% CI]"] = s9.apply(
        lambda row: (
            f"{row['Difference']:+.4f} "
            f"[{row['CI low']:+.4f}, {row['CI high']:+.4f}]"
        ),
        axis=1,
    )
    s9 = s9[
        ["Scenario", "Comparison", "Metric", "n", "Difference [95% CI]"]
    ]

    def _simultaneous_view(name: str, keep_mask) -> pd.DataFrame:
        raw = _read(name)
        view = raw.loc[
            keep_mask(raw),
            [
                "scenario",
                "left",
                "right",
                "metric",
                "n_cases",
                "mean_difference",
                "simultaneous_ci95_low",
                "simultaneous_ci95_high",
                "supports_benefit_fwer95",
            ],
        ].copy()
        view["Comparison"] = (
            view["left"].str.replace("_", " ")
            + " vs "
            + view["right"].str.replace("_", " ")
        )
        view["Difference [simultaneous 95% CI]"] = view.apply(
            lambda row: (
                f"{row['mean_difference']:+.4f} "
                f"[{row['simultaneous_ci95_low']:+.4f}, "
                f"{row['simultaneous_ci95_high']:+.4f}]"
            ),
            axis=1,
        )
        view["FWER support"] = view["supports_benefit_fwer95"].map(
            {True: "Yes", False: "No", "True": "Yes", "False": "No"}
        )
        return view.rename(
            columns={"scenario": "Scenario", "metric": "Metric", "n_cases": "n"}
        )[["Scenario", "Comparison", "Metric", "n", "Difference [simultaneous 95% CI]", "FWER support"]]

    s10 = _simultaneous_view(
        "tableS10_m7_4_multiplicity_adjusted.csv",
        lambda raw: (
            (raw["scenario"] == "all")
            & (raw["left"] == "affine_sheaf")
            & raw["right"].isin(["weighted_graph", "graph_laplacian", "affine_sheaf_permuted"])
            & raw["metric"].isin(["pr_auc", "brier", "log_loss"])
        )
        | (
            (raw["scenario"] == "incompatible_cycles")
            & (raw["left"] == "affine_sheaf")
            & (raw["right"] == "weighted_graph")
            & raw["metric"].isin(["pr_auc", "conflict_localisation_pr_auc"])
        ),
    )

    s11 = _simultaneous_view(
        "tableS11_m7_5_multiplicity_adjusted.csv",
        lambda raw: (
            (raw["calibration_regime"] == "separate_crossfit")
            & (raw["scenario"] == "all")
            & (
                (
                    (raw["left"] == "robust_hybrid")
                    & raw["right"].isin(["edge_local", "robust_hybrid_permuted"])
                    & raw["metric"].isin(["pr_auc", "brier", "log_loss", "conflict_localisation_pr_auc"])
                )
                | (
                    (raw["left"] == "robust_affine_global")
                    & (raw["right"] == "original_affine_global")
                    & raw["metric"].isin(["brier", "log_loss"])
                )
            )
        ),
    )

    s12_raw = _read("tableS12_precision_and_power.csv")
    s12 = s12_raw[
        [
            "design",
            "metric",
            "source_n_cases",
            "planned_n_cases",
            "planning_margin_benefit_scale",
            "probability_ci_excludes_zero_in_favourable_direction",
            "probability_ci_clears_planning_margin",
            "interpretation_status",
        ]
    ].rename(
        columns={
            "design": "Design",
            "metric": "Metric",
            "source_n_cases": "Source n",
            "planned_n_cases": "Planned n",
            "planning_margin_benefit_scale": "Margin",
            "probability_ci_excludes_zero_in_favourable_direction": "P(CI favourable)",
            "probability_ci_clears_planning_margin": "P(CI clears margin)",
            "interpretation_status": "Status",
        }
    )
    s12["Design"] = s12["Design"].replace(
        {
            "Process-based evidence-panel locked test": "Evidence panel",
            "Locked sheaf-versus-weighted-graph representation test": "Representation",
            "Locked local-first/global-fallback estimator diagnostic": "Estimator",
            "Topology-conditioned-age test (informative)": "Topology-age: two-tracer",
            "Topology-conditioned-age test (tritium_only)": "Topology-age: tritium",
            "Reaction-family recovery test": "Reaction recovery",
        }
    )
    s12["Metric"] = s12["Metric"].replace(
        {
            "age_mae_years": "age MAE",
            "mean_interval_width_years": "interval width",
            "modal_family_accuracy": "modal accuracy",
        }
    )
    s12["Status"] = s12["Status"].replace(
        {
            "DEVELOPMENT_ONLY_EMPIRICAL_PLANNING": "Development planning",
            "POST_TEST_REPLICATION_PLANNING_ONLY": "Replication planning",
        }
    )

    s13 = _read("tableS13_public_pipeline_selection.csv").rename(
        columns={
            "condition": "Condition",
            "candidate_edges": "Candidates",
            "selected_edges": "Selected",
            "minimum_selected_probability": "Min probability",
            "conditional_true_positive": "TP",
            "conditional_false_positive": "FP",
            "conditional_false_negative": "Conditional FN",
            "end_to_end_false_negative": "End-to-end FN",
            "end_to_end_selected_f1": "End-to-end F1",
            "candidate_recall": "Candidate recall",
        }
    )[["Condition", "Candidates", "Selected", "Min probability", "TP", "FP", "Conditional FN", "End-to-end FN", "End-to-end F1", "Candidate recall"]]

    s14 = _read("tableS14_m7_6_m3_mechanism.csv").rename(
        columns={
            "contrast": "Contrast",
            "mean_difference": "Difference",
            "ci95_low": "CI low",
            "ci95_high": "CI high",
            "n_cases": "Cases",
            "n_bootstrap": "Bootstrap",
            "decision": "Decision",
        }
    )[["Contrast", "Difference", "CI low", "CI high", "Cases", "Bootstrap", "Decision"]]
    s14["Contrast"] = s14["Contrast"].str.replace("_", " ", regex=False)

    parts = [
        "# Supplementary Figures and Tables",
        "",
        "**Conditional evidence integration and the incremental contribution of sheaf "
        "structure in controlled-synthetic groundwater benchmarks**",
        "",
        "The complete authoritative supplementary tables are the version-controlled CSV "
        "files in `tables/publication/`. The compact views below preserve the "
        "claim-bearing comparisons while keeping the Word and PDF tables legible; omitted "
        "columns and per-case or per-edge rows remain available without loss in the cited "
        "CSV. Tables S1-S6 derive from locked process-based integration outputs, Table S7 "
        "from the prospectively locked representation benchmark, Table S8 from the strict "
        "public-pipeline acceptance run, Table S9 from the prospectively locked estimator "
        "diagnostic, and Tables S10-S13 "
        "from the post-review family-wise, precision and selection audit. Table S14 "
        "is the auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic.",
        "",
        "## Figure S1",
        "",
        "![](figures/supporting_validation/figure_s1_model_domain_map.png)",
        "",
        "**Figure S1.** Locked synthetic model domain (representative realization 4101). "
        "The MODFLOW/MODPATH truth network, observation nodes and hydraulic heads are shown "
        "in synthetic model-space coordinates. This is not a geographic map.",
        "",
        "## Table S1",
        "",
        "**Table S1.** Evidence-panel performance in every locked process-based integration condition. Compact "
        "view of `tables/publication/tableS1_all_evidence_conditions.csv`.",
        "",
        _table(s1),
        "",
        "## Table S2",
        "",
        "**Table S2.** Predeclared evidence contrasts for discrimination and calibration "
        "(10,000 case-block bootstrap replicates). Complete five-metric record: "
        "`tables/publication/tableS2_case_block_contrasts.csv`.",
        "",
        _table(s2),
        "",
        "## Table S3",
        "",
        "**Table S3.** Topology-to-age sensitivity aggregated over twelve locked cases. "
        "Values are means except median ESS; the complete CSV also reports bias, entropy, "
        "order violations and importance weights. Complete 96-row case audit: "
        "`tables/publication/tableS3_topology_age_sensitivity.csv`.",
        "",
        _table(s3),
        "",
        "## Table S4",
        "",
        "**Table S4.** Reaction-family non-uniqueness aggregated by benchmark tier and true "
        "process (C, core; E, enhanced). Complete 216-row edgewise audit: "
        "`tables/publication/tableS4_reaction_edge_nonuniqueness.csv`.",
        "",
        _table(s4),
        "",
        "## Table S5",
        "",
        "**Table S5.** Conflict diagnostics under native and adverse evidence conditions. "
        "The complete CSV additionally reports error rates within flagged and concordant "
        "subsets: `tables/publication/tableS5_conflict_diagnostics.csv`.",
        "",
        _table(s5),
        "",
        "## Table S6",
        "",
        "**Table S6.** Exact paired-permutation tests with Benjamini-Hochberg correction "
        "for the predeclared process-based integration contrast family. Complete record: "
        "`tables/publication/tableS6_multiplicity_correction.csv`.",
        "",
        _table(s6),
        "",
        "## Table S7",
        "",
        "**Table S7.** Claim-bearing representation-benchmark sheaf-versus-graph contrasts, including the "
        "identity limit and incompatible-cycle diagnostics (10,000 case-block bootstrap "
        "replicates). Complete 120-row contrast matrix: "
        "`tables/publication/tableS7_sheaf_vs_graph_contrasts.csv`.",
        "",
        _table(s7),
        "",
        "## Table S8",
        "",
        "**Table S8.** Strict public-pipeline system acceptance. The execution criterion "
        "passed, whereas a general full-sheaf incremental-performance claim did not.",
        "",
        "### Condition summary",
        "",
        _table(s8_summary),
        "",
        "### Paired case-block contrasts",
        "",
        _table(s8_contrasts),
        "",
        "Complete record: `tables/publication/tableS8_public_pipeline_acceptance.csv`.",
        "",
        "## Table S9",
        "",
        "**Table S9.** Claim-bearing estimator-diagnostic robust/hybrid contrasts under the primary "
        "separately cross-fitted calibration regime. Differences are left minus "
        "comparator; lower Brier score and log loss are favourable. Estimator "
        "abbreviations are EL, edge-local; OG, original affine-global; OH, original "
        "hybrid; RG, robust affine-global; RH, robust hybrid; and PRH, permuted robust "
        "hybrid. PRAUC denotes precision-recall area under the curve. The full "
        "560-row matrix, including the shared-calibrator regime and secondary metrics, "
        "is `tables/publication/tableS9_robust_hybrid_contrasts.csv`.",
        "",
        _table(s9),
        "",
        "## Table S10",
        "",
        "**Table S10.** Selected representation-benchmark contrasts with simultaneous 95% intervals "
        "controlling all 120 published contrasts as one family. Complete record: "
        "`tables/publication/tableS10_m7_4_multiplicity_adjusted.csv`.",
        "",
        _table(s10),
        "",
        "## Table S11",
        "",
        "**Table S11.** Selected estimator-diagnostic contrasts with simultaneous 95% intervals "
        "controlling all 560 published contrasts as one family. The selected "
        "robust-hybrid arm had local weight 1.0 and is the local-first/global-"
        "fallback estimator in the main text. Complete record: "
        "`tables/publication/tableS11_m7_5_multiplicity_adjusted.csv`.",
        "",
        _table(s11),
        "",
        "## Table S12",
        "",
        "**Table S12.** Post-review empirical precision and future-replication "
        "planning audit (20,000 simulations). Margins were not prespecified or "
        "field validated; POST_TEST rows are not evidence for completed tests. "
        "Complete record: `tables/publication/tableS12_precision_and_power.csv`.",
        "",
        _table(s12),
        "",
        "## Table S13",
        "",
        "**Table S13.** Public-pipeline selection and confusion-count audit. All "
        "generated candidates were retained; no scalar probability threshold was "
        "applied. Complete record: "
        "`tables/publication/tableS13_public_pipeline_selection.csv`.",
        "",
        _table(s13),
        "",
        "## Table S14",
        "",
        "**Table S14.** Auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic. "
        "The complete result record is `tables/publication/tableS14_m7_6_m3_mechanism.csv`; "
        "this table is not field validation and does not identify the USGS cause.",
        "",
        _table(s14),
    ]
    return "\n".join(parts).rstrip() + "\n"


def main() -> int:
    methods = (SUPPLEMENT / "Supplementary-Methods.md").read_text(encoding="utf-8")
    figures_tables = _compact_figures_and_tables()
    FIGURES_TABLES.write_text(figures_tables, encoding="utf-8")
    assembled = methods.rstrip() + "\n\n" + figures_tables.lstrip()
    if "[[TAB:" in assembled or "[[FIG:" in assembled:
        raise RuntimeError("Unresolved artifact token remains in supplementary output")
    OUTPUT.write_text(assembled.rstrip() + "\n", encoding="utf-8")
    print(f"M7 supplementary information -> {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
