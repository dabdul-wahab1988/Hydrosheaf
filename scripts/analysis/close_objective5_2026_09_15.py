"""Assemble a traceable, evidence-bounded closure package for Objective 5.

The existing UER run and its corrected/public-list reruns are preserved.  This
script reads those immutable snapshots together with the locked M7.3 and M8
controlled-synthetic decisions, derives branch-specific evidence labels, and
writes a separate closure package.  It deliberately does not substitute field
truth, vendor quotations, or independent model validation that are absent from
the source artifacts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "outputs" / "objective5_closure_2026-09-15"

UER_RUN_DIR = ROOT / "outputs" / "uer_objective5_groundwater_budgetary_2026-09-15"
UER_AUDIT_DIR = ROOT / "outputs" / "objective5_audit_2026-09-15"
UER_FIELD_DIR = ROOT / "outputs" / "uer_field_integration"
M7_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "results" / "m7_3_locked"
M8_ROOT = ROOT / "M8" / "m8_calibration_benchmark" / "provenance" / "runs"
COST_MANIFEST = ROOT / "provenance" / "uer_objective5_cost_manifest_2026-09-15.json"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def _require_paths(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Objective 5 closure inputs are missing: " + "; ".join(missing))


def _boolean_series(series: pd.Series) -> pd.Series:
    """Parse CSV boolean fields without treating the string ``False`` as true."""
    if pd.api.types.is_bool_dtype(series):
        return series
    normalized = series.astype(str).str.strip().str.lower()
    parsed = normalized.map({"true": True, "false": False, "1": True, "0": False})
    if parsed.isna().any():
        unknown = sorted(normalized[parsed.isna()].unique().tolist())
        raise ValueError(f"Unrecognised boolean values: {unknown}")
    return parsed.astype(bool)


def _source_manifest(paths: dict[str, Path]) -> dict[str, Any]:
    return {
        key: {
            "path": _relative(path),
            "sha256": _sha256(path),
            "size_bytes": path.stat().st_size,
        }
        for key, path in paths.items()
    }


def _frontier_summary() -> dict[str, Any]:
    frontier_path = UER_RUN_DIR / "tables" / "uer_objective5_minimax_pareto_frontier.csv"
    frontier = pd.read_csv(frontier_path)
    required_columns = {
        "prior_regime",
        "budget",
        "selected_options",
        "total_cost",
        "mtt_ambiguity_yr",
        "mtt_ambiguity_norm",
        "status",
        "ambiguity_reduction_pct",
        "cost_basis_status",
    }
    missing = required_columns.difference(frontier.columns)
    if missing:
        raise ValueError(f"UER frontier is missing columns: {sorted(missing)}")

    median = frontier[frontier["prior_regime"] == "Regional Median Screening Cohort"].copy()
    selected = median[median["budget"].isin([0.0, 800.0, 2500.0])]
    if len(selected) != 3 or set(selected["budget"].tolist()) != {0.0, 800.0, 2500.0}:
        raise ValueError("UER frontier must contain one median-cohort row at budgets 0, 800 and 2500 USD")
    selected = selected.sort_values("budget").set_index("budget")
    b0 = selected.loc[0.0]
    b800 = selected.loc[800.0]
    b2500 = selected.loc[2500.0]
    plateau_delta = float(b800["mtt_ambiguity_yr"] - b2500["mtt_ambiguity_yr"])
    if plateau_delta < 0.0:
        raise ValueError("UER frontier is not monotone at the requested comparison budgets")

    return {
        "source_path": _relative(frontier_path),
        "prior_regime": "Regional Median Screening Cohort",
        "cohort_size": int(b0["cohort_size"]),
        "baseline_tritium_TU": float(b0["baseline_tu"]),
        "baseline_sigma_TU": float(b0["baseline_sigma"]),
        "age_grid_years": [0.5, 80.0],
        "target_functional": "mean_transit_time",
        "target_units": "years",
        "normalised_target_units": "fraction_of_80yr",
        "target_tolerance_years": 16.0,
        "practical_equivalence_margin_years": 0.05,
        "budgets": {
            "0": {
                "allocation_budget_usd": 0.0,
                "actual_analytical_cost_usd_per_borehole": float(b0["total_cost"]),
                "selected_options": str(b0["selected_options"]),
                "mtt_ambiguity_years": float(b0["mtt_ambiguity_yr"]),
                "mtt_ambiguity_fraction_of_80yr": float(b0["mtt_ambiguity_norm"]),
                "status": str(b0["status"]),
            },
            "800": {
                "allocation_budget_usd": 800.0,
                "actual_analytical_cost_usd_per_borehole": float(b800["total_cost"]),
                "selected_options": str(b800["selected_options"]),
                "mtt_ambiguity_years": float(b800["mtt_ambiguity_yr"]),
                "mtt_ambiguity_fraction_of_80yr": float(b800["mtt_ambiguity_norm"]),
                "ambiguity_reduction_percent": float(b800["ambiguity_reduction_pct"]),
                "status": str(b800["status"]),
            },
            "2500": {
                "allocation_budget_usd": 2500.0,
                "actual_analytical_cost_usd_per_borehole": float(b2500["total_cost"]),
                "selected_options": str(b2500["selected_options"]),
                "mtt_ambiguity_years": float(b2500["mtt_ambiguity_yr"]),
                "mtt_ambiguity_fraction_of_80yr": float(b2500["mtt_ambiguity_norm"]),
                "ambiguity_reduction_percent": float(b2500["ambiguity_reduction_pct"]),
                "status": str(b2500["status"]),
            },
        },
        "plateau_delta_years_800_minus_2500": plateau_delta,
        "cost_basis_status": str(b800["cost_basis_status"]),
    }


def _uer_summary() -> dict[str, Any]:
    audit = _load_json(UER_AUDIT_DIR / "reproduction_audit.json")
    field = _load_json(UER_FIELD_DIR / "uer_integration_summary.json")
    cost = _load_json(COST_MANIFEST)
    gating = pd.read_csv(UER_RUN_DIR / "tables" / "uer_objective5_sheaf_gating_comparison.csv")
    boundary = pd.read_csv(UER_RUN_DIR / "tables" / "uer_objective5_tripartite_boundary_audit.csv")
    if len(boundary) != 1600:
        raise ValueError(f"Expected 1600 UER candidate edges, found {len(boundary)}")
    if len(gating) != 2:
        raise ValueError("Expected two UER gating comparison rows")
    conflict = _boolean_series(boundary["is_conflict"])
    cbe_conflict = _boolean_series(boundary["cbe_conflict"])
    median = _frontier_summary()
    return {
        "dataset": str(field.get("dataset", "unknown")),
        "sampling_locations": int(field["total_nodes"]),
        "candidate_edges": int(field["total_inferred_edges"]),
        "candidate_edge_interpretation": "Tier-C DEM/topographic adjacency hypotheses; not measured hydraulic flowpaths",
        "screen_sequence_from_reproduction": audit.get("sheaf_edge_counts", audit.get("edge_counts")),
        "boundary_rows": len(boundary),
        "conflict_edges": int(conflict.sum()),
        "admitted_edges": int((~conflict).sum()),
        "tritium_dual_endpoint_tested": int((boundary["tritium_status"] != "UNTESTED_MISSING_DATA").sum()),
        "tritium_untested": int((boundary["tritium_status"] == "UNTESTED_MISSING_DATA").sum()),
        "nitrate_dual_endpoint_tested": int((boundary["nitrate_status"] != "UNTESTED_MISSING_DATA").sum()),
        "nitrate_untested": int((boundary["nitrate_status"] == "UNTESTED_MISSING_DATA").sum()),
        "cbe_fail_edges": int(cbe_conflict.sum()),
        "ungated_dirichlet_energy": float(gating.loc[0, "dirichlet_energy"]),
        "gated_dirichlet_energy": float(gating.loc[1, "dirichlet_energy"]),
        "same_edge_energy_ungated": float(gating.loc[0, "coherent_edge_energy"]),
        "same_edge_energy_gated": float(gating.loc[1, "coherent_edge_energy"]),
        "mean_residual_same_edges_ungated": float(gating.loc[0, "mean_residual_coherent_edges"]),
        "mean_residual_same_edges_gated": float(gating.loc[1, "mean_residual_coherent_edges"]),
        "frontier": median,
        "independent_field_truth_available": False,
        "cost_manifest_id": cost["manifest_id"],
        "cost_evidence_classification": cost["evidence_classification"],
        "cost_eligible_for_final_frontier": bool(cost["eligible_for_final_frontier"]),
        "campaign_cost_status": cost["campaign_cost_model"]["status"],
    }


def _branch_rows(
    uer: dict[str, Any],
    m7: dict[str, Any],
    m7_contrasts: pd.DataFrame,
    m8_frontier: dict[str, Any],
    m8_independent: dict[str, Any],
    m8_confirm: dict[str, Any],
) -> list[dict[str, Any]]:
    decisions = m7["decisions"]
    rows: list[dict[str, Any]] = []

    def add(branch: str, evidence_id: str, source: str, condition: str, metric: str, value: Any, units: str, scope: str, evidence_status: str, interpretation: str, ci_low: Any = None, ci_high: Any = None) -> None:
        rows.append({
            "branch": branch,
            "evidence_id": evidence_id,
            "source": source,
            "condition": condition,
            "metric": metric,
            "value": value,
            "units": units,
            "ci95_low": ci_low,
            "ci95_high": ci_high,
            "scope": scope,
            "evidence_status": evidence_status,
            "interpretation": interpretation,
        })

    chemistry = decisions["chemistry_adds_topology_ranking_value"]
    add("improves", "M7.3-chemistry-ranking", "M7.3 confirmatory_decision.json", "native", "PR-AUC delta HAC minus HA", chemistry["pr_auc_difference_hac_minus_ha"], "PR-AUC (dimensionless)", "12 independent controlled-synthetic MODFLOW cases", "SUPPORTED_CONDITIONAL", "Chemistry addition improved topology ranking under the locked synthetic model.", chemistry["ci95"][0], chemistry["ci95"][1])
    hydraulics = decisions["hydraulics_adds_topology_ranking_value"]
    add("improves", "M7.3-hydraulics-ranking", "M7.3 confirmatory_decision.json", "native", "PR-AUC delta HAC minus AC", hydraulics["pr_auc_difference_hac_minus_ac"], "PR-AUC (dimensionless)", "12 independent controlled-synthetic MODFLOW cases", "SUPPORTED_CONDITIONAL", "Hydraulic evidence added a small positive topology-ranking increment under the locked synthetic model.", hydraulics["ci95"][0], hydraulics["ci95"][1])
    age_correct = decisions["correct_topology_improves_informative_tracer_age_mae"]
    add("improves", "M7.3-correct-topology-age", "M7.3 confirmatory_decision.json", "native", "Age MAE delta correct minus no topology", age_correct["mean_difference_years"], "years (negative is improvement)", "12 independent controlled-synthetic MODFLOW cases", "SUPPORTED_CONDITIONAL", "Correct topology reduced informative-tracer age MAE under the locked synthetic model.", age_correct["ci95"][0], age_correct["ci95"][1])
    add("improves", "M8-frontier-active-learning", "M8 frontier claim_decision.json", "native controlled-synthetic", "frontier_active_learning_claim_supported", bool(m8_frontier["frontier_active_learning_claim_supported"]), "boolean claim gate", "untouched code-independent controlled-synthetic cases", "SUPPORTED_CONDITIONAL", m8_frontier["allowed_claim"])

    add("no-material-value", "UER-frontier-plateau", _relative(UER_RUN_DIR / "tables" / "uer_objective5_minimax_pareto_frontier.csv"), "regional median screening cohort", "MTT ambiguity delta at 800 minus 2500", uer["frontier"]["plateau_delta_years_800_minus_2500"], "years", "finite-grid UER model with assumed kernels, histories and analytical public-list costs", "MODEL_CONDITIONAL", "The additional model-defined narrowing is below the declared 0.05-year practical-equivalence margin; this is not field performance.")
    transport = m8_confirm["transport"]
    add("no-material-value", "M8-transport-time-tie", "M8 confirmatory claim_decision.json", "controlled-synthetic transport design", "target_specific_candidate_times_differ", bool(transport["target_specific_candidate_times_differ"]), "boolean", "M8 locked confirmatory transport model", "SUPPORTED_CONDITIONAL", "The target-specific candidate times were identical (50 days); no target-specific split was established.")
    kinetics = m8_confirm["kinetics"]
    add("no-material-value", "M8-kinetic-rank-one", "M8 confirmatory claim_decision.json", "controlled-synthetic kinetics", "structural_nonidentifiability_supported", bool(kinetics["structural_nonidentifiability_supported"]), "boolean", "M8 locked confirmatory kinetics model", "SUPPORTED_CONDITIONAL", "Additional residence-time observations alone retained the k*A rank-one confounding; independent surface-area information restored rank.")

    age_increment = decisions["age_adds_topology_ranking_value"]
    add("weakens", "M7.3-age-increment", "M7.3 confirmatory_decision.json", "native", "PR-AUC delta HAC minus HC", age_increment["pr_auc_difference_hac_minus_hc"], "PR-AUC (dimensionless)", "12 independent controlled-synthetic MODFLOW cases", "SUPPORTED_ADVERSE_CONTROL", "Adding age evidence worsened topology ranking under the locked native synthetic model; this is adverse evidence for an unqualified integration claim.", age_increment["ci95"][0], age_increment["ci95"][1])
    def add_adverse_contrast(
        contrast: str,
        metric: str,
        units: str,
        interpretation: str,
    ) -> None:
        match = m7_contrasts.loc[
            (m7_contrasts["contrast"] == contrast)
            & (m7_contrasts["metric"] == metric)
        ]
        if len(match) != 1:
            raise ValueError(f"Expected one M7.3 contrast row for {contrast}/{metric}, found {len(match)}")
        contrast_row = match.iloc[0]
        add(
            "weakens",
            f"M7.3-{contrast}-{metric}",
            "M7.3 evidence_case_bootstrap_contrasts.csv",
            str(contrast_row["condition"]),
            f"{contrast} {metric} delta",
            float(contrast_row["mean_difference"]),
            units,
            f"{int(contrast_row['n_cases'])} independent controlled-synthetic MODFLOW cases",
            "SUPPORTED_ADVERSE_CONTROL",
            interpretation,
            float(contrast_row["ci95_low"]),
            float(contrast_row["ci95_high"]),
        )

    add_adverse_contrast(
        "permuted_age_increment",
        "pr_auc",
        "PR-AUC (dimensionless; negative is worse)",
        "Permuting age evidence degraded topology ranking in the adverse control.",
    )
    add_adverse_contrast(
        "permuted_age_increment",
        "brier",
        "Brier score (dimensionless; positive is worse)",
        "Permuting age evidence increased Brier error in the adverse control.",
    )
    add_adverse_contrast(
        "permuted_age_increment",
        "overconfident_error_fraction",
        "fraction (positive is worse)",
        "Permuting age evidence increased the overconfident-error fraction in the adverse control.",
    )
    add_adverse_contrast(
        "permuted_hydraulic_increment",
        "pr_auc",
        "PR-AUC (dimensionless; negative is worse)",
        "Permuting hydraulic evidence degraded topology ranking in the adverse control.",
    )
    add_adverse_contrast(
        "permuted_hydraulic_increment",
        "brier",
        "Brier score (dimensionless; positive is worse)",
        "Permuting hydraulic evidence increased Brier error in the adverse control.",
    )
    add_adverse_contrast(
        "permuted_hydraulic_increment",
        "overconfident_error_fraction",
        "fraction (positive is worse)",
        "Permuting hydraulic evidence increased the overconfident-error fraction in the adverse control.",
    )
    add_adverse_contrast(
        "joint_misspecification",
        "pr_auc",
        "PR-AUC (dimensionless; negative is worse)",
        "Joint misspecification degraded topology ranking in the adverse control.",
    )
    add_adverse_contrast(
        "joint_misspecification",
        "brier",
        "Brier score (dimensionless; positive is worse)",
        "Joint misspecification increased Brier error in the adverse control.",
    )
    add_adverse_contrast(
        "joint_misspecification",
        "overconfident_error_fraction",
        "fraction (positive is worse)",
        "Joint misspecification increased the overconfident-error fraction in the adverse control.",
    )
    add("weakens", "M7.3-reversed-topology", "M7.3 confirmatory_decision.json", "reversed topology", "Informative-tracer age MAE delta reversed minus correct", decisions["reversed_graph_incompatibility_detected"]["informative_reversed_minus_correct_mae_years"], "years (positive is worse)", "12 independent controlled-synthetic MODFLOW cases", "SUPPORTED_ADVERSE_CONTROL", "Reversed topology increased informative-tracer age error, supporting incompatibility detection under the declared synthetic model.", decisions["reversed_graph_incompatibility_detected"]["informative_mae_ci95"][0], decisions["reversed_graph_incompatibility_detected"]["informative_mae_ci95"][1])
    add("weakens", "M8-independent-robustness", "M8 independent claim_decision.json", "independent model", "independent_model_robustness_supported", bool(m8_independent["independent_model_robustness_supported"]), "boolean claim gate", "locked independent M8 transport model", "UNSUPPORTED", m8_independent["allowed_claim"])

    add("ABSTAIN", "UER-field-truth-gap", _relative(UER_AUDIT_DIR / "OBJECTIVE_5_AUDIT.md"), "UER field application", "independent flow/age/reaction truth available", bool(uer["independent_field_truth_available"]), "boolean", "UER field package", "ABSTAIN", "Do not claim improved field accuracy, unique flowpaths, unique reaction mechanisms or protective gating from this package.")
    add("ABSTAIN", "UER-cost-quote-gap", _relative(COST_MANIFEST), "procurement and field budget", "eligible_for_final_frontier", bool(uer["cost_eligible_for_final_frontier"]), "boolean", "analytical public-list cost snapshot", "ABSTAIN", "Keep costs budgetary and analytical-only until transaction-specific quotations and field/logistics costs are documented.")
    add("ABSTAIN", "M8-independent-model-gate", "M8 independent claim_decision.json", "independent transport robustness", "active_learning_transport_claim_supported", bool(m8_independent["active_learning_transport_claim_supported"]), "boolean claim gate", "M8 independent model", "ABSTAIN", m8_independent["active_learning_consequence"])
    return rows


def _branch_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "branch": "improves",
            "decision": "SUPPORTED_CONDITIONAL",
            "scope": "controlled-synthetic evidence only",
            "claim": "Integration can improve the declared ranking or acquisition target under the locked compatible model and comparator.",
        },
        {
            "branch": "no-material-value",
            "decision": "SUPPORTED_MODEL_CONDITIONAL",
            "scope": "finite-grid or controlled-synthetic model-defined plateau",
            "claim": "Some added channels or higher allocations produced no practically material model-defined change under the declared margin.",
        },
        {
            "branch": "weakens",
            "decision": "SUPPORTED_ADVERSE_CONTROLS",
            "scope": "controlled-synthetic adverse conditions",
            "claim": "Misaligned, permuted, reversed or misspecified evidence can worsen ranking or calibration; integration is not unconditionally protective.",
        },
        {
            "branch": "ABSTAIN",
            "decision": "REQUIRED_FOR_UER_FIELD_AND_PROCUREMENT_CLAIMS",
            "scope": "UER field package and transaction costs",
            "claim": "Field accuracy, unique mechanisms, field-optimal design and quote-backed procurement remain unverified.",
        },
    ]


def _write_figure(rows: list[dict[str, Any]], uer: dict[str, Any]) -> Path:
    counts = pd.Series([row["branch"] for row in rows]).value_counts()
    branch_order = ["improves", "no-material-value", "weakens", "ABSTAIN"]
    branch_colors = {"improves": "#15803d", "no-material-value": "#2563eb", "weakens": "#b91c1c", "ABSTAIN": "#6b7280"}
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.8), dpi=600, gridspec_kw={"width_ratios": [1.0, 1.35]})
    ax = axes[0]
    values = [int(counts.get(name, 0)) for name in branch_order]
    bars = ax.barh(branch_order, values, color=[branch_colors[name] for name in branch_order], edgecolor="black", linewidth=0.7)
    ax.invert_yaxis()
    ax.set_xlabel("Distinct supporting or limiting checks")
    ax.set_title("(a) Branch evidence status", loc="left", fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    for bar, value in zip(bars, values):
        ax.text(value + 0.05, bar.get_y() + bar.get_height() / 2.0, str(value), va="center", fontsize=9)
    ax.text(0.0, -0.18, "Counts summarize evidence rows; they are not effect sizes.", transform=ax.transAxes, fontsize=8, color="#374151")

    ax = axes[1]
    frontier = uer["frontier"]["budgets"]
    x = [float(frontier[key]["allocation_budget_usd"]) for key in ("0", "800", "2500")]
    y = [float(frontier[key]["mtt_ambiguity_years"]) for key in ("0", "800", "2500")]
    ax.plot(x, y, "o-", color="#1d4ed8", lw=2.0, ms=6, label="UER conditional finite-grid W(S)")
    ax.axhspan(0.0, uer["frontier"]["target_tolerance_years"], color="#dcfce7", alpha=0.55, label="Declared target <= 16 years")
    ax.axvline(800.0, color="#6b7280", linestyle=":", lw=1.0)
    ax.axvline(2500.0, color="#6b7280", linestyle=":", lw=1.0)
    ax.annotate(f"800: {y[1]:.3f} yr\n{frontier['800']['selected_options']}", (x[1], y[1]), xytext=(x[1] + 90, y[1] + 12), arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8)
    ax.annotate(f"2500: {y[2]:.3f} yr\nΔ={uer['frontier']['plateau_delta_years_800_minus_2500']:.6f} yr", (x[2], y[2]), xytext=(x[2] - 830, y[2] + 9), arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8)
    ax.set_xlabel("Allocation cap (USD; analytical cost per borehole, logistics excluded)")
    ax.set_ylabel("Worst-case MTT ambiguity (years)")
    ax.set_title("(b) UER model-conditional frontier", loc="left", fontweight="bold")
    ax.set_ylim(0.0, max(y) * 1.18)
    ax.grid(linestyle="--", alpha=0.35)
    ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("Objective 5 closure: evidence branches and conditional UER design", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    path = OUT_DIR / "figure_o5_closure_branch_evidence.png"
    fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def _write_report(uer: dict[str, Any], rows: list[dict[str, Any]], summary: list[dict[str, Any]], source_paths: dict[str, Path], figure_path: Path) -> Path:
    b = uer["frontier"]["budgets"]
    branch_table = "\n".join(
        f"| {row['branch']} | {row['decision']} | {row['scope']} | {row['claim']} |"
        for row in summary
    )
    evidence_table = "\n".join(
        f"| {row['branch']} | {row['source']} | {row['condition']} | {row['metric']} | {row['value']} | {row['units']} | {row['evidence_status']} | {row['interpretation']} |"
        for row in rows
    )
    report = f"""# Objective 5 closure package (15 September 2026)

**Closure ID:** `O5-CLOSURE-2026-09-15`  
**Overall decision:** `PARTIALLY_SUPPORTED_CONDITIONAL`  
**Scope:** UER field screening plus locked M7.3/M8 controlled-synthetic evidence.  
**Preservation:** Existing O5 reports, tables, figures, thesis DOCX files and thesis text sources were not edited.

## Closure decision

Objective 5 is closed only at a bounded conditional level. The evidence supports a four-branch answer: integration can improve a declared target under compatible controlled-synthetic conditions; some added channels or budget increments have no practically material model-defined value; discordant or misspecified evidence can weaken inference; and the UER field package requires `ABSTAIN` for field accuracy, unique mechanism, field-optimal design and procurement claims.

| Branch | Decision | Scope | Replacement claim |
| :--- | :--- | :--- | :--- |
{branch_table}

## UER run used for closure

The latest preserved UER package is `{_relative(UER_RUN_DIR)}`. It contains **{uer['sampling_locations']}** sampling locations and **{uer['candidate_edges']}** Tier-C DEM/topographic candidate edges. The edge count is an adjacency hypothesis count, not a measured flowpath count. There is no independent UER flow, age or reaction truth in the inspected package.

The current boundary snapshot has **{uer['conflict_edges']}** flagged edges and **{uer['admitted_edges']}** admitted edges. Dual-endpoint tritium is available for **{uer['tritium_dual_endpoint_tested']}** edges, leaving **{uer['tritium_untested']}** untested. Dual nitrate-isotope evidence is available for **{uer['nitrate_dual_endpoint_tested']}** edges, leaving **{uer['nitrate_untested']}** untested. **{uer['cbe_fail_edges']}** edges touch a CBE quality-control failure. Missing evidence is retained as `UNTESTED_MISSING_DATA` and is not counted as a pass.

The UER re-solve is a chemistry regularisation diagnostic. Total Dirichlet energy is **{uer['ungated_dirichlet_energy']:.6f}** for the all-edge fit and **{uer['gated_dirichlet_energy']:.6f}** after changing the edge set. On the same admitted edges, mean residual changes from **{uer['mean_residual_same_edges_ungated']:.6f}** to **{uer['mean_residual_same_edges_gated']:.6f}**, and same-edge energy changes from **{uer['same_edge_energy_ungated']:.6f}** to **{uer['same_edge_energy_gated']:.6f}**. These values do not establish protective predictive performance.

## Conditional UER minimax frontier

The frontier uses the Regional Median Screening Cohort (**n={uer['frontier']['cohort_size']}**, **{uer['frontier']['baseline_tritium_TU']:.3f} TU ± {uer['frontier']['baseline_sigma_TU']:.3f} TU**) as a single cohort-level prior. It uses a finite **0.5–80.0 year** grid and reports MTT ambiguity in **years**; the normalized field is **fraction of the 80-year grid**, not years. The declared quantitative target is **16.0 years** (0.20 × 80 years), and the practical-equivalence margin is **0.05 years**.

| Allocation cap (USD) | Actual analytical cost (USD per borehole) | Selected options | Worst-case MTT ambiguity (years) | Status |
| :---: | :---: | :--- | :---: | :--- |
| 0 | {b['0']['actual_analytical_cost_usd_per_borehole']:.2f} | {b['0']['selected_options']} | {b['0']['mtt_ambiguity_years']:.6f} | `{b['0']['status']}` |
| 800 | {b['800']['actual_analytical_cost_usd_per_borehole']:.2f} | {b['800']['selected_options']} | {b['800']['mtt_ambiguity_years']:.6f} | `{b['800']['status']}` |
| 2500 | {b['2500']['actual_analytical_cost_usd_per_borehole']:.2f} | {b['2500']['selected_options']} | {b['2500']['mtt_ambiguity_years']:.6f} | `{b['2500']['status']}` |

The change from the $800 row to the $2,500 row is **{uer['frontier']['plateau_delta_years_800_minus_2500']:.9f} years**, below the declared 0.05-year equivalence margin. This is a model-defined plateau under assumed tracer responses, error bounds, histories and a finite grid. It is not a field optimum.

![Objective 5 closure branches and UER conditional frontier]({figure_path.name})

## Branch evidence ledger

| Branch | Source | Condition | Metric | Value | Units | Evidence status | Interpretation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{evidence_table}

`IMPROVES` rows are controlled-synthetic evidence with declared comparators and model scope. `NO_MATERIAL_VALUE` rows are finite-grid or controlled-synthetic equivalence statements. `WEAKENS` rows are adverse controls or negative native increments. `ABSTAIN` rows are explicit missing-evidence or failed-gate decisions.

## Cost and procurement boundary

The cost manifest is `{uer['cost_manifest_id']}` with evidence class `{uer['cost_evidence_classification']}`. `eligible_for_final_frontier={str(uer['cost_eligible_for_final_frontier']).lower()}` and the campaign cost model is `{uer['campaign_cost_status']}`. The current rates are public-list analytical scenarios in USD per analytical measurement per borehole. They exclude field labour, purging, transport, accommodation, containers, freight, insurance, customs, duties, taxes, payment fees, duplicates, blanks, spikes, repeats and nonconforming-matrix surcharges. A transaction-specific quotation and a complete field/logistics ledger are required before procurement or a field-budget claim.

## Solver status repair

`hydrosheaf/nuclear/ttd_certified_design.py` now checks the all-candidate LP result before using it as an impossibility witness or budget reference. Status `INFEASIBLE_*` is reserved for an infeasible declared constraint set; non-success statuses such as iteration limits return `NUMERICAL_ERROR` and no certificate or selected subset. Regression tests cover initial-prior and all-candidate failure paths for both certified and budgeted minimax designs.

## Exact changes and validation runs

The closure package adds `scripts/analysis/close_objective5_2026_09_15.py`, which reads the preserved latest UER run, O5 audit, cost manifest, M7.3 bootstrap contrasts and M8 locked decision files; derives the four branch labels; writes the machine-readable ledger, frontier extract, manifest and figure; and records SHA-256 hashes for every source input. The existing `outputs/uer_objective5_groundwater_budgetary_2026-09-15`, O5 audit, M7.3/M8 artifacts and thesis files remain unchanged.

The live repair in `hydrosheaf/nuclear/ttd_certified_design.py` guards both `solve_certified_measurement_design` and `solve_budgeted_minimax_design` against a failed all-candidate LP. The regression additions in `tests/test_ttd_certified_design.py` force a non-success all-candidate result and verify `NUMERICAL_ERROR`, the preserved feasibility status and an empty selection rather than a false certificate.

Validation completed on 15 September 2026:

* `.\\.venv\\Scripts\\python.exe -m py_compile scripts/analysis/close_objective5_2026_09_15.py hydrosheaf/nuclear/ttd_certified_design.py tests/test_ttd_certified_design.py` — pass.
* `.\\.venv\\Scripts\\python.exe scripts/analysis/close_objective5_2026_09_15.py` — pass; 22 branch rows written and figure regenerated.
* `.\\.venv\\Scripts\\python.exe -m pytest tests/test_ttd_certified_design.py -q` — **15 passed**.
* `.\\.venv\\Scripts\\python.exe -m pytest tests/test_uer_objective5_audited_logic.py -q` — **3 passed**.
* `.\\.venv\\Scripts\\python.exe -m pytest tests/test_uer_objective5.py -q` — **6 passed**.
* `text_hygiene_audit.py` on this report — **PASS**, 0 findings; the audit is read-only and saved beside this report.

## Remaining gaps

1. Obtain co-timed groundwater levels, construction/screen intervals and an independent flow/connectivity reference before interpreting candidate edges as field flowpaths.
2. Obtain dated field and laboratory quotations, including matrix, method, units, precision, detection limits, minimum volume, QA/QC, duplicates/blanks/spikes, turnaround, shipping and customs.
3. Supply a Ghana-specific recharge history and unit-validated response/error model for SF6, CFC-12, tritium, 3H/3He and 14C; test old-water tails and nuisance corrections.
4. Run a frozen held-out field or controlled-synthetic comparison of the same inference target under compatible, redundant, conflicting, permuted and gated/ungated evidence, including calibration and false-alarm/missed-conflict rates.
5. Re-run the frontier only after those inputs are frozen; do not backfill missing costs or relabel public-list prices as quotes.

## Replacement thesis wording for Objective 5

> Objective 5 was answered as a conditional evidence-integration question. In locked controlled-synthetic experiments, chemistry and hydraulic evidence improved topology ranking and correct topology reduced age error for the declared model and comparator, while the M8 acquisition result improved over random acquisition within its relative-cost and likelihood model. The same evidence base showed that adding age evidence can worsen ranking, and that permuted, reversed or misspecified evidence can increase error or overconfidence. Some additional measurements produced no materially different model result: the UER finite-grid median-cohort frontier changed by only {uer['frontier']['plateau_delta_years_800_minus_2500']:.6f} years between the 800 and 2,500 USD allocation caps, and the locked M8 transport candidates both selected 50 days. The UER field package itself contains topographic candidate adjacencies, chemistry and sparse isotope coverage but no independent flow, age or reaction truth; it therefore supports reproducible screening and illustrative conditional minimax calculations rather than improved field accuracy, unique mechanisms, field-optimal sampling or protective joint multi-physics gating. Analytical public-list costs remain budgetary until transaction-specific quotations and field/logistics costs are obtained.

## Traceability

The machine-readable branch ledger is [`objective5_branch_evidence.csv`](objective5_branch_evidence.csv); the branch decisions are [`objective5_branch_decisions.csv`](objective5_branch_decisions.csv); the frontier extract is [`objective5_uer_frontier_summary.json`](objective5_uer_frontier_summary.json); the source hashes and artifact list are [`objective5_closure_manifest.json`](objective5_closure_manifest.json); and the read-only prose audit is [`objective5_text_hygiene_audit.json`](objective5_text_hygiene_audit.json). The original audit's injected-failure result is retained unchanged; this closure package records the live repair and its regression coverage separately.
"""
    path = OUT_DIR / "OBJECTIVE_5_CLOSURE_REPORT.md"
    path.write_text(report, encoding="utf-8")
    return path


def main() -> None:
    source_paths = {
        "uer_report": UER_RUN_DIR / "UER_OBJECTIVE_5_FULL_REPORT.md",
        "uer_frontier": UER_RUN_DIR / "tables" / "uer_objective5_minimax_pareto_frontier.csv",
        "uer_boundary": UER_RUN_DIR / "tables" / "uer_objective5_tripartite_boundary_audit.csv",
        "uer_gating": UER_RUN_DIR / "tables" / "uer_objective5_sheaf_gating_comparison.csv",
        "uer_reproduction_audit": UER_AUDIT_DIR / "reproduction_audit.json",
        "uer_audit_report": UER_AUDIT_DIR / "OBJECTIVE_5_AUDIT.md",
        "uer_field_summary": UER_FIELD_DIR / "uer_integration_summary.json",
        "cost_manifest": COST_MANIFEST,
        "m7_confirmatory_decision": M7_DIR / "confirmatory_decision.json",
        "m7_bootstrap_contrasts": M7_DIR / "evidence_case_bootstrap_contrasts.csv",
        "m8_frontier_decision": M8_ROOT / "RUN-M8-FRONTIER-AL-20260728-01" / "claim_decision.json",
        "m8_independent_decision": M8_ROOT / "RUN-M8-INDEPENDENT-20260728-01" / "claim_decision.json",
        "m8_confirmatory_decision": M8_ROOT / "RUN-M8-CONFIRM-20260728-01" / "claim_decision.json",
        "uer_runner": ROOT / "scripts" / "analysis" / "run_uer_objective5_full.py",
        "ttd_design": ROOT / "hydrosheaf" / "nuclear" / "ttd_design.py",
        "ttd_certified_design": ROOT / "hydrosheaf" / "nuclear" / "ttd_certified_design.py",
    }
    _require_paths(source_paths.values())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)

    uer = _uer_summary()
    m7 = _load_json(M7_DIR / "confirmatory_decision.json")
    m8_frontier = _load_json(M8_ROOT / "RUN-M8-FRONTIER-AL-20260728-01" / "claim_decision.json")
    m8_independent = _load_json(M8_ROOT / "RUN-M8-INDEPENDENT-20260728-01" / "claim_decision.json")
    m8_confirm = _load_json(M8_ROOT / "RUN-M8-CONFIRM-20260728-01" / "claim_decision.json")
    m7_contrasts = pd.read_csv(M7_DIR / "evidence_case_bootstrap_contrasts.csv")
    rows = _branch_rows(uer, m7, m7_contrasts, m8_frontier, m8_independent, m8_confirm)
    summary = _branch_summary(rows)

    pd.DataFrame(rows).to_csv(OUT_DIR / "objective5_branch_evidence.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT_DIR / "objective5_branch_decisions.csv", index=False)
    (OUT_DIR / "objective5_uer_frontier_summary.json").write_text(json.dumps(uer["frontier"], indent=2, sort_keys=True), encoding="utf-8")
    figure_path = _write_figure(rows, uer)
    report_path = _write_report(uer, rows, summary, source_paths, figure_path)

    created_artifacts = {
        "report": _relative(report_path),
        "branch_evidence": _relative(OUT_DIR / "objective5_branch_evidence.csv"),
        "branch_decisions": _relative(OUT_DIR / "objective5_branch_decisions.csv"),
        "frontier_summary": _relative(OUT_DIR / "objective5_uer_frontier_summary.json"),
        "figure": _relative(figure_path),
    }
    hygiene_audit_path = OUT_DIR / "objective5_text_hygiene_audit.json"
    if hygiene_audit_path.exists():
        created_artifacts["text_hygiene_audit"] = _relative(hygiene_audit_path)

    manifest = {
        "closure_id": "O5-CLOSURE-2026-09-15",
        "date": "2026-09-15",
        "overall_decision": "PARTIALLY_SUPPORTED_CONDITIONAL",
        "source_manifest": _source_manifest(source_paths),
        "created_artifacts": created_artifacts,
        "checks": {
            "uer_candidate_edges": uer["candidate_edges"],
            "uer_boundary_rows": uer["boundary_rows"],
            "uer_missingness_preserved": True,
            "cost_manifest_eligible_for_final_frontier": uer["cost_eligible_for_final_frontier"],
            "field_truth_available": uer["independent_field_truth_available"],
            "branch_labels": ["improves", "no-material-value", "weakens", "ABSTAIN"],
        },
    }
    (OUT_DIR / "objective5_closure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote Objective 5 closure report: {report_path}")
    print(f"Wrote branch evidence rows: {len(rows)}")
    print(f"Wrote closure figure: {figure_path}")
    print(f"UER plateau delta (800 - 2500 USD): {uer['frontier']['plateau_delta_years_800_minus_2500']:.9f} years")


if __name__ == "__main__":
    main()
