"""Create M4 manuscript-ready analysis tables from benchmark outputs."""

from __future__ import annotations

import importlib.metadata
import platform
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_ROOT / "results"
TABLE_DIR = BENCHMARK_ROOT / "tables" / "Manuscript_Ready"
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _round_numeric(df: pd.DataFrame, digits: int = 3) -> pd.DataFrame:
    out = df.copy()
    for column in out.columns:
        if pd.api.types.is_numeric_dtype(out[column]):
            out[column] = out[column].round(digits)
    return out


def _join_edges(edges: Iterable[str]) -> str:
    values = [str(edge) for edge in edges if str(edge)]
    return ";".join(values) if values else "none"


def make_table1_benchmark_design() -> None:
    rows = [
        {
            "evidence_block": "controlled_topology_scenarios",
            "source_file": "results/independent_graph_vs_modpath.csv",
            "independent_validation": True,
            "current_status": "implemented",
            "main_use": "Primary controlled test of directed-edge agreement and scale mismatch.",
            "required_guardrail": "Supports reduced-order topology reproduction only for the analysed benchmark scenarios.",
        },
        {
            "evidence_block": "edge_failure_diagnostics",
            "source_file": "results/edge_classification.csv",
            "independent_validation": True,
            "current_status": "implemented",
            "main_use": "Trace true-positive, false-positive, false-negative, and true-negative edges.",
            "required_guardrail": "Do not convert edge labels into hydrogeologic causes without supporting domain metadata.",
        },
        {
            "evidence_block": "endpoint_based_modpath_validation",
            "source_file": "results/m4_topology_benchmark_summary.csv",
            "independent_validation": True,
            "current_status": "implemented when M2 MODPATH endpoint files are present",
            "main_use": "External endpoint-derived check against selected MODPATH output.",
            "required_guardrail": "Endpoint transitions are not full pathline-shape validation.",
        },
        {
            "evidence_block": "external_usgs_savage_modpath_archive",
            "source_file": "results/external_modpath_archive_summary.csv",
            "independent_validation": True,
            "current_status": "implemented from M2 external MODPATH validation outputs",
            "main_use": "End-to-end public archive validation of endpoint/pathline topology, direction, source-receptor overlap, and travel-time rank.",
            "required_guardrail": "Validates MODPATH topology ingestion and projection, not field geochemistry or capture-zone polygons.",
        },
        {
            "evidence_block": "sparse_node_sensitivity",
            "source_file": "results/m4_sparsity_sensitivity.csv",
            "independent_validation": True,
            "current_status": "implemented when endpoint validation has run",
            "main_use": "Quantify how graph skill changes as available nodes are reduced.",
            "required_guardrail": "Random subsampling sensitivity does not represent all field-data missingness mechanisms.",
        },
        {
            "evidence_block": "modpath_informed_priors",
            "source_file": "results/modpath_informed_priors.csv",
            "independent_validation": False,
            "current_status": "implemented",
            "main_use": "Document how MODPATH connectivity can enter Hydrosheaf as prior information.",
            "required_guardrail": "Never report prior-mode rows as independent validation evidence.",
        },
        {
            "evidence_block": "future_external_archives",
            "source_file": "not yet generated",
            "independent_validation": False,
            "current_status": "future work",
            "main_use": "Great Miami, Long Island, and coastal archives should remain external-validation candidates.",
            "required_guardrail": "Do not claim multi-archive validation until processed outputs exist.",
        },
    ]
    pd.DataFrame(rows).to_csv(TABLE_DIR / "Manuscript_Table1_M4_Benchmark_Design.csv", index=False)


def make_table2_independent_metrics() -> None:
    metrics = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    if metrics.empty:
        return
    columns = [
        "scenario",
        "validation_mode",
        "precision",
        "recall",
        "f1",
        "false_positive_rate",
        "false_negative_rate",
        "tp",
        "fp",
        "fn",
        "tn",
        "scale_mismatch",
        "median_reference_length",
        "median_inferred_length",
    ]
    table = metrics[[column for column in columns if column in metrics.columns]].copy()
    table["evidence_class"] = np.select(
        [
            (table["f1"] >= 0.95) & (table["fp"] == 0) & (table["fn"] == 0),
            table["f1"] >= 0.80,
            table["f1"] > 0,
        ],
        ["near_exact_topology", "partial_topology_skill", "diagnostic_failure_mode"],
        default="failed_topology_recovery",
    )
    table["manuscript_guardrail"] = np.where(
        table.get("scale_mismatch", False).astype(str).str.lower().eq("true"),
        "Report as scale-mismatch diagnostic, not only as edge-count performance.",
        "Report with independent-validation scope.",
    )
    _round_numeric(table).to_csv(TABLE_DIR / "Manuscript_Table2_Independent_Topology_Metrics.csv", index=False)


def make_table3_edge_failure_diagnostics() -> None:
    edges = _read_csv(RESULT_DIR / "edge_classification.csv")
    metrics = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    if edges.empty:
        return

    grouped: List[Dict[str, object]] = []
    for scenario, subset in edges.groupby("scenario", sort=False):
        labels = {
            label: _join_edges(subset.loc[subset["classification"] == label, "edge"])
            for label in ("TP", "FP", "FN", "TN")
        }
        metric_row = metrics.loc[metrics["scenario"] == scenario].head(1) if not metrics.empty else pd.DataFrame()
        scale = False
        f1 = np.nan
        if not metric_row.empty:
            scale = str(metric_row.iloc[0].get("scale_mismatch", False)).lower() == "true"
            f1 = float(metric_row.iloc[0].get("f1", np.nan))
        failure_notes = []
        if labels["FP"] != "none":
            failure_notes.append("false_positive_connections")
        if labels["FN"] != "none":
            failure_notes.append("false_negative_connections")
        if scale:
            failure_notes.append("scale_mismatch")
        if not failure_notes and np.isfinite(f1) and f1 >= 0.95:
            failure_notes.append("near_exact_reproduction")
        grouped.append(
            {
                "scenario": scenario,
                "true_positive_edges": labels["TP"],
                "false_positive_edges": labels["FP"],
                "false_negative_edges": labels["FN"],
                "true_negative_edges": labels["TN"],
                "failure_mode": ";".join(failure_notes),
                "interpretation_limit": "Hydrogeologic cause requires model metadata; this table only classifies graph evidence.",
            }
        )
    pd.DataFrame(grouped).to_csv(TABLE_DIR / "Manuscript_Table3_Edge_Failure_Diagnostics.csv", index=False)


def make_table4_prior_mode_audit() -> None:
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")
    if priors.empty:
        return
    table = priors.copy()
    table["allowed_use"] = "Prior-assisted Hydrosheaf configuration or sensitivity analysis."
    table["prohibited_claim"] = "Independent topology validation."
    table["independent_validation"] = False
    _round_numeric(table).to_csv(TABLE_DIR / "Manuscript_Table4_MODPATH_Prior_Mode_Audit.csv", index=False)


def make_table5_external_archive_validation() -> None:
    external = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    if external.empty:
        return
    columns = [
        "archive_name",
        "source_doi",
        "n_endpoint_records",
        "n_pathline_points",
        "n_particles",
        "n_endpoint_edges",
        "n_pathline_edges",
        "true_positive_edges",
        "false_positive_edges",
        "false_negative_edges",
        "edge_precision",
        "edge_recall",
        "edge_f1",
        "direction_agreement_rate",
        "mean_source_receptor_overlap",
        "endpoint_projection_preservation_rate",
        "median_compressed_path_cells",
        "travel_time_spearman_rho",
        "travel_time_kendall_tau",
        "travel_time_rank_supported",
        "capture_envelope_iou",
        "source_cell_jaccard",
        "harmonized_travel_time_spearman_rho",
        "harmonized_travel_time_kendall_tau",
        "harmonized_travel_time_median_abs_log10_difference",
        "travel_time_rank_interpretation",
        "claim_guardrail",
    ]
    table = external[[col for col in columns if col in external.columns]].copy()
    _round_numeric(table).to_csv(TABLE_DIR / "Manuscript_Table5_External_MODPATH_Archive_Validation.csv", index=False)


def make_table6_external_pathline_time_diagnostics() -> None:
    archive = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    rank = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv")
    boot = _read_csv(RESULT_DIR / "external_modpath_bootstrap_ci.csv")
    rows: list[dict[str, object]] = []
    if not archive.empty:
        top = archive.iloc[0]
        rows.extend(
            [
                {
                    "diagnostic": "compressed_path_cells_median",
                    "estimate": top.get("median_compressed_path_cells"),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n": top.get("n_particles"),
                    "evidence_scope": "full raw pathline sequence complexity",
                    "guardrail": "Compressed cell sequence is not a capture-zone polygon.",
                },
                {
                    "diagnostic": "endpoint_projection_preservation_rate",
                    "estimate": top.get("endpoint_projection_preservation_rate"),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n": top.get("n_particles"),
                    "evidence_scope": "particle-level start/end preservation",
                    "guardrail": "Projection preserves source and receptor cells, not every intermediate cell.",
                },
            ]
        )
    if not rank.empty:
        top = rank.iloc[0]
        rows.extend(
            [
                {
                    "diagnostic": "travel_time_spearman_rho",
                    "estimate": top.get("spearman_rho"),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n": top.get("n_edges"),
                    "evidence_scope": "edge-level travel-time rank diagnostic",
                    "guardrail": top.get("guardrail"),
                },
                {
                    "diagnostic": "travel_time_kendall_tau",
                    "estimate": top.get("kendall_tau"),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "n": top.get("n_edges"),
                    "evidence_scope": "edge-level travel-time rank diagnostic",
                    "guardrail": top.get("guardrail"),
                },
            ]
        )
    if not boot.empty:
        for _, row in boot.iterrows():
            rows.append(
                {
                    "diagnostic": row.get("metric"),
                    "estimate": row.get("estimate"),
                    "ci_low": row.get("ci_low"),
                    "ci_high": row.get("ci_high"),
                    "n": row.get("n"),
                    "evidence_scope": row.get("method"),
                    "guardrail": "Bootstrap interval describes the analysed archive only.",
                }
            )
    if rows:
        _round_numeric(pd.DataFrame(rows)).to_csv(
            TABLE_DIR / "Manuscript_Table6_External_Pathline_Time_Diagnostics.csv",
            index=False,
        )


def make_table7_capture_travel_time_validation() -> None:
    capture = _read_csv(RESULT_DIR / "external_modpath_capture_envelope_summary.csv")
    time = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv")
    rows: list[dict[str, object]] = []
    if not capture.empty:
        top = capture.iloc[0]
        rows.extend(
            [
                {
                    "metric": "particle_weighted_capture_envelope_iou",
                    "estimate": top.get("particle_weighted_capture_envelope_iou"),
                    "n": top.get("n_receptor_targets"),
                    "evidence_scope": "endpoint/pathline point-cloud capture envelope",
                    "guardrail": top.get("guardrail"),
                },
                {
                    "metric": "median_source_cell_jaccard",
                    "estimate": top.get("median_source_cell_jaccard"),
                    "n": top.get("n_receptor_targets"),
                    "evidence_scope": "source-cell set overlap by receptor target",
                    "guardrail": top.get("guardrail"),
                },
            ]
        )
    if not time.empty:
        top = time.iloc[0]
        rows.extend(
            [
                {
                    "metric": "harmonized_travel_time_spearman_rho",
                    "estimate": top.get("spearman_rho"),
                    "n": top.get("n_edges"),
                    "evidence_scope": "endpoint total time versus MODPATH-derived Hydrosheaf edge weights",
                    "guardrail": top.get("guardrail"),
                },
                {
                    "metric": "harmonized_travel_time_kendall_tau",
                    "estimate": top.get("kendall_tau"),
                    "n": top.get("n_edges"),
                    "evidence_scope": "endpoint total time versus MODPATH-derived Hydrosheaf edge weights",
                    "guardrail": top.get("guardrail"),
                },
                {
                    "metric": "harmonized_travel_time_median_abs_log10_difference",
                    "estimate": top.get("median_abs_log10_difference"),
                    "n": top.get("n_edges"),
                    "evidence_scope": "scale agreement for endpoint total time and MODPATH-derived Hydrosheaf edge weights",
                    "guardrail": top.get("guardrail"),
                },
            ]
        )
    if rows:
        _round_numeric(pd.DataFrame(rows)).to_csv(
            TABLE_DIR / "Manuscript_Table7_Capture_Travel_Time_Validation.csv",
            index=False,
        )


def make_table8_endpoint_scenarios() -> None:
    endpoint = _read_csv(RESULT_DIR / "m4_topology_benchmark_summary.csv")
    if endpoint.empty:
        return
    table = endpoint.copy()
    if "false_positive_rate" in table.columns:
        table["false_positive_rate_note"] = np.where(
            table["false_positive_rate"].isna(),
            "not estimated because a full candidate-edge universe was not defined",
            "computed",
        )
    table["evidence_scope"] = "endpoint-derived MODPATH connectivity"
    table["required_guardrail"] = "This is endpoint transition agreement, not pathline-shape or capture-zone validation."
    _round_numeric(table).to_csv(TABLE_DIR / "Manuscript_Table8_MODPATH_Endpoint_Scenarios.csv", index=False)


def make_table9_sparsity_sensitivity() -> None:
    sensitivity = _read_csv(RESULT_DIR / "m4_sparsity_sensitivity.csv")
    if sensitivity.empty:
        return
    table = sensitivity.copy()
    table["evidence_scope"] = "node-sparsity robustness"
    table["required_guardrail"] = "Subsampling sensitivity is a controlled robustness diagnostic, not a field uncertainty model."
    _round_numeric(table).to_csv(TABLE_DIR / "Manuscript_Table9_Sparsity_Sensitivity.csv", index=False)


def _package_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def make_table7_reproducibility() -> None:
    rows = [
        {"component": "Python", "version": platform.python_version(), "role": "analysis runtime"},
        {"component": "hydrosheaf", "version": _package_version("hydrosheaf"), "role": "validation API"},
        {"component": "numpy", "version": _package_version("numpy"), "role": "numeric summaries"},
        {"component": "pandas", "version": _package_version("pandas"), "role": "table generation"},
        {"component": "matplotlib", "version": _package_version("matplotlib"), "role": "figure generation"},
        {"component": "networkx", "version": _package_version("networkx"), "role": "graph support"},
        {"component": "flopy", "version": _package_version("flopy"), "role": "MODPATH readers"},
    ]
    pd.DataFrame(rows).to_csv(TABLE_DIR / "Manuscript_Table10_Reproducibility.csv", index=False)


def main() -> None:
    make_table1_benchmark_design()
    make_table2_independent_metrics()
    make_table3_edge_failure_diagnostics()
    make_table4_prior_mode_audit()
    make_table5_external_archive_validation()
    make_table6_external_pathline_time_diagnostics()
    make_table7_capture_travel_time_validation()
    make_table8_endpoint_scenarios()
    make_table9_sparsity_sensitivity()
    make_table7_reproducibility()
    print(f"Wrote M4 manuscript-ready tables to {TABLE_DIR}")


if __name__ == "__main__":
    main()
