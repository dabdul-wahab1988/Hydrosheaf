"""Create Q1 manuscript tables matching figure_Table_guide.txt specification.

Produces 2 main-text tables + 10 supplementary tables from M4 benchmark results.
"""

from __future__ import annotations

import importlib.metadata
import platform
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np
import pandas as pd

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_ROOT / "results"
TABLE_DIR = BENCHMARK_ROOT / "tables" / "Manuscript_Ready"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

PUBLIC_DIR = RESULT_DIR / "public_archives" / "savage"

EVIDENCE_LADDER = {
    "negative_random": {
        "evidence_level": 0,
        "graph_scenario": "Negative control (random edges)",
        "evidence_source": "Random edge generation",
        "independent_validation": False,
        "main_use": "Negative-control baseline for topology metrics",
        "allowed_claim": "Random edges serve as a negative-control baseline establishing the floor for precision.",
        "required_guardrail": "Random edges are not a proposed graph model; they establish the floor for precision.",
    },
    "negative_wrong_direction": {
        "evidence_level": 0,
        "graph_scenario": "Negative control (wrong direction)",
        "evidence_source": "Reversed reference edges",
        "independent_validation": False,
        "main_use": "Directionality sensitivity diagnostic",
        "allowed_claim": "Reversing reference edge direction quantifies directionality sensitivity.",
        "required_guardrail": "Wrong-direction graph is a diagnostic negative control, not a Hydrosheaf inference.",
    },
    "negative_shortcut": {
        "evidence_level": 0,
        "graph_scenario": "Negative control (shortcut edges) — INAPPLICABLE on this reference",
        "evidence_source": "Two-hop paths in the reference graph (none exist here)",
        "independent_validation": False,
        "main_use": "Would test sensitivity to skipped intermediate nodes",
        "allowed_claim": "NONE. The Savage reference is a fan-in graph whose sources and receptors are disjoint, so no two-hop path exists and the shortcut edge set is empty. F1 = 0.000 with 0 inferred edges is an inapplicable control, not a rejected topology hypothesis.",
        "required_guardrail": "Must not be reported as a failed or rejected control. Use negative_misrouted_sink, which is well-posed on fan-in reference topologies.",
    },
    "sink_aware_baseline": {
        "evidence_level": 0,
        "graph_scenario": "Sink-aware structural baseline (informed control)",
        "evidence_source": "The 3-element reference receptor set only; no hydraulic, spatial or MODPATH connectivity information",
        "independent_validation": False,
        "main_use": "Establishes the informed performance floor for a fan-in reference topology",
        "allowed_claim": "Knowledge of the receptor set alone recovers the reference topology at precision 0.382, recall 1.000, F1 0.552. Any inference scenario must be judged against this floor, not only against the uninformed spatial/random/wrong-direction controls.",
        "required_guardrail": "This is a structural baseline, not a Hydrosheaf inference. It is informed by the reference receptor set and is therefore not an independent method; it bounds how much credit hydraulic directionality can claim.",
    },
    "negative_misrouted_sink": {
        "evidence_level": 0,
        "graph_scenario": "Negative control (misrouted receptor)",
        "evidence_source": "Reference edges with each source reassigned to a different reference receptor",
        "independent_validation": False,
        "main_use": "Receptor-attribution sensitivity at constant graph size and fan-in shape",
        "allowed_claim": "Reassigning sources to the wrong receptor collapses F1 to 0.138, showing that the benchmark penalises wrong receptor attribution and is not merely rewarding fan-in shape or graph size.",
        "required_guardrail": "Misrouted-sink graph is a diagnostic negative control, not a Hydrosheaf inference. It replaces the two-hop shortcut control, which is inapplicable to fan-in reference topologies.",
    },
    "spatial_only": {
        "evidence_level": 0,
        "graph_scenario": "Spatial proximity only (geometry-only control)",
        "evidence_source": "Haversine distance, nearest-neighbour graph; no hydraulic, stratigraphic, or MODPATH data",
        "independent_validation": False,
        "main_use": "Geometry-only control demonstrating that spatial proximity alone cannot recover directed groundwater topology",
        "allowed_claim": "Spatial proximity alone completely fails to recover directed groundwater topology (F1 = 0.0), establishing the geometry-only floor and demonstrating that hydraulic evidence is necessary.",
        "required_guardrail": "Spatial-only graph uses no hydraulic head, hydrostratigraphic, or MODPATH connectivity information. Its F1 = 0.0 result is a control outcome, not evidence of any inference skill.",
    },
    "sparse_node": {
        "evidence_level": "S",
        "graph_scenario": "Sparse-node robustness (sensitivity analysis — not a scenario)",
        "evidence_source": "Random node subsampling (10–100% of nodes) using head-gradient inference rule; mean F1 over 20 trials",
        "independent_validation": False,
        "main_use": "Quantify robustness of head-gradient topology recovery to reduced node density; this is a sensitivity diagnostic, not an independent inference scenario",
        "allowed_claim": "Node-sparsity sensitivity analysis quantifies how head-gradient topology recovery degrades as available nodes decrease from 100% to 10%; reported as mean F1 across 20 random subsampling trials.",
        "required_guardrail": "Sparse-node is a sensitivity diagnostic, NOT an independent inference scenario. It applies the same hydraulic inference rule as head_gradient but with fewer nodes. Its aggregate F1 cannot be directly compared with single-run scenario metrics. See Supp Table S9 for full sensitivity curve.",
    },
    "head_gradient": {
        "evidence_level": 3,
        "graph_scenario": "Head-gradient constrained",
        "evidence_source": "Elevation-as-head proxy, downhill flow direction",
        "independent_validation": True,
        "main_use": "Primary independent topology validation",
        "allowed_claim": "Elevation-as-head proxy constrains flow direction and improves topology recovery.",
        "required_guardrail": "Elevation is a proxy for hydraulic head, not a MODFLOW-simulated head value.",
    },
    "head_gradient_bayesian_hodge": {
        "evidence_level": 3,
        "graph_scenario": "Head-gradient Hodge-pruned",
        "evidence_source": "Elevation-as-head proxy with Hodge decomposition edge pruning",
        "independent_validation": True,
        "main_use": "Diagnostic: Hodge decomposition applied post-inference to remove harmonic/curl artefacts",
        "allowed_claim": "Hodge decomposition prunes spurious edges after head-gradient constrained inference, diagnosing harmonic-flow artefacts.",
        "required_guardrail": "Hodge pruning is a diagnostic post-processing step, not an independent edge-selection mechanism; it produces identical topology metrics to head_gradient on this benchmark.",
    },
    "real_head_projected_gradient": {
        "evidence_level": 3,
        "graph_scenario": "Real-head projected gradient",
        "evidence_source": "MODFLOW-simulated head projected onto steepest-descent gradient",
        "independent_validation": True,
        "main_use": "Diagnostic: MODFLOW head substituted into steepest-descent framework to isolate head-source effect",
        "allowed_claim": "Substituting MODFLOW-simulated head into the steepest-descent framework isolates the contribution of head-source fidelity to topology recovery.",
        "required_guardrail": "Real-head projected gradient uses MODFLOW head that is part of the same model generating the MODPATH reference; it is diagnostic for head-source sensitivity, not independent validation.",
    },
    "head_depth": {
        "evidence_level": 4,
        "graph_scenario": "Head-depth constrained",
        "evidence_source": "Depth-tiered elevation, stratified candidate edges",
        "independent_validation": True,
        "main_use": "Depth-stratified edge filtering",
        "allowed_claim": "Depth-tiered elevation stratifies candidate edges and permits more neighbors.",
        "required_guardrail": "Depth tiers are based on elevation terciles, not formal hydrostratigraphic units.",
    },
    "hydrostratigraphic": {
        "evidence_level": 5,
        "graph_scenario": "Hydrostratigraphic constrained",
        "evidence_source": "Depth-based aquifer unit filtering",
        "independent_validation": True,
        "main_use": "Same-aquifer edge filtering",
        "allowed_claim": "Same-aquifer edge filtering removes cross-unit false positives.",
        "required_guardrail": "Aquifer units are depth-based splits from elevation, not lithostratigraphic or model-layer assignments.",
    },
    "modpath_prior_override": {
        "evidence_level": 6,
        "graph_scenario": "MODPATH prior override",
        "evidence_source": "MODPATH reference edges as priors",
        "independent_validation": False,
        "main_use": "Prior-assisted graph construction (override mode)",
        "allowed_claim": "MODPATH-informed priors are useful for prior-assisted topology construction.",
        "required_guardrail": "MODPATH-informed prior rows must never be reported as independent validation evidence.",
    },
    "modpath_prior_merge": {
        "evidence_level": 6,
        "graph_scenario": "MODPATH prior merge",
        "evidence_source": "MODPATH reference edges merged with Hydrosheaf edges",
        "independent_validation": False,
        "main_use": "Prior-assisted graph construction (merge mode)",
        "allowed_claim": "MODPATH-informed priors are useful for prior-assisted topology construction.",
        "required_guardrail": "MODPATH-informed prior rows must never be reported as independent validation evidence.",
    },
    "modpath_prior_only": {
        "evidence_level": 6,
        "graph_scenario": "MODPATH prior only",
        "evidence_source": "MODPATH reference edges only",
        "independent_validation": False,
        "main_use": "Prior-assisted graph construction (only-prior mode)",
        "allowed_claim": "MODPATH-informed priors are useful for prior-assisted topology construction.",
        "required_guardrail": "MODPATH-informed prior rows must never be reported as independent validation evidence.",
    },
    "external_archive_savage": {
        "evidence_level": 5,
        "graph_scenario": "Public archive endpoint-pathline projection (USGS Savage)",
        "evidence_source": "Public USGS Savage MODFLOW/MODPATH archive outputs",
        "independent_validation": False,
        "main_use": "Public-archive ingestion and endpoint-pathline projection diagnostic",
        "allowed_claim": "The workflow can ingest public MODPATH archive outputs and project endpoint/pathline evidence into graph-validation form.",
        "required_guardrail": "External archive rows report endpoint-pathline projection consistency from the same MODPATH archive; they are not independent Hydrosheaf validation, capture-zone polygon validation, field truth, or independent travel-time prediction.",
    },
}


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


# ---------------------------------------------------------------------------
# MAIN TABLE 1: Evidence ladder and validation contract (concise)
# ---------------------------------------------------------------------------

def make_main_table1_evidence_ladder() -> None:
    perf = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")

    rows: List[Dict[str, Any]] = []
    for scenario_key, info in EVIDENCE_LADDER.items():
        row = {
            "evidence_level": info["evidence_level"],
            "graph_scenario": info["graph_scenario"],
            "evidence_source": info["evidence_source"],
            "independent_validation": info["independent_validation"],
            "main_use": info["main_use"],
            "allowed_claim": info["allowed_claim"],
            "required_guardrail": info["required_guardrail"],
        }
        # Add F1 from results if available
        if "modpath_prior" in scenario_key and not priors.empty:
            mode_key = scenario_key.replace("modpath_prior_", "")
            prior_row = priors[priors["prior_mode"] == mode_key]
            if not prior_row.empty:
                row["f1"] = "prior-assisted (see Supp S4)"
            else:
                row["f1"] = np.nan
        elif scenario_key == "external_archive_savage":
            row["f1"] = "self-consistency (see Supp S5)"
        elif not perf.empty:
            perf_row = perf[perf["scenario"] == scenario_key]
            if not perf_row.empty:
                f1 = perf_row["f1"].iloc[0]
                row["f1"] = round(float(f1), 3) if pd.notna(f1) else np.nan
            else:
                row["f1"] = np.nan
        else:
            row["f1"] = np.nan
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "Main_Table1_Evidence_Ladder_Validation_Contract.csv", index=False)


# ---------------------------------------------------------------------------
# MAIN TABLE 2: Topology-validation performance and failure summary (concise)
# ---------------------------------------------------------------------------

def make_main_table2_performance_summary() -> None:
    perf = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")

    if perf.empty:
        return

    ind = perf.copy()
    # sparse_node is a sensitivity analysis — excluded from the scenario comparison rows.
    # spatial_only is a geometry-only control (level 0, F1 = 0.0).
    scenario_ordered = ["spatial_only", "sink_aware_baseline",
                        "head_gradient", "head_gradient_bayesian_hodge",
                        "real_head_projected_gradient", "head_depth", "hydrostratigraphic",
                        "negative_random", "negative_wrong_direction",
                        "negative_misrouted_sink", "negative_shortcut"]
    ind["_order"] = ind["scenario"].map({k: i for i, k in enumerate(scenario_ordered)})
    ind = ind.sort_values("_order")
    ind = ind[ind["scenario"].isin(scenario_ordered)]

    rows = []
    for _, row in ind.iterrows():
        scenario = str(row["scenario"])
        level = EVIDENCE_LADDER.get(scenario, {}).get("evidence_level", "—")
        f1 = float(row["f1"]) if pd.notna(row["f1"]) else np.nan

        if f1 >= 0.95:
            interpretation = "Near-exact topology reproduction"
        elif f1 >= 0.5:
            interpretation = "Partial topology skill"
        elif f1 > 0.0:
            interpretation = "Diagnostic failure mode"
        else:
            interpretation = "No topology recovery (control)"

        if scenario in ("negative_random", "negative_wrong_direction", "negative_shortcut",
                        "negative_misrouted_sink", "spatial_only"):
            interpretation = "No topology recovery (control)"

        # A control that emitted no edges is inapplicable to this reference
        # topology, not a rejected hypothesis; say so instead of reporting a
        # bare F1 = 0.000 that reads as a passed falsification test.
        if str(row.get("control_applicable", "")).lower() == "false":
            interpretation = (
                "INAPPLICABLE — no edges could be constructed on this reference "
                "topology (sources and receptors are disjoint, so no two-hop path "
                "exists); not evidence of rejection. See 'Negative: misrouted sink'"
            )
        if scenario == "sink_aware_baseline":
            interpretation = (
                "Informed structural floor — uses only the 3-receptor set, no "
                "hydraulic information; every inference scenario must be judged "
                "against this row"
            )
        if scenario == "negative_misrouted_sink":
            interpretation = (
                "Wrong receptor attribution at constant graph size (well-posed "
                "replacement for the inapplicable shortcut control)"
            )

        rows.append({
            "scenario": _display_scenario(scenario),
            "evidence_level": level,
            "validation_mode": row.get("validation_mode", ""),
            "independent_validation": row.get("independent_validation", ""),
            "precision": round(float(row["precision"]), 3) if pd.notna(row.get("precision")) else np.nan,
            "recall": round(float(row["recall"]), 3) if pd.notna(row.get("recall")) else np.nan,
            "f1": round(f1, 3) if np.isfinite(f1) else np.nan,
            "false_positives": int(row["fp"]) if pd.notna(row.get("fp")) else np.nan,
            "false_negatives": int(row["fn"]) if pd.notna(row.get("fn")) else np.nan,
            "scale_mismatch": str(row.get("scale_mismatch", "")).lower() == "true",
            "interpretation": interpretation,
        })

    # Append sparse_node as a clearly labelled sensitivity-analysis row (not a scenario peer)
    sparse_rows = perf[perf["scenario"] == "sparse_node"]
    if not sparse_rows.empty:
        sr = sparse_rows.iloc[0]
        rows.append({
            "scenario": "Sparse node (sensitivity analysis — see Supp Table S9)",
            "evidence_level": "S",
            "validation_mode": "sensitivity_analysis",
            "independent_validation": False,
            "precision": "—",
            "recall": "—",
            "f1": round(float(sr["f1"]), 3) if pd.notna(sr.get("f1")) else np.nan,
            "false_positives": "—",
            "false_negatives": "—",
            "scale_mismatch": False,
            "interpretation": "Sensitivity diagnostic: mean F1 across 20 node-subsampling trials (10–100% nodes); not comparable with single-run scenario metrics",
        })

    # Add prior-assisted rows
    if not priors.empty:
        for _, row in priors.iterrows():
            scenario = f"modpath_prior_{row['prior_mode']}"
            rows.append({
                "scenario": _display_scenario(scenario),
                "evidence_level": 6,
                "validation_mode": "prior_assisted",
                "independent_validation": False,
                "precision": np.nan,
                "recall": np.nan,
                "f1": "prior-assisted",
                "false_positives": np.nan,
                "false_negatives": np.nan,
                "scale_mismatch": False,
                "interpretation": "Prior-assisted (NOT independent validation)",
            })

    pd.DataFrame(rows).to_csv(
        TABLE_DIR / "Main_Table2_Topology_Performance_Failure_Summary.csv", index=False
    )


def _display_scenario(value: str) -> str:
    labels = {
        "spatial_only": "Spatial only",
        "head_gradient": "Head gradient",
        "head_gradient_bayesian_hodge": "Hodge pruned",
        "real_head_projected_gradient": "Proj. gradient",
        "head_depth": "Head depth",
        "hydrostratigraphic": "Hydrostratigraphic",
        "sparse_node": "Sparse node",
        "negative_random": "Negative: random",
        "negative_wrong_direction": "Negative: wrong direction",
        "negative_shortcut": "Negative: shortcut",
        "sink_aware_baseline": "Sink-aware baseline",
        "negative_misrouted_sink": "Negative: misrouted sink",
        "modpath_prior_override": "Prior: override",
        "modpath_prior_merge": "Prior: merge",
        "modpath_prior_only": "Prior: only",
    }
    return labels.get(value, value.replace("_", " "))


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S1: Full benchmark design and evidence contract
# ---------------------------------------------------------------------------

def make_supp_table_s1_benchmark_design() -> None:
    rows = []
    for scenario_key, info in EVIDENCE_LADDER.items():
        rows.append({
            "evidence_level": info["evidence_level"],
            "graph_scenario": info["graph_scenario"],
            "evidence_source": info["evidence_source"],
            "independent_validation": info["independent_validation"],
            "main_use": info["main_use"],
            "allowed_claim": info["allowed_claim"],
            "required_guardrail": info["required_guardrail"],
        })
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / "Supp_TableS1_Full_Benchmark_Design_Evidence_Contract.csv", index=False
    )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S2: Independent topology metrics
# ---------------------------------------------------------------------------

def make_supp_table_s2_independent_metrics() -> None:
    metrics = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    if metrics.empty:
        return
    metrics = metrics[metrics["validation_mode"] == "independent_graph_inference"].copy()
    metrics = metrics[metrics["scenario"].isin(EVIDENCE_LADDER.keys())].copy()
    columns = [
        "scenario", "validation_mode", "independent_validation", "result_class", "precision", "recall", "f1",
        "false_positive_rate", "false_negative_rate", "tp", "fp", "fn", "tn",
        "scale_mismatch", "median_reference_length", "median_inferred_length",
        "allowed_claim", "required_guardrail",
    ]
    table = metrics[[c for c in columns if c in metrics.columns]].copy()
    for c in table.columns:
        if pd.api.types.is_numeric_dtype(table[c]):
            table[c] = table[c].round(4)
    table.to_csv(TABLE_DIR / "Supp_TableS2_Independent_Topology_Metrics.csv", index=False)


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S3: Edge failure diagnostics
# ---------------------------------------------------------------------------

def make_supp_table_s3_edge_failure_diagnostics() -> None:
    edges = _read_csv(RESULT_DIR / "edge_classification.csv")
    metrics = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    if edges.empty:
        return

    grouped: List[Dict[str, object]] = []
    for scenario, subset in edges.groupby("scenario", sort=False):
        edge_col = "edge_id" if "edge_id" in subset.columns else "edge"
        labels = {
            label: _join_edges(subset.loc[subset["classification"] == label, edge_col])
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
        grouped.append({
            "scenario": scenario,
            "true_positive_edges": labels["TP"],
            "false_positive_edges": labels["FP"],
            "false_negative_edges": labels["FN"],
            "true_negative_edges": labels["TN"],
            "failure_mode": ";".join(failure_notes),
            "interpretation_limit": "Hydrogeologic cause requires model metadata; this table only classifies graph evidence.",
        })
    pd.DataFrame(grouped).to_csv(
        TABLE_DIR / "Supp_TableS3_Edge_Failure_Diagnostics.csv", index=False
    )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S4: MODPATH prior-mode audit
# ---------------------------------------------------------------------------

def make_supp_table_s4_prior_mode_audit() -> None:
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")
    if priors.empty:
        return
    table = priors.copy()
    table["allowed_use"] = "Prior-assisted Hydrosheaf configuration or sensitivity analysis."
    table["prohibited_claim"] = "Independent topology validation."
    table["independent_validation"] = False
    _round_numeric(table).to_csv(TABLE_DIR / "Supp_TableS4_MODPATH_Prior_Mode_Audit.csv", index=False)


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S5: External MODPATH archive validation
# ---------------------------------------------------------------------------

def make_supp_table_s5_external_archive_validation() -> None:
    external = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    if external.empty:
        return
    columns = [
        "validation_tier", "archive_name", "source_doi",
        "n_endpoint_records", "n_pathline_points", "n_particles",
        "n_endpoint_edges", "n_pathline_edges",
        "true_positive_edges", "false_positive_edges", "false_negative_edges",
        "edge_precision", "edge_recall", "edge_f1",
        "direction_agreement_rate", "mean_source_receptor_overlap",
        "endpoint_projection_preservation_rate",
        "median_compressed_path_cells",
        "travel_time_spearman_rho", "travel_time_kendall_tau",
        "travel_time_rank_supported",
        "capture_envelope_iou", "source_cell_jaccard",
        "harmonized_travel_time_spearman_rho", "harmonized_travel_time_kendall_tau",
        "harmonized_travel_time_median_abs_log10_difference",
        "travel_time_rank_interpretation", "claim_guardrail",
    ]
    table = external[[c for c in columns if c in external.columns]].copy()
    table["evidence_scope"] = "Public MODPATH endpoint-pathline projection diagnostic; not independent Hydrosheaf validation."
    # A tier whose archive was never ingested must not sit in a validation table
    # looking like a processed result merely because it carries a live DOI.
    # Flag processing status explicitly as the first column after the tier id.
    processed = pd.to_numeric(
        table.get("n_particles", pd.Series(0, index=table.index)), errors="coerce"
    ).fillna(0) > 0
    table["processing_status"] = np.where(
        processed,
        "PROCESSED",
        "NOT PROCESSED - fallback stub; archive was not ingested and contributes "
        "no validation evidence. Reported for scope transparency only.",
    )
    cols = list(table.columns)
    cols.insert(1, cols.pop(cols.index("processing_status")))
    table = table[cols]
    _round_numeric(table).to_csv(
        TABLE_DIR / "Supp_TableS5_External_MODPATH_Archive_Validation.csv", index=False
    )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S6: External pathline and time diagnostics
# ---------------------------------------------------------------------------

def make_supp_table_s6_external_pathline_time_diagnostics() -> None:
    archive = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    rank = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv")
    boot = _read_csv(RESULT_DIR / "external_modpath_bootstrap_ci.csv")
    rows: List[Dict[str, object]] = []

    if not archive.empty:
        top = archive.iloc[0]
        rows.extend([
            {
                "diagnostic": "compressed_path_cells_median",
                "estimate": top.get("median_compressed_path_cells"),
                "ci_low": np.nan, "ci_high": np.nan,
                "n": top.get("n_particles"),
                "evidence_scope": "full raw pathline sequence complexity",
                "guardrail": "Compressed cell sequence is not a capture-zone polygon.",
            },
            {
                "diagnostic": "endpoint_projection_preservation_rate",
                "estimate": top.get("endpoint_projection_preservation_rate"),
                "ci_low": np.nan, "ci_high": np.nan,
                "n": top.get("n_particles"),
                "evidence_scope": "particle-level start/end preservation",
                "guardrail": "Projection preserves source and receptor cells, not every intermediate cell.",
            },
        ])
    if not rank.empty:
        top = rank.iloc[0]
        rows.extend([
            {
                "diagnostic": "travel_time_spearman_rho",
                "estimate": top.get("spearman_rho"),
                "ci_low": np.nan, "ci_high": np.nan,
                "n": top.get("n_edges"),
                "evidence_scope": "edge-level travel-time rank diagnostic",
                "guardrail": top.get("guardrail"),
            },
            {
                "diagnostic": "travel_time_kendall_tau",
                "estimate": top.get("kendall_tau"),
                "ci_low": np.nan, "ci_high": np.nan,
                "n": top.get("n_edges"),
                "evidence_scope": "edge-level travel-time rank diagnostic",
                "guardrail": top.get("guardrail"),
            },
        ])
    if not boot.empty:
        for _, row in boot.iterrows():
            rows.append({
                "diagnostic": row.get("metric"),
                "estimate": row.get("estimate"),
                "ci_low": row.get("ci_low"),
                "ci_high": row.get("ci_high"),
                "n": row.get("n"),
                "evidence_scope": row.get("method"),
                "guardrail": "Bootstrap interval describes the analysed archive only.",
            })

    if rows:
        _round_numeric(pd.DataFrame(rows)).to_csv(
            TABLE_DIR / "Supp_TableS6_External_Pathline_Time_Diagnostics.csv", index=False
        )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S7: Capture-envelope and harmonised travel-time
# ---------------------------------------------------------------------------

def make_supp_table_s7_capture_travel_time_validation() -> None:
    capture = _read_csv(RESULT_DIR / "external_modpath_capture_envelope_summary.csv")
    time_summary = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv")
    rows: List[Dict[str, object]] = []

    if not capture.empty:
        top = capture.iloc[0]
        rows.extend([
            {
                "metric": "particle_weighted_source_start_envelope_iou",
                "estimate": top.get("particle_weighted_capture_envelope_iou"),
                "n": top.get("n_receptor_targets"),
                "evidence_scope": "endpoint/pathline source-start point-cloud envelope",
                "guardrail": top.get("guardrail"),
            },
            {
                "metric": "median_source_cell_jaccard",
                "estimate": top.get("median_source_cell_jaccard"),
                "n": top.get("n_receptor_targets"),
                "evidence_scope": "source-cell set overlap by receptor target",
                "guardrail": top.get("guardrail"),
            },
        ])
    if not time_summary.empty:
        top = time_summary.iloc[0]
        rows.extend([
            {
                "metric": "harmonized_travel_time_spearman_rho",
                "estimate": top.get("spearman_rho"),
                "n": top.get("n_edges"),
                "evidence_scope": "endpoint total time vs MODPATH-derived Hydrosheaf edge weights",
                "guardrail": top.get("guardrail"),
            },
            {
                "metric": "harmonized_travel_time_kendall_tau",
                "estimate": top.get("kendall_tau"),
                "n": top.get("n_edges"),
                "evidence_scope": "endpoint total time vs MODPATH-derived Hydrosheaf edge weights",
                "guardrail": top.get("guardrail"),
            },
            {
                "metric": "harmonized_travel_time_median_abs_log10_difference",
                "estimate": top.get("median_abs_log10_difference"),
                "n": top.get("n_edges"),
                "evidence_scope": "scale agreement for endpoint total time and MODPATH-derived edge weights",
                "guardrail": top.get("guardrail"),
            },
        ])

    if rows:
        _round_numeric(pd.DataFrame(rows)).to_csv(
            TABLE_DIR / "Supp_TableS7_Capture_Travel_Time_Validation.csv", index=False
        )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S8: MODPATH endpoint scenario metrics
# ---------------------------------------------------------------------------

def make_supp_table_s8_endpoint_scenarios() -> None:
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
    table["required_guardrail"] = (
        "This is endpoint transition agreement, not pathline-shape or capture-zone validation."
    )
    _round_numeric(table).to_csv(
        TABLE_DIR / "Supp_TableS8_MODPATH_Endpoint_Scenarios.csv", index=False
    )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S9: Node-sparsity sensitivity
# ---------------------------------------------------------------------------

def make_supp_table_s9_sparsity_sensitivity() -> None:
    sensitivity = _read_csv(RESULT_DIR / "m4_sparsity_sensitivity.csv")
    if sensitivity.empty:
        return
    table = sensitivity.copy()
    table["evidence_scope"] = "node-sparsity robustness"
    table["required_guardrail"] = (
        "Subsampling sensitivity is a controlled robustness diagnostic, not a field uncertainty model."
    )

    def _survivorship_note(row: pd.Series) -> str:
        successful = int(row["successful_trials"])
        planned = int(row["planned_trials"])
        if successful < planned:
            return (
                f"Only {successful}/{planned} trials met the minimum reference-edge threshold "
                f"(≥5 edges); mean F1 may be upward-biased."
            )
        return "All trials succeeded."

    table["survivorship_note"] = table.apply(_survivorship_note, axis=1)
    table["bias_warning"] = table["successful_trials"] < table["planned_trials"]

    _round_numeric(table).to_csv(
        TABLE_DIR / "Supp_TableS9_Sparsity_Sensitivity.csv", index=False
    )


# ---------------------------------------------------------------------------
# SUPPLEMENTARY TABLE S10: Reproducibility
# ---------------------------------------------------------------------------

def _package_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def make_supp_table_s10_reproducibility() -> None:
    rows = [
        {"component": "Python", "version": platform.python_version(), "role": "analysis runtime"},
        {"component": "hydrosheaf", "version": _package_version("hydrosheaf"), "role": "validation API"},
        {"component": "numpy", "version": _package_version("numpy"), "role": "numeric summaries"},
        {"component": "pandas", "version": _package_version("pandas"), "role": "table generation"},
        {"component": "matplotlib", "version": _package_version("matplotlib"), "role": "figure generation"},
        {"component": "scipy", "version": _package_version("scipy"), "role": "convex hull, statistics"},
        {"component": "networkx", "version": _package_version("networkx"), "role": "graph support"},
        {"component": "flopy", "version": _package_version("flopy"), "role": "MODPATH readers"},
    ]
    pd.DataFrame(rows).to_csv(TABLE_DIR / "Supp_TableS10_Reproducibility.csv", index=False)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print("Generating main-text tables...")
    make_main_table1_evidence_ladder()
    make_main_table2_performance_summary()

    print("Generating supplementary tables...")
    make_supp_table_s1_benchmark_design()
    make_supp_table_s2_independent_metrics()
    make_supp_table_s3_edge_failure_diagnostics()
    make_supp_table_s4_prior_mode_audit()
    make_supp_table_s5_external_archive_validation()
    make_supp_table_s6_external_pathline_time_diagnostics()
    make_supp_table_s7_capture_travel_time_validation()
    make_supp_table_s8_endpoint_scenarios()
    make_supp_table_s9_sparsity_sensitivity()
    make_supp_table_s10_reproducibility()

    print(f"All tables written to {TABLE_DIR}")


if __name__ == "__main__":
    main()
