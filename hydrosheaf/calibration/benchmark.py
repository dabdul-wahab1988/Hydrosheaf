"""
Assumption-calibration benchmark workflow.

Compares baseline Hydrosheaf topology against assumption-calibrated topology
on the same held-out validation labels, producing defensible manuscript claims.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..config import Config as HConfig
from ..graph.types import Edge
from ..log import get_logger
from ..sheaf.topology_refine import refine_edges_with_sheaf
from .adapters import TopologyCalibrationAdapter, TopologyCalibrationObservation
from .benchmark_bootstrap import bootstrap_benchmark_metrics
from .factory import setup_topology_adapter
from .validation_workflow import (
    _independence_status,
    _load_topology_observations,
    _extract_assumption_params_from_optimal,
    _extract_config_only_thresholds,
    _resolve_path,
)

logger = get_logger("calibration.benchmark")


# ── variant evaluation ───────────────────────────────────────────────

def _evaluate_variant(
    samples: Optional[List[Dict[str, Any]]],
    candidates: List[Edge],
    cfg: HConfig,
    val_observations: List[TopologyCalibrationObservation],
    variant_name: str,
) -> Dict[str, Any]:
    """Run sheaf refinement for one config variant and compute metrics
    against the validation label set only."""
    if samples is not None:
        selected = refine_edges_with_sheaf(samples, candidates, cfg)
    else:
        # Fallback: threshold on edge_confidence without sheaf refinement
        selected = [
            e for e in candidates
            if float((e.attrs or {}).get("edge_confidence", 0.5)) >= 0.5
        ]
    selected_ids = {str(e.edge_id) for e in selected}

    tp = fp = tn = fn = 0
    for obs in val_observations:
        is_selected = obs.edge_id in selected_ids
        if obs.observed_present >= 0.5:
            if is_selected:
                tp += 1
            else:
                fn += 1
        else:
            if is_selected:
                fp += 1
            else:
                tn += 1

    total = tp + fp + tn + fn
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if precision + recall else 0.0)
    accuracy = (tp + tn) / total if total > 0 else float("nan")

    return {
        "variant": variant_name,
        "n_selected_edges": len(selected),
        "n_evaluated_validation_labels": len(val_observations),
        "selected_edge_ids": sorted(list(selected_ids)),
        "confusion_matrix": {
            "true_positives": tp,
            "false_positives": fp,
            "true_negatives": tn,
            "false_negatives": fn,
        },
        "metrics": {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "accuracy": accuracy,
        },
    }


# ── main benchmark runner ────────────────────────────────────────────

def run_assumption_benchmark(
    config: Any,
    calibration_result: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compare baseline, null-model-defaults, and assumption-calibrated
    topology variants against held-out validation labels.

    Parameters
    ----------
    config : CalibrationConfig
        Topology calibration config with ``validation_observations_file``
        in ``adapter_settings``.
    calibration_result : dict, optional
        Pre-computed calibration result. When provided, calibration is
        skipped and its optimal parameters are used for the
        assumption-calibrated variant.

    Returns
    -------
    dict
        Benchmark report (also written to
        ``assumption_benchmark_results.json`` and
        ``assumption_benchmark_summary.csv`` in ``output_dir``).

    Raises
    ------
    ValueError
        If validation_observations_file is missing, or calibration and
        validation labels resolve to the same file.
    """
    settings = config.adapter_settings

    # ── 1. Load data and enforce independence ──────────────────────────
    cal_obs_file = (
        settings.get("observations_file") or config.observations_file
    )
    val_obs_file = settings.get("validation_observations_file")
    if not val_obs_file:
        raise ValueError(
            "Benchmark requires model.validation_observations_file."
        )
    if not cal_obs_file:
        raise ValueError(
            "Benchmark requires observations_file (calibration labels)."
        )

    cal_path = _resolve_path(cal_obs_file)
    val_path = _resolve_path(val_obs_file)
    if cal_path == val_path:
        raise ValueError(
            "Calibration and validation label files must be different "
            "for an independent benchmark. "
            f"Both resolved to: {cal_path}"
        )

    cal_observations = _load_topology_observations(cal_obs_file)
    val_observations = _load_topology_observations(val_obs_file)
    cal_edge_ids = {obs.edge_id for obs in cal_observations}
    val_edge_ids = {obs.edge_id for obs in val_observations}
    overlapping = cal_edge_ids & val_edge_ids
    independent, independence_reason = _independence_status(
        cal_obs_file,
        val_obs_file,
        settings,
        overlapping,
    )

    logger.info(
        "Benchmark setup — cal labels: %d, val labels: %d, overlap: %d",
        len(cal_observations), len(val_observations), len(overlapping),
    )

    if not val_observations:
        raise ValueError("No validation observations loaded.")

    # ── 2. Build adapter to get candidates, samples, and base config ───
    problem = setup_topology_adapter(config)
    candidates = list(problem.candidate_edges)
    samples = problem.samples
    assumption_params = settings.get("assumption_params") or []

    # Build base config from the adapter's config
    if problem.config is not None and isinstance(problem.config, HConfig):
        base_cfg = replace(problem.config)
    else:
        base_cfg = HConfig()

    # ── 3. Resolve optimal params if calibration result is provided ────
    optimal_params: Dict[str, float] = {}
    has_calibration = False
    if calibration_result is not None and calibration_result.get(
        "optimal_parameters"
    ):
        optimal_params = calibration_result["optimal_parameters"]
        has_calibration = True
        logger.info(
            "Reusing calibration result for assumption-calibrated variant."
        )
    elif calibration_result is not None:
        logger.warning(
            "calibration_result provided but optimal_parameters is empty "
            "— assumption_calibrated variant will use defaults."
        )
    else:
        logger.warning(
            "No calibration_result provided — assumption_calibrated variant "
            "uses default settings. manuscript_claim_allowed will be false."
        )

    calibrated_assumption_params = _extract_assumption_params_from_optimal(
        optimal_params, assumption_params,
    )
    config_only_thresholds = _extract_config_only_thresholds(base_cfg)

    # ── 4. Build three config variants and evaluate ────────────────────

    # Variant A: baseline — default config, no null/evidence ladder
    baseline_cfg = replace(base_cfg)
    baseline_cfg.null_model_enabled = False
    baseline_cfg.evidence_ladder_enabled = False
    baseline_cfg.assumption_calibration_enabled = False
    baseline_candidates = list(candidates)

    # Variant B: null-model defaults — null/evidence enabled, default params
    null_default_cfg = replace(base_cfg)
    null_default_cfg.null_model_enabled = True
    null_default_cfg.evidence_ladder_enabled = True
    null_default_cfg.assumption_calibration_enabled = True
    null_default_candidates = list(candidates)

    # Variant C: assumption-calibrated — null/evidence enabled,
    # calibrated assumption params, calibrated edge confidences
    calibrated_cfg = replace(base_cfg)
    calibrated_cfg.null_model_enabled = True
    calibrated_cfg.evidence_ladder_enabled = True
    calibrated_cfg.assumption_calibration_enabled = True
    for a_name in assumption_params:
        if a_name in optimal_params and hasattr(calibrated_cfg, a_name):
            setattr(calibrated_cfg, a_name, float(optimal_params[a_name]))

    # Apply calibrated edge confidences to variant C candidates
    calibrated_candidates: List[Edge] = []
    for edge in candidates:
        p_name = TopologyCalibrationAdapter._param_name(edge.edge_id)
        if p_name in optimal_params:
            prob = TopologyCalibrationAdapter._sigmoid(
                float(optimal_params[p_name])
            )
            new_attrs = dict(edge.attrs or {})
            new_attrs["edge_confidence"] = prob
            calibrated_candidates.append(
                Edge(
                    edge_id=edge.edge_id,
                    u=edge.u,
                    v=edge.v,
                    attrs=new_attrs,
                )
            )
        else:
            calibrated_candidates.append(edge)

    variant_results = [
        _evaluate_variant(
            samples, baseline_candidates, baseline_cfg,
            val_observations, "baseline",
        ),
        _evaluate_variant(
            samples, null_default_candidates, null_default_cfg,
            val_observations, "null_model_defaults",
        ),
        _evaluate_variant(
            samples, calibrated_candidates, calibrated_cfg,
            val_observations, "assumption_calibrated",
        ),
    ]

    # ── 5. Improvement summary ─────────────────────────────────────────
    baseline_metrics = variant_results[0]["metrics"]
    calibrated_metrics = variant_results[2]["metrics"]

    def _delta(key: str) -> float:
        return round(float(calibrated_metrics.get(key, 0.0))
                     - float(baseline_metrics.get(key, 0.0)), 6)

    improvement_summary = {
        "delta_precision": _delta("precision"),
        "delta_recall": _delta("recall"),
        "delta_f1": _delta("f1"),
        "delta_accuracy": _delta("accuracy"),
    }

    # ── 5a. Bootstrap uncertainty (optional) ──────────────────────────
    uncertainty: Optional[Dict[str, Any]] = None
    bootstrap_enabled = bool(settings.get("benchmark_bootstrap", False))
    if bootstrap_enabled:
        n_boot = int(settings.get("benchmark_bootstrap_n", 1000))
        if n_boot <= 0:
            raise ValueError(
                f"benchmark_bootstrap_n must be > 0, got {n_boot}"
            )
        boot_seed = int(settings.get("benchmark_bootstrap_seed", 123))
        selected_ids_by_variant = {
            vr["variant"]: set(vr["selected_edge_ids"])
            for vr in variant_results
        }
        boot_result = bootstrap_benchmark_metrics(
            selected_ids_by_variant=selected_ids_by_variant,
            val_observations=val_observations,
            n_boot=n_boot,
            seed=boot_seed,
        )
        uncertainty = {
            "n_boot": boot_result["n_boot"],
            "seed": boot_result["seed"],
            "method": "percentile_bootstrap_on_validation_resamples",
            "variant_cis": {
                vname: {
                    "precision": list(ci.precision),
                    "recall": list(ci.recall),
                    "f1": list(ci.f1),
                    "accuracy": list(ci.accuracy),
                }
                for vname, ci in boot_result["variant_cis"].items()
            },
            "improvement_summary_ci": {
                "delta_precision": list(
                    boot_result["delta_cis"].delta_precision
                ),
                "delta_recall": list(
                    boot_result["delta_cis"].delta_recall
                ),
                "delta_f1": list(
                    boot_result["delta_cis"].delta_f1
                ),
                "delta_accuracy": list(
                    boot_result["delta_cis"].delta_accuracy
                ),
                "probability_delta_f1_gt_0": (
                    boot_result["delta_cis"].probability_delta_f1_gt_0
                ),
                "probability_delta_precision_gt_0": (
                    boot_result["delta_cis"].probability_delta_precision_gt_0
                ),
            },
        }
        if "bootstrap_warning" in boot_result:
            uncertainty["bootstrap_warning"] = boot_result["bootstrap_warning"]
        logger.info(
            "Bootstrap complete — n_boot=%d, delta_f1_ci95=[%.4f, %.4f]",
            n_boot,
            uncertainty["improvement_summary_ci"]["delta_f1"][0],
            uncertainty["improvement_summary_ci"]["delta_f1"][1],
        )

    # ── 6. Build report ────────────────────────────────────────────────
    report = {
        "calibration_label_file": cal_obs_file,
        "validation_label_file": val_obs_file,
        "independent_validation": independent,
        "independence_reason": independence_reason,
        "manuscript_claim_allowed": independent and has_calibration,
        "n_calibration_labels": len(cal_observations),
        "n_validation_labels": len(val_observations),
        "n_overlapping_edge_ids": len(overlapping),
        "calibrated_assumption_parameters": calibrated_assumption_params,
        "fixed_config_only_thresholds": config_only_thresholds,
        "variants": {
            vr["variant"]: {
                "n_selected_edges": vr["n_selected_edges"],
                "n_evaluated_validation_labels": vr["n_evaluated_validation_labels"],
                "confusion_matrix": vr["confusion_matrix"],
                "metrics": vr["metrics"],
                "selected_edge_ids": vr["selected_edge_ids"],
            }
            for vr in variant_results
        },
        "improvement_summary": improvement_summary,
    }
    if uncertainty is not None:
        report["uncertainty"] = uncertainty

    # ── 7. Write output files ──────────────────────────────────────────
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = out_dir / "assumption_benchmark_results.json"
    import numpy as np

    def _convert(obj: Any) -> Any:
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        raise TypeError

    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=_convert)
    logger.info("Benchmark JSON written to %s", json_path)

    # CSV summary
    csv_path = out_dir / "assumption_benchmark_summary.csv"
    _write_benchmark_csv(csv_path, variant_results, improvement_summary, uncertainty)
    logger.info("Benchmark CSV written to %s", csv_path)

    # Optional Markdown report
    md_path = out_dir / "assumption_benchmark_report.md"
    _write_benchmark_md(md_path, report)
    logger.info("Benchmark Markdown report written to %s", md_path)

    return report


# ── output helpers ────────────────────────────────────────────────────

def _write_benchmark_csv(
    path: Path,
    variant_results: List[Dict[str, Any]],
    improvement_summary: Dict[str, float],
    uncertainty: Optional[Dict[str, Any]] = None,
) -> None:
    """Write a flat CSV summary of benchmark results.

    When ``uncertainty`` is provided, CI columns are appended after the
    point-estimate metric columns.
    """
    ci_columns = [
        "precision_low95", "precision_high95",
        "recall_low95", "recall_high95",
        "f1_low95", "f1_high95",
        "accuracy_low95", "accuracy_high95",
    ]
    has_ci = uncertainty is not None

    rows: List[Dict[str, Any]] = []
    for vr in variant_results:
        m = vr["metrics"]
        cm = vr["confusion_matrix"]
        row: Dict[str, Any] = {
            "variant": vr["variant"],
            "n_selected_edges": vr["n_selected_edges"],
            "n_evaluated_validation_labels": vr["n_evaluated_validation_labels"],
            "precision": m["precision"],
            "recall": m["recall"],
            "f1": m["f1"],
            "accuracy": m["accuracy"],
            "tp": cm["true_positives"],
            "fp": cm["false_positives"],
            "tn": cm["true_negatives"],
            "fn": cm["false_negatives"],
        }
        if has_ci:
            vcis = uncertainty["variant_cis"].get(vr["variant"], {})
            row["precision_low95"] = vcis.get("precision", [None, None])[0]
            row["precision_high95"] = vcis.get("precision", [None, None])[1]
            row["recall_low95"] = vcis.get("recall", [None, None])[0]
            row["recall_high95"] = vcis.get("recall", [None, None])[1]
            row["f1_low95"] = vcis.get("f1", [None, None])[0]
            row["f1_high95"] = vcis.get("f1", [None, None])[1]
            row["accuracy_low95"] = vcis.get("accuracy", [None, None])[0]
            row["accuracy_high95"] = vcis.get("accuracy", [None, None])[1]
        else:
            for c in ci_columns:
                row[c] = ""
        rows.append(row)

    # Add improvement row
    imp_row: Dict[str, Any] = {
        "variant": "improvement_vs_baseline",
        "n_selected_edges": "",
        "n_evaluated_validation_labels": "",
        "precision": improvement_summary["delta_precision"],
        "recall": improvement_summary["delta_recall"],
        "f1": improvement_summary["delta_f1"],
        "accuracy": improvement_summary["delta_accuracy"],
        "tp": "", "fp": "", "tn": "", "fn": "",
    }
    if has_ci:
        delta_ci = uncertainty.get("improvement_summary_ci", {})
        imp_row["precision_low95"] = delta_ci.get("delta_precision", [None, None])[0]
        imp_row["precision_high95"] = delta_ci.get("delta_precision", [None, None])[1]
        imp_row["recall_low95"] = delta_ci.get("delta_recall", [None, None])[0]
        imp_row["recall_high95"] = delta_ci.get("delta_recall", [None, None])[1]
        imp_row["f1_low95"] = delta_ci.get("delta_f1", [None, None])[0]
        imp_row["f1_high95"] = delta_ci.get("delta_f1", [None, None])[1]
        imp_row["accuracy_low95"] = delta_ci.get("delta_accuracy", [None, None])[0]
        imp_row["accuracy_high95"] = delta_ci.get("delta_accuracy", [None, None])[1]
    else:
        for c in ci_columns:
            imp_row[c] = ""
    rows.append(imp_row)

    fieldnames = [
        "variant", "n_selected_edges", "n_evaluated_validation_labels",
        "precision", "recall", "f1", "accuracy",
    ] + ci_columns + [
        "tp", "fp", "tn", "fn",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_benchmark_md(path: Path, report: Dict[str, Any]) -> None:
    """Write a human-readable Markdown benchmark report."""
    variants = report["variants"]
    imp = report["improvement_summary"]

    def _fmt(v: float) -> str:
        return f"{v:.4f}"

    lines = [
        "# Assumption Calibration Benchmark Report",
        "",
        f"- **Independent validation**: {report['independent_validation']}",
        f"- **Manuscript claim allowed**: {report['manuscript_claim_allowed']}",
        f"- Calibration labels: {report['n_calibration_labels']}",
        f"- Validation labels: {report['n_validation_labels']}",
        f"- Overlapping edge IDs: {report['n_overlapping_edge_ids']}",
        "",
        "## Variant Metrics",
        "",
        "| Variant | Edges | Precision | Recall | F1 | Accuracy | TP | FP | TN | FN |",
        "|---------|-------|-----------|--------|-----|----------|----|----|----|----|",
    ]

    for vname, vdata in variants.items():
        m = vdata["metrics"]
        cm = vdata["confusion_matrix"]
        lines.append(
            f"| {vname} | {vdata['n_selected_edges']} "
            f"| {_fmt(m['precision'])} | {_fmt(m['recall'])} "
            f"| {_fmt(m['f1'])} | {_fmt(m['accuracy'])} "
            f"| {cm['true_positives']} | {cm['false_positives']} "
            f"| {cm['true_negatives']} | {cm['false_negatives']} |"
        )

    lines += [
        "",
        "## Improvement (Assumption-Calibrated vs. Baseline)",
        "",
        f"| Metric | Delta |",
        f"|--------|-------|",
        f"| Precision | {_fmt(imp['delta_precision'])} |",
        f"| Recall | {_fmt(imp['delta_recall'])} |",
        f"| F1 | {_fmt(imp['delta_f1'])} |",
        f"| Accuracy | {_fmt(imp['delta_accuracy'])} |",
        "",
    ]

    if report.get("calibrated_assumption_parameters"):
        lines += [
            "## Calibrated Assumption Parameters",
            "",
        ]
        for k, v in report["calibrated_assumption_parameters"].items():
            lines.append(f"- **{k}**: {v:.6f}")

    # ── Uncertainty section (when bootstrap is enabled) ────────────
    uncertainty = report.get("uncertainty")
    if uncertainty:
        lines += [
            "",
            "## Uncertainty",
            "",
            f"Percentile bootstrap (n={uncertainty['n_boot']}, "
            f"seed={uncertainty['seed']}) on validation edge resamples "
            f"with pre-computed selected edge sets.",
            "",
        ]
        if uncertainty.get("bootstrap_warning"):
            lines.append(f"> **Warning**: {uncertainty['bootstrap_warning']}")
            lines.append("")

        # Per-variant CI table
        lines += [
            "### Variant 95% Bootstrap CIs",
            "",
            "| Variant | Precision 95% CI | Recall 95% CI | F1 95% CI | Accuracy 95% CI |",
            "|---------|-----------------|---------------|-----------|-----------------|",
        ]
        for vname, vcis in uncertainty["variant_cis"].items():
            lines.append(
                f"| {vname} "
                f"| [{_fmt(vcis['precision'][0])}, {_fmt(vcis['precision'][1])}] "
                f"| [{_fmt(vcis['recall'][0])}, {_fmt(vcis['recall'][1])}] "
                f"| [{_fmt(vcis['f1'][0])}, {_fmt(vcis['f1'][1])}] "
                f"| [{_fmt(vcis['accuracy'][0])}, {_fmt(vcis['accuracy'][1])}] |"
            )

        # Delta CI table
        dci = uncertainty["improvement_summary_ci"]
        lines += [
            "",
            "### Paired Delta 95% Bootstrap CIs (Assumption-Calibrated - Baseline)",
            "",
            "| Delta Metric | 95% CI | P(delta > 0) |",
            "|-------------|--------|--------------|",
            f"| Precision | [{_fmt(dci['delta_precision'][0])}, "
            f"{_fmt(dci['delta_precision'][1])}] "
            f"| {dci.get('probability_delta_precision_gt_0', 0.0):.4f} |",
            f"| Recall | [{_fmt(dci['delta_recall'][0])}, "
            f"{_fmt(dci['delta_recall'][1])}] | — |",
            f"| F1 | [{_fmt(dci['delta_f1'][0])}, "
            f"{_fmt(dci['delta_f1'][1])}] "
            f"| {dci.get('probability_delta_f1_gt_0', 0.0):.4f} |",
            f"| Accuracy | [{_fmt(dci['delta_accuracy'][0])}, "
            f"{_fmt(dci['delta_accuracy'][1])}] | — |",
            "",
            f"Bootstrap interval for delta F1 was "
            f"[{_fmt(dci['delta_f1'][0])}, {_fmt(dci['delta_f1'][1])}].",
            "",
        ]

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
