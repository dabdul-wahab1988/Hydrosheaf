"""
Active-learning measurement recommendations for topology uncertainty reduction.

Ranks edges by priority for the next field/lab campaign using five signals:
  1. Variant disagreement (baseline vs null-default vs assumption-calibrated)
  2. Bootstrap instability (report-level CI width / probability of delta > 0)
  3. Evidence ambiguity (flags such as age_reversal, missing_isotopes, etc.)
  4. Validation error status (FP / FN / selected_unlabeled)
  5. Missing data gaps in sample node measurements

This is a decision-support workflow — it never claims an edge is validated.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..config import Config as HConfig
from ..graph.types import Edge

# ── measurement recommendation mapping ──────────────────────────────

FLAG_TO_MEASUREMENTS: Dict[str, List[str]] = {
    "age_reversal": [
        "groundwater age tracer", "tritium/14C/SF6/CFC",
    ],
    "evap_candidate": [
        "d18O/d2H sampling", "evaporation-line check",
    ],
    "null_chemistry_similar": [
        "major ion chemistry", "chloride/bromide",
    ],
    "null_common_lithology": [
        "lithologic log / aquifer unit assignment",
    ],
    "null_common_lithology_explicit": [
        "lithologic log / aquifer unit assignment",
    ],
    "null_shared_recharge": [
        "d18O/d2H sampling", "evaporation-line check",
    ],
    "null_isotope_proximity": [
        "d18O/d2H sampling",
    ],
    "missing_evidence": [
        "independent connectivity evidence",
    ],
    "iso_missing_u": [
        "d18O/d2H sampling at upstream node",
    ],
    "iso_missing_v": [
        "d18O/d2H sampling at downstream node",
    ],
    "cl_missing": [
        "chloride sampling",
    ],
    "missing_head": [
        "repeat hydraulic head", "datum/elevation survey",
    ],
    "missing_elevation": [
        "datum/elevation survey",
    ],
    "missing_isotopes": [
        "d18O/d2H sampling", "evaporation-line check",
    ],
    "missing_age": [
        "groundwater age tracer", "tritium/14C/SF6/CFC",
    ],
    "missing_lithology": [
        "lithologic log / aquifer unit assignment",
    ],
    "missing_screen_interval": [
        "well construction / screen interval survey",
    ],
    "missing_aquifer_layer": [
        "lithologic log / aquifer unit assignment",
    ],
    "missing_major_ions": [
        "major ion chemistry", "chloride/bromide",
    ],
    "selected_unlabeled": [
        "independent connectivity evidence",
        "MODPATH/pathline check", "tracer test",
    ],
}

# Default field keys to check for missing data (configurable via HConfig)
_DEFAULT_SAMPLE_FIELDS = [
    "hydraulic_head", "elevation", "d18O", "d2H", "Cl",
    "age", "lithology", "screen_interval", "aquifer_layer",
]


# ── internal helpers ─────────────────────────────────────────────────

def _parse_evidence_flags(edge: Edge) -> List[str]:
    """Parse comma-separated evidence_flags from edge attrs."""
    flags_str = (edge.attrs or {}).get("evidence_flags", "")
    if not flags_str or not isinstance(flags_str, str):
        return []
    return [f.strip() for f in flags_str.split(",") if f.strip()]


def _compute_disagreement_score(
    edge_id: str,
    variant_selected: Dict[str, bool],
) -> tuple:
    """Score how much the three variants disagree on this edge.

    Returns (score, reason_string).
    """
    bl = variant_selected.get("baseline", False)
    nd = variant_selected.get("null_model_defaults", False)
    ac = variant_selected.get("assumption_calibrated", False)

    if ac and not bl:
        return (1.0, "assumption-sensitive (calibrated selects, baseline does not)")
    if bl and not ac:
        return (0.8, "calibration-removed (baseline selected, calibrated removed)")
    if nd and not bl and not ac:
        return (0.4, "null-model-ambiguous (only null-model-defaults selects)")
    if bl and ac and not nd:
        return (0.3, "null-model-rejects (baseline+calibrated agree, null-model-default disagrees)")
    if bl == nd == ac:
        return (0.0, "all-variants-agree")
    return (0.2, "mixed-disagreement")


def _compute_bootstrap_instability_modulation(
    benchmark_report: dict,
) -> tuple:
    """Extract report-level bootstrap instability boost.

    Returns (modulation_float, reason_string_or_empty).
    """
    uncertainty = benchmark_report.get("uncertainty")
    if not uncertainty:
        return (0.0, "")

    dci = uncertainty.get("improvement_summary_ci", {})
    reasons = []
    modulation = 0.0

    delta_f1_ci = dci.get("delta_f1", [0.0, 0.0])
    ci_width = delta_f1_ci[1] - delta_f1_ci[0]
    if ci_width > 0.5:
        modulation += 0.15
        reasons.append("wide delta_f1 CI (>{:.2f})".format(ci_width))

    prob_f1 = dci.get("probability_delta_f1_gt_0", 1.0)
    if prob_f1 < 0.95:
        modulation += 0.10
        reasons.append("P(delta_f1>0) = {:.2f} < 0.95".format(prob_f1))

    return (modulation, "; ".join(reasons))


def _score_evidence_ambiguity(edge: Edge) -> tuple:
    """Score evidence ambiguity from edge attrs flags.

    Returns (score, reason_string).
    """
    flags = _parse_evidence_flags(edge) or []
    attrs = edge.attrs or {}

    # Also check sheaf_flags as fallback
    sheaf_str = attrs.get("sheaf_flags", "")
    if sheaf_str and isinstance(sheaf_str, str):
        sheaf_flags = [f.strip() for f in sheaf_str.split(",") if f.strip()]
    else:
        sheaf_flags = []

    all_flags = set(flags + sheaf_flags)
    if not all_flags:
        evidence_class = attrs.get("evidence_class", "")
        if evidence_class == "VALIDATED":
            return (0.0, "")
        return (0.0, "")

    # Tiered priority
    tier_map = [
        ({"age_reversal"}, 0.7, "age_reversal"),
        ({"evap_candidate"}, 0.6, "isotope_evaporation"),
        ({"null_chemistry_similar", "null_chemistry_error"}, 0.5, "null_chemistry_similar"),
        ({"missing_evidence"}, 0.4, "missing_evidence"),
        ({"iso_missing_u", "iso_missing_v"}, 0.35, "missing_isotope_data"),
        ({"cl_missing"}, 0.3, "missing_chloride_data"),
        ({"null_common_lithology", "null_common_lithology_explicit"}, 0.2, "null_common_lithology"),
        ({"null_shared_recharge", "null_isotope_proximity"}, 0.15, "null_model_flag"),
        ({"null_spatial_autocorr", "null_common_anthropogenic"}, 0.1, "null_model_flag"),
    ]

    best_score = 0.0
    best_reasons = []
    for flag_set, score, reason in tier_map:
        if all_flags & flag_set:
            matched = sorted(all_flags & flag_set)
            if score > best_score:
                best_score = score
                best_reasons = matched
            elif score == best_score:
                best_reasons.extend(matched)

    if best_score == 0.0 and all_flags:
        best_score = 0.05
        best_reasons = sorted(all_flags)

    return (best_score, ", ".join(best_reasons))


def _score_validation_error(
    edge_id: str,
    validation_report: Optional[dict],
    val_observations_by_edge: Dict[str, float],
    cal_selected_set: set,
) -> tuple:
    """Score based on validation label status.

    Returns (score, validation_status_string, reason_string).
    """
    in_selected = edge_id in cal_selected_set
    in_val_labels = edge_id in val_observations_by_edge

    if not validation_report and not val_observations_by_edge:
        return (0.0, "unlabeled", "")

    if in_val_labels:
        observed_present = val_observations_by_edge[edge_id]
        if observed_present >= 0.5:
            if in_selected:
                return (0.1, "observed_present", "correctly selected (TP)")
            else:
                return (0.9, "false_negative", "missed detection (FN)")
        else:
            if in_selected:
                return (1.0, "false_positive", "incorrectly selected (FP)")
            else:
                return (0.0, "true_negative", "correctly excluded (TN)")
    else:
        if in_selected:
            return (0.5, "selected_unlabeled", "model selected; validation label missing")
        else:
            return (0.0, "unlabeled", "")


def _score_missing_data(
    edge: Edge,
    samples: Optional[List[dict]],
    config: Optional[HConfig],
) -> tuple:
    """Score missing measurement data for the two nodes of an edge.

    Returns (score, list_of_reason_strings).
    """
    if not samples or not config:
        return (0.0, [])

    # Build sample lookup
    sample_by_node: Dict[str, dict] = {}
    for s in samples:
        for key in ("site_id", "node_id", "sample_id", "well_id"):
            if key in s:
                sample_by_node[str(s[key])] = s
                break

    su = sample_by_node.get(edge.u, {})
    sv = sample_by_node.get(edge.v, {})

    # Field name mapping from config
    field_keys = {
        "hydraulic_head": getattr(config, "edge_head_key", "hydraulic_head") or "hydraulic_head",
        "elevation": getattr(config, "edge_elevation_key", "elevation") or "elevation",
        "d18O": getattr(config, "isotope_d18o_key", "d18O") or "d18O",
        "d2H": getattr(config, "isotope_d2h_key", "d2H") or "d2H",
        "Cl": "Cl",
        "age": "age",
        "lithology": "lithology",
        "screen_interval": "screen_interval",
        "aquifer_layer": getattr(config, "layer_key", "aquifer_layer") or "aquifer_layer",
    }

    reason_to_field = {
        "missing_head": "hydraulic_head",
        "missing_elevation": "elevation",
        "missing_isotopes": "d18O",  # or d2H
        "missing_age": "age",
        "missing_lithology": "lithology",
        "missing_screen_interval": "screen_interval",
        "missing_aquifer_layer": "aquifer_layer",
        "missing_major_ions": "Cl",
    }

    missing_count = 0
    total_fields = 0
    reasons_set: set = set()

    for display_name, csv_key in field_keys.items():
        total_fields += 2  # one per node
        for node_name, sample in [(edge.u, su), (edge.v, sv)]:
            val = sample.get(csv_key)
            if val is None or (isinstance(val, float) and np.isnan(val)) or val == "":
                missing_count += 1
                # Map to reason
                for reason, field in reason_to_field.items():
                    if field in display_name or field == csv_key:
                        reasons_set.add(f"{reason} (node {node_name})")
                        break
                else:
                    reasons_set.add(f"missing_data (node {node_name})")

    if total_fields == 0:
        return (0.0, [])

    fraction = missing_count / total_fields
    return (fraction * 0.5, sorted(reasons_set))


def _recommended_measurements(uncertainty_reasons: List[str]) -> List[str]:
    """Map uncertainty reasons to recommended field/lab measurements.

    Reasons may include a node suffix (e.g. ``missing_age (node A)``)
    emitted by ``_score_missing_data``.  The suffix is stripped before
    looking up the matching measurement recommendation.
    """
    import re
    measurements: set = set()
    for reason in uncertainty_reasons:
        # Strip "(node ...)" suffix added by _score_missing_data
        base_reason = re.sub(r"\s*\(node\s+\S+\)\s*$", "", reason).strip()
        if base_reason in FLAG_TO_MEASUREMENTS:
            for m in FLAG_TO_MEASUREMENTS[base_reason]:
                measurements.add(m)
    return sorted(measurements)


def _score_posterior_uncertainty(edge: Edge) -> Tuple[float, str]:
    """Signal F: Posterior uncertainty from Bayesian topology sampling.

    Edges with posterior probability near 0.5 are maximally uncertain
    and get highest priority. Also factors in overall posterior entropy.
    """
    attrs = edge.attrs or {}
    prob = attrs.get("posterior_edge_probability")
    entropy = attrs.get("posterior_topology_entropy")

    if prob is None:
        return 0.0, ""

    try:
        prob = float(prob)
    except (ValueError, TypeError):
        return 0.0, ""

    # Maximum uncertainty at p=0.5, zero at p=0 or p=1
    # Score = 1 - |2p - 1| which peaks at 1 when p=0.5
    prob_uncertainty = 1.0 - abs(2.0 * prob - 1.0)

    # Boost by system-level entropy if available
    entropy_boost = 0.0
    if entropy is not None:
        try:
            entropy = float(entropy)
            # Normalise entropy boost (max per-edge entropy is ln(2) ~ 0.69)
            entropy_boost = min(1.0, entropy / 0.69)
        except (ValueError, TypeError):
            pass

    score = prob_uncertainty * (1.0 + entropy_boost) / 2.0
    score = min(1.0, score)

    reason = ""
    if prob_uncertainty > 0.5:
        reason = (
            f"Posterior edge probability near 0.5 (p={prob:.3f}); "
            "additional measurements would reduce topology uncertainty"
        )

    return score, reason


def _compute_expected_benefit(
    priority_score: float,
    uncertainty_reasons: List[str],
    evidence_class: str,
    validation_status: str,
) -> str:
    """Generate a human-readable expected benefit summary."""
    parts = []
    if evidence_class == "FALSIFIED":
        parts.append("Resolving flags could reclassify edge from FALSIFIED to PROBABLE")
    elif evidence_class == "AMBIGUOUS":
        parts.append("Additional data would reduce ambiguity and allow definitive classification")
    elif evidence_class == "PRIOR_ASSISTED":
        parts.append("Independent evidence would strengthen or refute prior-based selection")

    if validation_status == "false_positive":
        parts.append("measurements could explain why the edge is selected despite validation evidence")
    elif validation_status == "selected_unlabeled":
        parts.append("measurements could provide independent confirmation before field validation")

    if not parts:
        if priority_score > 0.5:
            parts.append("High-priority edge for uncertainty reduction")
        else:
            parts.append("Additional data would improve topology confidence")

    return " ".join(parts)


def _manuscript_safe_note(evidence_class: str, validation_status: str) -> str:
    """Generate a manuscript-safe note that never claims an edge is validated."""
    if evidence_class == "VALIDATED":
        return "Edge is already VALIDATED; new measurements provide independent confirmation only."
    if evidence_class == "FALSIFIED":
        return "Evidence currently falsifies this edge; measurements would test alternative hypotheses."
    if evidence_class == "AMBIGUOUS":
        return "Edge classification is currently ambiguous; measurements would add independent constraints."
    if validation_status == "observed_present":
        return "This measurement recommendation does not alter the edge's validation status."
    return "Recommendations are for planning purposes; they do not validate the edge."


# ── main entry point ────────────────────────────────────────────────

def rank_next_measurements(
    benchmark_report: dict,
    validation_report: Optional[dict] = None,
    samples: Optional[List[dict]] = None,
    candidate_edges: Optional[List[Edge]] = None,
    config: Optional[HConfig] = None,
    top_k: int = 20,
    output_dir: Optional[str] = None,
) -> dict:
    """Rank edges by priority for the next field/lab measurement campaign.

    Parameters
    ----------
    benchmark_report : dict
        Full benchmark report from ``run_assumption_benchmark``
        (typically loaded from ``assumption_benchmark_results.json``).
    validation_report : dict, optional
        Validation report from ``run_assumption_calibration_validation``
        (typically loaded from ``assumption_validation_results.json``).
    samples : list of dict, optional
        Sample/well measurement data used for missing-data scoring.
    candidate_edges : list of Edge, optional
        Full candidate edge list with attrs (evidence_flags, etc.).
    config : Config, optional
        Hydrosheaf config for field-name resolution.
    top_k : int
        Maximum number of recommendations to return.
    output_dir : str, optional
        When provided, writes JSON/CSV/MD output files.

    Returns
    -------
    dict
        ``recommendations`` list, ``summary``, and ``inputs_used``.
    """
    variants = benchmark_report.get("variants", {})
    if not variants:
        result = {
            "recommendations": [],
            "summary": {
                "n_recommendations": 0,
                "n_edges_scored": 0,
                "top_priority_score": 0.0,
                "bootstrap_instability_detected": False,
                "warning": "No variant data in benchmark report.",
            },
            "inputs_used": {
                "benchmark_report": True,
                "validation_report": validation_report is not None,
                "samples": samples is not None,
                "candidate_edges": candidate_edges is not None,
                "config": config is not None,
            },
        }
        if output_dir:
            _write_outputs(result, output_dir)
        return result

    # ── Build per-variant selected-edge-id sets ─────────────────────
    variant_selected_sets: Dict[str, set] = {}
    for vname, vdata in variants.items():
        variant_selected_sets[vname] = set(vdata.get("selected_edge_ids", []))

    # ── Build edge lookup table ─────────────────────────────────────
    edge_lookup: Dict[str, Edge] = {}
    if candidate_edges:
        for e in candidate_edges:
            edge_lookup[e.edge_id] = e

    # ── Validation observations ─────────────────────────────────────
    val_obs_by_edge: Dict[str, float] = {}
    if validation_report:
        val_label_file = validation_report.get("validation_label_file")
        if val_label_file and os.path.isfile(val_label_file):
            try:
                from .validation_workflow import _load_topology_observations
                val_obs = _load_topology_observations(val_label_file)
                for o in val_obs:
                    val_obs_by_edge[o.edge_id] = o.observed_present
            except Exception:
                pass
        # Also check validation_report for direct observations
        if not val_obs_by_edge:
            eval_ids = validation_report.get("evaluated_validation_edge_ids", [])
            for eid in eval_ids:
                val_obs_by_edge.setdefault(eid, 0.5)

    # ── Calibrated-variant selected set (for validation scoring) ─────
    ac_selected = variant_selected_sets.get("assumption_calibrated", set())

    # ── Report-level bootstrap instability ──────────────────────────
    boot_modulation, boot_reason = _compute_bootstrap_instability_modulation(
        benchmark_report,
    )

    # ── Collect all edge IDs to rank ────────────────────────────────
    all_edge_ids: set = set()
    for vset in variant_selected_sets.values():
        all_edge_ids |= vset
    if candidate_edges:
        all_edge_ids |= {e.edge_id for e in candidate_edges}

    # ── Score each edge ─────────────────────────────────────────────
    scored: List[dict] = []
    for edge_id in all_edge_ids:
        edge = edge_lookup.get(edge_id, Edge(edge_id=edge_id, u="", v=""))

        # Which variants select this edge?
        vs = {
            vname: edge_id in vset
            for vname, vset in variant_selected_sets.items()
        }

        # Signal A: disagreement
        disagree_score, disagree_reason = _compute_disagreement_score(edge_id, vs)

        # Signal B: bootstrap instability (report-level, same for all edges)
        bootstrap_score = boot_modulation

        # Signal C: evidence ambiguity
        ev_score, ev_reason = _score_evidence_ambiguity(edge)

        # Signal D: validation error
        val_score, val_status, val_reason = _score_validation_error(
            edge_id, validation_report, val_obs_by_edge, ac_selected,
        )

        # Signal E: missing data
        md_score, md_reasons = _score_missing_data(edge, samples, config)

        # Signal F: posterior uncertainty (Bayesian topology posterior)
        post_score, post_reason = _score_posterior_uncertainty(edge)

        # Aggregate priority
        priority = (
            disagree_score * 0.25
            + bootstrap_score * 0.10
            + ev_score * 0.20
            + val_score * 0.20
            + md_score * 0.10
            + post_score * 0.15
        )
        priority = round(min(max(priority, 0.0), 1.0), 6)

        # Collect reasons
        all_reasons = []
        if disagree_reason:
            all_reasons.append(disagree_reason)
        if boot_reason:
            all_reasons.append(boot_reason)
        if ev_reason:
            all_reasons.append(ev_reason)
        if val_reason:
            all_reasons.append(val_reason)
        if md_reasons:
            all_reasons.extend(md_reasons)
        if post_reason:
            all_reasons.append(post_reason)

        evidence_class = (edge.attrs or {}).get("evidence_class", "")
        if not evidence_class:
            evidence_class = ""

        scored.append({
            "edge_id": edge_id,
            "u": edge.u or "",
            "v": edge.v or "",
            "priority_score": priority,
            "uncertainty_reasons": all_reasons,
            "evidence_flags": _parse_evidence_flags(edge),
            "evidence_reason": (edge.attrs or {}).get("evidence_reason", ""),
            "evidence_class": evidence_class,
            "baseline_selected": vs.get("baseline", False),
            "null_model_default_selected": vs.get("null_model_defaults", False),
            "assumption_calibrated_selected": vs.get("assumption_calibrated", False),
            "validation_status": val_status,
        })

    # ── Sort by priority descending ─────────────────────────────────
    scored.sort(key=lambda x: x["priority_score"], reverse=True)
    top = scored[:top_k]

    # ── Build final recommendation dicts ────────────────────────────
    recommendations = []
    for i, item in enumerate(top):
        reasons = item["uncertainty_reasons"]
        recs = _recommended_measurements(reasons)
        status_item = {
            "baseline_selected": item["baseline_selected"],
            "null_model_default_selected": item["null_model_default_selected"],
            "assumption_calibrated_selected": item["assumption_calibrated_selected"],
        }
        recommendations.append({
            "rank": i + 1,
            "edge_id": item["edge_id"],
            "u": item["u"],
            "v": item["v"],
            "priority_score": item["priority_score"],
            "uncertainty_reasons": reasons,
            "recommended_measurements": recs if recs else ["review edge evidence"],
            "expected_benefit": _compute_expected_benefit(
                item["priority_score"], reasons, item["evidence_class"],
                item["validation_status"],
            ),
            "evidence_flags": item["evidence_flags"],
            "evidence_reason": item["evidence_reason"],
            "evidence_class": item["evidence_class"],
            "current_selection_status": status_item,
            "validation_status": item["validation_status"],
            "manuscript_safe_note": _manuscript_safe_note(
                item["evidence_class"], item["validation_status"],
            ),
        })

    summary = {
        "n_recommendations": len(recommendations),
        "top_priority_score": recommendations[0]["priority_score"] if recommendations else 0.0,
        "n_edges_scored": len(scored),
        "bootstrap_instability_detected": bool(boot_reason),
    }

    result = {
        "recommendations": recommendations,
        "summary": summary,
        "inputs_used": {
            "benchmark_report": True,
            "validation_report": validation_report is not None,
            "samples": samples is not None,
            "candidate_edges": candidate_edges is not None,
            "config": config is not None,
        },
    }

    if output_dir:
        _write_outputs(result, output_dir)

    return result


# ── output writers ───────────────────────────────────────────────────

def _write_outputs(result: dict, output_dir: str) -> None:
    """Write JSON, CSV, and Markdown output files."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = out_dir / "active_learning_recommendations.json"

    def _convert(obj: Any) -> Any:
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        raise TypeError

    with open(json_path, "w") as f:
        json.dump(result, f, indent=2, default=_convert)

    # CSV
    csv_path = out_dir / "active_learning_recommendations.csv"
    _write_active_learning_csv(csv_path, result["recommendations"])

    # Markdown
    md_path = out_dir / "active_learning_report.md"
    _write_active_learning_md(md_path, result)


def _write_active_learning_csv(
    path: Path,
    recommendations: List[dict],
) -> None:
    """Write a flat CSV of recommendations."""
    fieldnames = [
        "rank", "edge_id", "u", "v", "priority_score",
        "uncertainty_reasons", "recommended_measurements",
        "expected_benefit", "validation_status",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for rec in recommendations:
            row = dict(rec)
            row["uncertainty_reasons"] = "; ".join(rec.get("uncertainty_reasons", []))
            row["recommended_measurements"] = "; ".join(rec.get("recommended_measurements", []))
            writer.writerow(row)


def _write_active_learning_md(path: Path, result: dict) -> None:
    """Write a human-readable Markdown report."""
    recs = result["recommendations"]
    summary = result["summary"]

    lines = [
        "# Active Learning: Next-Measurement Recommendations",
        "",
        "This report ranks edges for the next field/lab measurement campaign.",
        "Recommendations are decision-support only — they do not validate any edge.",
        "",
        f"- **Recommendations**: {summary['n_recommendations']}",
        f"- **Edges scored**: {summary['n_edges_scored']}",
        f"- **Top priority score**: {summary['top_priority_score']:.4f}",
    ]

    if summary.get("bootstrap_instability_detected"):
        lines.append("- **Bootstrap instability**: detected in benchmark delta CIs")
    lines.append("")

    if not recs:
        lines.append("No recommendations generated.")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        return

    lines += [
        "## Top Recommendations",
        "",
        "| Rank | Edge ID | Priority | Variant Selection | Validation | Key Reasons | Recommended Measurements |",
        "|------|---------|----------|-------------------|------------|-------------|--------------------------|",
    ]

    for rec in recs[:20]:
        sel_parts = []
        if rec["current_selection_status"]["baseline_selected"]:
            sel_parts.append("B")
        else:
            sel_parts.append("-")
        if rec["current_selection_status"]["null_model_default_selected"]:
            sel_parts.append("N")
        else:
            sel_parts.append("-")
        if rec["current_selection_status"]["assumption_calibrated_selected"]:
            sel_parts.append("C")
        else:
            sel_parts.append("-")
        sel_str = "".join(sel_parts)

        reasons_str = "; ".join(rec["uncertainty_reasons"][:2])
        if len(rec["uncertainty_reasons"]) > 2:
            reasons_str += " ..."
        if not reasons_str:
            reasons_str = "—"

        meas_str = "; ".join(rec["recommended_measurements"][:2])
        if len(rec["recommended_measurements"]) > 2:
            meas_str += " ..."

        lines.append(
            f"| {rec['rank']} | {rec['edge_id']} | {rec['priority_score']:.4f} "
            f"| {sel_str} | {rec['validation_status']} "
            f"| {reasons_str} | {meas_str} |"
        )

    # Detail section
    lines += [
        "",
        "## Recommendation Details",
        "",
    ]
    for rec in recs[:10]:
        lines += [
            f"### {rec['rank']}. {rec['edge_id']} (priority: {rec['priority_score']:.4f})",
            "",
            f"- **Edge**: {rec['u']} → {rec['v']}",
            f"- **Evidence class**: {rec['evidence_class'] or '—'}",
            f"- **Validation status**: {rec['validation_status']}",
            f"- **Selection**: baseline={rec['current_selection_status']['baseline_selected']}, "
            f"null={rec['current_selection_status']['null_model_default_selected']}, "
            f"calibrated={rec['current_selection_status']['assumption_calibrated_selected']}",
            "",
        ]
        if rec["uncertainty_reasons"]:
            lines.append("- **Uncertainty reasons**:")
            for r in rec["uncertainty_reasons"]:
                lines.append(f"  - {r}")
            lines.append("")
        if rec["evidence_flags"]:
            lines.append(f"- **Evidence flags**: {', '.join(rec['evidence_flags'])}")
            lines.append("")
        if rec.get("evidence_reason"):
            lines.append(f"- **Evidence detail**: {rec['evidence_reason']}")
            lines.append("")
        lines.append(f"- **Recommended**: {', '.join(rec['recommended_measurements'])}")
        lines.append("")
        lines.append(f"- **Expected benefit**: {rec['expected_benefit']}")
        lines.append("")
        lines.append(f"- **Note**: {rec['manuscript_safe_note']}")
        lines.append("")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
