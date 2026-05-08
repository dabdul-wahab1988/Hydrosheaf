"""Directed graph age-coherence audits for groundwater age networks."""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _node_record(raw: Any) -> Dict[str, Any]:
    if isinstance(raw, Mapping):
        age = (
            _finite_float(raw.get("age_years"))
            or _finite_float(raw.get("mean_age_years"))
            or _finite_float(raw.get("best_mean_age_years"))
        )
        low = (
            _finite_float(raw.get("ci_low_years"))
            or _finite_float(raw.get("age_low_years"))
            or _finite_float(raw.get("age_95_low_years"))
        )
        high = (
            _finite_float(raw.get("ci_high_years"))
            or _finite_float(raw.get("age_high_years"))
            or _finite_float(raw.get("age_95_high_years"))
        )
        flags = raw.get("flags") or raw.get("diagnostic_flags") or {}
        disagreement = raw.get("disagreement") or {}
        if isinstance(disagreement, Mapping):
            flags = {**dict(disagreement.get("flags", {})), **dict(flags)}
        return {
            "age_years": age,
            "ci_low_years": low if low is not None else age,
            "ci_high_years": high if high is not None else age,
            "flags": flags if isinstance(flags, Mapping) else {},
        }
    age = _finite_float(raw)
    return {"age_years": age, "ci_low_years": age, "ci_high_years": age, "flags": {}}


def _edge_pair(edge: Any) -> Optional[Tuple[str, str, Dict[str, Any]]]:
    if isinstance(edge, Mapping):
        upstream = edge.get("upstream") or edge.get("source") or edge.get("from") or edge.get("u")
        downstream = edge.get("downstream") or edge.get("target") or edge.get("to") or edge.get("v")
        if upstream is None or downstream is None:
            return None
        return str(upstream), str(downstream), dict(edge)
    if isinstance(edge, Sequence) and not isinstance(edge, (str, bytes)) and len(edge) >= 2:
        meta = dict(edge[2]) if len(edge) >= 3 and isinstance(edge[2], Mapping) else {}
        return str(edge[0]), str(edge[1]), meta
    return None


def _intervals_overlap(a_low: float, a_high: float, b_low: float, b_high: float) -> bool:
    if a_high < a_low:
        a_low, a_high = a_high, a_low
    if b_high < b_low:
        b_low, b_high = b_high, b_low
    return max(a_low, b_low) <= min(a_high, b_high)


def _explanation(
    upstream_flags: Mapping[str, Any],
    downstream_flags: Mapping[str, Any],
    *,
    overlaps: bool,
) -> str:
    combined = {**dict(upstream_flags), **dict(downstream_flags)}
    if overlaps:
        return "age_intervals_overlap; violation is not resolved at stated uncertainty"
    if combined.get("possible_gas_contamination") or combined.get("uncorroborated_young_gas_with_old_14c"):
        return "possible tracer contamination or young-fraction bias"
    if combined.get("mixed_age_ambiguity") or combined.get("not_identifiable_from_available_tracers"):
        return "possible mixing or non-identifiable tracer age"
    return "possible bad graph direction, recharge input error, or unmodelled mixing"


def audit_graph_age_coherence(
    edges: Iterable[Any],
    node_ages: Mapping[Any, Any],
    *,
    min_downstream_increase_years: float = 0.0,
    severe_log10_threshold: float = 0.3,
) -> Dict[str, Any]:
    """Audit directed edges for physically implausible age reversals.

    A directed edge is interpreted as upstream -> downstream. The downstream
    node is expected to be at least as old as the upstream node unless local
    recharge, mixing, contamination, or graph-direction uncertainty explains
    the reversal.
    """

    records = {str(node): _node_record(value) for node, value in node_ages.items()}
    edge_rows: List[Dict[str, Any]] = []
    for edge in edges:
        parsed = _edge_pair(edge)
        if parsed is None:
            continue
        upstream, downstream, meta = parsed
        up = records.get(upstream)
        down = records.get(downstream)
        if up is None or down is None:
            edge_rows.append(
                {
                    "upstream": upstream,
                    "downstream": downstream,
                    "checked": False,
                    "violation": False,
                    "reason": "missing node age",
                    "metadata": meta,
                }
            )
            continue
        up_age = _finite_float(up.get("age_years"))
        down_age = _finite_float(down.get("age_years"))
        if up_age is None or down_age is None:
            edge_rows.append(
                {
                    "upstream": upstream,
                    "downstream": downstream,
                    "checked": False,
                    "violation": False,
                    "reason": "non-finite node age",
                    "metadata": meta,
                }
            )
            continue

        expected_min = up_age + float(min_downstream_increase_years)
        severity_years = max(0.0, expected_min - down_age)
        severity_log10 = max(0.0, math.log10(max(expected_min, 0.1)) - math.log10(max(down_age, 0.1)))
        violation = severity_years > 0.0
        overlaps = _intervals_overlap(
            float(up.get("ci_low_years") or up_age),
            float(up.get("ci_high_years") or up_age),
            float(down.get("ci_low_years") or down_age),
            float(down.get("ci_high_years") or down_age),
        )
        edge_rows.append(
            {
                "upstream": upstream,
                "downstream": downstream,
                "checked": True,
                "upstream_age_years": float(up_age),
                "downstream_age_years": float(down_age),
                "age_difference_down_minus_up_years": float(down_age - up_age),
                "uncertainty_overlap": bool(overlaps),
                "violation": bool(violation),
                "severity_years": float(severity_years),
                "severity_log10": float(severity_log10),
                "severe": bool(violation and severity_log10 >= severe_log10_threshold and not overlaps),
                "explanation": _explanation(
                    up.get("flags", {}),
                    down.get("flags", {}),
                    overlaps=overlaps,
                )
                if violation
                else "age order coherent",
                "metadata": meta,
            }
        )

    checked = [row for row in edge_rows if row.get("checked")]
    violations = [row for row in checked if row.get("violation")]
    severe = [row for row in checked if row.get("severe")]
    return {
        "n_edges": len(edge_rows),
        "n_checked": len(checked),
        "n_violations": len(violations),
        "n_severe_violations": len(severe),
        "violation_fraction": (len(violations) / len(checked)) if checked else float("nan"),
        "edges": edge_rows,
        "interpretation": (
            "Directed age reversals detected; inspect severe edges before using graph regularization."
            if violations
            else "No downstream-younger-than-upstream age reversals detected."
        ),
    }
