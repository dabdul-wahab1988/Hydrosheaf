"""Reusable helpers for the synthetic programme validation gate.

The programme runner is deliberately thin: benchmark-specific code prepares a
blind case and calls HydroSheaf, while this module owns the shared truth-blind
boundary, stage interpretation, and basic scoring definitions.  This keeps
the M2-M8 runners from growing another set of subtly different validators.
"""

from __future__ import annotations

import math
import random
from typing import Iterable, Mapping, Sequence

from .metrics import (
    classification_metrics,
    interval_coverage,
    regression_metrics,
    topology_metrics,
)
from .programme_contract import ProgrammeStage, StageStatus, assert_truth_blind


def prepare_truth_blind_rows(
    observations: Iterable[Mapping[str, object]],
    *,
    tracer_source: str | None = None,
    tracer_target: str | None = None,
    permute_tracer_seed: int | None = None,
    forbidden_fields: Iterable[str] = (),
) -> list[dict[str, object]]:
    """Copy observations, optionally remap/permutate one tracer, and audit them.

    Truth is never accepted as an input argument.  The caller supplies only
    the observed rows and, when needed, a named observed tracer to remap into
    the public API's expected field name.
    """

    rows = [dict(row) for row in observations]
    assert_truth_blind(rows, forbidden_fields=forbidden_fields)

    if (tracer_source is None) != (tracer_target is None):
        raise ValueError("tracer_source and tracer_target must be supplied together.")
    if tracer_source is not None and tracer_target is not None:
        values = []
        for row in rows:
            if tracer_source not in row:
                raise KeyError(f"Observed tracer field is missing: {tracer_source}")
            values.append(row[tracer_source])
        if permute_tracer_seed is not None:
            random.Random(int(permute_tracer_seed)).shuffle(values)
        for row, value in zip(rows, values):
            row[tracer_target] = value

    # The remapping itself must not have introduced a truth-looking field.
    assert_truth_blind(rows, forbidden_fields=forbidden_fields)
    return rows


def programme_stages_from_status(
    stage_status: Mapping[str, Mapping[str, object]],
    *,
    name_prefix: str = "",
) -> tuple[ProgrammeStage, ...]:
    """Convert the public pipeline's stage records into programme records."""

    stages: list[ProgrammeStage] = []
    for stage_name in sorted(stage_status):
        record = stage_status[stage_name]
        status = StageStatus(str(record.get("status", "failed")))
        detail = str(record.get("detail", ""))
        truth_fields = record.get("truth_fields_seen", ())
        if isinstance(truth_fields, str):
            truth_fields = (truth_fields,)
        stages.append(
            ProgrammeStage(
                name=f"{name_prefix}{stage_name}",
                status=status,
                # Missing instrumentation is not evidence of truth blindness;
                # fail closed so an incomplete stage cannot pass the gate.
                truth_blind=bool(record.get("truth_blind", False)),
                truth_fields_seen=tuple(str(item) for item in truth_fields),
                detail=detail,
            )
        )
    return tuple(stages)


def required_stages_completed(
    stage_status: Mapping[str, Mapping[str, object]],
    required_stages: Iterable[str],
) -> bool:
    """Return whether every requested stage completed successfully."""

    return all(
        str(stage_status.get(name, {}).get("status")) == StageStatus.COMPLETED.value
        for name in required_stages
    )


def _edge_pair(edge: object) -> tuple[str, str]:
    if hasattr(edge, "u") and hasattr(edge, "v"):
        return str(getattr(edge, "u")), str(getattr(edge, "v"))
    if isinstance(edge, Mapping):
        return str(edge["u"]), str(edge["v"])
    values = tuple(edge)  # type: ignore[arg-type]
    if len(values) < 2:
        raise ValueError(f"Edge must contain at least two values: {edge!r}")
    return str(values[0]), str(values[1])


def score_topology(
    reference_edges: Iterable[Sequence[object]],
    candidate_edges: Iterable[object],
    selected_edges: Iterable[object],
) -> dict[str, float]:
    """Score the candidate and selected directed edge sets consistently."""

    reference = [tuple(edge[:2]) for edge in reference_edges]
    candidates = [_edge_pair(edge) for edge in candidate_edges]
    selected = [_edge_pair(edge) for edge in selected_edges]
    candidate_metrics = topology_metrics(
        reference,
        candidates,
        candidate_edges=candidates,
    )
    selected_metrics = topology_metrics(
        reference,
        selected,
        candidate_edges=candidates,
    )
    return {
        "candidate_precision": float(candidate_metrics["precision"]),
        "candidate_recall": float(candidate_metrics["recall"]),
        "candidate_f1": float(candidate_metrics["f1"]),
        "candidate_false_positive_rate": float(
            candidate_metrics["false_positive_rate"]
        ),
        "candidate_false_discovery_rate": float(
            candidate_metrics["false_discovery_rate"]
        ),
        "selected_precision": float(selected_metrics["precision"]),
        "selected_recall": float(selected_metrics["recall"]),
        "selected_f1": float(selected_metrics["f1"]),
        "selected_false_positive_rate": float(
            selected_metrics["false_positive_rate"]
        ),
        "selected_false_discovery_rate": float(
            selected_metrics["false_discovery_rate"]
        ),
        "n_reference_edges": float(candidate_metrics["reference_edges"]),
        "n_candidate_edges": float(candidate_metrics["inferred_edges"]),
        "n_selected_edges": float(selected_metrics["inferred_edges"]),
    }


def score_age_posteriors(
    true_ages_years: Mapping[str, object],
    posterior_results: Mapping[str, Mapping[str, object]] | None,
) -> dict[str, object]:
    """Score point and 95% interval age outputs without pooling other units."""

    if not posterior_results:
        return {"status": "not_available", "n": 0}

    observed: list[float] = []
    predicted: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    for node, true_age in true_ages_years.items():
        result = posterior_results.get(node)
        if not isinstance(result, Mapping):
            continue
        required = ("mean_age_years", "age_95_low", "age_95_high")
        if not all(key in result for key in required):
            continue
        observed.append(float(true_age))
        predicted.append(float(result["mean_age_years"]))
        lower.append(float(result["age_95_low"]))
        upper.append(float(result["age_95_high"]))

    if not observed:
        return {"status": "no_comparable_outputs", "n": 0}
    point = dict(regression_metrics(observed, predicted))
    interval = dict(interval_coverage(observed, lower, upper))
    finite = all(
        math.isfinite(float(value))
        for metric in (point, interval)
        for value in metric.values()
        if isinstance(value, (int, float))
    )
    return {
        "status": "scored",
        "n": len(observed),
        "point": point,
        "interval": interval,
        "all_finite": finite,
    }


def _reaction_family(label: object) -> str:
    text = str(label or "").lower()
    if any(token in text for token in ("calcite", "dolomite", "carbonate")):
        return "carbonate"
    if any(
        token in text
        for token in ("albite", "anorthite", "feldspar", "silicate", "exch")
    ):
        return "silicate_exchange"
    if "sulfate_reduction" in text:
        return "sulfate_reduction"
    if "iron_reduction" in text:
        return "iron_reduction"
    if "denit" in text:
        return "denitrification"
    if "pyrite" in text:
        return "other_redox"
    if "gypsum" in text or "so4" in text:
        return "sulfate_source"
    return "other"


def _predicted_reaction_family(result: object) -> str:
    labels = list(getattr(result, "z_labels", ()) or ())
    extents = list(getattr(result, "z_extents", ()) or ())
    if not labels or not extents:
        return "none"
    pairs = []
    for label, extent in zip(labels, extents):
        try:
            magnitude = abs(float(extent))
        except (TypeError, ValueError):
            continue
        if math.isfinite(magnitude):
            pairs.append((magnitude, label))
    if not pairs:
        return "none"
    magnitude, label = max(pairs, key=lambda item: item[0])
    return _reaction_family(label) if magnitude > 1.0e-8 else "none"


def score_reaction_families(
    true_processes: Mapping[str, object],
    edge_results: Iterable[object] | None,
) -> dict[str, object]:
    """Score the highest-magnitude inferred reaction family per true edge."""

    if not edge_results:
        return {
            "status": "not_available",
            "n": 0,
            "n_truth": len(true_processes),
            "n_missing_outputs": len(true_processes),
        }
    result_by_edge = {
        str(getattr(result, "edge_id")): result for result in edge_results
    }
    expected: list[str] = []
    predicted: list[str] = []
    for edge_id, process in true_processes.items():
        result = result_by_edge.get(str(edge_id))
        if result is None:
            continue
        expected.append(_reaction_family(process))
        predicted.append(_predicted_reaction_family(result))
    if not expected:
        return {
            "status": "no_comparable_outputs",
            "n": 0,
            "n_truth": len(true_processes),
            "n_missing_outputs": len(true_processes),
        }
    metrics = dict(classification_metrics(expected, predicted))
    metrics["unresolved_rate"] = predicted.count("none") / len(predicted)
    return {
        "status": "scored",
        "n": len(expected),
        "n_truth": len(true_processes),
        "n_missing_outputs": len(true_processes) - len(expected),
        "metrics": metrics,
        "expected_families": sorted(set(expected)),
        "predicted_families": sorted(set(predicted)),
    }


def finite_numeric_mapping(mapping: Mapping[str, object]) -> bool:
    """Check all numeric values in a flat metric mapping for finiteness."""

    for value in mapping.values():
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)) and not math.isfinite(float(value)):
            return False
    return True


__all__ = [
    "finite_numeric_mapping",
    "prepare_truth_blind_rows",
    "programme_stages_from_status",
    "required_stages_completed",
    "score_age_posteriors",
    "score_reaction_families",
    "score_topology",
]
