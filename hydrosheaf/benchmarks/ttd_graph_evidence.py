"""Evidence contracts and scoring utilities for graph-TTD virtual benchmarks.

The functions here are intentionally independent of a particular TTD solver.
They make the boundary between a truth-blind inference stage and a truth-aware
scoring stage executable, and keep a completed smoke run from being reported
as either field validation or a general-superiority result.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Iterable, Mapping, Sequence

from hydrosheaf.validation.programme_contract import assert_truth_blind


CONTROLLED_SYNTHETIC_CLAIM_BOUNDARY = (
    "Controlled-synthetic recovery and calibrated abstention under declared "
    "generator, observation, and discrepancy scenarios only; not field "
    "validation or universal groundwater-model superiority."
)
FIELD_VALIDATION_STATUS = "DEFERRED"


def _text(value: object, *, name: str) -> str:
    result = str(value).strip()
    if not result:
        raise ValueError(f"{name} must be non-empty.")
    return result


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    return number


def _normalise_edges(edges: Iterable[Sequence[object]]) -> frozenset[tuple[str, str]]:
    normalised: set[tuple[str, str]] = set()
    for index, edge in enumerate(edges):
        if isinstance(edge, (str, bytes)):
            raise ValueError(f"edge {index} must contain two endpoint IDs.")
        pair = tuple(edge)
        if len(pair) != 2:
            raise ValueError(f"edge {index} must contain exactly two endpoint IDs.")
        u = _text(pair[0], name=f"edge {index} source")
        v = _text(pair[1], name=f"edge {index} target")
        if u == v:
            raise ValueError(f"edge {index} is a self-loop.")
        normalised.add((u, v))
    return frozenset(normalised)


def _stable_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str)


def payload_sha256(value: object) -> str:
    """Return a deterministic hash for a JSON-compatible benchmark artifact."""

    return hashlib.sha256(_stable_json(value).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class IntervalScore:
    """Truth-aware score for one interval-valued benchmark estimand.

    ``appropriate_abstention_rate`` is defined only when the virtual scenario
    declares which cases should be unidentifiable.  It measures agreement with
    that declared stress label, not a claim that field non-identifiability is
    known exactly.
    """

    n_total: int
    n_scored: int
    n_abstained: int
    coverage_nonabstained: float | None
    coverage_including_abstention: float
    mean_interval_width: float | None
    abstention_rate: float
    selective_mae: float | None
    appropriate_abstention_rate: float | None

    def to_dict(self) -> dict[str, object]:
        return {
            "n_total": self.n_total,
            "n_scored": self.n_scored,
            "n_abstained": self.n_abstained,
            "coverage_nonabstained": self.coverage_nonabstained,
            "coverage_including_abstention": self.coverage_including_abstention,
            "mean_interval_width": self.mean_interval_width,
            "abstention_rate": self.abstention_rate,
            "selective_mae": self.selective_mae,
            "appropriate_abstention_rate": self.appropriate_abstention_rate,
        }


def score_intervals(
    truth: Sequence[object],
    lower: Sequence[object | None],
    upper: Sequence[object | None],
    abstained: Sequence[bool],
    *,
    point_estimates: Sequence[object | None] | None = None,
    expected_abstention: Sequence[bool] | None = None,
) -> IntervalScore:
    """Score intervals after inference has completed.

    All truth values are consumed only in this scorer.  An abstained target is
    not counted as covered in ``coverage_including_abstention``; this avoids
    rewarding a method merely for withholding all predictions.
    """

    sizes = {len(truth), len(lower), len(upper), len(abstained)}
    if point_estimates is not None:
        sizes.add(len(point_estimates))
    if expected_abstention is not None:
        sizes.add(len(expected_abstention))
    if len(sizes) != 1:
        raise ValueError("truth, intervals, and decision arrays must have equal length.")
    n_total = len(truth)
    if n_total == 0:
        raise ValueError("At least one benchmark target is required.")

    covered = 0
    widths: list[float] = []
    absolute_errors: list[float] = []
    n_abstained = 0
    agreement: list[bool] = []

    for index, value in enumerate(truth):
        actual = _finite(value, name=f"truth[{index}]")
        did_abstain = abstained[index]
        if not isinstance(did_abstain, bool):
            raise ValueError(f"abstained[{index}] must be boolean.")
        if expected_abstention is not None:
            expected = expected_abstention[index]
            if not isinstance(expected, bool):
                raise ValueError(f"expected_abstention[{index}] must be boolean.")
            agreement.append(did_abstain is expected)
        if did_abstain:
            n_abstained += 1
            if lower[index] is not None or upper[index] is not None:
                raise ValueError("Abstained targets must not contain an interval.")
            if point_estimates is not None and point_estimates[index] is not None:
                raise ValueError("Abstained targets must not contain a point estimate.")
            continue

        lo = _finite(lower[index], name=f"lower[{index}]")
        hi = _finite(upper[index], name=f"upper[{index}]")
        if lo > hi:
            raise ValueError(f"lower[{index}] cannot exceed upper[{index}].")
        widths.append(hi - lo)
        if lo <= actual <= hi:
            covered += 1
        if point_estimates is not None:
            estimate = _finite(point_estimates[index], name=f"point_estimates[{index}]")
            absolute_errors.append(abs(estimate - actual))

    n_scored = n_total - n_abstained
    return IntervalScore(
        n_total=n_total,
        n_scored=n_scored,
        n_abstained=n_abstained,
        coverage_nonabstained=(covered / n_scored if n_scored else None),
        coverage_including_abstention=covered / n_total,
        mean_interval_width=(sum(widths) / len(widths) if widths else None),
        abstention_rate=n_abstained / n_total,
        selective_mae=(sum(absolute_errors) / len(absolute_errors) if absolute_errors else None),
        appropriate_abstention_rate=(sum(agreement) / len(agreement) if agreement else None),
    )


def score_prediction_rmse(
    observed: Sequence[object], prediction: Sequence[object]
) -> float:
    """Return RMSE for a held-out tracer series without silently dropping data."""

    if len(observed) != len(prediction):
        raise ValueError("observed and prediction must have equal length.")
    if not observed:
        raise ValueError("At least one held-out observation is required.")
    squared = [
        (_finite(prediction[index], name=f"prediction[{index}]") - _finite(value, name=f"observed[{index}]") ) ** 2
        for index, value in enumerate(observed)
    ]
    return math.sqrt(sum(squared) / len(squared))


def score_topology(
    truth_edges: Iterable[Sequence[object]], selected_edges: Iterable[Sequence[object]]
) -> dict[str, float | int | None]:
    """Score directed topology after a truth-blind selection stage."""

    truth = _normalise_edges(truth_edges)
    selected = _normalise_edges(selected_edges)
    tp = len(truth & selected)
    fp = len(selected - truth)
    fn = len(truth - selected)
    precision = tp / (tp + fp) if (tp + fp) else None
    recall = tp / (tp + fn) if (tp + fn) else None
    return {
        "true_positive": tp,
        "false_positive": fp,
        "false_negative": fn,
        "topology_precision": precision,
        "topology_recall": recall,
        "false_edge_acceptance_rate": fp / len(selected) if selected else 0.0,
    }


@dataclass(frozen=True)
class BenchmarkRecord:
    """One persisted inference result, with no generator truth embedded."""

    case_id: str
    generator_family: str
    regime: str
    scenario: str
    method: str
    held_out: bool
    truth_blind: bool
    metrics: Mapping[str, float | int | None]
    abstention_reason: str | None = None
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("case_id", "generator_family", "regime", "scenario", "method"):
            object.__setattr__(self, name, _text(getattr(self, name), name=name))
        if not isinstance(self.held_out, bool) or not isinstance(self.truth_blind, bool):
            raise ValueError("held_out and truth_blind must be boolean.")
        normalised: dict[str, float | int | None] = {}
        for key, value in self.metrics.items():
            metric = _text(key, name="metric name")
            if value is None:
                normalised[metric] = None
            else:
                number = _finite(value, name=f"metrics[{metric}]")
                normalised[metric] = int(number) if isinstance(value, int) and not isinstance(value, bool) else number
        object.__setattr__(self, "metrics", dict(sorted(normalised.items())))
        reason = None if self.abstention_reason is None else str(self.abstention_reason).strip()
        if reason == "":
            reason = None
        object.__setattr__(self, "abstention_reason", reason)
        object.__setattr__(self, "notes", tuple(str(item) for item in self.notes))

    def to_dict(self) -> dict[str, object]:
        return {
            "case_id": self.case_id,
            "generator_family": self.generator_family,
            "regime": self.regime,
            "scenario": self.scenario,
            "method": self.method,
            "held_out": self.held_out,
            "truth_blind": self.truth_blind,
            "metrics": dict(self.metrics),
            "abstention_reason": self.abstention_reason,
            "notes": list(self.notes),
        }


def assess_claim_readiness(
    records: Iterable[BenchmarkRecord],
    *,
    required_generator_families: Iterable[object],
    required_comparators: Iterable[object],
    required_scenarios: Iterable[object],
) -> dict[str, object]:
    """Audit whether results are complete enough to *adjudicate* a claim.

    This function never returns a performance pass.  It reports either an
    execution pass with a claim still awaiting preregistered adjudication, or
    an explicit abstention describing the missing evidence.
    """

    rows = list(records)
    required_generators = {_text(value, name="required generator family") for value in required_generator_families}
    required_methods = {_text(value, name="required comparator") for value in required_comparators}
    required_conditions = {_text(value, name="required scenario") for value in required_scenarios}
    if not required_generators or not required_methods or not required_conditions:
        raise ValueError("Required generators, comparators, and scenarios must be non-empty.")

    observed_generators = {row.generator_family for row in rows if row.held_out}
    observed_methods = {row.method for row in rows}
    observed_conditions = {row.scenario for row in rows}
    missing_generators = sorted(required_generators - observed_generators)
    missing_methods = sorted(required_methods - observed_methods)
    missing_conditions = sorted(required_conditions - observed_conditions)
    truth_blind_failures = sorted(
        f"{row.case_id}:{row.method}" for row in rows if not row.truth_blind
    )
    no_rows = not rows
    complete_execution = not (
        no_rows
        or missing_generators
        or missing_methods
        or missing_conditions
        or truth_blind_failures
    )
    missing: list[str] = []
    if no_rows:
        missing.append("no_benchmark_records")
    missing.extend(f"generator:{item}" for item in missing_generators)
    missing.extend(f"comparator:{item}" for item in missing_methods)
    missing.extend(f"scenario:{item}" for item in missing_conditions)
    missing.extend(f"truth_blind:{item}" for item in truth_blind_failures)
    return {
        "execution_status": "PASS" if complete_execution else "FAIL",
        "controlled_synthetic_claim_status": (
            "READY_FOR_PREREGISTERED_ADJUDICATION" if complete_execution else "ABSTAIN"
        ),
        "field_validation_status": FIELD_VALIDATION_STATUS,
        "claim_boundary": CONTROLLED_SYNTHETIC_CLAIM_BOUNDARY,
        "n_records": len(rows),
        "n_held_out_records": sum(row.held_out for row in rows),
        "observed_generator_families": sorted(observed_generators),
        "observed_methods": sorted(observed_methods),
        "observed_scenarios": sorted(observed_conditions),
        "missing_requirements": missing,
    }


def assert_observation_view_is_truth_blind(rows: Iterable[Mapping[str, object]]) -> None:
    """Validate an inference-visible observation table without mutating it."""

    assert_truth_blind(rows)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True, default=str) + "\n", encoding="utf-8")


def write_virtual_benchmark_artifacts(
    output_dir: Path | str,
    *,
    run_id: str,
    protocol_path: Path | str,
    config: Mapping[str, object],
    observations: Sequence[Mapping[str, object]],
    records: Sequence[BenchmarkRecord],
    truth_for_scoring: Mapping[str, object],
    generator_provenance: Mapping[str, object],
    readiness: Mapping[str, object],
    overwrite: bool = False,
) -> dict[str, str]:
    """Write an auditable virtual-benchmark artifact set.

    ``truth_for_scoring`` is deliberately written only after the inference
    records and truth-blind observation table have been checked.  The function
    refuses to overwrite a populated output directory unless explicitly asked.
    """

    target = Path(output_dir)
    if target.exists() and any(target.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite populated benchmark directory: {target}"
        )
    observations_copy = [dict(row) for row in observations]
    assert_observation_view_is_truth_blind(observations_copy)
    if any(not row.truth_blind for row in records):
        raise ValueError("Cannot persist benchmark records marked truth_blind=False.")
    protocol = Path(protocol_path)
    if not protocol.exists():
        raise FileNotFoundError(f"Protocol is missing: {protocol}")

    target.mkdir(parents=True, exist_ok=True)
    observation_path = target / "inference_observations.json"
    record_path = target / "inference_records.json"
    truth_path = target / "truth_scoring_only.json"
    provenance_path = target / "generator_provenance.json"
    readiness_path = target / "claim_readiness.json"
    _write_json(observation_path, observations_copy)
    _write_json(record_path, [record.to_dict() for record in records])
    _write_json(truth_path, dict(truth_for_scoring))
    _write_json(provenance_path, dict(generator_provenance))
    _write_json(readiness_path, dict(readiness))

    files = {
        path.name: payload_sha256(path.read_text(encoding="utf-8"))
        for path in (observation_path, record_path, truth_path, provenance_path, readiness_path)
    }
    manifest = {
        "schema": "hydrosheaf-ttd-graph-virtual-run-v1",
        "run_id": _text(run_id, name="run_id"),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_path": str(protocol.resolve()),
        "protocol_sha256": hashlib.sha256(protocol.read_bytes()).hexdigest(),
        "config": dict(config),
        "config_sha256": payload_sha256(config),
        "generator_independent": bool(
            generator_provenance.get("independent_from_hydrosheaf_inference") is True
            or generator_provenance.get("imports_hydrosheaf") is False
        ),
        "truth_blind_contract": "observation and inference-record artifacts were checked before scoring truth was written",
        "claim_boundary": CONTROLLED_SYNTHETIC_CLAIM_BOUNDARY,
        "field_validation_status": FIELD_VALIDATION_STATUS,
        "python": sys.version,
        "platform": platform.platform(),
        "artifacts": files,
    }
    manifest_path = target / "run_manifest.json"
    _write_json(manifest_path, manifest)
    return {
        "observations": str(observation_path),
        "records": str(record_path),
        "truth_scoring_only": str(truth_path),
        "generator_provenance": str(provenance_path),
        "claim_readiness": str(readiness_path),
        "manifest": str(manifest_path),
    }


__all__ = [
    "BenchmarkRecord",
    "CONTROLLED_SYNTHETIC_CLAIM_BOUNDARY",
    "FIELD_VALIDATION_STATUS",
    "IntervalScore",
    "assert_observation_view_is_truth_blind",
    "assess_claim_readiness",
    "payload_sha256",
    "score_intervals",
    "score_prediction_rmse",
    "score_topology",
    "write_virtual_benchmark_artifacts",
]
