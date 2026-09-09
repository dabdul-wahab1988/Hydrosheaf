"""Fail-closed QA audit for the two-tier age-adjacency protocol.

The benchmark has two deliberately different information tiers:

``T1_temporal_order``
    Endpoint ages may test downstream-older compatibility.  Direct adjacency
    is not identifiable and every directness result must be ``ABSTAIN``.

``T2_segment_transport``
    Endpoint ages are paired with independently supplied, edge-specific direct
    and indirect travel-time hypotheses.  A directness score is permitted only
    when both hypotheses, their uncertainty, and their provenance are present.

This module audits a run package, not a model.  It does not fit, calibrate,
threshold, or infer a graph.  Truth is an evaluation-only sidecar and is
never accepted inside prediction records.  The implementation uses only the
Python standard library so it can be run before the full benchmark
environment is installed.

Package contract (``age-adjacency-two-tier-v1``)
-------------------------------------------------
The JSON package contains ``metadata``, ``units``, ``provenance``, ``ages``,
``transport_hypotheses``, ``predictions``, and a ``truth_artifact`` reference.
The latter points to a separately sealed JSON truth file for independent
synthetic runs.  An Aiken/calibrated-model emulation must instead declare that
the truth artifact is unavailable and cannot receive a directness score.

Example
-------
    python audit_age_adjacency_tiers.py \
        --package results/RUN-AGE-ADJACENCY-01/package.json \
        --output results/RUN-AGE-ADJACENCY-01/qa.json
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence


PROTOCOL_ID = "age-adjacency-two-tier-v1"
TRUTH_SCHEMA = "age-adjacency-truth-v1"

TIER_TEMPORAL = "T1_temporal_order"
TIER_SEGMENT_TRANSPORT = "T2_segment_transport"
TIERS = frozenset({TIER_TEMPORAL, TIER_SEGMENT_TRANSPORT})

SOURCE_KINDS = frozenset(
    {
        "observed_field",
        "calibrated_model_reference",
        "independent_synthetic_truth",
        "conceptual_soft_path",
    }
)
CLAIM_TIERS = frozenset(
    {
        "direction_diagnostic",
        "controlled_synthetic_component",
        "calibrated_model_reference",
        "field_transfer_screening",
    }
)
IDENTIFIABILITY_STRATA = frozenset(
    {
        "endpoint_age_order_only",
        "edge_transport_comparison",
        "transport_censored",
        "age_observation_censored",
        "unidentifiable",
    }
)
PREDICTION_STATUSES = frozenset({"scored", "ABSTAIN", "invalid"})
ADJACENCY_STATUSES = frozenset(
    {"scored", "ambiguous", "insufficient_information", "ABSTAIN", "invalid"}
)
AGE_STATUSES = frozenset({"observed", "interval", "censored", "missing", "invalid"})
HYPOTHESES = frozenset({"direct", "indirect"})

# These are evaluation labels or truth-bearing aliases.  They are forbidden
# in model-input declarations and prediction records, but are allowed in the
# separate truth sidecar.  ``direct_probability`` is intentionally not here:
# it is a prediction output, not a truth label.
TRUTH_TOKENS = (
    "truth",
    "true",
    "relation",
    "closure",
    "reachable",
    "label",
    "is_edge",
    "is_direct",
)

REQUIRED_TOP_LEVEL = frozenset(
    {
        "schema",
        "protocol_id",
        "run_id",
        "tier",
        "claim_tier",
        "metadata",
        "units",
        "provenance",
        "ages",
        "transport_hypotheses",
        "predictions",
        "truth_artifact",
    }
)
REQUIRED_METADATA = frozenset(
    {
        "source_kind",
        "truth_sealed",
        "truth_access_mode",
        "candidate_set_frozen",
        "integrated_scoring_allowed",
        "aiken_emulation",
        "prediction_input_columns",
        "evaluation_only_columns",
    }
)
REQUIRED_UNITS = {
    "age_years": "years",
    "age_sigma_years": "years",
    "age_covariance_years2": "years^2",
    "travel_time_years": "years",
    "travel_sigma_years": "years",
    "probability": "1",
}
REQUIRED_AGE_FIELDS = frozenset(
    {
        "case_id",
        "node_id",
        "age_years",
        "age_sigma_years",
        "age_status",
        "source_id",
    }
)
REQUIRED_TRANSPORT_FIELDS = frozenset(
    {
        "case_id",
        "edge_id",
        "u",
        "v",
        "hypothesis",
        "travel_time_years",
        "travel_sigma_years",
        "evidence_source_id",
        "independent_of_endpoint_age",
        "path_basis",
        "intermediate_nodes",
    }
)
REQUIRED_PREDICTION_FIELDS = frozenset(
    {
        "case_id",
        "edge_id",
        "u",
        "v",
        "direction_probability",
        "prediction_status",
        "adjacency_status",
        "identifiability_stratum",
    }
)
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class TierAuditError(ValueError):
    """A package violates the machine-checkable scientific contract."""


class _Audit:
    """Collect deterministic errors/warnings without throwing early."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(str(message))

    def warning(self, message: str) -> None:
        self.warnings.append(str(message))

    def finish(self) -> None:
        self.errors[:] = sorted(set(self.errors))
        self.warnings[:] = sorted(set(self.warnings))


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _nonempty(value: Any) -> bool:
    return value is not None and bool(str(value).strip())


def _finite(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _nonnegative(value: Any) -> bool:
    number = _finite(value)
    return number is not None and number >= 0.0


def _boolean(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "yes", "1"}:
            return True
        if lowered in {"false", "no", "0"}:
            return False
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    return None


def _required_keys(
    value: Any,
    required: Sequence[str] | set[str] | frozenset[str],
    location: str,
    audit: _Audit,
) -> None:
    if not _is_mapping(value):
        audit.error(f"{location} must be an object")
        return
    missing = sorted(set(required) - set(value))
    for field in missing:
        audit.error(f"{location}.{field} is required")


def _forbidden_truth_keys(value: Any, location: str, audit: _Audit) -> None:
    """Reject truth-bearing keys in any prediction/input declaration."""

    if isinstance(value, Mapping):
        for key, nested in value.items():
            key_text = str(key).strip().lower()
            if any(token in key_text for token in TRUTH_TOKENS):
                audit.error(
                    f"{location}.{key}: truth-bearing field is not allowed in model input/prediction"
                )
            _forbidden_truth_keys(nested, f"{location}.{key}", audit)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for index, nested in enumerate(value):
            _forbidden_truth_keys(nested, f"{location}[{index}]", audit)


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_relative_artifact(package_path: Path, relative_path: Any) -> Path | None:
    if not isinstance(relative_path, str) or not relative_path.strip():
        return None
    candidate = Path(relative_path)
    if candidate.is_absolute():
        return None
    root = package_path.parent.resolve()
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        return None
    return resolved


def _validate_provenance(package: Mapping[str, Any], audit: _Audit) -> None:
    provenance = package.get("provenance")
    _required_keys(
        provenance,
        {"created_utc", "protocol_hash", "generator", "inference", "environment", "inputs"},
        "provenance",
        audit,
    )
    if not _is_mapping(provenance):
        return
    created = provenance.get("created_utc")
    if not isinstance(created, str) or not created.strip():
        audit.error("provenance.created_utc must be a non-empty ISO-8601 string")
    else:
        try:
            datetime.fromisoformat(created.replace("Z", "+00:00"))
        except ValueError:
            audit.error("provenance.created_utc must be parseable ISO-8601")
    protocol_hash = provenance.get("protocol_hash")
    if not isinstance(protocol_hash, str) or not SHA256_RE.fullmatch(protocol_hash):
        audit.error("provenance.protocol_hash must be a 64-character SHA-256")
    for section in ("generator", "inference", "environment"):
        value = provenance.get(section)
        if not _is_mapping(value):
            audit.error(f"provenance.{section} must be an object")
            continue
        if section == "environment":
            required = {"runtime", "dependency_lock_hash"}
        else:
            required = {"name", "version", "revision"}
        _required_keys(value, required, f"provenance.{section}", audit)
        for key in required:
            if key.endswith("hash"):
                if not isinstance(value.get(key), str) or not SHA256_RE.fullmatch(value[key]):
                    audit.error(f"provenance.{section}.{key} must be a 64-character SHA-256")
            elif key in {"name", "version", "revision", "runtime"} and not _nonempty(value.get(key)):
                audit.error(f"provenance.{section}.{key} must be non-empty")
    inputs = provenance.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        audit.error("provenance.inputs must be a non-empty list")
        return
    source_ids: set[str] = set()
    for index, item in enumerate(inputs):
        location = f"provenance.inputs[{index}]"
        _required_keys(item, {"source_id", "role", "path", "sha256", "size_bytes"}, location, audit)
        if not _is_mapping(item):
            continue
        source_id = str(item.get("source_id", "")).strip()
        if not source_id:
            audit.error(f"{location}.source_id must be non-empty")
        elif source_id in source_ids:
            audit.error(f"{location}.source_id is duplicated: {source_id}")
        source_ids.add(source_id)
        for key in ("role", "path"):
            if not _nonempty(item.get(key)):
                audit.error(f"{location}.{key} must be non-empty")
        if not isinstance(item.get("sha256"), str) or not SHA256_RE.fullmatch(item["sha256"]):
            audit.error(f"{location}.sha256 must be a 64-character SHA-256")
        size = item.get("size_bytes")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            audit.error(f"{location}.size_bytes must be a non-negative integer")


def _validate_metadata(package: Mapping[str, Any], audit: _Audit) -> None:
    metadata = package.get("metadata")
    _required_keys(metadata, REQUIRED_METADATA, "metadata", audit)
    if not _is_mapping(metadata):
        return
    source_kind = metadata.get("source_kind")
    if source_kind not in SOURCE_KINDS:
        audit.error(f"metadata.source_kind must be one of {sorted(SOURCE_KINDS)}")
    for key in ("truth_sealed", "candidate_set_frozen", "integrated_scoring_allowed", "aiken_emulation"):
        if _boolean(metadata.get(key)) is None:
            audit.error(f"metadata.{key} must be boolean")
    if metadata.get("truth_access_mode") != "evaluation_only":
        audit.error("metadata.truth_access_mode must be 'evaluation_only'")
    input_columns = metadata.get("prediction_input_columns")
    eval_columns = metadata.get("evaluation_only_columns")
    if not isinstance(input_columns, list) or not all(_nonempty(item) for item in input_columns):
        audit.error("metadata.prediction_input_columns must be a non-empty list")
    if not isinstance(eval_columns, list) or not all(_nonempty(item) for item in eval_columns):
        audit.error("metadata.evaluation_only_columns must be a non-empty list")
    if isinstance(input_columns, list) and isinstance(eval_columns, list):
        overlap = sorted(set(map(str, input_columns)) & set(map(str, eval_columns)))
        if overlap:
            audit.error("prediction input and evaluation-only columns overlap: " + ", ".join(overlap))
    _forbidden_truth_keys(input_columns, "metadata.prediction_input_columns", audit)
    # The declaration itself is a model-input surface; reject aliases even if
    # they are hidden in a nested list/object.
    _forbidden_truth_keys(metadata.get("prediction_features", {}), "metadata.prediction_features", audit)


def _validate_units(package: Mapping[str, Any], audit: _Audit) -> None:
    units = package.get("units")
    _required_keys(units, REQUIRED_UNITS, "units", audit)
    if not _is_mapping(units):
        return
    for key, expected in REQUIRED_UNITS.items():
        if units.get(key) != expected:
            audit.error(f"units.{key} must be exactly {expected!r}")


def _validate_age_rows(package: Mapping[str, Any], audit: _Audit) -> tuple[set[tuple[str, str]], dict[tuple[str, str], Mapping[str, Any]]]:
    rows = package.get("ages")
    if not isinstance(rows, list) or not rows:
        audit.error("ages must be a non-empty list")
        return set(), {}
    provenance = package.get("provenance")
    input_source_ids = {
        str(item.get("source_id", "")).strip()
        for item in provenance.get("inputs", [])
        if _is_mapping(item)
    } if _is_mapping(provenance) and isinstance(provenance.get("inputs"), list) else set()
    keys: set[tuple[str, str]] = set()
    by_key: dict[tuple[str, str], Mapping[str, Any]] = {}
    for index, row in enumerate(rows):
        location = f"ages[{index}]"
        _required_keys(row, REQUIRED_AGE_FIELDS, location, audit)
        if not _is_mapping(row):
            continue
        case_id = str(row.get("case_id", "")).strip()
        node_id = str(row.get("node_id", "")).strip()
        key = (case_id, node_id)
        if not case_id or not node_id:
            audit.error(f"{location}.case_id and node_id must be non-empty")
        elif key in keys:
            audit.error(f"{location} duplicates age record {key!r}")
        keys.add(key)
        by_key[key] = row
        status = row.get("age_status")
        if status not in AGE_STATUSES:
            audit.error(f"{location}.age_status must be one of {sorted(AGE_STATUSES)}")
        age = _finite(row.get("age_years"))
        sigma = _finite(row.get("age_sigma_years"))
        if status not in {"missing", "invalid"} and age is None:
            audit.error(f"{location}.age_years must be finite for status {status!r}")
        if sigma is None or sigma < 0.0:
            audit.error(f"{location}.age_sigma_years must be finite and non-negative")
        if status in {"interval", "censored"}:
            lower = _finite(row.get("age_lower_years"))
            upper = _finite(row.get("age_upper_years"))
            if lower is None or upper is None or lower > upper:
                audit.error(f"{location} interval/censoring requires finite age_lower_years <= age_upper_years")
        source_id = str(row.get("source_id", "")).strip()
        if not source_id:
            audit.error(f"{location}.source_id must be non-empty")
        elif source_id not in input_source_ids:
            audit.error(f"{location}.source_id is not declared in provenance.inputs")
        _forbidden_truth_keys(row, location, audit)
    return keys, by_key


def _validate_transport_rows(
    package: Mapping[str, Any],
    audit: _Audit,
) -> dict[tuple[str, str, str, str], Mapping[str, Any]]:
    rows = package.get("transport_hypotheses")
    if not isinstance(rows, list):
        audit.error("transport_hypotheses must be a list")
        return {}
    by_key: dict[tuple[str, str, str, str], Mapping[str, Any]] = {}
    provenance = package.get("provenance")
    input_source_ids = {
        str(item.get("source_id", "")).strip()
        for item in provenance.get("inputs", [])
        if _is_mapping(item)
    } if _is_mapping(provenance) and isinstance(provenance.get("inputs"), list) else set()
    for index, row in enumerate(rows):
        location = f"transport_hypotheses[{index}]"
        _required_keys(row, REQUIRED_TRANSPORT_FIELDS, location, audit)
        if not _is_mapping(row):
            continue
        case_id = str(row.get("case_id", "")).strip()
        edge_id = str(row.get("edge_id", "")).strip()
        u = str(row.get("u", "")).strip()
        v = str(row.get("v", "")).strip()
        hypothesis = row.get("hypothesis")
        key = (case_id, u, v, str(hypothesis))
        if not case_id or not edge_id or not u or not v or u == v:
            audit.error(f"{location} requires non-empty case/edge/endpoints with u != v")
        if hypothesis not in HYPOTHESES:
            audit.error(f"{location}.hypothesis must be 'direct' or 'indirect'")
        if key in by_key:
            audit.error(f"{location} duplicates transport hypothesis {key!r}")
        by_key[key] = row
        if not _nonnegative(row.get("travel_time_years")):
            audit.error(f"{location}.travel_time_years must be finite and non-negative")
        sigma = _finite(row.get("travel_sigma_years"))
        if sigma is None or sigma < 0.0:
            audit.error(f"{location}.travel_sigma_years must be finite and non-negative")
        if not _nonempty(row.get("evidence_source_id")):
            audit.error(f"{location}.evidence_source_id must be non-empty")
        elif str(row.get("evidence_source_id")).strip() not in input_source_ids:
            audit.error(f"{location}.evidence_source_id is not declared in provenance.inputs")
        if _boolean(row.get("independent_of_endpoint_age")) is None:
            audit.error(f"{location}.independent_of_endpoint_age must be boolean")
        if not _nonempty(row.get("path_basis")):
            audit.error(f"{location}.path_basis must be non-empty")
        intermediate = row.get("intermediate_nodes")
        if not isinstance(intermediate, list) or any(not _nonempty(item) for item in intermediate):
            audit.error(f"{location}.intermediate_nodes must be a list of non-empty node IDs")
        elif hypothesis == "direct" and intermediate:
            audit.error(f"{location}: direct hypothesis cannot contain intermediate_nodes")
        elif hypothesis == "indirect" and not intermediate:
            audit.error(f"{location}: indirect hypothesis requires intermediate_nodes")
        _forbidden_truth_keys(row, location, audit)
    return by_key


def _validate_prediction_rows(
    package: Mapping[str, Any],
    audit: _Audit,
) -> tuple[set[tuple[str, str, str]], Counter[str]]:
    rows = package.get("predictions")
    if not isinstance(rows, list) or not rows:
        audit.error("predictions must be a non-empty list")
        return set(), Counter()
    keys: set[tuple[str, str, str]] = set()
    statuses: Counter[str] = Counter()
    for index, row in enumerate(rows):
        location = f"predictions[{index}]"
        _required_keys(row, REQUIRED_PREDICTION_FIELDS, location, audit)
        if not _is_mapping(row):
            continue
        # Reject relation labels and truth aliases even if they are not in the
        # declared input columns.  This prevents accidental post-inference
        # joins from masquerading as model outputs.
        _forbidden_truth_keys(row, location, audit)
        case_id = str(row.get("case_id", "")).strip()
        u = str(row.get("u", "")).strip()
        v = str(row.get("v", "")).strip()
        key = (case_id, u, v)
        if not case_id or not u or not v or u == v:
            audit.error(f"{location} requires non-empty case/endpoints with u != v")
        if key in keys:
            audit.error(f"{location} duplicates candidate pair {key!r}")
        keys.add(key)
        status = row.get("prediction_status")
        adjacency_status = row.get("adjacency_status")
        strata = row.get("identifiability_stratum")
        statuses[str(status)] += 1
        if status not in PREDICTION_STATUSES:
            audit.error(f"{location}.prediction_status must be one of {sorted(PREDICTION_STATUSES)}")
        if adjacency_status not in ADJACENCY_STATUSES:
            audit.error(f"{location}.adjacency_status must be one of {sorted(ADJACENCY_STATUSES)}")
        if strata not in IDENTIFIABILITY_STRATA:
            audit.error(f"{location}.identifiability_stratum must be one of {sorted(IDENTIFIABILITY_STRATA)}")
        direction_probability = _finite(row.get("direction_probability"))
        if direction_probability is None or not 0.0 <= direction_probability <= 1.0:
            audit.error(f"{location}.direction_probability must be finite and in [0, 1]")
        direct_probability = row.get("direct_probability")
        if direct_probability is not None:
            direct_number = _finite(direct_probability)
            if direct_number is None or not 0.0 <= direct_number <= 1.0:
                audit.error(f"{location}.direct_probability must be null or finite in [0, 1]")
        log_bf = row.get("log_bayes_factor_direct_vs_indirect")
        if log_bf is not None and _finite(log_bf) is None:
            audit.error(f"{location}.log_bayes_factor_direct_vs_indirect must be null or finite")
        if status == "scored" and adjacency_status != "scored":
            audit.error(f"{location}: scored prediction must have adjacency_status='scored'")
        if status == "ABSTAIN" and direct_probability is not None:
            audit.error(f"{location}: ABSTAIN prediction cannot carry direct_probability")
        if status == "ABSTAIN" and log_bf is not None:
            audit.error(f"{location}: ABSTAIN prediction cannot carry a Bayes factor")
        if status == "invalid" and adjacency_status != "invalid":
            audit.error(f"{location}: invalid prediction must have adjacency_status='invalid'")
    return keys, statuses


def _validate_truth_artifact(
    package_path: Path,
    package: Mapping[str, Any],
    audit: _Audit,
) -> tuple[bool, dict[str, Any]]:
    reference = package.get("truth_artifact")
    _required_keys(reference, {"available", "sealed", "role"}, "truth_artifact", audit)
    if not _is_mapping(reference):
        return False, {}
    available = _boolean(reference.get("available"))
    sealed = _boolean(reference.get("sealed"))
    if available is None:
        audit.error("truth_artifact.available must be boolean")
        available = False
    if sealed is not True:
        audit.error("truth_artifact.sealed must be true")
    if available:
        if reference.get("role") != "evaluation_only":
            audit.error("available truth_artifact.role must be 'evaluation_only'")
        _required_keys(reference, {"path", "sha256"}, "truth_artifact", audit)
        truth_path = _safe_relative_artifact(package_path, reference.get("path"))
        if truth_path is None:
            audit.error("truth_artifact.path must be a safe path relative to the package")
            return True, {}
        expected_hash = reference.get("sha256")
        if not isinstance(expected_hash, str) or not SHA256_RE.fullmatch(expected_hash):
            audit.error("truth_artifact.sha256 must be a 64-character SHA-256")
        elif not truth_path.is_file():
            audit.error(f"truth_artifact.path does not exist: {reference.get('path')}")
        else:
            actual_hash = _hash_file(truth_path)
            if actual_hash != expected_hash:
                audit.error("truth_artifact.sha256 does not match the sidecar")
        if truth_path.is_file():
            try:
                truth = json.loads(truth_path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                audit.error(f"truth_artifact.sidecar is not valid UTF-8 JSON: {exc}")
                return True, {}
            if not isinstance(truth, Mapping):
                audit.error("truth_artifact.sidecar must contain an object")
                return True, {}
            _required_keys(
                truth,
                {"schema", "run_id", "sealed", "direct_edge_count", "reachable_ordered_pair_count"},
                "truth_artifact.sidecar",
                audit,
            )
            if truth.get("schema") != TRUTH_SCHEMA:
                audit.error(f"truth_artifact.sidecar.schema must be {TRUTH_SCHEMA!r}")
            if truth.get("run_id") != package.get("run_id"):
                audit.error("truth_artifact.sidecar.run_id must match package.run_id")
            if _boolean(truth.get("sealed")) is not True:
                audit.error("truth_artifact.sidecar.sealed must be true")
            for key in ("direct_edge_count", "reachable_ordered_pair_count"):
                value = truth.get(key)
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    audit.error(f"truth_artifact.sidecar.{key} must be a non-negative integer")
            return True, dict(truth)
        return True, {}
    if reference.get("role") != "unavailable":
        audit.error("unavailable truth_artifact.role must be 'unavailable'")
    if not _nonempty(reference.get("reason")):
        audit.error("unavailable truth_artifact requires a non-empty reason")
    if reference.get("path") is not None or reference.get("sha256") is not None:
        audit.error("unavailable truth_artifact must not carry a path or hash")
    return False, {}


def _validate_tier_rules(
    package: Mapping[str, Any],
    package_path: Path,
    age_rows: Mapping[tuple[str, str], Mapping[str, Any]],
    transport_rows: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    prediction_keys: set[tuple[str, str, str]],
    truth_available: bool,
    audit: _Audit,
) -> dict[str, Any]:
    tier = package.get("tier")
    claim_tier = package.get("claim_tier")
    metadata = package.get("metadata", {})
    source_kind = metadata.get("source_kind") if _is_mapping(metadata) else None
    aiken = metadata.get("aiken_emulation") is True if _is_mapping(metadata) else False
    integrated = metadata.get("integrated_scoring_allowed") is True if _is_mapping(metadata) else False
    predictions = package.get("predictions", [])
    scored_count = 0
    abstain_count = 0
    for index, row in enumerate(predictions if isinstance(predictions, list) else []):
        if not _is_mapping(row):
            continue
        location = f"predictions[{index}]"
        case_id = str(row.get("case_id", "")).strip()
        u = str(row.get("u", "")).strip()
        v = str(row.get("v", "")).strip()
        if (case_id, u) not in age_rows or (case_id, v) not in age_rows:
            audit.error(f"{location}: both endpoint age records are required")
        status = row.get("prediction_status")
        if status == "scored":
            scored_count += 1
        elif status == "ABSTAIN":
            abstain_count += 1
        if tier == TIER_TEMPORAL:
            if row.get("identifiability_stratum") not in {"endpoint_age_order_only", "unidentifiable", "age_observation_censored"}:
                audit.error(f"{location}: T1 requires an endpoint-age-only identifiability stratum")
            if status == "scored" or row.get("adjacency_status") == "scored":
                audit.error(f"{location}: T1 cannot score direct adjacency")
            if row.get("direct_probability") is not None or row.get("log_bayes_factor_direct_vs_indirect") is not None:
                audit.error(f"{location}: T1 must not carry a directness score")
        elif tier == TIER_SEGMENT_TRANSPORT:
            direct_key = (case_id, u, v, "direct")
            indirect_key = (case_id, u, v, "indirect")
            direct = transport_rows.get(direct_key)
            indirect = transport_rows.get(indirect_key)
            if status == "scored":
                if direct is None or indirect is None:
                    audit.error(f"{location}: scored T2 prediction requires direct and indirect hypotheses")
                else:
                    for hypothesis, item in (("direct", direct), ("indirect", indirect)):
                        if item.get("hypothesis") != hypothesis:
                            audit.error(f"{location}: malformed {hypothesis} transport hypothesis")
                        if item.get("independent_of_endpoint_age") is not True:
                            audit.error(f"{location}: scored T2 transport must be independent_of_endpoint_age=true")
                        sigma = _finite(item.get("travel_sigma_years"))
                        if sigma is None or sigma <= 0.0:
                            audit.error(f"{location}: scored T2 {hypothesis} uncertainty must be > 0")
                    if (
                        _finite(direct.get("travel_time_years")) == _finite(indirect.get("travel_time_years"))
                        and _finite(direct.get("travel_sigma_years")) == _finite(indirect.get("travel_sigma_years"))
                    ):
                        audit.error(f"{location}: equal direct/indirect hypotheses are non-discriminating")
                if row.get("direct_probability") is None or row.get("log_bayes_factor_direct_vs_indirect") is None:
                    audit.error(f"{location}: scored T2 prediction requires direct_probability and log Bayes factor")
                if row.get("identifiability_stratum") != "edge_transport_comparison":
                    audit.error(f"{location}: scored T2 prediction requires edge_transport_comparison stratum")
            elif status == "ABSTAIN" and row.get("identifiability_stratum") == "edge_transport_comparison":
                # It is valid to abstain after transport metadata are present
                # (for example due to censoring), but the row should say why.
                flags = row.get("flags", [])
                if not isinstance(flags, list) or not flags:
                    audit.error(f"{location}: T2 ABSTAIN requires a non-empty flags list")
        else:
            audit.error(f"tier must be one of {sorted(TIERS)}")

    if tier == TIER_TEMPORAL and transport_rows:
        audit.error("T1_temporal_order must have an empty transport_hypotheses list")
    if tier == TIER_SEGMENT_TRANSPORT and not transport_rows and not aiken:
        audit.error("T2_segment_transport requires transport_hypotheses")

    # Claim and source restrictions are deliberately stricter than the
    # numerical checks.  They define what a passing run is allowed to say.
    if claim_tier not in CLAIM_TIERS:
        audit.error(f"claim_tier must be one of {sorted(CLAIM_TIERS)}")
    if source_kind == "independent_synthetic_truth":
        if claim_tier != "controlled_synthetic_component":
            audit.error("independent_synthetic_truth requires controlled_synthetic_component claim_tier")
        if not truth_available:
            audit.error("independent_synthetic_truth requires a sealed truth sidecar")
        if integrated is not True:
            audit.error("independent_synthetic_truth requires integrated_scoring_allowed=true")
    if source_kind == "calibrated_model_reference":
        if claim_tier != "calibrated_model_reference":
            audit.error("calibrated_model_reference requires calibrated_model_reference claim_tier")
        if not aiken:
            audit.warning("calibrated_model_reference is not marked aiken_emulation; verify the archive-specific restriction")
    if source_kind in {"observed_field", "conceptual_soft_path"} and integrated:
        audit.error(f"{source_kind} cannot set integrated_scoring_allowed=true without independent truth")

    if aiken:
        if source_kind != "calibrated_model_reference":
            audit.error("aiken_emulation requires source_kind='calibrated_model_reference'")
        if integrated:
            audit.error("Aiken emulation must set integrated_scoring_allowed=false")
        if truth_available:
            audit.error("Aiken emulation must not attach independent truth")
        if scored_count:
            audit.error("Aiken emulation cannot score direct adjacency")
        if claim_tier != "calibrated_model_reference":
            audit.error("Aiken emulation can only use calibrated_model_reference claim_tier")

    return {
        "tier": tier,
        "claim_tier": claim_tier,
        "source_kind": source_kind,
        "aiken_emulation": bool(aiken),
        "n_predictions": len(predictions) if isinstance(predictions, list) else 0,
        "n_scored_directness": scored_count,
        "n_abstain_directness": abstain_count,
        "directness_estimand": "ABSTAIN" if tier == TIER_TEMPORAL or aiken else "direct_vs_indirect",
        "prediction_pair_count": len(prediction_keys),
    }


def audit_package(path: str | Path) -> dict[str, Any]:
    """Audit one age-adjacency package and return a deterministic report."""

    package_path = Path(path)
    audit = _Audit()
    package: Mapping[str, Any]
    try:
        value = json.loads(package_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return {
            "protocol": PROTOCOL_ID,
            "package": str(package_path),
            "valid": False,
            "decision": "FAIL",
            "errors": [f"package is not valid UTF-8 JSON: {exc}"],
            "warnings": [],
        }
    if not isinstance(value, Mapping):
        return {
            "protocol": PROTOCOL_ID,
            "package": str(package_path),
            "valid": False,
            "decision": "FAIL",
            "errors": ["package root must be an object"],
            "warnings": [],
        }
    package = value
    missing = sorted(REQUIRED_TOP_LEVEL - set(package))
    for field in missing:
        audit.error(f"package.{field} is required")
    if package.get("schema") != PROTOCOL_ID:
        audit.error(f"package.schema must be {PROTOCOL_ID!r}")
    if package.get("protocol_id") != PROTOCOL_ID:
        audit.error(f"package.protocol_id must be {PROTOCOL_ID!r}")
    if not _nonempty(package.get("run_id")):
        audit.error("package.run_id must be non-empty")
    if package.get("tier") not in TIERS:
        audit.error(f"package.tier must be one of {sorted(TIERS)}")
    _validate_metadata(package, audit)
    _validate_units(package, audit)
    _validate_provenance(package, audit)
    age_keys, age_rows = _validate_age_rows(package, audit)
    transport_rows = _validate_transport_rows(package, audit)
    prediction_keys, prediction_statuses = _validate_prediction_rows(package, audit)
    truth_available, truth = _validate_truth_artifact(package_path, package, audit)
    tier_summary = _validate_tier_rules(
        package,
        package_path,
        age_rows,
        transport_rows,
        prediction_keys,
        truth_available,
        audit,
    )
    # Keep this check after row validation so it reports every bad candidate,
    # not only the first one.
    if isinstance(package.get("predictions"), list):
        for index, row in enumerate(package["predictions"]):
            if not _is_mapping(row):
                continue
            case_id = str(row.get("case_id", "")).strip()
            u = str(row.get("u", "")).strip()
            v = str(row.get("v", "")).strip()
            if (case_id, u) not in age_keys or (case_id, v) not in age_keys:
                audit.error(f"predictions[{index}]: endpoint is absent from ages")

    # Truth sidecar is never joined into prediction rows by this audit.  It is
    # summarized only as a denominator for a QA report.
    audit.finish()
    valid = not audit.errors
    if not valid:
        decision = "FAIL"
    elif tier_summary.get("aiken_emulation"):
        decision = "PASS_AIKEN_EMULATION_ONLY"
    elif tier_summary.get("tier") == TIER_TEMPORAL:
        decision = "PASS_T1_TEMPORAL_ONLY"
    elif tier_summary.get("n_scored_directness", 0) > 0:
        decision = "PASS_T2_DIRECTNESS_SCORING"
    else:
        decision = "PASS_T2_NO_DIRECTNESS_SCORE"
    return {
        "protocol": PROTOCOL_ID,
        "package": str(package_path),
        "run_id": package.get("run_id"),
        "valid": valid,
        "decision": decision,
        "errors": audit.errors,
        "warnings": audit.warnings,
        "summary": {
            **tier_summary,
            "n_age_records": len(age_keys),
            "n_transport_hypotheses": len(transport_rows),
            "prediction_status_counts": dict(sorted(prediction_statuses.items())),
            "truth_available": truth_available,
            "truth_direct_edge_count": truth.get("direct_edge_count") if truth else None,
            "truth_reachable_ordered_pair_count": truth.get("reachable_ordered_pair_count") if truth else None,
            "units": dict(package.get("units", {})) if _is_mapping(package.get("units")) else {},
        },
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", required=True, type=Path, help="age-adjacency package JSON")
    parser.add_argument("--output", type=Path, default=None, help="optional JSON audit output")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = audit_package(args.package)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["valid"] else 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["PROTOCOL_ID", "TierAuditError", "audit_package", "main"]
