"""Contracts for integrated groundwater benchmark evidence.

The project uses several kinds of reference information that must not be
pooled as if they were interchangeable labels.  This module provides small,
dependency-light audits for those panels and a deterministic provenance
manifest helper.  It deliberately does not download data or infer a missing
edge, age, or reaction label.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


REFERENCE_TYPES = frozenset(
    {
        "observed_field",
        "calibrated_model_reference",
        "independent_synthetic_truth",
        "conceptual_soft_path",
    }
)

CLAIM_GUARDRAILS: dict[str, str] = {
    "observed_field": (
        "Observed field records support transfer, prediction, and hold-forward "
        "diagnostics only unless an independent age, flow-path, adjacency, or "
        "reaction reference is supplied."
    ),
    "calibrated_model_reference": (
        "A calibrated MODFLOW/MODPATH or tracer-LPM product is a model-conditioned "
        "reference. It can test implementation consistency and model-reference "
        "reproduction, not absolute field truth or universal predictive accuracy."
    ),
    "independent_synthetic_truth": (
        "Independent generated truth supports controlled end-to-end scoring under "
        "the declared generator and does not establish field performance."
    ),
    "conceptual_soft_path": (
        "A conceptual path or geology association is a soft prior/context feature; "
        "it is not a direct-adjacency, transport-time, or reaction-truth label."
    ),
}

_COMMON_NODE_FIELDS = ("node_id",)
_COMMON_EDGE_FIELDS = ("u", "v")
_FIELD_NODE_FIELDS = ("sample_date", "lat", "lon")
_MODEL_FIELDS = ("model_name", "source_doi")
_SYNTHETIC_FIELDS = ("truth_edges", "truth_ages_years", "truth_processes")


def sha256_file(path: Path | str) -> str:
    """Hash one file in streaming blocks."""

    file_path = Path(path)
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _columns(value: Any) -> set[str]:
    """Return column/key names for a dataframe-like or record collection."""

    if value is None:
        return set()
    columns = getattr(value, "columns", None)
    if columns is not None:
        return {str(column) for column in columns}
    if isinstance(value, Mapping):
        return {str(key) for key in value}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for item in value:
            if isinstance(item, Mapping):
                return {str(key) for key in item}
        return set()
    return set()


def _records(value: Any) -> list[Mapping[str, Any]]:
    if value is None:
        return []
    if isinstance(value, Mapping):
        return [value]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [item for item in value if isinstance(item, Mapping)]
    # pandas.DataFrame without importing pandas in this core module
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        try:
            rows = to_dict(orient="records")
        except TypeError:
            rows = []
        return [item for item in rows if isinstance(item, Mapping)]
    return []


def _present(value: Any) -> bool:
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return bool(value)


def _nonempty_ids(value: Any, key: str) -> set[str]:
    out: set[str] = set()
    for row in _records(value):
        item = row.get(key)
        if item is not None and str(item).strip():
            out.add(str(item).strip())
    return out


def _status(missing: Sequence[str], n_rows: int) -> str:
    if n_rows <= 0:
        return "MISSING"
    return "COMPLETE" if not missing else "PARTIAL"


@dataclass(frozen=True)
class ReferenceAudit:
    """Machine-readable audit of one reference panel."""

    panel_id: str
    reference_type: str
    status: str
    n_nodes: int
    n_edges: int
    n_observations: int
    missing_requirements: tuple[str, ...]
    capabilities: Mapping[str, bool]
    claim_guardrail: str
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "panel_id": self.panel_id,
            "reference_type": self.reference_type,
            "status": self.status,
            "n_nodes": self.n_nodes,
            "n_edges": self.n_edges,
            "n_observations": self.n_observations,
            "missing_requirements": list(self.missing_requirements),
            "capabilities": dict(self.capabilities),
            "claim_guardrail": self.claim_guardrail,
            "notes": list(self.notes),
        }


def audit_reference_panel(
    panel_id: str,
    reference_type: str,
    *,
    nodes: Any = None,
    edges: Any = None,
    observations: Any = None,
    metadata: Mapping[str, Any] | None = None,
) -> ReferenceAudit:
    """Audit a panel without filling or guessing any missing field.

    ``metadata`` contains declarations supplied by the archive or caller.  A
    field present in a table is not treated as independent truth unless the
    corresponding metadata flag explicitly says so.
    """

    panel_id = str(panel_id).strip()
    reference_type = str(reference_type).strip()
    if not panel_id:
        raise ValueError("panel_id must be non-empty")
    if reference_type not in REFERENCE_TYPES:
        raise ValueError(
            f"reference_type must be one of {sorted(REFERENCE_TYPES)}, "
            f"got {reference_type!r}"
        )
    metadata = dict(metadata or {})
    node_columns = _columns(nodes)
    edge_columns = _columns(edges)
    observation_columns = _columns(observations)
    n_nodes = len(_records(nodes))
    n_edges = len(_records(edges))
    n_observations = len(_records(observations))

    required: list[str] = []
    components = metadata.get("required_components")
    if components is None:
        components = ("nodes", "edges")
    components = {str(component) for component in components}
    if "nodes" in components:
        required.extend(f"nodes.{field}" for field in _COMMON_NODE_FIELDS)
    if "edges" in components:
        required.extend(f"edges.{field}" for field in _COMMON_EDGE_FIELDS)
    if "observations" in components and n_observations <= 0:
        required.append("observations.rows")
    if reference_type == "observed_field":
        required.extend(f"observations.{field}" for field in _FIELD_NODE_FIELDS)
        required.append("metadata.independent_reference_available")
    elif reference_type == "calibrated_model_reference":
        required.extend(f"metadata.{field}" for field in _MODEL_FIELDS)
    elif reference_type == "independent_synthetic_truth":
        required.extend(f"metadata.{field}" for field in _SYNTHETIC_FIELDS)
        required.append("metadata.generator_independent")
    elif reference_type == "conceptual_soft_path":
        required.append("metadata.source_description")

    missing: list[str] = []
    for item in required:
        scope, field = item.split(".", 1)
        if scope == "nodes":
            ok = field in node_columns
        elif scope == "edges":
            ok = field in edge_columns
        else:
            ok = _present(metadata.get(field))
        if not ok:
            missing.append(item)

    node_ids = _nonempty_ids(nodes, "node_id")
    edge_rows = _records(edges)
    edge_endpoints = {
        str(endpoint).strip()
        for row in edge_rows
        for endpoint in (row.get("u"), row.get("v"))
        if endpoint is not None and str(endpoint).strip()
    }
    orphan_edges = sorted(edge_endpoints - node_ids) if node_ids else []
    if orphan_edges:
        missing.append("edges.endpoint_nodes_present")

    independent_edge_truth = bool(
        metadata.get(
            "independent_edge_truth",
            metadata.get("independent_direct_adjacency_truth", False),
        )
    )
    independent_age_truth = bool(metadata.get("independent_age_truth", False))
    independent_reaction_truth = bool(metadata.get("independent_reaction_truth", False))
    capabilities = {
        "node_identity": "node_id" in node_columns and n_nodes > 0,
        "directed_edge_table": all(field in edge_columns for field in _COMMON_EDGE_FIELDS)
        and n_edges > 0,
        "independent_direct_adjacency_truth": independent_edge_truth,
        "independent_age_truth": independent_age_truth,
        "independent_reaction_truth": independent_reaction_truth,
        "coordinates": all(field in observation_columns for field in ("lat", "lon")),
        "sample_time": "sample_date" in observation_columns,
        "transport_time": bool(metadata.get("transport_time_available", False)),
        "screen_intervals": bool(metadata.get("screen_intervals_available", False)),
    }
    notes: list[str] = []
    if orphan_edges:
        notes.append(f"{len(orphan_edges)} edge endpoints are absent from nodes")
    if reference_type == "calibrated_model_reference":
        notes.append("model-derived reference; do not report as field truth")
    if reference_type == "observed_field" and not independent_edge_truth:
        notes.append("field observations have no independent direct-edge labels")

    return ReferenceAudit(
        panel_id=panel_id,
        reference_type=reference_type,
        status=_status(missing, n_nodes + n_edges + n_observations),
        n_nodes=n_nodes,
        n_edges=n_edges,
        n_observations=n_observations,
        missing_requirements=tuple(sorted(set(missing))),
        capabilities=capabilities,
        claim_guardrail=CLAIM_GUARDRAILS[reference_type],
        notes=tuple(notes),
    )


def audit_integrated_panels(audits: Iterable[ReferenceAudit]) -> dict[str, Any]:
    """Summarise whether panels can be joined for a genuine integrated score."""

    values = list(audits)
    complete = [item for item in values if item.status == "COMPLETE"]
    # A single complete independent-generator panel is internally crosswalked
    # by construction: its node, edge, age, and process truth share the same
    # generated case identifiers.  This is controlled-synthetic integration,
    # not field validation.
    if (
        len(values) == 1
        and values[0].status == "COMPLETE"
        and values[0].reference_type == "independent_synthetic_truth"
    ):
        return {
            "status": "COMPLETE",
            "n_panels": 1,
            "complete_panels": [values[0].panel_id],
            "crosswalk_declared": True,
            "linked_common_node_count": values[0].n_nodes,
            "integrated_scoring_allowed": True,
            "reason": "Node, edge, age, and reaction truth share the independent generator's declared case identifiers.",
            "claim_guardrail": CLAIM_GUARDRAILS["independent_synthetic_truth"],
        }
    # The audit intentionally does not assume that similarly named IDs match;
    # callers supply the overlap count after an explicit crosswalk.  In this
    # generic summary, separate panels are therefore not considered linked.
    return {
        # A collection of complete component panels is still only a partial
        # integrated benchmark until an explicit crosswalk is supplied.
        "status": "PARTIAL",
        "n_panels": len(values),
        "complete_panels": [item.panel_id for item in complete],
        "crosswalk_declared": False,
        "linked_common_node_count": 0,
        "integrated_scoring_allowed": False,
        "reason": (
            "No explicit panel crosswalk was supplied; panel-specific metrics "
            "must remain separate until identifiers, timing, and provenance are aligned."
        ),
        "claim_guardrail": (
            "Do not combine model-derived age, model-derived topology, synthetic "
            "truth, and Ghana observations into one field-accuracy score without "
            "an explicit crosswalk and independent labels."
        ),
    }


def source_manifest_entry(
    path: Path | str,
    *,
    source_id: str,
    role: str,
    doi: str | None = None,
    url: str | None = None,
    license_name: str | None = None,
) -> dict[str, Any]:
    """Return an auditable source entry, failing closed for missing files."""

    file_path = Path(path).resolve()
    entry: dict[str, Any] = {
        "source_id": str(source_id),
        "role": str(role),
        "path": str(file_path),
        "exists": file_path.exists(),
        "doi": doi,
        "url": url,
        "license": license_name,
    }
    if file_path.exists() and file_path.is_file():
        entry["size_bytes"] = file_path.stat().st_size
        entry["sha256"] = sha256_file(file_path)
    else:
        entry["size_bytes"] = None
        entry["sha256"] = None
    return entry


def write_json(path: Path | str, value: Mapping[str, Any]) -> None:
    """Write stable UTF-8 JSON, creating only the requested parent directory."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True, default=str)
        + "\n",
        encoding="utf-8",
    )


def build_run_manifest(
    *,
    run_id: str,
    protocol: str,
    reference_audits: Sequence[ReferenceAudit],
    sources: Sequence[Mapping[str, Any]],
    claim_boundary: str,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a run manifest with explicit panel and integration boundaries."""

    manifest: dict[str, Any] = {
        "schema": "integrated-benchmark-manifest-v1",
        "run_id": str(run_id),
        "protocol": str(protocol),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [dict(source) for source in sources],
        "reference_panels": [audit.to_dict() for audit in reference_audits],
        "integration": audit_integrated_panels(reference_audits),
        "claim_boundary": str(claim_boundary),
    }
    if extra:
        manifest.update(dict(extra))
    return manifest


__all__ = [
    "CLAIM_GUARDRAILS",
    "REFERENCE_TYPES",
    "ReferenceAudit",
    "audit_integrated_panels",
    "audit_reference_panel",
    "build_run_manifest",
    "sha256_file",
    "source_manifest_entry",
    "write_json",
]
