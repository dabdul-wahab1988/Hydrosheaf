"""Explicit geology metadata contracts for soft, evidence-aware priors.

Mapped geology is contextual metadata.  It is not screened-interval
lithology, hydraulic truth, or a reaction observation.  This module provides
the small normalization boundary needed before such metadata is allowed to
influence a soft reaction prior.

The boundary is intentionally conservative:

* a context is ``matched`` only when a non-empty unit is present and the join
  status says that the assignment is valid;
* ``unmatched`` and ``unknown`` contexts remain usable as audit states but
  never imply a geological preference;
* malformed or conflicting metadata is ``invalid`` and also receives a
  neutral prior downstream;
* no nearest-unit, lithology, or geological interpretation is inferred.

The accepted mapping keys are explicit aliases used by the field-data
contracts.  Callers with another schema should normalize it before invoking
the prior layer instead of relying on an accidental column name.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Literal, Mapping


GeologyStatus = Literal["matched", "unmatched", "unknown", "invalid"]
GeologyRelation = Literal["same_unit", "different_unit", "unknown"]

_MATCHED = frozenset({"matched", "match", "ok", "valid", "known", "joined"})
_UNMATCHED = frozenset(
    {
        "unmatched",
        "no_match",
        "no-match",
        "no match",
        "boundary_review",
        "boundary-review",
        "overlap_review",
        "overlap-review",
    }
)
_UNKNOWN = frozenset(
    {
        "unknown",
        "unavailable",
        "missing",
        "not_joined",
        "not-joined",
        "not joined",
        "not_available",
        "not-available",
        "not available",
        "",
    }
)

# These are the only aliases interpreted by ``normalise_geology_context``.
# Keeping the list public makes the schema boundary reviewable.
GEOLOGY_UNIT_KEYS: tuple[str, ...] = (
    # The field-loader contract's ordered mapped-geology fields.  A
    # stratigraphic unit is the most specific stable key; symbol/code and
    # legend text are progressively weaker but still explicit metadata.
    "geology_stratigraphic_unit",
    "geology_symbol",
    "geology_code_1000",
    "geology_legend_text",
    "geology_unit",
    "geology_group",
    "geology",
    "mapped_geology",
    "lithology",
    "aquifer_unit",
    "aquifer_layer",
)
GEOLOGY_STATUS_KEYS: tuple[str, ...] = (
    "geology_join_status",
    "geology_status",
    "join_status",
)
GEOLOGY_SOURCE_KEYS: tuple[str, ...] = (
    "geology_source",
    "geology_source_hash",
    "source",
)


def _clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.casefold() in {
        "nan",
        "none",
        "null",
        "unknown",
        "not_available",
        "not available",
    }:
        return None
    return text


def _finite(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _status(value: Any) -> GeologyStatus:
    text = (_clean_text(value) or "").lower().replace("__", "_")
    if text in _MATCHED:
        return "matched"
    if text in _UNMATCHED:
        return "unmatched"
    if text in _UNKNOWN:
        return "unknown"
    return "invalid"


def _first_value(mapping: Mapping[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return None


@dataclass(frozen=True)
class GeologyContext:
    """Normalized geology metadata and its explicit join state.

    ``unit`` is retained as supplied (trimmed) and is meaningful for prior
    lookup only when ``status == "matched"``.  ``confidence`` is optional
    metadata; the prior layer never invents it and never uses it unless a
    caller explicitly requests confidence scaling.
    """

    unit: str | None = None
    status: GeologyStatus = "unknown"
    source: str | None = None
    confidence: float | None = None
    reason: str = ""

    def __post_init__(self) -> None:
        unit = _clean_text(self.unit)
        source = _clean_text(self.source)
        status = _status(self.status)
        confidence = self.confidence
        if confidence is not None:
            confidence = _finite(confidence)
            if confidence is None or not 0.0 <= confidence <= 1.0:
                status = "invalid"
                confidence = None
        if status == "matched" and unit is None:
            status = "invalid"
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "confidence", confidence)
        object.__setattr__(self, "reason", str(self.reason or "").strip())

    @property
    def has_informative_unit(self) -> bool:
        """Whether this context can support a non-neutral geological prior."""

        return self.status == "matched" and self.unit is not None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GeologyPairEvidence:
    """Comparison of two normalized contexts without a flow interpretation."""

    upstream: GeologyContext
    downstream: GeologyContext
    relation: GeologyRelation
    status: Literal["known", "ambiguous", "unknown"]
    same_unit: bool | None
    neutral_prior_required: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "upstream": self.upstream.to_dict(),
            "downstream": self.downstream.to_dict(),
            "relation": self.relation,
            "status": self.status,
            "same_unit": self.same_unit,
            "neutral_prior_required": self.neutral_prior_required,
            "reason": self.reason,
        }


def normalise_geology_context(
    value: Any,
    *,
    preferred_key: str | None = None,
) -> GeologyContext:
    """Normalize a declared geology context without inferring missing fields.

    Accepted inputs are an existing :class:`GeologyContext`, a mapping using
    the public alias keys, a non-empty string (treated as a matched unit), or
    ``None`` (unknown).  ``preferred_key`` lets a run's user-supplied geology
    dictionary choose one declared field before the conservative fallback
    precedence is applied.  A mapping with an unrecognized status is invalid.
    An explicit ``NO_MATCH``/``UNMATCHED`` status remains unmatched even if a
    stale unit value is present; it cannot activate a prior.
    """

    if isinstance(value, GeologyContext):
        return value
    if value is None:
        return GeologyContext(status="unknown", reason="geology context not supplied")
    if isinstance(value, str):
        unit = _clean_text(value)
        if unit is None:
            return GeologyContext(status="unknown", reason="geology unit is empty")
        status = _status(unit)
        if status == "unknown":
            return GeologyContext(status="unknown", reason="geology unit is an unknown token")
        if status == "unmatched":
            return GeologyContext(status="unmatched", reason="geology unit is an unmatched token")
        if status == "matched":
            return GeologyContext(
                status="unknown",
                reason="matched status token was supplied without a unit",
            )
        return GeologyContext(unit=unit, status="matched", reason="unit supplied directly")
    if not isinstance(value, Mapping):
        return GeologyContext(status="invalid", reason="unsupported geology context type")

    unit = None
    if preferred_key:
        preferred = str(preferred_key).strip()
        if preferred:
            unit = _clean_text(value.get(preferred))
    if unit is None:
        unit = _clean_text(_first_value(value, GEOLOGY_UNIT_KEYS))
    raw_status = _first_value(value, GEOLOGY_STATUS_KEYS)
    status = (
        "matched"
        if raw_status is None and unit is not None
        else "unknown"
        if raw_status is None
        else _status(raw_status)
    )
    source = _clean_text(_first_value(value, GEOLOGY_SOURCE_KEYS))
    confidence = value.get("geology_confidence", value.get("join_confidence"))
    if status == "matched" and unit is None:
        return GeologyContext(
            status="invalid",
            source=source,
            confidence=confidence if confidence is not None else None,
            reason="matched status lacks a geology unit",
        )
    if status == "unmatched":
        return GeologyContext(
            unit=None,
            status="unmatched",
            source=source,
            confidence=confidence if confidence is not None else None,
            reason="geology join did not produce a valid match",
        )
    if status == "unknown":
        return GeologyContext(
            unit=None,
            status="unknown",
            source=source,
            confidence=confidence if confidence is not None else None,
            reason="geology assignment is unavailable",
        )
    if status == "invalid":
        return GeologyContext(
            unit=unit,
            status="invalid",
            source=source,
            confidence=confidence if confidence is not None else None,
            reason="geology join status is not recognized",
        )
    return GeologyContext(
        unit=unit,
        status="matched",
        source=source,
        confidence=confidence if confidence is not None else None,
        reason="explicit matched geology assignment",
    )


def compare_geology_contexts(upstream: Any, downstream: Any) -> GeologyPairEvidence:
    """Compare two contexts; unknown/unmatched pairs remain neutral."""

    first = normalise_geology_context(upstream)
    second = normalise_geology_context(downstream)
    if first.has_informative_unit and second.has_informative_unit:
        assert first.unit is not None and second.unit is not None
        same = first.unit.casefold() == second.unit.casefold()
        return GeologyPairEvidence(
            upstream=first,
            downstream=second,
            relation="same_unit" if same else "different_unit",
            status="known",
            same_unit=same,
            neutral_prior_required=False,
            reason="both endpoints have explicitly matched geology units",
        )
    if first.status in {"unmatched", "unknown"} or second.status in {"unmatched", "unknown"}:
        return GeologyPairEvidence(
            upstream=first,
            downstream=second,
            relation="unknown",
            status="ambiguous",
            same_unit=None,
            neutral_prior_required=True,
            reason="at least one endpoint has no matched geology unit",
        )
    return GeologyPairEvidence(
        upstream=first,
        downstream=second,
        relation="unknown",
        status="unknown",
        same_unit=None,
        neutral_prior_required=True,
        reason="at least one geology context is invalid",
    )


def add_boundary_confidence(
    value: Any,
    *,
    sigma_m: float,
) -> Any:
    """Attach a bounded map-boundary confidence when it is explicitly available.

    The confidence is a descriptive spatial-quality factor, not a probability
    that the mapped surface unit equals screened-interval lithology.  For a
    matched sidecar record with boundary distance ``d`` it is
    ``1 - exp(-d/sigma_m)``; records without a declared distance are returned
    unchanged.  Existing caller-supplied confidence is never overwritten.
    """

    if not isinstance(value, Mapping):
        return value
    if value.get("geology_confidence") is not None or value.get("join_confidence") is not None:
        return value
    status = _status(_first_value(value, GEOLOGY_STATUS_KEYS))
    distance = _finite(value.get("geology_boundary_distance_m"))
    sigma = _finite(sigma_m)
    if status != "matched" or distance is None or distance < 0.0 or sigma is None or sigma <= 0.0:
        return value
    enriched = dict(value)
    enriched["geology_confidence"] = float(1.0 - math.exp(-distance / sigma))
    return enriched


__all__ = [
    "GEOLOGY_SOURCE_KEYS",
    "GEOLOGY_STATUS_KEYS",
    "GEOLOGY_UNIT_KEYS",
    "GeologyContext",
    "GeologyPairEvidence",
    "GeologyRelation",
    "GeologyStatus",
    "compare_geology_contexts",
    "add_boundary_confidence",
    "normalise_geology_context",
]
