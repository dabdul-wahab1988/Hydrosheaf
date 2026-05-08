"""Claim-to-evidence guardrails for Hydrosheaf manuscripts and benchmarks."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Iterable, List, Mapping, Optional


class EvidenceLevel(IntEnum):
    """Ordered evidence levels used to prevent overclaiming.

    The order is deliberate: a stronger level can support weaker claims, but
    weaker evidence should not be used to support stronger manuscript language.
    """

    EXPERIMENTAL = 0
    DEMONSTRATION = 1
    SCREENING = 2
    VALIDATED = 3

    @classmethod
    def parse(cls, value: object) -> "EvidenceLevel":
        if isinstance(value, EvidenceLevel):
            return value
        text = str(value).strip().lower().replace("-", "_")
        aliases = {
            "experimental": cls.EXPERIMENTAL,
            "prototype": cls.EXPERIMENTAL,
            "demonstration": cls.DEMONSTRATION,
            "field_demonstration": cls.DEMONSTRATION,
            "screening": cls.SCREENING,
            "screening_level": cls.SCREENING,
            "validated": cls.VALIDATED,
            "validation": cls.VALIDATED,
        }
        if text not in aliases:
            raise ValueError(f"Unknown evidence level: {value!r}.")
        return aliases[text]


@dataclass(frozen=True)
class ClaimRecord:
    """Trace one manuscript claim back to code, data, metrics, and limits."""

    claim_id: str
    manuscript: str
    claim: str
    supported_level: EvidenceLevel
    requested_level: EvidenceLevel
    code_modules: List[str] = field(default_factory=list)
    datasets: List[str] = field(default_factory=list)
    scripts: List[str] = field(default_factory=list)
    metrics: Mapping[str, object] = field(default_factory=dict)
    guardrail: str = ""
    figure_or_table: str = ""

    @classmethod
    def from_mapping(cls, row: Mapping[str, object]) -> "ClaimRecord":
        def _list(key: str) -> List[str]:
            value = row.get(key, [])
            if value is None:
                return []
            if isinstance(value, str):
                return [item.strip() for item in value.split(";") if item.strip()]
            return [str(item) for item in value]

        return cls(
            claim_id=str(row.get("claim_id", "")),
            manuscript=str(row.get("manuscript", "")),
            claim=str(row.get("claim", "")),
            supported_level=EvidenceLevel.parse(row.get("supported_level", "experimental")),
            requested_level=EvidenceLevel.parse(row.get("requested_level", "experimental")),
            code_modules=_list("code_modules"),
            datasets=_list("datasets"),
            scripts=_list("scripts"),
            metrics=row.get("metrics", {}) if isinstance(row.get("metrics"), Mapping) else {},
            guardrail=str(row.get("guardrail", "")),
            figure_or_table=str(row.get("figure_or_table", "")),
        )

    def is_overclaim(self) -> bool:
        return not evidence_level_allows(self.supported_level, self.requested_level)

    def assessment(self) -> Mapping[str, object]:
        return {
            "claim_id": self.claim_id,
            "manuscript": self.manuscript,
            "supported_level": self.supported_level.name.lower(),
            "requested_level": self.requested_level.name.lower(),
            "overclaim": self.is_overclaim(),
            "guardrail": self.guardrail,
            "figure_or_table": self.figure_or_table,
        }


def evidence_level_allows(
    supported_level: EvidenceLevel,
    requested_level: EvidenceLevel,
) -> bool:
    """Return True when evidence is strong enough for the requested wording."""

    return EvidenceLevel.parse(supported_level) >= EvidenceLevel.parse(requested_level)


def assess_claim_records(
    claims: Iterable[ClaimRecord],
    *,
    manuscript: Optional[str] = None,
) -> List[Mapping[str, object]]:
    """Assess a collection of claim records for unsupported overclaims."""

    out: List[Mapping[str, object]] = []
    for claim in claims:
        if manuscript is not None and claim.manuscript != manuscript:
            continue
        out.append(claim.assessment())
    return out
