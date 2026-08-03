"""Executable contracts for the programme-level validation workflow.

These small data contracts make the synthetic validation boundary explicit. In
particular, an inference stage must not see generator truth, and every decision
must report either an action, a compatible set, or an abstention reason.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import math
from typing import Iterable, Mapping, Sequence


class DecisionKind(str, Enum):
    """Allowed outputs from a next-measurement or inference decision."""

    ACTION = "ACTION"
    SET_REPORT = "SET_REPORT"
    ABSTAIN = "ABSTAIN"


class IdentifiabilityStatus(str, Enum):
    """Status of the latent quantity being reported."""

    IDENTIFIED = "IDENTIFIED"
    PARTIALLY_IDENTIFIED = "PARTIALLY_IDENTIFIED"
    NON_IDENTIFIABLE = "NON_IDENTIFIABLE"
    MODEL_DISAGREEMENT = "MODEL_DISAGREEMENT"
    UNKNOWN = "UNKNOWN"


class StageStatus(str, Enum):
    """Status values for a programme validation stage."""

    NOT_REQUESTED = "not_requested"
    PENDING = "pending"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


DEFAULT_TRUTH_FIELD_PREFIXES = ("true_", "truth_")
DEFAULT_TRUTH_FIELD_NAMES = frozenset(
    {
        "age_years",
        "travel_age_years",
        "true_age",
        "true_age_years",
        "true_edges",
        "true_ages_years",
        "true_processes",
        "truth_edges",
        "truth_ages",
        "truth_processes",
        "reference_edges",
        "reference_ages",
        "reference_processes",
        "pathline_rows",
        "particle_id",
        "particle_ids",
    }
)


def _finite_number(value: object, *, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return number


def assert_truth_blind(
    rows: Iterable[Mapping[str, object]],
    *,
    forbidden_fields: Iterable[str] = (),
) -> None:
    """Raise when inference rows contain generator truth fields.

    Explicit field names cover known truth payloads. The ``true_`` and
    ``truth_`` prefixes catch newly added truth fields without requiring every
    benchmark runner to maintain a second leak-prone list.
    """

    explicit = {str(field) for field in forbidden_fields}
    leaked: set[str] = set()
    exact_names = {item.lower() for item in DEFAULT_TRUTH_FIELD_NAMES}

    def visit(value: object) -> None:
        if isinstance(value, Mapping):
            for key, nested in value.items():
                key_text = str(key)
                lowered = key_text.lower()
                if (
                    key_text in explicit
                    or lowered in exact_names
                    or any(
                        lowered.startswith(prefix)
                        for prefix in DEFAULT_TRUTH_FIELD_PREFIXES
                    )
                ):
                    leaked.add(key_text)
                visit(nested)
        elif isinstance(value, (list, tuple)):
            for nested in value:
                visit(nested)

    for row in rows:
        visit(row)
    if leaked:
        raise ValueError(
            "Generator truth fields leaked into inference rows: "
            f"{sorted(leaked)}"
        )


@dataclass(frozen=True)
class ProgrammeDecision:
    """One contract-compliant inference or measurement decision."""

    decision_kind: DecisionKind | str
    identifiability: IdentifiabilityStatus | str = IdentifiabilityStatus.UNKNOWN
    reason: str = ""
    measurement: str | None = None
    target: str | None = None
    cost: float | None = None
    expected_utility: float | None = None
    compatible_hypotheses: Sequence[str] = ()
    scenario_count: int = 0
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        try:
            decision_kind = DecisionKind(self.decision_kind)
        except ValueError as exc:
            raise ValueError(
                f"Unknown decision kind: {self.decision_kind!r}."
            ) from exc
        try:
            identifiability = IdentifiabilityStatus(self.identifiability)
        except ValueError as exc:
            raise ValueError(
                f"Unknown identifiability status: {self.identifiability!r}."
            ) from exc

        reason = str(self.reason).strip()
        if not reason:
            raise ValueError("Every programme decision requires a reason.")

        hypotheses = tuple(str(item) for item in self.compatible_hypotheses)
        if decision_kind is DecisionKind.ACTION:
            if not str(self.measurement or "").strip():
                raise ValueError("ACTION requires a measurement.")
            if not str(self.target or "").strip():
                raise ValueError("ACTION requires a target.")
            cost = _finite_number(self.cost, name="cost")
            utility = _finite_number(self.expected_utility, name="expected_utility")
            if cost < 0:
                raise ValueError("cost must be non-negative.")
            object.__setattr__(self, "cost", cost)
            object.__setattr__(self, "expected_utility", utility)
        elif decision_kind is DecisionKind.SET_REPORT:
            if not hypotheses:
                raise ValueError("SET_REPORT requires compatible_hypotheses.")
        elif decision_kind is DecisionKind.ABSTAIN:
            if self.measurement is not None or self.target is not None:
                raise ValueError("ABSTAIN cannot contain an actionable target.")

        if int(self.scenario_count) < 0:
            raise ValueError("scenario_count must be non-negative.")

        object.__setattr__(self, "decision_kind", decision_kind)
        object.__setattr__(self, "identifiability", identifiability)
        object.__setattr__(self, "reason", reason)
        object.__setattr__(self, "compatible_hypotheses", hypotheses)
        object.__setattr__(self, "scenario_count", int(self.scenario_count))

    def to_dict(self) -> dict[str, object]:
        """Return a stable JSON-compatible representation."""

        return {
            "decision_kind": self.decision_kind.value,
            "identifiability": self.identifiability.value,
            "reason": self.reason,
            "measurement": self.measurement,
            "target": self.target,
            "cost": self.cost,
            "expected_utility": self.expected_utility,
            "compatible_hypotheses": list(self.compatible_hypotheses),
            "scenario_count": self.scenario_count,
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class ProgrammeStage:
    """Truth-blind stage record for a synthetic validation run."""

    name: str
    status: StageStatus | str
    truth_blind: bool = True
    truth_fields_seen: Sequence[str] = ()
    detail: str = ""

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        if not name:
            raise ValueError("Programme stages require a name.")
        try:
            status = StageStatus(self.status)
        except ValueError as exc:
            raise ValueError(f"Unknown stage status: {self.status!r}.") from exc
        leaked = tuple(sorted(str(item) for item in self.truth_fields_seen))
        if self.truth_blind and leaked:
            raise ValueError(
                f"Truth-blind stage {name!r} saw truth fields: {list(leaked)}"
            )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "truth_fields_seen", leaked)
        object.__setattr__(self, "detail", str(self.detail))

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "status": self.status.value,
            "truth_blind": self.truth_blind,
            "truth_fields_seen": list(self.truth_fields_seen),
            "detail": self.detail,
        }


@dataclass(frozen=True)
class ProgrammeRun:
    """Top-level contract for one independent synthetic validation run."""

    run_id: str
    generator: str
    generator_independent: bool
    stages: Sequence[ProgrammeStage] = ()
    decisions: Sequence[ProgrammeDecision] = ()
    claim_boundary: str = "controlled-synthetic and conditional"
    truth_released_for_scoring: bool = False

    def __post_init__(self) -> None:
        run_id = str(self.run_id).strip()
        generator = str(self.generator).strip()
        boundary = str(self.claim_boundary).strip()
        if not run_id or not generator:
            raise ValueError("Programme runs require run_id and generator.")
        if not self.generator_independent:
            raise ValueError("Programme validation requires an independent generator.")
        if not boundary:
            raise ValueError("Programme runs require a claim boundary.")
        stage_names = [stage.name for stage in self.stages]
        if len(stage_names) != len(set(stage_names)):
            raise ValueError("Programme stage names must be unique.")
        if self.truth_released_for_scoring:
            raise ValueError(
                "Truth must remain sealed until all inference stages complete."
            )
        object.__setattr__(self, "run_id", run_id)
        object.__setattr__(self, "generator", generator)
        object.__setattr__(self, "claim_boundary", boundary)
        object.__setattr__(self, "stages", tuple(self.stages))
        object.__setattr__(self, "decisions", tuple(self.decisions))

    def to_dict(self) -> dict[str, object]:
        """Return a stable JSON-compatible representation."""

        return {
            "run_id": self.run_id,
            "generator": self.generator,
            "generator_independent": self.generator_independent,
            "claim_boundary": self.claim_boundary,
            "truth_released_for_scoring": self.truth_released_for_scoring,
            "stages": [stage.to_dict() for stage in self.stages],
            "decisions": [decision.to_dict() for decision in self.decisions],
        }
