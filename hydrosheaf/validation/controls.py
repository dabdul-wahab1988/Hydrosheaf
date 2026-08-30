"""Truth-blind negative and permutation controls for programme validation.

The controls in this module alter only the observed channels.  They preserve
row identities and the marginal distribution of every permuted field, while
destroying the declared association between an observation and its channel
value.  Generator truth is never an input to a control.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import random
from typing import Iterable, Mapping, Sequence

from .programme_contract import assert_truth_blind


@dataclass(frozen=True)
class ControlResult:
    """One reproducible truth-blind observation-control result."""

    scenario: str
    seed: int
    fields: tuple[str, ...]
    observations: tuple[dict[str, object], ...]
    changed_values: int
    preserved_marginal_counts: Mapping[str, int]

    def to_dict(self) -> dict[str, object]:
        return {
            "scenario": self.scenario,
            "seed": self.seed,
            "fields": list(self.fields),
            "n_observations": len(self.observations),
            "changed_values": self.changed_values,
            "preserved_marginal_counts": dict(self.preserved_marginal_counts),
        }


def _normalise_fields(fields: Iterable[object]) -> tuple[str, ...]:
    normalised = tuple(dict.fromkeys(str(field) for field in fields))
    if not normalised or any(not field.strip() for field in normalised):
        raise ValueError("At least one non-empty observed field is required.")
    return normalised


def _annotate(row: dict[str, object], *, scenario: str, seed: int, fields: Sequence[str]) -> None:
    """Attach auditable control metadata without adding latent truth."""

    existing = row.get("observation_flags")
    flags = dict(existing) if isinstance(existing, Mapping) else {}
    flags["control_scenario"] = str(scenario)
    flags["control_seed"] = int(seed)
    flags["permuted_fields"] = list(fields)
    row["observation_flags"] = flags


def apply_field_permutation(
    observations: Iterable[Mapping[str, object]],
    fields: Iterable[object],
    *,
    seed: int,
    scenario: str = "field_permutation",
) -> ControlResult:
    """Permute declared observed fields across rows using a local RNG.

    The operation is row-safe: identifiers, coordinates, metadata, and
    undeclared channels remain attached to their original row.  A field is
    permuted only among rows where it is present; missing/censored values are
    therefore preserved as part of that field's marginal distribution.
    """

    scenario = str(scenario).strip()
    if not scenario:
        raise ValueError("scenario must be non-empty.")
    selected_fields = _normalise_fields(fields)
    rows = [deepcopy(dict(row)) for row in observations]
    assert_truth_blind(rows)
    rng = random.Random(int(seed))
    changed_values = 0
    marginal_counts: dict[str, int] = {}

    for field in selected_fields:
        present_indices = [index for index, row in enumerate(rows) if field in row]
        values = [rows[index][field] for index in present_indices]
        marginal_counts[field] = len(values)
        shuffled = list(values)
        rng.shuffle(shuffled)
        for index, value in zip(present_indices, shuffled):
            if rows[index].get(field) != value:
                changed_values += 1
            rows[index][field] = value

    for row in rows:
        _annotate(row, scenario=scenario, seed=int(seed), fields=selected_fields)
    assert_truth_blind(rows)
    return ControlResult(
        scenario=scenario,
        seed=int(seed),
        fields=selected_fields,
        observations=tuple(rows),
        changed_values=changed_values,
        preserved_marginal_counts=dict(sorted(marginal_counts.items())),
    )


def make_no_flow_control(
    observations: Iterable[Mapping[str, object]],
    *,
    seed: int,
    head_fields: Sequence[str] = ("head_meas", "hydraulic_head"),
) -> ControlResult:
    """Destroy spatial head association while preserving head marginals.

    This is a negative control for topology inference.  It does not assert
    that a real aquifer has no flow; it tests whether a method reports an
    informative graph when the observed head channels have been detached from
    their original locations.
    """

    return apply_field_permutation(
        observations,
        head_fields,
        seed=int(seed),
        scenario="no_flow_head_permutation",
    )


def make_permutation_control(
    observations: Iterable[Mapping[str, object]],
    fields: Sequence[str],
    *,
    seed: int,
    scenario: str = "declared_field_permutation",
) -> ControlResult:
    """Create a named negative control for any observed channel family."""

    return apply_field_permutation(
        observations,
        fields,
        seed=int(seed),
        scenario=scenario,
    )


__all__ = [
    "ControlResult",
    "apply_field_permutation",
    "make_no_flow_control",
    "make_permutation_control",
]
