"""Truth-blind missingness and censoring scenarios for synthetic validation."""

from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Iterable, Mapping


OBSERVATION_STRESS_SCENARIOS = (
    "complete",
    "structured_missingness",
    "left_censoring",
    "combined_stress",
)

_MISSING_GROUPS = {
    "age_tracers": ("tritium_TU", "argon39_pmc", "cfc11_pptv", "sf6_pptv"),
    "source_isotopes": ("d18O_permil", "d2H_permil"),
    "reaction_chemistry": ("SO4", "NO3", "Fe"),
    "head_channel": ("head_meas",),
}
_CENSOR_FIELDS = (
    "tritium_TU",
    "sf6_pptv",
    "cfc11_pptv",
    "cfc113_pptv",
    "NO3",
    "Fe",
)


@dataclass(frozen=True)
class ObservationStressResult:
    """One reproducible observation-process variant."""

    scenario: str
    seed: int
    observations: tuple[dict[str, object], ...]
    missing_count: int
    censored_count: int
    field_counts: Mapping[str, int]

    def to_dict(self) -> dict[str, object]:
        return {
            "scenario": self.scenario,
            "seed": self.seed,
            "n_observations": len(self.observations),
            "missing_count": self.missing_count,
            "censored_count": self.censored_count,
            "field_counts": dict(self.field_counts),
        }


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a censoring limit without observations.")
    position = (len(ordered) - 1) * float(fraction)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return float(ordered[lower] + weight * (ordered[upper] - ordered[lower]))


def _merge_flags(
    row: dict[str, object],
    *,
    scenario: str,
    missing_fields: list[str],
    censored_fields: list[str],
    censor_limits: Mapping[str, float],
) -> None:
    row["observation_flags"] = {
        "scenario": scenario,
        "missing_fields": sorted(set(missing_fields)),
        "censored_fields": sorted(set(censored_fields)),
        "censor_limits": dict(censor_limits),
    }


def _apply_structured_missingness(
    rows: list[dict[str, object]],
    *,
    seed: int,
    scenario: str,
) -> None:
    rng = random.Random(int(seed))
    n_rows = len(rows)
    if n_rows == 0:
        return

    for group_index, (group_name, fields) in enumerate(_MISSING_GROUPS.items()):
        count = max(1, int(round(n_rows * (0.18 + 0.03 * group_index))))
        selected = rng.sample(range(n_rows), min(count, n_rows))
        for row_index in selected:
            row = rows[row_index]
            missing_fields = [
                str(field)
                for field in row.get("observation_flags", {}).get("missing_fields", [])
            ] if isinstance(row.get("observation_flags"), Mapping) else []
            for field in fields:
                if field in row:
                    row[field] = None
                    missing_fields.append(field)
            existing_censored = (
                row.get("observation_flags", {}).get("censored_fields", [])
                if isinstance(row.get("observation_flags"), Mapping)
                else []
            )
            limits = (
                row.get("observation_flags", {}).get("censor_limits", {})
                if isinstance(row.get("observation_flags"), Mapping)
                else {}
            )
            _merge_flags(
                row,
                scenario=scenario,
                missing_fields=missing_fields,
                censored_fields=[str(item) for item in existing_censored],
                censor_limits=limits if isinstance(limits, Mapping) else {},
            )


def _apply_left_censoring(
    rows: list[dict[str, object]],
    *,
    scenario: str,
) -> None:
    limits: dict[str, float] = {}
    for field in _CENSOR_FIELDS:
        values = [
            value
            for row in rows
            for value in [_finite(row.get(field))]
            if value is not None
        ]
        if len(values) < 2:
            continue
        limits[field] = max(0.0, _quantile(values, 0.35))

    for row in rows:
        flags = row.get("observation_flags")
        missing_fields = (
            [str(item) for item in flags.get("missing_fields", [])]
            if isinstance(flags, Mapping)
            else []
        )
        censored_fields = (
            [str(item) for item in flags.get("censored_fields", [])]
            if isinstance(flags, Mapping)
            else []
        )
        row_limits = (
            dict(flags.get("censor_limits", {}))
            if isinstance(flags, Mapping)
            else {}
        )
        for field, limit in limits.items():
            value = _finite(row.get(field))
            if value is not None and value <= limit:
                row[field] = f"<{limit:.8g}"
                censored_fields.append(field)
                row_limits[field] = float(limit)
        _merge_flags(
            row,
            scenario=scenario,
            missing_fields=missing_fields,
            censored_fields=censored_fields,
            censor_limits=row_limits,
        )


def apply_observation_stress(
    observations: Iterable[Mapping[str, object]],
    scenario: str,
    *,
    seed: int,
) -> ObservationStressResult:
    """Apply a preregistered missingness/censoring transform to blind rows.

    The transform sees only the observed values.  It never receives ages,
    routes, process labels, or any other generator truth.  Missing values are
    represented as ``None`` and detection-limit observations as strings such
    as ``"<0.25"`` so HydroSheaf's existing parsing policy is exercised.
    """

    scenario = str(scenario)
    if scenario not in OBSERVATION_STRESS_SCENARIOS:
        raise ValueError(
            f"Unknown observation stress scenario {scenario!r}; "
            f"choose from {OBSERVATION_STRESS_SCENARIOS}."
        )
    rows = [dict(row) for row in observations]
    for row in rows:
        _merge_flags(
            row,
            scenario=scenario,
            missing_fields=[],
            censored_fields=[],
            censor_limits={},
        )

    if scenario in {"structured_missingness", "combined_stress"}:
        _apply_structured_missingness(rows, seed=int(seed), scenario=scenario)
    if scenario in {"left_censoring", "combined_stress"}:
        _apply_left_censoring(rows, scenario=scenario)

    field_counts: dict[str, int] = {}
    missing_count = 0
    censored_count = 0
    for row in rows:
        flags = row.get("observation_flags")
        if not isinstance(flags, Mapping):
            continue
        missing_fields = [str(item) for item in flags.get("missing_fields", [])]
        censored_fields = [str(item) for item in flags.get("censored_fields", [])]
        missing_count += len(set(missing_fields))
        censored_count += len(set(censored_fields))
        for field in set(missing_fields + censored_fields):
            field_counts[field] = field_counts.get(field, 0) + 1

    return ObservationStressResult(
        scenario=scenario,
        seed=int(seed),
        observations=tuple(rows),
        missing_count=missing_count,
        censored_count=censored_count,
        field_counts=dict(sorted(field_counts.items())),
    )


__all__ = [
    "OBSERVATION_STRESS_SCENARIOS",
    "ObservationStressResult",
    "apply_observation_stress",
]
