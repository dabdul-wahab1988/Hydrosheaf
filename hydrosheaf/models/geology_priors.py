"""Soft geology-to-reaction priors with explicit neutrality for missing joins.

Geology can supply a contextual prior over candidate reaction families, but a
mapped surface unit is not a direct reaction observation.  This module keeps
that prior separate from the reaction matrix and makes its information state
explicit.  A caller supplies the unit-to-family multipliers; the module does
not invent a lithology-to-reaction mapping.

For a matched unit and an explicitly supplied positive family multiplier
``m``, the returned penalty scale is ``1 / m`` (after optional strength and
explicit confidence scaling).  The scale is intended for an L1 reaction penalty: values below one
favour a family softly, values above one discourage it softly.  Unknown,
unmatched, invalid, or unmapped contexts always return a multiplier and
penalty scale of one.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from .geology_context import GeologyContext, normalise_geology_context


GeologyPriorStatus = Literal[
    "matched_prior",
    "neutral_unknown",
    "neutral_unmatched",
    "neutral_invalid",
    "neutral_unmapped",
]

_MIN_MULTIPLIER = 1.0e-6
_MAX_MULTIPLIER = 1.0e6


def load_geology_family_multipliers(path: str | Path) -> dict[str, dict[str, float]]:
    """Load and validate a user-supplied geology dictionary from JSON.

    The accepted JSON forms are either a direct unit mapping or a mapping
    under ``family_multipliers``/``units``/``entries``.  A unit profile may
    likewise be direct or nested under ``multipliers``.  No geological
    association is invented when a unit or family is absent: the prior layer
    remains neutral for that entry.
    """

    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"geology dictionary is not valid JSON: {source}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("geology dictionary must be a JSON object")
    profiles: Any = payload
    for key in ("family_multipliers", "units", "entries"):
        candidate = payload.get(key)
        if isinstance(candidate, Mapping):
            profiles = candidate
            break
    if not isinstance(profiles, Mapping):
        raise ValueError("geology dictionary profiles must be a JSON object")

    result: dict[str, dict[str, float]] = {}
    for raw_unit, raw_profile in profiles.items():
        unit = str(raw_unit).strip()
        if not unit:
            raise ValueError("geology dictionary contains an empty unit key")
        if not isinstance(raw_profile, Mapping):
            raise ValueError(f"geology profile for {unit!r} must be an object")
        profile_payload: Any = raw_profile.get("multipliers", raw_profile)
        if not isinstance(profile_payload, Mapping):
            raise ValueError(f"geology profile for {unit!r} must map families to multipliers")
        profile: dict[str, float] = {}
        for raw_family, raw_multiplier in profile_payload.items():
            family = str(raw_family).strip()
            if not family:
                raise ValueError(f"geology profile {unit!r} contains an empty family")
            multiplier = _positive_multiplier(
                f"geology multiplier {unit!r}/{family!r}", raw_multiplier
            )
            profile[family] = multiplier
        result[unit] = profile
    return result


def reaction_family(reaction: Any) -> str:
    """Return a conservative family token for a declared reaction label.

    The classifier only groups well-known labels used by HydroSheaf's
    reaction dictionary.  Unknown labels are returned as ``"other"`` rather
    than assigned a geological interpretation.
    """

    label = str(reaction or "").strip().casefold().replace("-", "_").replace(" ", "_")
    if not label:
        return "other"
    if any(token in label for token in ("calcite", "dolomite", "magnesite", "aragonite", "carbonate")):
        return "carbonate"
    if any(
        token in label
        for token in (
            "albite",
            "anorthite",
            "feldspar",
            "biotite",
            "chlorite",
            "pyroxene",
            "enstatite",
            "diopside",
            "olivine",
            "silicate",
        )
    ):
        return "silicate"
    if any(token in label for token in ("gypsum", "anhydrite", "halite", "sylvite", "fluorite", "evaporite")):
        return "evaporite"
    if any(token in label for token in ("pyrite", "sulfate_reduction", "sulphate_reduction", "iron_reduction", "denit", "nitrate")):
        return "redox"
    if "exch" in label or "exchange" in label:
        return "exchange"
    if label in {"mix", "mixing", "conservative", "no_reaction", "none"}:
        return "conservative"
    return "other"


def _finite(name: str, value: Any) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_multiplier(name: str, value: Any) -> float:
    result = _finite(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _unit_profile(
    context: GeologyContext,
    family_multipliers: Mapping[str, Mapping[str, Any]],
) -> Mapping[str, Any] | None:
    if not context.has_informative_unit or context.unit is None:
        return None
    for unit, profile in family_multipliers.items():
        if str(unit).strip().casefold() == context.unit.casefold():
            if isinstance(profile, Mapping):
                return profile
            return None
    return None


@dataclass(frozen=True)
class GeologyReactionPrior:
    """Soft contextual prior for one reaction label."""

    reaction: str
    family: str
    geology_unit: str | None
    context_status: str
    multiplier: float
    penalty_scale: float
    status: GeologyPriorStatus
    source: str | None
    reason: str

    @property
    def is_neutral(self) -> bool:
        return self.status != "matched_prior"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def geology_reaction_prior(
    reaction: Any,
    context: Any,
    *,
    family_multipliers: Mapping[str, Mapping[str, Any]] | None = None,
    strength: float = 1.0,
    preferred_key: str | None = None,
) -> GeologyReactionPrior:
    """Construct a soft reaction prior from explicitly mapped geology.

    ``family_multipliers`` maps a matched geology unit to family multipliers,
    for example ``{"granite": {"silicate": 2.0, "carbonate": 0.5}}``.
    Every multiplier is interpreted as relative support, not a calibrated
    probability.  Omitting the mapping leaves even a matched context neutral,
    which prevents an implicit geological lookup from masquerading as data.
    """

    label = str(reaction or "").strip()
    if not label:
        raise ValueError("reaction label must not be empty")
    prior_strength = _finite("strength", strength)
    if prior_strength < 0.0:
        raise ValueError("strength must be non-negative")
    normalised = normalise_geology_context(context, preferred_key=preferred_key)
    family = reaction_family(label)
    unit = normalised.unit
    source = normalised.source

    if normalised.status == "unknown":
        return GeologyReactionPrior(
            reaction=label,
            family=family,
            geology_unit=None,
            context_status=normalised.status,
            multiplier=1.0,
            penalty_scale=1.0,
            status="neutral_unknown",
            source=source,
            reason="geology assignment is unavailable; prior is neutral",
        )
    if normalised.status == "unmatched":
        return GeologyReactionPrior(
            reaction=label,
            family=family,
            geology_unit=None,
            context_status=normalised.status,
            multiplier=1.0,
            penalty_scale=1.0,
            status="neutral_unmatched",
            source=source,
            reason="geology join is unmatched; prior is neutral",
        )
    if normalised.status == "invalid":
        return GeologyReactionPrior(
            reaction=label,
            family=family,
            geology_unit=unit,
            context_status=normalised.status,
            multiplier=1.0,
            penalty_scale=1.0,
            status="neutral_invalid",
            source=source,
            reason="geology metadata is invalid; prior is neutral",
        )

    profiles = family_multipliers or {}
    profile = _unit_profile(normalised, profiles)
    if profile is None:
        return GeologyReactionPrior(
            reaction=label,
            family=family,
            geology_unit=unit,
            context_status=normalised.status,
            multiplier=1.0,
            penalty_scale=1.0,
            status="neutral_unmapped",
            source=source,
            reason="matched geology unit has no declared family multipliers",
        )

    raw_multiplier: Any = None
    wanted_keys = {family.casefold(), label.casefold()}
    for profile_key, profile_value in profile.items():
        if str(profile_key).strip().casefold() in wanted_keys:
            raw_multiplier = profile_value
            break
    if raw_multiplier is None:
        return GeologyReactionPrior(
            reaction=label,
            family=family,
            geology_unit=unit,
            context_status=normalised.status,
            multiplier=1.0,
            penalty_scale=1.0,
            status="neutral_unmapped",
            source=source,
            reason=f"matched geology unit has no declared multiplier for family {family!r}",
        )

    base_multiplier = _positive_multiplier(
        f"family multiplier for {label!r}", raw_multiplier
    )
    # Exponentiation keeps strength=0 exactly neutral and treats reciprocal
    # multipliers symmetrically in log-prior space.  Clipping protects a
    # caller-provided but extreme contextual weight from destabilising the L1
    # solver while leaving the choice explicit in the returned record.
    confidence = normalised.confidence
    effective_strength = prior_strength
    if confidence is not None:
        effective_strength *= confidence
    log_multiplier = max(
        math.log(_MIN_MULTIPLIER),
        min(math.log(_MAX_MULTIPLIER), effective_strength * math.log(base_multiplier)),
    )
    multiplier = float(math.exp(log_multiplier))
    penalty_scale = float(1.0 / multiplier)
    return GeologyReactionPrior(
        reaction=label,
        family=family,
        geology_unit=unit,
        context_status=normalised.status,
        multiplier=multiplier,
        penalty_scale=penalty_scale,
        status="matched_prior",
        source=source,
        reason="explicitly mapped geology supplied a soft family multiplier",
    )


def build_geology_reaction_priors(
    reaction_labels: Sequence[Any],
    context: Any,
    *,
    family_multipliers: Mapping[str, Mapping[str, Any]] | None = None,
    strength: float = 1.0,
    preferred_key: str | None = None,
) -> dict[str, GeologyReactionPrior]:
    """Return deterministic priors keyed by reaction label."""

    result: dict[str, GeologyReactionPrior] = {}
    for raw_label in reaction_labels:
        label = str(raw_label).strip()
        if label in result:
            raise ValueError(f"reaction labels must be unique: {label!r}")
        result[label] = geology_reaction_prior(
            label,
            context,
            family_multipliers=family_multipliers,
            strength=strength,
            preferred_key=preferred_key,
        )
    return result


def apply_geology_penalty_scales(
    reaction_labels: Sequence[Any],
    base_penalty_scales: Sequence[Any],
    context: Any,
    *,
    family_multipliers: Mapping[str, Mapping[str, Any]] | None = None,
    strength: float = 1.0,
    preferred_key: str | None = None,
) -> tuple[list[float], dict[str, GeologyReactionPrior]]:
    """Apply soft geology multipliers while preserving a supplied base scale."""

    if len(reaction_labels) != len(base_penalty_scales):
        raise ValueError("reaction_labels and base_penalty_scales must have equal length")
    base = [_positive_multiplier("base penalty scale", value) for value in base_penalty_scales]
    priors = build_geology_reaction_priors(
        reaction_labels,
        context,
        family_multipliers=family_multipliers,
        strength=strength,
        preferred_key=preferred_key,
    )
    scales = [
        float(base[index] * priors[str(raw_label).strip()].penalty_scale)
        for index, raw_label in enumerate(reaction_labels)
    ]
    return scales, priors


__all__ = [
    "GeologyPriorStatus",
    "GeologyReactionPrior",
    "apply_geology_penalty_scales",
    "build_geology_reaction_priors",
    "geology_reaction_prior",
    "load_geology_family_multipliers",
    "reaction_family",
]
