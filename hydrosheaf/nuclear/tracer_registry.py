"""Pre-registered environmental tracer definitions and observational metadata.

This module provides the TracerSpec typed data contract and registry for
multi-tracer graph inversion across stable isotopes, tritium, SF6, 14C,
and geochemical constraints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
from typing import Any, Mapping, Optional


def _alias_lookup_key(value: str) -> str:
    """Return a separator-insensitive lookup key for a tracer alias."""
    text = str(value).strip()
    text = text.replace("δ", "delta").replace("Δ", "delta")
    text = text.translate(str.maketrans({"¹": "1", "²": "2", "³": "3", "⁴": "4"}))
    text = text.upper()
    return re.sub(r"[\s/_-]+", "", text)


# Keep the canonical IDs stable because they are used as serialized keys by
# history loaders, temporal samples, and the nuclear inference APIs.  The
# groups include the aliases that were already accepted by tracer_inputs plus
# the stable-water-isotope spellings used in field data exports.
TRACER_ALIAS_GROUPS: dict[str, tuple[str, ...]] = {
    "d18O": ("d18O", "18O", "D18O", "delta18O", "Δ18O", "δ18O"),
    "d2H": ("d2H", "2H", "D2H", "delta2H", "Δ2H", "δ2H"),
    "3H": ("3H", "H3", "H-3", "tritium", "Tritium"),
    "SF6": ("SF6", "SF_6", "sulfur hexafluoride", "Sulfur hexafluoride"),
    "14C": ("14C", "C14", "carbon14", "Carbon-14", "radiocarbon", "Radiocarbon"),
    "39Ar": ("39Ar", "Ar39", "argon39", "Argon-39"),
    "85Kr": ("85Kr", "Kr85", "krypton85", "Krypton-85"),
    "4He": ("4He", "He4", "helium4", "Helium-4"),
    "3H/3He": ("3H/3He", "3H3He", "H3/He3", "tritium-helium3"),
    "CFC11": ("CFC11", "CFC_11", "CFC-11"),
    "CFC12": ("CFC12", "CFC_12", "CFC-12"),
    "CFC113": ("CFC113", "CFC_113", "CFC-113"),
}


# Retain a flat public alias map for callers that imported the old constant
# from tracer_inputs.  Lookup itself is performed by canonicalize_tracer_alias
# so that Unicode delta symbols and separator variants share one implementation.
TRACER_ALIASES: dict[str, str] = {}
for _canonical, _aliases in TRACER_ALIAS_GROUPS.items():
    for _alias in _aliases:
        TRACER_ALIASES[_alias] = _canonical
        TRACER_ALIASES[_alias_lookup_key(_alias)] = _canonical

_TRACER_ALIAS_INDEX = {
    _alias_lookup_key(alias): canonical
    for canonical, aliases in TRACER_ALIAS_GROUPS.items()
    for alias in aliases
}


def canonicalize_tracer_alias(tracer: str) -> str:
    """Return HydroSheaf's canonical ID for a known tracer alias.

    Unknown tracers are returned in their stripped, original spelling.  That
    fallback preserves the historical behaviour for chemistry keys such as
    ``Cl`` and ``Na`` while making stable-isotope and nuclear aliases resolve
    consistently across the temporal and nuclear loaders.
    """
    text = str(tracer).strip()
    if not text:
        return text
    return _TRACER_ALIAS_INDEX.get(_alias_lookup_key(text), text)


@dataclass(frozen=True)
class TracerSpec:
    """Specification of an environmental tracer for graph transport and inversion.

    Parameters
    ----------
    tracer_id : str
        Unique standard identifier (e.g. "d18O", "d2H", "3H", "SF6", "14C", "Cl").
    display_name : str
        Human-readable scientific label (e.g. "δ¹⁸O", "Tritium (³H)").
    kind : str
        One of 'stable_isotope', 'radioactive', 'dissolved_gas', 'chemistry'.
    units : str
        Measurement units (e.g. "permil", "TU", "pptv", "pmC", "mg/L").
    decay_constant_per_day : Optional[float]
        Radioactive decay rate lambda [1/day]. None for stable tracers.
    source_history_required : bool
        Whether a time-varying atmospheric/recharge input curve is mandatory.
    response_model : str
        Observation response type: 'conservative', 'fractionated', 'gas_exchange',
        'carbonate_corrected', 'unbalanced_ot'.
    measurement_sd : float
        Default analytical 1-sigma uncertainty.
    detection_limit : Optional[float]
        Lower analytical detection limit. Observations below this are left-censored.
    fractionation_model : Optional[str]
        Fractionation model name for isotopes (e.g. 'craig_gordon', 'rayleigh').
    exchange_model : Optional[str]
        Dissolved gas exchange model (e.g. 'closed_system', 'excess_air').
    response_family : str
        Forward-response family used by the inference layer.  This is
        intentionally distinct from ``response_model``: stable water
        isotopes use a recharge-history response, whereas radioactive and
        dissolved-gas tracers can use scalar age-grid kernels when supported.
    supports_scalar_age_grid_kernel : bool
        Whether the scalar age-grid ``tracer_response_kernel`` can represent
        this tracer without a time-indexed stable-isotope source model.
    enabled : bool
        Whether this tracer is included in the active inversion panel.
    aliases : tuple[str, ...]
        Accepted input spellings for this canonical tracer ID.
    metadata : Mapping[str, Any]
        Additional physical parameters (e.g. half-life, Henry's constants).
    """

    tracer_id: str
    display_name: str
    kind: str
    units: str
    decay_constant_per_day: Optional[float] = None
    source_history_required: bool = False
    response_model: str = "conservative"
    measurement_sd: float = 0.05
    detection_limit: Optional[float] = None
    fractionation_model: Optional[str] = None
    exchange_model: Optional[str] = None
    enabled: bool = True
    metadata: Mapping[str, Any] = field(default_factory=dict)
    response_family: str = "generic"
    supports_scalar_age_grid_kernel: bool = False
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        valid_kinds = {"stable_isotope", "radioactive", "dissolved_gas", "chemistry"}
        valid_response_families = {
            "generic",
            "stable_isotope_time_history",
            "radioactive_decay",
            "gas_input_history",
            "chemistry",
        }
        if self.kind not in valid_kinds:
            raise ValueError(f"kind must be one of {valid_kinds}, got {self.kind!r}")
        if self.response_family not in valid_response_families:
            raise ValueError(
                "response_family must be one of "
                f"{valid_response_families}, got {self.response_family!r}"
            )
        if self.decay_constant_per_day is not None and self.decay_constant_per_day < 0.0:
            raise ValueError("decay_constant_per_day must be non-negative.")
        if self.measurement_sd <= 0.0:
            raise ValueError("measurement_sd must be strictly positive.")


# Half-life constants in solar years
TRITIUM_HALF_LIFE_YEARS = 12.32
CARBON14_HALF_LIFE_YEARS = 5730.0
DAYS_PER_YEAR = 365.25

LAMBDA_3H_PER_DAY = math.log(2.0) / (TRITIUM_HALF_LIFE_YEARS * DAYS_PER_YEAR)
LAMBDA_14C_PER_DAY = math.log(2.0) / (CARBON14_HALF_LIFE_YEARS * DAYS_PER_YEAR)


def build_default_tracer_registry() -> dict[str, TracerSpec]:
    """Return dictionary of pre-registered standard groundwater tracers."""
    return {
        "d18O": TracerSpec(
            tracer_id="d18O",
            display_name="δ¹⁸O",
            kind="stable_isotope",
            units="permil VSMOW",
            source_history_required=True,
            response_model="conservative",
            measurement_sd=0.10,
            detection_limit=None,
            response_family="stable_isotope_time_history",
            supports_scalar_age_grid_kernel=False,
            aliases=TRACER_ALIAS_GROUPS["d18O"],
            metadata={
                "isotope_system": "water_molecule",
                "canonical_id": "d18O",
                "aliases": TRACER_ALIAS_GROUPS["d18O"],
            },
        ),
        "d2H": TracerSpec(
            tracer_id="d2H",
            display_name="δ²H",
            kind="stable_isotope",
            units="permil VSMOW",
            source_history_required=True,
            response_model="conservative",
            measurement_sd=0.80,
            detection_limit=None,
            response_family="stable_isotope_time_history",
            supports_scalar_age_grid_kernel=False,
            aliases=TRACER_ALIAS_GROUPS["d2H"],
            metadata={
                "isotope_system": "water_molecule",
                "lmwl_slope": 8.0,
                "lmwl_intercept": 10.0,
                "canonical_id": "d2H",
                "aliases": TRACER_ALIAS_GROUPS["d2H"],
            },
        ),
        "3H": TracerSpec(
            tracer_id="3H",
            display_name="Tritium (³H)",
            kind="radioactive",
            units="TU",
            decay_constant_per_day=LAMBDA_3H_PER_DAY,
            source_history_required=True,
            response_model="radioactive_decay",
            measurement_sd=0.15,
            detection_limit=0.05,
            response_family="radioactive_decay",
            supports_scalar_age_grid_kernel=True,
            aliases=TRACER_ALIAS_GROUPS["3H"],
            metadata={
                "half_life_years": TRITIUM_HALF_LIFE_YEARS,
                "canonical_id": "3H",
                "aliases": TRACER_ALIAS_GROUPS["3H"],
            },
        ),
        "SF6": TracerSpec(
            tracer_id="SF6",
            display_name="Sulfur Hexafluoride (SF₆)",
            kind="dissolved_gas",
            units="pptv",
            decay_constant_per_day=None,
            source_history_required=True,
            response_model="gas_exchange",
            measurement_sd=0.20,
            detection_limit=0.01,
            exchange_model="excess_air_unf",
            response_family="gas_input_history",
            supports_scalar_age_grid_kernel=True,
            aliases=TRACER_ALIAS_GROUPS["SF6"],
            metadata={
                "excess_air_default_cc_g": 3.0e-3,
                "recharge_temp_c": 25.0,
                "canonical_id": "SF6",
                "aliases": TRACER_ALIAS_GROUPS["SF6"],
            },
        ),
        "14C": TracerSpec(
            tracer_id="14C",
            display_name="Radiocarbon (¹⁴C)",
            kind="radioactive",
            units="pmC",
            decay_constant_per_day=LAMBDA_14C_PER_DAY,
            source_history_required=True,
            response_model="carbonate_corrected",
            measurement_sd=1.0,
            detection_limit=0.5,
            response_family="radioactive_decay",
            supports_scalar_age_grid_kernel=True,
            aliases=TRACER_ALIAS_GROUPS["14C"],
            metadata={
                "half_life_years": CARBON14_HALF_LIFE_YEARS,
                "default_q_correction": 0.85,
                "q_uncertainty": 0.10,
                "canonical_id": "14C",
                "aliases": TRACER_ALIAS_GROUPS["14C"],
            },
        ),
        "chemistry": TracerSpec(
            tracer_id="chemistry",
            display_name="Reaction-Aware Chemistry",
            kind="chemistry",
            units="dimensionless_penalty",
            source_history_required=False,
            response_model="unbalanced_ot",
            measurement_sd=1.0,
            enabled=False,  # Auxiliary by default
            response_family="chemistry",
            supports_scalar_age_grid_kernel=False,
            metadata={"conservative_ions": ["Cl", "Br"]},
        ),
    }


__all__ = [
    "TracerSpec",
    "TRACER_ALIASES",
    "TRACER_ALIAS_GROUPS",
    "canonicalize_tracer_alias",
    "build_default_tracer_registry",
]
