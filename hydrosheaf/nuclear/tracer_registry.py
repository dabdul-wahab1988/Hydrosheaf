"""Pre-registered environmental tracer definitions and observational metadata.

This module provides the TracerSpec typed data contract and registry for
multi-tracer graph inversion across stable isotopes, tritium, SF6, 14C,
and geochemical constraints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Optional, Sequence


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
    enabled : bool
        Whether this tracer is included in the active inversion panel.
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

    def __post_init__(self) -> None:
        valid_kinds = {"stable_isotope", "radioactive", "dissolved_gas", "chemistry"}
        if self.kind not in valid_kinds:
            raise ValueError(f"kind must be one of {valid_kinds}, got {self.kind!r}")
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
            metadata={"isotope_system": "water_molecule"},
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
            metadata={"isotope_system": "water_molecule", "lmwl_slope": 8.0, "lmwl_intercept": 10.0},
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
            metadata={"half_life_years": TRITIUM_HALF_LIFE_YEARS},
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
            metadata={"excess_air_default_cc_g": 3.0e-3, "recharge_temp_c": 25.0},
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
            metadata={
                "half_life_years": CARBON14_HALF_LIFE_YEARS,
                "default_q_correction": 0.85,
                "q_uncertainty": 0.10,
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
            metadata={"conservative_ions": ["Cl", "Br"]},
        ),
    }
