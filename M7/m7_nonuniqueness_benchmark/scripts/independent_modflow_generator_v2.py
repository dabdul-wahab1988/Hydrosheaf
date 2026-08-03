"""Environmental-isotope and nuclear-tracer extension of the M7 generator.

This module is an additive layer on top of
``independent_modflow_generator``.  The v1 module is byte-locked by
``m7_system_acceptance_protocol.lock.json`` and by the M8 frontier lock, so it
is imported and reused here rather than edited.

Like v1, this module intentionally imports no ``hydrosheaf`` code.  MODFLOW 6
supplies heads, MODPATH 7 supplies routes and travel times, the v1 nonlinear
process simulator supplies major-ion chemistry and the 3H/39Ar pair, and the
formulations below supply everything else.  Every atmospheric input history,
solubility relationship, mixing model and decay law here is a closed-form
analytic surrogate written from published endmember values; no HydroSheaf
input-history table, correction routine or solver is consulted.

Three things are added:

1. Environmental isotopes.  ``d18O_permil`` and ``d2H_permil`` are set at the
   recharge location and then carried conservatively, so they are informative
   about recharge source and mixing and carry no age information at all.
   ``d13C_permil`` and ``sr87_sr86`` are *not* conservative: they evolve along
   the realised reaction path.
2. A nuclear-tracer panel: 14C, CFC-11, CFC-12, CFC-113, SF6, terrigenic 4He
   and tritiogenic 3He.
3. An explicit per-case, per-family shared-nuisance mechanism whose true value
   is returned as ground truth so a benchmark can score its recoverability.

Reproducibility guarantee
-------------------------
All v2 randomness is drawn from auxiliary generators seeded independently of
the v1 ``rng`` stream (``seed + 900_000`` and ``seed + 910_000``).  The v1
``rng`` is consumed in exactly the same order, with exactly the same calls, as
in :func:`independent_modflow_generator.generate_independent_aquifer`.  Every
value produced by v1 is therefore bit-identical under v2 for every seed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import flopy
import numpy as np

from independent_modflow_generator import (
    ION_ORDER,
    IndependentAquifer,
    _build_flow_model,
    _chemistry_step,
    _run_pathlines,
    _sha256,
    _tracer_observations,
    _version,
)

__all__ = [
    "ION_ORDER",
    "IndependentAquifer",
    "IndependentAquiferV2",
    "generate_independent_aquifer_v2",
]

# ---------------------------------------------------------------------------
# Random-stream offsets.  Large and distinct so the isotope, tracer and
# nuisance streams cannot alias the locked v1 stream.
# ---------------------------------------------------------------------------
ISOTOPE_STREAM_SEED_OFFSET = 900_000
NUISANCE_STREAM_SEED_OFFSET = 910_000

# --- Stable water isotopes -------------------------------------------------
# Craig (1961) global meteoric water line, d2H = 8 * d18O + 10.
GMWL_SLOPE = 8.0
GMWL_INTERCEPT_PERMIL = 10.0
# Weighted-mean tropical West African recharge composition.
D18O_RECHARGE_BASE_PERMIL = -3.20
# The spatial gradients are amplified relative to continental-scale field
# values because the synthetic domain is only 3.0 km x 1.5 km.  Without
# amplification the recharge-source signal would be far smaller than analytical
# precision and the benchmark could not test whether d18O/d2H carry usable
# source information.
D18O_CONTINENTAL_GRADIENT_PERMIL_PER_KM = -0.35
D18O_LATITUDE_GRADIENT_PERMIL_PER_KM = 0.18
# Altitude effect, -0.20 permil per 100 m of recharge elevation.
D18O_ALTITUDE_GRADIENT_PERMIL_PER_M = -0.0020
D18O_ANALYTICAL_SIGMA_PERMIL = 0.08
DEUTERIUM_EXCESS_SIGMA_PERMIL = 1.5
# Declared fraction of particles whose recharge was evaporatively enriched.
EVAPORATED_PARTICLE_FRACTION = 0.18
# Local evaporation lines in semi-arid settings have slopes near 4-6.
EVAPORATION_LINE_SLOPE = 5.0
EVAPORATION_ENRICHMENT_RANGE_PERMIL = (0.5, 3.0)

# --- Carbon and strontium isotopes -----------------------------------------
# C3/C4 soil-zone CO2 endmember and marine carbonate endmember.
SOIL_CO2_D13C_PERMIL = -23.0
MARINE_CARBONATE_D13C_PERMIL = 0.0
# Fractional approach of DIC toward the carbonate endmember per weathering
# milestone.  Closed-system congruent calcite dissolution converges near 50 %
# carbonate-derived carbon; open-system exchange pushes it higher.
CARBONATE_CARBON_GAIN_PER_STEP = 0.42
# Calcite is roughly +2 permil relative to HCO3-, so precipitation leaves the
# residual DIC slightly lighter.
D13C_SHIFT_CARBONATE_PRECIPITATION_PERMIL = -0.70
# Organic-matter oxidation adds isotopically light, radiocarbon-dead carbon.
D13C_SHIFT_SULFATE_REDUCTION_PERMIL = -1.60
D13C_SHIFT_IRON_REDUCTION_PERMIL = -0.90
D13C_ANALYTICAL_SIGMA_PERMIL = 0.15
# Two-endmember 87Sr/86Sr mixing: marine carbonate versus crystalline basement.
SR87_SR86_CARBONATE = 0.7080
SR87_SR86_SILICATE = 0.7200
SR87_SR86_SILICATE_GAIN_PER_STEP = 0.30
SR87_SR86_ANALYTICAL_SIGMA = 4.0e-5

# --- Radiocarbon -----------------------------------------------------------
CARBON14_HALF_LIFE_YEARS = 5730.0
CARBON14_INITIAL_PMC = 100.0
# Radiocarbon-dead carbon added by each redox milestone, on top of the d13C
# mixing dilution.
C14_REDOX_DEAD_CARBON_PER_STEP = 0.06
C14_DILUTION_BOUNDS = (0.15, 1.0)
C14_ANALYTICAL_SIGMA_PMC = 0.6

# --- Atmospheric halocarbon and SF6 surrogates -----------------------------
# (plateau_pptv, rise_midpoint_year, rise_width_years, peak_year, decline_rate)
# A smooth logistic rise with an exponential post-peak decline.  The parameters
# reproduce the published Northern-Hemisphere shape -- CFC-11 near 272 pptv
# peaking in 1993-1994, CFC-12 near 545 pptv peaking near 2002, CFC-113 near
# 84 pptv peaking near 1996, and SF6 still rising through 2025 -- without
# importing or reading any tabulated history.
ATMOSPHERIC_TRACER_HISTORY: Dict[str, Tuple[float, float, float, float, float]] = {
    "cfc11_pptv": (272.0, 1972.0, 6.5, 1993.5, 0.0080),
    "cfc12_pptv": (545.0, 1974.0, 7.5, 2002.0, 0.0032),
    "cfc113_pptv": (84.0, 1982.0, 5.0, 1996.0, 0.0068),
    "sf6_pptv": (11.5, 2003.0, 12.0, 1.0e6, 0.0),
}
HALOCARBON_GASES = ("cfc11_pptv", "cfc12_pptv", "cfc113_pptv", "sf6_pptv")
# Recharge temperature the reporting laboratory assumes when it converts the
# measured dissolved concentration back into an equivalent mixing ratio.
REFERENCE_RECHARGE_TEMPERATURE_C = 27.0
# -d(ln K_H)/dT near 25 C.  A 1 C recharge-temperature error therefore biases
# the reported mixing ratio by roughly 3-4.5 % depending on the gas, which is
# the classic common-mode CFC/SF6 dating error.
GAS_SOLUBILITY_TEMPERATURE_SENSITIVITY = {
    "cfc11_pptv": 0.042,
    "cfc12_pptv": 0.036,
    "cfc113_pptv": 0.045,
    "sf6_pptv": 0.028,
}
# Relative susceptibility to reductive dehalogenation.  CFC-11 is removed
# rapidly under sulfate-reducing and methanogenic conditions, CFC-113 more
# slowly, CFC-12 is nearly recalcitrant, and SF6 is unaffected.
CFC_DEGRADATION_SUSCEPTIBILITY = {
    "cfc11_pptv": 1.00,
    "cfc113_pptv": 0.35,
    "cfc12_pptv": 0.03,
    "sf6_pptv": 0.00,
}
CFC_DEGRADATION_STRENGTH = 0.85
REDUCING_INTENSITY_WEIGHTS = {
    "sulfate_reduction": 1.00,
    "iron_reduction": 0.55,
}
HALOCARBON_ANALYTICAL_CV = 0.03

# --- Helium ----------------------------------------------------------------
# Air-saturated-water 4He near 27 C plus a combined in-situ and external
# crustal accumulation rate typical of a thin sedimentary/basement aquifer.
HE4_AIR_SATURATED_CCPG = 4.5e-8
HE4_ACCUMULATION_CCPG_PER_YEAR = 1.0e-10
HE4_ANALYTICAL_CV = 0.05
# Used only to invert the *reported* tritium into its tritiogenic 3He partner.
# The tritium value itself is produced by v1 and is never recomputed.
TRITIUM_HALF_LIFE_YEARS = 12.32
TRITIOGENIC_HE3_RETENTION = 0.90
# The closed-system 3H/3He clock loses meaning once essentially all tritium has
# decayed, so the ingrowth exponent saturates beyond this age.
TRITIOGENIC_MAX_AGE_YEARS = 70.0

# --- Shared nuisance magnitudes --------------------------------------------
NUISANCE_RECHARGE_TEMPERATURE_SIGMA_C = 2.0
NUISANCE_C14_A0_SIGMA_PERMIL = 2.5

REDOX_PROCESSES = ("sulfate_reduction", "iron_reduction")


@dataclass(frozen=True)
class IndependentAquiferV2:
    """v1 case content plus the v2 isotope/tracer ground truth."""

    seed: int
    observations: Tuple[Dict[str, object], ...]
    true_edges: Tuple[Tuple[str, str], ...]
    true_ages_years: Dict[str, float]
    true_processes: Dict[str, str]
    pathline_rows: Tuple[Dict[str, object], ...]
    provenance: Dict[str, object]
    true_nuisance: Mapping[str, float] = field(default_factory=dict)
    true_recharge_isotopes: Mapping[str, Mapping[str, float]] = field(
        default_factory=dict
    )
    true_reaction_isotopes: Mapping[str, Mapping[str, float]] = field(
        default_factory=dict
    )

    def as_v1(self) -> IndependentAquifer:
        """Return the plain v1 record for consumers that expect it."""

        return IndependentAquifer(
            seed=self.seed,
            observations=self.observations,
            true_edges=self.true_edges,
            true_ages_years=self.true_ages_years,
            true_processes=self.true_processes,
            pathline_rows=self.pathline_rows,
            provenance=self.provenance,
        )


# ---------------------------------------------------------------------------
# Shared nuisance
# ---------------------------------------------------------------------------
def draw_shared_nuisance(
    seed: int,
    *,
    shared_nuisance: bool = False,
    nuisance_scale: float = 0.0,
) -> Dict[str, float]:
    """Draw the per-case, per-family systematic errors.

    These values are latent: they bias a whole tracer family at once but are
    never written into ``observations``.  They are returned as ground truth so
    a benchmark can score whether an estimator recovers them.  With the default
    arguments the result is exactly zero, so the default v2 case is unperturbed.
    """

    if not shared_nuisance or float(nuisance_scale) == 0.0:
        return {
            "recharge_temperature_error_c": 0.0,
            "c14_a0_error_permil": 0.0,
        }
    nuisance_rng = np.random.default_rng(int(seed) + NUISANCE_STREAM_SEED_OFFSET)
    scale = float(nuisance_scale)
    return {
        # Common-mode across the whole CFC/SF6 family through the solubility
        # conversion.
        "recharge_temperature_error_c": float(
            scale * NUISANCE_RECHARGE_TEMPERATURE_SIGMA_C * nuisance_rng.normal()
        ),
        # A bias on the assumed soil-CO2 d13C endmember, which propagates into
        # the 14C initial activity and therefore couples 14C to the chemistry
        # stream instead of acting as independent noise.
        "c14_a0_error_permil": float(
            scale * NUISANCE_C14_A0_SIGMA_PERMIL * nuisance_rng.normal()
        ),
    }


# ---------------------------------------------------------------------------
# Environmental isotopes
# ---------------------------------------------------------------------------
def draw_recharge_evaporation(
    iso_rng: Optional[np.random.Generator],
    *,
    evaporated_fraction: float = EVAPORATED_PARTICLE_FRACTION,
) -> Dict[str, object]:
    """Decide whether one particle's recharge was evaporatively enriched."""

    if iso_rng is None:
        return {"evaporated": False, "evaporation_extent_permil": 0.0}
    evaporated = bool(iso_rng.random() < float(evaporated_fraction))
    low, high = EVAPORATION_ENRICHMENT_RANGE_PERMIL
    extent = float(iso_rng.uniform(low, high)) if evaporated else 0.0
    return {"evaporated": evaporated, "evaporation_extent_permil": extent}


def environmental_isotopes(
    x_m: float,
    y_m: float,
    elevation_m: float,
    *,
    evaporated: bool = False,
    evaporation_extent_permil: float = 0.0,
    iso_rng: Optional[np.random.Generator] = None,
) -> Dict[str, object]:
    """Conservative stable water isotopes, fixed at the recharge location.

    ``d18O`` and ``d2H`` depend only on where the water entered the system, not
    on how long it has been underground.  Passing ``iso_rng=None`` gives the
    noise-free deterministic value, which is what the unit tests exercise.
    """

    d18o = (
        D18O_RECHARGE_BASE_PERMIL
        + D18O_CONTINENTAL_GRADIENT_PERMIL_PER_KM * (float(x_m) / 1000.0)
        + D18O_LATITUDE_GRADIENT_PERMIL_PER_KM * (float(y_m) / 1000.0)
        + D18O_ALTITUDE_GRADIENT_PERMIL_PER_M * float(elevation_m)
    )
    deuterium_excess = 0.0
    if iso_rng is not None:
        d18o += float(iso_rng.normal(0.0, D18O_ANALYTICAL_SIGMA_PERMIL))
        deuterium_excess = float(iso_rng.normal(0.0, DEUTERIUM_EXCESS_SIGMA_PERMIL))
    d2h = GMWL_SLOPE * d18o + GMWL_INTERCEPT_PERMIL + deuterium_excess
    if evaporated and float(evaporation_extent_permil) > 0.0:
        # Evaporated water leaves the meteoric line along a local evaporation
        # line of slope ~5, which lowers deuterium excess.
        extent = float(evaporation_extent_permil)
        d2h = d2h + EVAPORATION_LINE_SLOPE * extent
        d18o = d18o + extent
    return {
        "d18O_permil": float(d18o),
        "d2H_permil": float(d2h),
        "deuterium_excess_permil": float(d2h - GMWL_SLOPE * d18o),
        "evaporated": bool(evaporated),
    }


def reaction_path_isotopes(
    process_history: Sequence[str],
    *,
    iso_rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """Reactive d13C and 87Sr/86Sr driven by the realised reaction path.

    Unlike the water isotopes these are not conservative: they record which
    reactions the parcel has passed through, which is what lets them share a
    latent nuisance with the radiocarbon stream.
    """

    carbonate_fraction = 0.0
    silicate_fraction = 0.0
    shift_permil = 0.0
    for process in process_history:
        if process == "carbonate_weathering":
            carbonate_fraction += CARBONATE_CARBON_GAIN_PER_STEP * (
                1.0 - carbonate_fraction
            )
        elif process == "silicate_weathering":
            silicate_fraction += SR87_SR86_SILICATE_GAIN_PER_STEP * (
                1.0 - silicate_fraction
            )
        elif process == "carbonate_precipitation":
            shift_permil += D13C_SHIFT_CARBONATE_PRECIPITATION_PERMIL
        elif process == "sulfate_reduction":
            shift_permil += D13C_SHIFT_SULFATE_REDUCTION_PERMIL
        elif process == "iron_reduction":
            shift_permil += D13C_SHIFT_IRON_REDUCTION_PERMIL

    d13c = (
        (1.0 - carbonate_fraction) * SOIL_CO2_D13C_PERMIL
        + carbonate_fraction * MARINE_CARBONATE_D13C_PERMIL
        + shift_permil
    )
    strontium = (
        1.0 - silicate_fraction
    ) * SR87_SR86_CARBONATE + silicate_fraction * SR87_SR86_SILICATE
    if iso_rng is not None:
        d13c += float(iso_rng.normal(0.0, D13C_ANALYTICAL_SIGMA_PERMIL))
        strontium += float(iso_rng.normal(0.0, SR87_SR86_ANALYTICAL_SIGMA))
    return {
        "d13C_permil": float(d13c),
        "sr87_sr86": float(strontium),
        "carbonate_carbon_fraction": float(carbonate_fraction),
        "silicate_strontium_fraction": float(silicate_fraction),
    }


# ---------------------------------------------------------------------------
# Nuclear tracer panel
# ---------------------------------------------------------------------------
def atmospheric_mixing_ratio(gas: str, recharge_year: float) -> float:
    """Analytic Northern-Hemisphere input-history surrogate, in pptv."""

    plateau, midpoint, width, peak_year, decline = ATMOSPHERIC_TRACER_HISTORY[gas]
    exponent = float(np.clip(-(float(recharge_year) - midpoint) / width, -600.0, 600.0))
    rise = plateau / (1.0 + np.exp(exponent))
    fall = np.exp(-decline * max(0.0, float(recharge_year) - peak_year))
    return float(max(0.0, rise * fall))


def reducing_intensity(process_history: Sequence[str]) -> float:
    """Dimensionless cumulative reducing exposure along the path."""

    return float(
        sum(REDUCING_INTENSITY_WEIGHTS.get(process, 0.0) for process in process_history)
    )


def halocarbon_panel(
    age_years: float,
    *,
    sample_date: float = 2025.5,
    process_history: Sequence[str] = (),
    recharge_temperature_error_c: float = 0.0,
    cfc_degradation_strength: float = CFC_DEGRADATION_STRENGTH,
    iso_rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """CFC-11/12/113 and SF6 as reported equivalent mixing ratios."""

    recharge_year = float(sample_date) - float(age_years)
    intensity = reducing_intensity(process_history)
    panel: Dict[str, float] = {}
    for gas in HALOCARBON_GASES:
        atmospheric = atmospheric_mixing_ratio(gas, recharge_year)
        # The reported mixing ratio is the dissolved concentration divided by
        # the solubility at the *assumed* recharge temperature, so a shared
        # error in that assumed temperature biases every gas in the family in
        # the same direction.
        sensitivity = GAS_SOLUBILITY_TEMPERATURE_SENSITIVITY[gas]
        solubility_bias = float(
            np.exp(-sensitivity * float(recharge_temperature_error_c))
        )
        survival = float(
            np.exp(
                -float(cfc_degradation_strength)
                * CFC_DEGRADATION_SUSCEPTIBILITY[gas]
                * intensity
            )
        )
        value = atmospheric * solubility_bias * survival
        if iso_rng is not None:
            value *= float(iso_rng.lognormal(0.0, HALOCARBON_ANALYTICAL_CV))
        panel[gas] = float(max(0.0, value))
    return panel


def radiocarbon_pmc(
    age_years: float,
    d13c_permil: float,
    *,
    process_history: Sequence[str] = (),
    c14_a0_error_permil: float = 0.0,
    iso_rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """Radiocarbon activity diluted through the shared d13C mixing relation.

    The dilution factor is the two-endmember (Ingerson-Pearson style) ratio
    computed from the sample's own d13C, so carbonate dissolution lowers 14C
    and shifts d13C together.  That shared latent nuisance is the scientific
    point of the panel.

    ``c14_a0_error_permil`` is a bias on the assumed soil-CO2 d13C endmember.
    It rescales the initial activity through the same endmember relationship,
    which is what couples the 14C stream to the chemistry stream.  Zero leaves
    the initial activity at exactly ``CARBON14_INITIAL_PMC``.
    """

    denominator = SOIL_CO2_D13C_PERMIL - MARINE_CARBONATE_D13C_PERMIL
    dilution = (float(d13c_permil) - MARINE_CARBONATE_D13C_PERMIL) / denominator
    low, high = C14_DILUTION_BOUNDS
    dilution = float(np.clip(dilution, low, high))
    dead_carbon_steps = sum(
        1 for process in process_history if process in REDOX_PROCESSES
    )
    dilution *= float(np.exp(-C14_REDOX_DEAD_CARBON_PER_STEP * dead_carbon_steps))
    biased_endmember = SOIL_CO2_D13C_PERMIL + float(c14_a0_error_permil)
    if abs(biased_endmember) < 1e-9:
        initial_activity = float(CARBON14_INITIAL_PMC)
    else:
        initial_activity = float(
            CARBON14_INITIAL_PMC * (SOIL_CO2_D13C_PERMIL / biased_endmember)
        )
    activity = (
        initial_activity
        * float(np.exp(-np.log(2.0) * float(age_years) / CARBON14_HALF_LIFE_YEARS))
        * dilution
    )
    if iso_rng is not None:
        activity += float(iso_rng.normal(0.0, C14_ANALYTICAL_SIGMA_PMC))
    return {
        "c14_pmc": float(max(0.05, activity)),
        "c14_dilution_factor": float(dilution),
        "c14_initial_activity_pmc": float(initial_activity),
    }


def helium_panel(
    age_years: float,
    tritium_tu: float,
    *,
    iso_rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """Terrigenic 4He ingrowth and the tritiogenic 3He partner of the tritium.

    ``tritium_tu`` is the value already produced by the v1 tracer equations.  It
    is inverted, never recomputed, so the 3H/3He pair stays self-consistent with
    the locked tritium stream.
    """

    helium4 = HE4_AIR_SATURATED_CCPG + HE4_ACCUMULATION_CCPG_PER_YEAR * max(
        0.0, float(age_years)
    )
    clock_age = float(np.clip(float(age_years), 0.0, TRITIOGENIC_MAX_AGE_YEARS))
    ingrowth = float(np.exp(np.log(2.0) * clock_age / TRITIUM_HALF_LIFE_YEARS) - 1.0)
    helium3 = TRITIOGENIC_HE3_RETENTION * float(tritium_tu) * ingrowth
    if iso_rng is not None:
        helium4 *= float(iso_rng.lognormal(0.0, HE4_ANALYTICAL_CV))
        helium3 *= float(iso_rng.lognormal(0.0, HE4_ANALYTICAL_CV))
    return {
        "he4_ccpg": float(max(0.0, helium4)),
        "h3_he3_TU": float(max(0.0, helium3)),
    }


def tracer_observations_v2(
    age_years: float,
    rng: np.random.Generator,
    *,
    iso_rng: Optional[np.random.Generator] = None,
    d13c_permil: float,
    sample_date: float = 2025.5,
    process_history: Sequence[str] = (),
    recharge_temperature_error_c: float = 0.0,
    c14_a0_error_permil: float = 0.0,
    cfc_degradation_strength: float = CFC_DEGRADATION_STRENGTH,
) -> Dict[str, float]:
    """v1 3H/39Ar plus the v2 nuclear panel.

    The v1 function is called first and unmodified, so ``tritium_TU`` and
    ``argon39_pmc`` are bit-identical and ``rng`` is advanced by exactly the two
    draws it made before.  Everything added below consumes ``iso_rng`` only.
    """

    observations: Dict[str, float] = dict(_tracer_observations(age_years, rng))
    radiocarbon = radiocarbon_pmc(
        age_years,
        float(d13c_permil),
        process_history=process_history,
        c14_a0_error_permil=c14_a0_error_permil,
        iso_rng=iso_rng,
    )
    # The dilution factor and initial activity are latent, so only the
    # measurable activity is exposed to inference.
    observations["c14_pmc"] = radiocarbon["c14_pmc"]
    observations.update(
        halocarbon_panel(
            age_years,
            sample_date=sample_date,
            process_history=process_history,
            recharge_temperature_error_c=recharge_temperature_error_c,
            cfc_degradation_strength=cfc_degradation_strength,
            iso_rng=iso_rng,
        )
    )
    observations.update(
        helium_panel(
            age_years,
            observations["tritium_TU"],
            iso_rng=iso_rng,
        )
    )
    return observations


# ---------------------------------------------------------------------------
# Case generator
# ---------------------------------------------------------------------------
def generate_independent_aquifer_v2(
    seed: int,
    workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    *,
    shared_nuisance: bool = False,
    nuisance_scale: float = 0.0,
    cfc_degradation_strength: float = CFC_DEGRADATION_STRENGTH,
    evaporated_fraction: float = EVAPORATED_PARTICLE_FRACTION,
) -> IndependentAquiferV2:
    """Run the external simulators and return a blind HydroSheaf case.

    Derived from ``independent_modflow_generator.generate_independent_aquifer``.
    The v1 driver is a single monolithic function with no extension point, and
    its source file is byte-locked by the M7 system-acceptance and M8 frontier
    protocols, so the milestone loop is reproduced here rather than edited in
    place.  Every ``rng`` call below -- the per-particle ``lognormal``, the
    ``_chemistry_step`` calls, the v1 tracer draws inside
    :func:`tracer_observations_v2`, the two head perturbations, the pH draw and
    the temperature draw -- appears in exactly the v1 order with exactly the v1
    arguments.  All v2 additions consume ``iso_rng`` instead.
    """

    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    mf6_executable = Path(mf6_executable).resolve()
    mp7_executable = Path(mp7_executable).resolve()
    if not mf6_executable.exists() or not mp7_executable.exists():
        raise FileNotFoundError("Both official mf6 and mp7 executables are required.")

    rng = np.random.default_rng(int(seed))
    # Independent stream: never interleaved with ``rng``.
    iso_rng = np.random.default_rng(int(seed) + ISOTOPE_STREAM_SEED_OFFSET)
    true_nuisance = draw_shared_nuisance(
        int(seed),
        shared_nuisance=bool(shared_nuisance),
        nuisance_scale=float(nuisance_scale),
    )
    gwf, heads, model_metadata = _build_flow_model(int(seed), workspace, mf6_executable)
    pathlines = _run_pathlines(int(seed), workspace, gwf, mp7_executable)
    nrow = int(model_metadata["nrow"])
    ncol = int(model_metadata["ncol"])
    delr = float(model_metadata["delr_m"])
    delc = float(model_metadata["delc_m"])
    milestone_columns = (1, 6, 11, 16)
    # Each case exercises every preregistered reaction family.  The reducing
    # sequence follows the chemically defensible electron-acceptor progression
    # nitrate -> sulfate -> ferric iron.
    process_paths = (
        (
            "carbonate_weathering",
            "carbonate_precipitation",
            "silicate_weathering",
        ),
        (
            "denitrification",
            "sulfate_reduction",
            "iron_reduction",
        ),
        (
            "silicate_weathering",
            "carbonate_weathering",
            "denitrification",
        ),
    )

    observations: List[Dict[str, object]] = []
    true_edges: List[Tuple[str, str]] = []
    true_ages: Dict[str, float] = {}
    true_processes: Dict[str, str] = {}
    true_recharge_isotopes: Dict[str, Dict[str, float]] = {}
    true_reaction_isotopes: Dict[str, Dict[str, float]] = {}
    pathline_rows: List[Dict[str, object]] = []
    base = np.asarray(
        [0.62, 0.27, 0.48, 0.055, 1.62, 0.24, 0.42, 0.30, 0.007, 0.003, 0.004, 0.55],
        dtype=float,
    )

    for particle_index, path in enumerate(pathlines):
        previous_node: str | None = None
        previous_time = 0.0
        chemistry = base * rng.lognormal(0.0, 0.055, base.size)
        particle_id = f"MF{seed}_P{particle_index}"
        evaporation = draw_recharge_evaporation(
            iso_rng,
            evaporated_fraction=float(evaporated_fraction),
        )
        process_history: List[str] = []
        recharge_isotopes: Dict[str, object] = {}
        for milestone_index, column in enumerate(milestone_columns):
            records = path[
                np.floor(np.maximum(path["x"], 0.0) / delr).astype(int) >= column
            ]
            if records.size == 0:
                raise RuntimeError(
                    f"Particle {particle_index} did not reach column {column}."
                )
            record = records[0]
            row = min(
                nrow - 1,
                max(
                    0,
                    nrow - 1 - int(np.floor(max(float(record["y"]), 0.0) / delc)),
                ),
            )
            actual_column = min(
                ncol - 1,
                max(0, int(np.floor(max(float(record["x"]), 0.0) / delr))),
            )
            elapsed_days = float(record["time"])
            age_years = 1.5 + elapsed_days / 365.25
            node_id = f"MF{seed}_P{particle_index}_M{milestone_index}"
            if previous_node is not None:
                process = process_paths[particle_index][milestone_index - 1]
                chemistry = _chemistry_step(
                    chemistry,
                    (elapsed_days - previous_time) / 365.25,
                    process,
                    rng,
                )
                edge_id = f"{previous_node}->{node_id}"
                true_edges.append((previous_node, node_id))
                true_processes[edge_id] = process
                process_history.append(process)
            else:
                process = "source"

            x = float(record["x"])
            y = float(record["y"])
            head = float(heads[0, row, actual_column])
            elevation = head + 12.0
            if milestone_index == 0:
                # Conservative water isotopes are fixed once, at the recharge
                # location, and then carried unchanged down the flow path.
                recharge_isotopes = environmental_isotopes(
                    x,
                    y,
                    elevation,
                    evaporated=bool(evaporation["evaporated"]),
                    evaporation_extent_permil=float(
                        evaporation["evaporation_extent_permil"]
                    ),
                    iso_rng=iso_rng,
                )
                true_recharge_isotopes[particle_id] = {
                    "true_recharge_d18O_permil": float(
                        recharge_isotopes["d18O_permil"]
                    ),
                    "true_recharge_d2H_permil": float(recharge_isotopes["d2H_permil"]),
                    "true_recharge_deuterium_excess_permil": float(
                        recharge_isotopes["deuterium_excess_permil"]
                    ),
                    "true_recharge_evaporated": bool(recharge_isotopes["evaporated"]),
                    "true_recharge_x_m": x,
                    "true_recharge_y_m": y,
                    "true_recharge_elevation_m": elevation,
                }
            reaction_isotopes = reaction_path_isotopes(
                process_history,
                iso_rng=iso_rng,
            )
            tracer = tracer_observations_v2(
                age_years,
                rng,
                iso_rng=iso_rng,
                d13c_permil=float(reaction_isotopes["d13C_permil"]),
                sample_date=2025.5,
                process_history=process_history,
                recharge_temperature_error_c=float(
                    true_nuisance["recharge_temperature_error_c"]
                ),
                c14_a0_error_permil=float(true_nuisance["c14_a0_error_permil"]),
                cfc_degradation_strength=float(cfc_degradation_strength),
            )
            p_h = 6.85 + 0.018 * age_years
            if process == "carbonate_precipitation":
                p_h += 0.55
            sample: Dict[str, object] = {
                "site_id": node_id,
                "sample_id": node_id,
                "x_m": x,
                "y_m": y,
                "lat": 7.0 + y / 111_000.0,
                "lon": -1.5 + x / 111_000.0,
                "head_meas": head + rng.normal(0.0, 0.10),
                "hydraulic_head": head + rng.normal(0.0, 0.10),
                "elevation": head + 12.0,
                "screen_depth": 18.0,
                "well_depth": 40.0,
                "aquifer_unit": "external_modflow_synthetic",
                "aquifer_layer": 1,
                "pH": float(np.clip(p_h + rng.normal(0.0, 0.04), 6.3, 8.8)),
                "temp_c": float(25.0 + rng.normal(0.0, 0.35)),
                "sample_date": 2025.5,
                # Conservative recharge-source isotopes: identical at every
                # milestone of a particle, so they carry no age information.
                "d18O_permil": float(recharge_isotopes["d18O_permil"]),
                "d2H_permil": float(recharge_isotopes["d2H_permil"]),
                # Reaction-path isotopes; the latent mixing fractions are
                # withheld from the observation record.
                "d13C_permil": float(reaction_isotopes["d13C_permil"]),
                "sr87_sr86": float(reaction_isotopes["sr87_sr86"]),
                **tracer,
            }
            for ion_index, ion in enumerate(ION_ORDER):
                sample[ion] = float(chemistry[ion_index])
            observations.append(sample)
            true_ages[node_id] = age_years
            true_reaction_isotopes[node_id] = {
                "true_carbonate_carbon_fraction": float(
                    reaction_isotopes["carbonate_carbon_fraction"]
                ),
                "true_silicate_strontium_fraction": float(
                    reaction_isotopes["silicate_strontium_fraction"]
                ),
            }
            pathline_rows.append(
                {
                    "particle": particle_index,
                    "milestone": milestone_index,
                    "node_id": node_id,
                    "modpath_node_zero_based": int(record["node"]),
                    "row_zero_based": row,
                    "column_zero_based": actual_column,
                    "x_m": x,
                    "y_m": y,
                    "travel_time_days": elapsed_days,
                    "age_years": age_years,
                    "true_recharge_d18O_permil": float(
                        recharge_isotopes["d18O_permil"]
                    ),
                    "true_recharge_d2H_permil": float(recharge_isotopes["d2H_permil"]),
                    "true_recharge_evaporated": bool(recharge_isotopes["evaporated"]),
                    "true_carbonate_carbon_fraction": float(
                        reaction_isotopes["carbonate_carbon_fraction"]
                    ),
                    "true_silicate_strontium_fraction": float(
                        reaction_isotopes["silicate_strontium_fraction"]
                    ),
                }
            )
            previous_node = node_id
            previous_time = elapsed_days

    provenance: Dict[str, object] = {
        "generator": "independent_modflow_modpath_nonlinear_reactive_v2",
        "imports_hydrosheaf": False,
        "generator_streams": [
            "modflow6_heads",
            "modpath7_pathlines",
            "independent_nonlinear_chemistry",
            "independent_analytic_3H_39Ar",
            "independent_conservative_water_isotopes",
            "independent_reaction_path_d13C_87Sr86Sr",
            "independent_analytic_14C",
            "independent_analytic_CFC11_CFC12_CFC113_SF6",
            "independent_analytic_4He_tritiogenic_3He",
            "shared_family_nuisance",
        ],
        "mf6_executable": str(mf6_executable),
        "mf6_sha256": _sha256(mf6_executable),
        "mf6_version": _version(mf6_executable),
        "mp7_executable": str(mp7_executable),
        "mp7_sha256": _sha256(mp7_executable),
        "mp7_version": _version(mp7_executable),
        "flopy_version": flopy.__version__,
        "chemistry_generator": "independent nonlinear process equations",
        "tracer_generator": "independent analytic 3H/39Ar surrogate",
        "isotope_generator": (
            "independent GMWL/evaporation-line water isotopes with "
            "reaction-path d13C and 87Sr/86Sr"
        ),
        "nuclear_panel_generator": (
            "independent analytic 14C/CFC/SF6/He surrogates; no external "
            "input-history table is read"
        ),
        "base_generator_module": "independent_modflow_generator",
        "isotope_random_stream_seed": int(seed) + ISOTOPE_STREAM_SEED_OFFSET,
        "nuisance_random_stream_seed": int(seed) + NUISANCE_STREAM_SEED_OFFSET,
        "shared_nuisance_enabled": bool(shared_nuisance),
        "nuisance_scale": float(nuisance_scale),
        "cfc_degradation_strength": float(cfc_degradation_strength),
        "evaporated_fraction": float(evaporated_fraction),
        "reference_recharge_temperature_c": REFERENCE_RECHARGE_TEMPERATURE_C,
        **model_metadata,
    }
    return IndependentAquiferV2(
        seed=int(seed),
        observations=tuple(observations),
        true_edges=tuple(true_edges),
        true_ages_years=true_ages,
        true_processes=true_processes,
        pathline_rows=tuple(pathline_rows),
        provenance=provenance,
        true_nuisance=dict(true_nuisance),
        true_recharge_isotopes=true_recharge_isotopes,
        true_reaction_isotopes=true_reaction_isotopes,
    )
