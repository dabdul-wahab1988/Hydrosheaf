"""Unit tests for the v2 isotope / nuclear-tracer / shared-nuisance layer.

The MODFLOW 6 and MODPATH 7 executables are not assumed to be present, so these
tests exercise the pure functions directly rather than
``generate_independent_aquifer_v2``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import independent_modflow_generator as v1  # noqa: E402
import independent_modflow_generator_v2 as v2  # noqa: E402


# ---------------------------------------------------------------------------
# A. Conservative water isotopes
# ---------------------------------------------------------------------------
def test_water_isotopes_are_invariant_to_age():
    """d18O/d2H take no age argument and cannot vary along the flow path."""

    import inspect

    parameters = set(inspect.signature(v2.environmental_isotopes).parameters)
    assert "age_years" not in parameters
    assert "process_history" not in parameters

    recharge = v2.environmental_isotopes(150.0, 1200.0, 105.0)
    # The generator fixes the value once at recharge and reuses it verbatim for
    # every downstream milestone, so repeated evaluation is exactly equal.
    for _ in range(4):
        again = v2.environmental_isotopes(150.0, 1200.0, 105.0)
        assert again["d18O_permil"] == recharge["d18O_permil"]
        assert again["d2H_permil"] == recharge["d2H_permil"]


def test_water_isotopes_vary_with_recharge_location():
    near = v2.environmental_isotopes(150.0, 1200.0, 105.0)
    far = v2.environmental_isotopes(2400.0, 1200.0, 105.0)
    lateral = v2.environmental_isotopes(150.0, 150.0, 105.0)
    high = v2.environmental_isotopes(150.0, 1200.0, 305.0)

    # Continental effect: depleted away from the recharge boundary.
    assert far["d18O_permil"] < near["d18O_permil"] - 0.5
    # Latitude effect.
    assert lateral["d18O_permil"] != pytest.approx(near["d18O_permil"])
    # Altitude effect: depleted at higher recharge elevation.
    assert high["d18O_permil"] < near["d18O_permil"]


def test_non_evaporated_samples_sit_on_the_meteoric_line():
    iso_rng = np.random.default_rng(11)
    residuals = []
    for index in range(200):
        result = v2.environmental_isotopes(
            50.0 * (index % 40),
            30.0 * (index % 50),
            100.0 + (index % 25),
            iso_rng=iso_rng,
        )
        residuals.append(
            result["d2H_permil"]
            - (v2.GMWL_SLOPE * result["d18O_permil"] + v2.GMWL_INTERCEPT_PERMIL)
        )
    residuals = np.asarray(residuals)
    # Only the declared deuterium-excess perturbation separates them from the
    # global meteoric water line.
    assert abs(float(np.mean(residuals))) < 0.5
    assert float(np.max(np.abs(residuals))) < 6.0 * v2.DEUTERIUM_EXCESS_SIGMA_PERMIL


def test_evaporated_samples_leave_the_meteoric_line_on_a_slope_five_trend():
    base = v2.environmental_isotopes(600.0, 600.0, 110.0)
    evaporated = v2.environmental_isotopes(
        600.0,
        600.0,
        110.0,
        evaporated=True,
        evaporation_extent_permil=2.0,
    )
    slope = (evaporated["d2H_permil"] - base["d2H_permil"]) / (
        evaporated["d18O_permil"] - base["d18O_permil"]
    )
    assert slope == pytest.approx(v2.EVAPORATION_LINE_SLOPE)
    # Enriched in 18O and lower in deuterium excess.
    assert evaporated["d18O_permil"] > base["d18O_permil"]
    assert (
        evaporated["deuterium_excess_permil"] < base["deuterium_excess_permil"] - 1.0
    )


def test_declared_evaporated_fraction_is_respected():
    iso_rng = np.random.default_rng(3)
    draws = [v2.draw_recharge_evaporation(iso_rng) for _ in range(4000)]
    fraction = float(np.mean([bool(draw["evaporated"]) for draw in draws]))
    assert fraction == pytest.approx(v2.EVAPORATED_PARTICLE_FRACTION, abs=0.03)
    assert all(
        draw["evaporation_extent_permil"] == 0.0
        for draw in draws
        if not draw["evaporated"]
    )


# ---------------------------------------------------------------------------
# B. Reaction-path isotopes
# ---------------------------------------------------------------------------
def test_d13c_moves_toward_heavier_values_under_carbonate_weathering():
    source = v2.reaction_path_isotopes([])
    one_step = v2.reaction_path_isotopes(["carbonate_weathering"])
    two_steps = v2.reaction_path_isotopes(
        ["carbonate_weathering", "carbonate_weathering"]
    )

    assert source["d13C_permil"] == pytest.approx(v2.SOIL_CO2_D13C_PERMIL)
    assert one_step["d13C_permil"] > source["d13C_permil"]
    assert two_steps["d13C_permil"] > one_step["d13C_permil"]
    # It never overshoots the marine carbonate endmember.
    assert two_steps["d13C_permil"] < v2.MARINE_CARBONATE_D13C_PERMIL
    assert 0.0 < one_step["carbonate_carbon_fraction"] < 1.0


def test_redox_and_precipitation_shift_d13c_lighter():
    reference = v2.reaction_path_isotopes(["carbonate_weathering"])
    for process in (
        "carbonate_precipitation",
        "sulfate_reduction",
        "iron_reduction",
    ):
        shifted = v2.reaction_path_isotopes(["carbonate_weathering", process])
        assert shifted["d13C_permil"] < reference["d13C_permil"]


def test_strontium_moves_toward_the_silicate_endmember():
    source = v2.reaction_path_isotopes([])
    one_step = v2.reaction_path_isotopes(["silicate_weathering"])
    three_steps = v2.reaction_path_isotopes(["silicate_weathering"] * 3)

    assert source["sr87_sr86"] == pytest.approx(v2.SR87_SR86_CARBONATE)
    assert v2.SR87_SR86_CARBONATE < one_step["sr87_sr86"] < three_steps["sr87_sr86"]
    assert three_steps["sr87_sr86"] < v2.SR87_SR86_SILICATE
    # Carbonate weathering must not move the strontium ratio.
    assert v2.reaction_path_isotopes(["carbonate_weathering"])[
        "sr87_sr86"
    ] == pytest.approx(v2.SR87_SR86_CARBONATE)


# ---------------------------------------------------------------------------
# C. Radiocarbon
# ---------------------------------------------------------------------------
def test_c14_decreases_monotonically_with_age():
    ages = [0.0, 100.0, 1000.0, 5730.0, 12000.0, 30000.0]
    activities = [
        v2.radiocarbon_pmc(age, v2.SOIL_CO2_D13C_PERMIL)["c14_pmc"] for age in ages
    ]
    assert all(
        later < earlier for earlier, later in zip(activities, activities[1:])
    )
    # One half-life reproduces the Libby/Cambridge decay law exactly.
    assert activities[3] == pytest.approx(0.5 * v2.CARBON14_INITIAL_PMC)


def test_c14_is_diluted_when_d13c_indicates_carbonate_dissolution():
    soil_only = v2.radiocarbon_pmc(1000.0, v2.SOIL_CO2_D13C_PERMIL)
    half_carbonate = v2.radiocarbon_pmc(1000.0, 0.5 * v2.SOIL_CO2_D13C_PERMIL)

    assert soil_only["c14_dilution_factor"] == pytest.approx(1.0)
    assert half_carbonate["c14_dilution_factor"] == pytest.approx(0.5)
    assert half_carbonate["c14_pmc"] == pytest.approx(0.5 * soil_only["c14_pmc"])


def test_c14_and_d13c_share_the_reaction_path_nuisance():
    """The same carbonate weathering shifts d13C up and 14C down together."""

    age = 2500.0
    d13c_source = v2.reaction_path_isotopes([])["d13C_permil"]
    d13c_weathered = v2.reaction_path_isotopes(["carbonate_weathering"])["d13C_permil"]
    c14_source = v2.radiocarbon_pmc(age, d13c_source)["c14_pmc"]
    c14_weathered = v2.radiocarbon_pmc(
        age, d13c_weathered, process_history=["carbonate_weathering"]
    )["c14_pmc"]

    assert d13c_weathered > d13c_source
    assert c14_weathered < c14_source


def test_c14_a0_nuisance_scales_the_initial_activity_and_defaults_to_neutral():
    unbiased = v2.radiocarbon_pmc(1000.0, -15.0)
    biased = v2.radiocarbon_pmc(1000.0, -15.0, c14_a0_error_permil=2.5)

    assert unbiased["c14_initial_activity_pmc"] == pytest.approx(
        v2.CARBON14_INITIAL_PMC
    )
    assert biased["c14_initial_activity_pmc"] > v2.CARBON14_INITIAL_PMC
    assert biased["c14_pmc"] > unbiased["c14_pmc"]
    # The dilution factor, which depends only on the observed d13C, is
    # untouched: the nuisance acts purely on A0.
    assert biased["c14_dilution_factor"] == pytest.approx(
        unbiased["c14_dilution_factor"]
    )


# ---------------------------------------------------------------------------
# D. Halocarbons and SF6
# ---------------------------------------------------------------------------
def test_atmospheric_histories_have_the_expected_shape():
    # Effectively absent before industrial release.
    for gas in v2.HALOCARBON_GASES:
        assert v2.atmospheric_mixing_ratio(gas, 1900.0) < 0.5
    # CFC peaks are near their published years and magnitudes.
    assert v2.atmospheric_mixing_ratio("cfc11_pptv", 1993.5) > v2.atmospheric_mixing_ratio(
        "cfc11_pptv", 2025.0
    )
    assert v2.atmospheric_mixing_ratio("cfc12_pptv", 2002.0) > v2.atmospheric_mixing_ratio(
        "cfc12_pptv", 2025.0
    )
    assert v2.atmospheric_mixing_ratio("cfc113_pptv", 1996.0) > v2.atmospheric_mixing_ratio(
        "cfc113_pptv", 2025.0
    )
    # SF6 is still rising.
    assert v2.atmospheric_mixing_ratio("sf6_pptv", 2025.0) > v2.atmospheric_mixing_ratio(
        "sf6_pptv", 2000.0
    )


def test_cfc11_degrades_under_sulfate_reduction_but_cfc12_does_not():
    oxic = v2.halocarbon_panel(15.0)
    reduced = v2.halocarbon_panel(15.0, process_history=["sulfate_reduction"])

    cfc11_ratio = reduced["cfc11_pptv"] / oxic["cfc11_pptv"]
    cfc12_ratio = reduced["cfc12_pptv"] / oxic["cfc12_pptv"]
    cfc113_ratio = reduced["cfc113_pptv"] / oxic["cfc113_pptv"]
    sf6_ratio = reduced["sf6_pptv"] / oxic["sf6_pptv"]

    assert cfc11_ratio < 0.6
    assert cfc12_ratio > 0.97
    # CFC-113 degrades, but less than CFC-11.
    assert cfc11_ratio < cfc113_ratio < 1.0
    # SF6 is unaffected.
    assert sf6_ratio == pytest.approx(1.0)


def test_iron_reduction_also_degrades_cfc11_and_strength_is_a_parameter():
    oxic = v2.halocarbon_panel(15.0)
    iron = v2.halocarbon_panel(15.0, process_history=["iron_reduction"])
    assert iron["cfc11_pptv"] < oxic["cfc11_pptv"]

    weak = v2.halocarbon_panel(
        15.0, process_history=["sulfate_reduction"], cfc_degradation_strength=0.0
    )
    assert weak["cfc11_pptv"] == pytest.approx(oxic["cfc11_pptv"])


def test_recharge_temperature_error_is_common_mode_across_the_family():
    unbiased = v2.halocarbon_panel(15.0)
    biased = v2.halocarbon_panel(15.0, recharge_temperature_error_c=2.0)

    ratios = {gas: biased[gas] / unbiased[gas] for gas in v2.HALOCARBON_GASES}
    # Every gas moves in the same direction ...
    assert all(ratio < 1.0 for ratio in ratios.values())
    # ... by a similar, solubility-determined amount.
    assert max(ratios.values()) - min(ratios.values()) < 0.05

    cooled = v2.halocarbon_panel(15.0, recharge_temperature_error_c=-2.0)
    assert all(cooled[gas] > unbiased[gas] for gas in v2.HALOCARBON_GASES)


# ---------------------------------------------------------------------------
# E. Helium
# ---------------------------------------------------------------------------
def test_he4_accumulates_with_age_and_he3_follows_the_reported_tritium():
    young = v2.helium_panel(2.0, 5.0)
    old = v2.helium_panel(150.0, 5.0)
    assert old["he4_ccpg"] > young["he4_ccpg"] > v2.HE4_AIR_SATURATED_CCPG * 0.99

    # Tritiogenic 3He is the decay-consistent partner of the *given* tritium.
    tritium = 4.0
    age = 24.64  # exactly two tritium half-lives
    expected = (
        v2.TRITIOGENIC_HE3_RETENTION
        * tritium
        * (2.0 ** (age / v2.TRITIUM_HALF_LIFE_YEARS) - 1.0)
    )
    assert v2.helium_panel(age, tritium)["h3_he3_TU"] == pytest.approx(expected)
    assert v2.helium_panel(0.0, tritium)["h3_he3_TU"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# F. Shared nuisance
# ---------------------------------------------------------------------------
def test_shared_nuisance_is_exactly_zero_by_default():
    for seed in (5201, 5301, 6301):
        assert v2.draw_shared_nuisance(seed) == {
            "recharge_temperature_error_c": 0.0,
            "c14_a0_error_permil": 0.0,
        }
        # Explicitly disabled, or enabled with zero scale, is also exact zero.
        assert v2.draw_shared_nuisance(
            seed, shared_nuisance=False, nuisance_scale=1.0
        ) == {"recharge_temperature_error_c": 0.0, "c14_a0_error_permil": 0.0}
        assert v2.draw_shared_nuisance(
            seed, shared_nuisance=True, nuisance_scale=0.0
        ) == {"recharge_temperature_error_c": 0.0, "c14_a0_error_permil": 0.0}


def test_shared_nuisance_is_per_case_deterministic_and_nonzero_when_enabled():
    first = v2.draw_shared_nuisance(5201, shared_nuisance=True, nuisance_scale=1.0)
    repeat = v2.draw_shared_nuisance(5201, shared_nuisance=True, nuisance_scale=1.0)
    other = v2.draw_shared_nuisance(5202, shared_nuisance=True, nuisance_scale=1.0)

    assert first == repeat
    assert first != other
    assert first["recharge_temperature_error_c"] != 0.0
    assert first["c14_a0_error_permil"] != 0.0

    scaled = v2.draw_shared_nuisance(5201, shared_nuisance=True, nuisance_scale=2.0)
    assert scaled["recharge_temperature_error_c"] == pytest.approx(
        2.0 * first["recharge_temperature_error_c"]
    )


def test_generate_v2_exposes_zero_default_nuisance_keyword_arguments():
    import inspect

    parameters = inspect.signature(v2.generate_independent_aquifer_v2).parameters
    assert parameters["shared_nuisance"].default is False
    assert parameters["nuisance_scale"].default == 0.0
    assert parameters["shared_nuisance"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["nuisance_scale"].kind is inspect.Parameter.KEYWORD_ONLY


# ---------------------------------------------------------------------------
# G. Reproducibility guarantees for the locked v1 stream
# ---------------------------------------------------------------------------
def _rng_state(rng: np.random.Generator):
    return rng.bit_generator.state


def test_new_functions_do_not_advance_the_main_rng():
    rng = np.random.default_rng(5201)
    iso_rng = np.random.default_rng(5201 + v2.ISOTOPE_STREAM_SEED_OFFSET)
    before = _rng_state(rng)

    v2.draw_shared_nuisance(5201, shared_nuisance=True, nuisance_scale=1.0)
    v2.draw_recharge_evaporation(iso_rng)
    v2.environmental_isotopes(120.0, 900.0, 104.0, iso_rng=iso_rng)
    v2.reaction_path_isotopes(["carbonate_weathering"], iso_rng=iso_rng)
    v2.radiocarbon_pmc(500.0, -15.0, iso_rng=iso_rng)
    v2.halocarbon_panel(20.0, iso_rng=iso_rng)
    v2.helium_panel(20.0, 3.0, iso_rng=iso_rng)
    v2.atmospheric_mixing_ratio("cfc11_pptv", 1990.0)
    v2.reducing_intensity(["sulfate_reduction"])

    assert _rng_state(rng) == before


def test_tracer_observations_v2_consumes_the_v1_rng_identically():
    """The v1 3H/39Ar values and rng consumption are bit-identical."""

    for age in (1.5, 12.0, 45.0, 130.0):
        rng_v1 = np.random.default_rng(4242)
        rng_v2 = np.random.default_rng(4242)
        iso_rng = np.random.default_rng(4242 + v2.ISOTOPE_STREAM_SEED_OFFSET)

        baseline = v1._tracer_observations(age, rng_v1)
        extended = v2.tracer_observations_v2(
            age,
            rng_v2,
            iso_rng=iso_rng,
            d13c_permil=-18.0,
            process_history=["carbonate_weathering"],
        )

        assert extended["tritium_TU"] == baseline["tritium_TU"]
        assert extended["argon39_pmc"] == baseline["argon39_pmc"]
        assert _rng_state(rng_v2) == _rng_state(rng_v1)
        # And the new keys are all present.
        assert {
            "c14_pmc",
            "cfc11_pptv",
            "cfc12_pptv",
            "cfc113_pptv",
            "sf6_pptv",
            "he4_ccpg",
            "h3_he3_TU",
        } <= set(extended)


def test_v2_import_does_not_alter_v1_output_for_a_fixed_seed():
    """Importing and using v2 must leave v1 behaviour untouched."""

    def v1_sequence() -> list:
        rng = np.random.default_rng(7331)
        chemistry = np.asarray(
            [0.62, 0.27, 0.48, 0.055, 1.62, 0.24, 0.42, 0.30, 0.007, 0.003, 0.004, 0.55]
        ) * rng.lognormal(0.0, 0.055, len(v1.ION_ORDER))
        out = []
        for process in ("carbonate_weathering", "sulfate_reduction", "silicate_weathering"):
            chemistry = v1._chemistry_step(chemistry, 7.5, process, rng)
            out.append(chemistry.copy())
            out.append(v1._tracer_observations(20.0, rng))
            out.append(float(rng.normal(0.0, 0.10)))
        return out

    expected = v1_sequence()

    # Exercise the entire v2 surface, including its own generators.
    iso_rng = np.random.default_rng(7331 + v2.ISOTOPE_STREAM_SEED_OFFSET)
    v2.draw_shared_nuisance(7331, shared_nuisance=True, nuisance_scale=1.5)
    v2.draw_recharge_evaporation(iso_rng)
    v2.environmental_isotopes(300.0, 300.0, 101.0, iso_rng=iso_rng)
    v2.reaction_path_isotopes(["silicate_weathering"], iso_rng=iso_rng)
    v2.tracer_observations_v2(
        20.0,
        np.random.default_rng(999),
        iso_rng=iso_rng,
        d13c_permil=-20.0,
    )

    actual = v1_sequence()
    assert len(actual) == len(expected)
    for left, right in zip(actual, expected):
        if isinstance(left, np.ndarray):
            assert np.array_equal(left, right)
        else:
            assert left == right

    # v2 re-exports v1's objects rather than redefining them.
    assert v2.ION_ORDER is v1.ION_ORDER
    assert v2.IndependentAquifer is v1.IndependentAquifer


def test_v1_module_source_is_unmodified_by_this_extension():
    """The v1 file is byte-locked by the M7/M8 protocol locks."""

    import hashlib

    locked = "3b06a4bd5c6807e4701c0b92114586f895a9f87ce2445c1cb01e9274564c46ce"
    source = SCRIPT_DIR / "independent_modflow_generator.py"
    raw = source.read_bytes()
    # The lock digest is over the LF form committed to git.  Normalise before
    # hashing so the test reports a genuine content change rather than a local
    # ``core.autocrlf`` checkout artefact; note that the on-disk file must also
    # be LF for ``run_m7_system_acceptance._verify_lock`` itself to pass.
    normalised = raw.replace(b"\r\n", b"\n")
    assert hashlib.sha256(normalised).hexdigest() == locked, (
        "independent_modflow_generator.py must stay byte-identical to the "
        "M7 system-acceptance and M8 frontier protocol locks"
    )
