import math

from hydrosheaf.nuclear.multi_tracer import (
    build_atmospheric_tracer_input,
    infer_3h_3he_age,
    infer_atmospheric_tracer_age,
    infer_multi_tracer_age,
)


def test_3h_3he_closed_system_age():
    age = infer_3h_3he_age(2.0, 2.0, 0.1, 0.1)
    assert age is not None
    expected = math.log(2.0) / (math.log(2.0) / 12.32)
    assert abs(age.age_years - expected) < 0.1
    assert age.ci_low_years < age.age_years < age.ci_high_years


def test_sf6_atmospheric_history_age():
    history = build_atmospheric_tracer_input("SF6")
    sample_year = 2020.0
    recharge_year = 2000.0
    obs = history.interpolate([recharge_year])[0][0]
    age = infer_atmospheric_tracer_age(obs, sample_year, "SF6", observed_sigma_pptv=0.1)
    assert age is not None
    assert abs(age.age_years - 20.0) < 1.0


def test_multi_tracer_uses_richer_inputs():
    estimate = infer_multi_tracer_age(
        {
            "tritium_TU": 1.1,
            "tritium_sigma_TU": 0.1,
            "he3_trit_TU": 2.1,
            "he3_trit_sigma_TU": 0.4,
            "sf6_pptv": 3.7,
            "sf6_sigma_pptv": 0.7,
            "c14_pmc": 67.3,
            "c14_sigma_pmc": 13.5,
        },
        sample_year=2016.64,
    )
    assert estimate["n_tracers"] >= 4
    assert "3H/3He" in estimate["tracers"]
    assert "SF6" in estimate["tracers"]
    assert estimate["age_years"] > 0


def test_uncorroborated_gas_does_not_override_old_14c():
    estimate = infer_multi_tracer_age(
        {
            "sf6_pptv": 8.0,
            "sf6_sigma_pptv": 0.4,
            "c14_pmc": 25.0,
            "c14_sigma_pmc": 1.2,
        },
        sample_year=2016.0,
    )
    assert estimate["n_tracers"] == 1
    assert estimate["tracers"] == ["14C"]
    assert estimate["skipped_estimates"]
