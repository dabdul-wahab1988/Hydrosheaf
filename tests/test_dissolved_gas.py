from hydrosheaf.nuclear.dissolved_gas import (
    atmospheric_equivalent_from_dissolved,
    correct_environmental_tracers,
    equilibrium_concentration_ccstp_g,
    excess_air_concentration_ccstp_g,
    fit_and_correct_dissolved_gases,
    fit_dissolved_gas_model,
    monte_carlo_dissolved_gas_fit,
    predict_dissolved_gas_ccstp_g,
    to_ccstp_g,
)


def _synthetic_noble_gases(temp_c=12.0, excess_air=0.006):
    observations = {}
    for gas in ["NE", "AR", "KR", "XE", "N2"]:
        value = predict_dissolved_gas_ccstp_g(
            gas,
            recharge_temperature_c=temp_c,
            excess_air_ccstp_g=excess_air,
            model="UA",
            salinity_ppt=0.0,
            pressure_atm=1.0,
            excess_n2_ccstp_g=0.001 if gas == "N2" else 0.0,
        )
        lower = gas.lower()
        observations[f"{lower}_ccstp_g"] = value
        observations[f"{lower}_sigma_ccstp_g"] = max(value * 0.02, 1.0e-12)
    return observations


def _synthetic_ce_noble_gases(temp_c=14.0, excess_air=0.010, fractionation=0.5):
    observations = {}
    for gas in ["NE", "AR", "KR", "XE"]:
        value = predict_dissolved_gas_ccstp_g(
            gas,
            recharge_temperature_c=temp_c,
            excess_air_ccstp_g=excess_air,
            model="CE",
            fractionation=fractionation,
            salinity_ppt=0.0,
            pressure_atm=1.0,
        )
        lower = gas.lower()
        observations[f"{lower}_ccstp_g"] = value
        observations[f"{lower}_sigma_ccstp_g"] = max(value * 0.02, 1.0e-12)
    return observations


def test_fit_dissolved_gas_model_recovers_recharge_conditions():
    observations = _synthetic_noble_gases()
    fits = fit_dissolved_gas_model(
        observations,
        models=["UA", "CE"],
        temperature_range_c=(8.0, 16.0),
        temperature_step_c=1.0,
        excess_air_range_ccstp_g=(0.0, 0.012),
        excess_air_steps=25,
        excess_n2_range_ccstp_g=(0.0, 0.002),
        excess_n2_steps=5,
    )

    best = fits[0]
    assert best.converged
    assert abs(best.recharge_temperature_c - 12.0) <= 1.0
    assert abs(best.excess_air_ccstp_g - 0.006) <= 0.001
    assert best.n_gases == 5
    assert best.rmse_standardized < 1.0


def test_ce_model_recovers_fractionated_closed_system_solution():
    fits = fit_dissolved_gas_model(
        _synthetic_ce_noble_gases(),
        models=["CE"],
        temperature_range_c=(12.0, 16.0),
        temperature_step_c=0.5,
        excess_air_range_ccstp_g=(0.006, 0.014),
        excess_air_steps=33,
        fractionation_grid=[0.0, 0.25, 0.5, 0.75, 0.9],
    )

    best = fits[0]
    assert best.converged
    assert abs(best.recharge_temperature_c - 14.0) <= 0.5
    assert abs(best.excess_air_ccstp_g - 0.010) <= 0.001
    assert abs(best.fractionation - 0.5) <= 0.001
    assert best.rmse_standardized < 0.1


def test_correct_environmental_tracer_from_dissolved_sf6():
    observations = _synthetic_noble_gases()
    fit = fit_dissolved_gas_model(
        observations,
        models=["UA"],
        temperature_range_c=(12.0, 12.0),
        temperature_step_c=1.0,
        excess_air_range_ccstp_g=(0.006, 0.006),
        excess_air_steps=2,
        excess_n2_range_ccstp_g=(0.001, 0.001),
        excess_n2_steps=1,
    )[0]
    sf6_pptv = 7.5
    sf6_eq = equilibrium_concentration_ccstp_g(
        "SF6",
        fit.recharge_temperature_c,
        salinity_ppt=fit.salinity_ppt,
        pressure_atm=fit.recharge_pressure_atm,
    )
    dissolved_sf6 = sf6_pptv * (
        sf6_eq + excess_air_concentration_ccstp_g(
            "SF6",
            fit.excess_air_ccstp_g,
            model=fit.model,
            fractionation=fit.fractionation,
        )
    )
    corrected = correct_environmental_tracers(
        {**observations, "sf6_dissolved": dissolved_sf6, "sf6_unit": "ccSTP/g"},
        fit,
    )

    assert abs(corrected["sf6_pptv"] - sf6_pptv) < 0.2
    assert abs(atmospheric_equivalent_from_dissolved("SF6", dissolved_sf6, fit) - sf6_pptv) < 0.2


def test_fit_and_correct_returns_model_and_corrected_observations():
    observations = _synthetic_noble_gases()
    sf6_pptv = 5.0
    sf6_eq = equilibrium_concentration_ccstp_g("SF6", 12.0)
    best, corrected = fit_and_correct_dissolved_gases(
        {
            **observations,
            "sf6_dissolved": sf6_pptv * (sf6_eq + excess_air_concentration_ccstp_g("SF6", 0.006)),
            "sf6_unit": "ccSTP/g",
        },
        models=["UA"],
        temperature_range_c=(10.0, 14.0),
        temperature_step_c=1.0,
        excess_air_range_ccstp_g=(0.004, 0.008),
        excess_air_steps=9,
        excess_n2_range_ccstp_g=(0.0, 0.002),
        excess_n2_steps=5,
    )

    assert best.converged
    assert "sf6_pptv" in corrected
    assert corrected["dissolved_gas_model"]["corrected_tracers"]["SF6"] == corrected["sf6_pptv"]


def test_monte_carlo_dissolved_gas_fit_summary():
    result = monte_carlo_dissolved_gas_fit(
        _synthetic_noble_gases(),
        n=8,
        models=["UA"],
        temperature_range_c=(10.0, 14.0),
        temperature_step_c=1.0,
        excess_air_range_ccstp_g=(0.004, 0.008),
        excess_air_steps=9,
        excess_n2_range_ccstp_g=(0.0, 0.002),
        excess_n2_steps=5,
    )

    assert result["n"] == 8
    assert "recharge_temperature_c" in result["summary"]
    assert result["models"]["UA"] == 8


def test_unit_conversion_mg_l_to_ccstp_g():
    cc = to_ccstp_g(20.0, "N2", "mg/L")
    assert 0.01 < cc < 0.02
