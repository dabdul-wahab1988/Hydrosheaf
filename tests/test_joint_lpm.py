from hydrosheaf.nuclear.joint_lpm import (
    build_lpm_tracer_observations,
    fit_lpm_model,
    fit_lpm_models,
    predict_lpm_tracers,
)


def test_joint_dm_fit_recovers_synthetic_multi_tracer_age():
    sample_year = 2020.0
    parameters = {"mean_age_years": 25.0, "dispersion": 0.1}
    predicted = predict_lpm_tracers(
        "DM",
        parameters,
        sample_year,
        ["3H", "3H/3He", "SF6", "CFC12"],
        max_age_years=150.0,
    )
    observations = {
        "tritium_TU": predicted["3H"],
        "tritium_sigma_TU": 0.3,
        "he3_trit_TU": predicted["3H/3He"],
        "he3_trit_sigma_TU": 0.5,
        "sf6_pptv": predicted["SF6"],
        "sf6_sigma_pptv": 0.2,
        "cfc12_pptv": predicted["CFC12"],
        "cfc12_sigma_pptv": 5.0,
    }

    fit = fit_lpm_model(
        observations,
        sample_year=sample_year,
        model="DM",
        max_age_years=150.0,
        age_steps=60,
    )

    assert fit.converged
    assert fit.n_tracers == 4
    assert fit.objective < 5.0
    assert abs(fit.parameters["mean_age_years"] - 25.0) < 8.0
    assert fit.parameters["dispersion"] in {0.01, 0.03, 0.1, 0.3, 1.0}


def test_joint_model_ranking_and_binary_mixture_run():
    sample_year = 2020.0
    predicted = predict_lpm_tracers(
        "EPM",
        {"mean_age_years": 18.0, "piston_fraction": 0.5},
        sample_year,
        ["3H", "SF6", "CFC11"],
        max_age_years=100.0,
    )
    observations = {
        "tritium_TU": predicted["3H"],
        "tritium_sigma_TU": 0.3,
        "sf6_pptv": predicted["SF6"],
        "sf6_sigma_pptv": 0.2,
        "cfc11_pptv": predicted["CFC11"],
        "cfc11_sigma_pptv": 3.0,
    }

    ranked = fit_lpm_models(
        observations,
        sample_year=sample_year,
        models=["PFM", "EM", "EPM", "BMM-DM-DM"],
        max_age_years=100.0,
        age_steps=45,
    )

    assert len(ranked) == 4
    assert ranked[0].converged
    assert ranked[0].aic <= ranked[-1].aic
    assert any(fit.model == "BMM-DM-DM" and fit.converged for fit in ranked)


def test_observation_builder_requires_supported_values():
    observations = build_lpm_tracer_observations(
        {
            "tritium_TU": "",
            "sf6_pptv": -1.0,
            "c14_pmc": 80.0,
        }
    )
    assert len(observations) == 1
    assert observations[0].tracer == "14C"
