import math

import pytest

from hydrosheaf.nuclear.joint_lpm import (
    TracerFitObservation,
    age_fraction_predictions,
    build_lpm_tracer_observations,
    compute_gated_bma_age,
    fit_lpm_model,
    fit_lpm_models,
    predict_lpm_tracers,
    _standardized_observation_loss,
    _aicc,
    JointLpmFit,
    _evaluate_lpm_candidate,
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


def test_prediction_scale_factors_scale_model_output_not_observations():
    parameters = {"mean_age_years": 12.0, "dispersion": 0.1}
    unscaled = predict_lpm_tracers(
        "DM", parameters, 2020.0, ["3H", "3H/3He", "14C"], max_age_years=200.0
    )
    scaled = predict_lpm_tracers(
        "DM",
        parameters,
        2020.0,
        ["3H", "3H/3He", "14C"],
        max_age_years=200.0,
        prediction_scale_factors={"3H": 1.4, "3H/3He": 1.4, "14C": 0.3},
    )

    assert scaled["3H"] == pytest.approx(1.4 * unscaled["3H"])
    assert scaled["3H/3He"] == pytest.approx(1.4 * unscaled["3H/3He"])
    assert scaled["14C"] == pytest.approx(0.3 * unscaled["14C"])


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


def test_reported_parameter_configuration_keeps_declared_bmm_parameters_fixed():
    sample_year = 2020.0
    parameters = {
        "mean_age_1_years": 8.0,
        "mean_age_2_years": 75.0,
        "binary_fraction": 0.65,
        "dispersion": 0.01,
        "dispersion_2": 0.01,
    }
    predicted = predict_lpm_tracers(
        "BMM-DM-DM", parameters, sample_year, ["3H", "3H/3He", "14C"]
    )
    observations = {
        "tritium_TU": predicted["3H"],
        "tritium_sigma_TU": 0.2,
        "he3_trit_TU": predicted["3H/3He"],
        "he3_trit_sigma_TU": 0.3,
        "c14_pmc": predicted["14C"],
        "c14_sigma_pmc": 2.0,
    }
    fit = fit_lpm_model(
        observations,
        sample_year=sample_year,
        model="BMM-DM-DM",
        age_steps=30,
        refine=True,
        parameter_template={
            "mean_age_1_years": 5.0,
            "mean_age_2_years": 75.0,
            "binary_fraction": 0.5,
            "dispersion": 0.01,
            "dispersion_2": 0.01,
        },
        free_parameters=("mean_age_1_years", "binary_fraction"),
    )

    assert fit.converged
    assert fit.effective_n_params == 2
    assert fit.parameters["mean_age_2_years"] == pytest.approx(75.0)
    assert fit.parameters["dispersion"] == pytest.approx(0.01)
    assert fit.parameters["dispersion_2"] == pytest.approx(0.01)


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


def test_old_tracer_window_preserves_subannual_pfm_resolution():
    predictions = [
        predict_lpm_tracers(
            "PFM",
            {"mean_age_years": age},
            2020.0,
            ["4He"],
            max_age_years=50000.0,
        )["4He"]
        for age in (0.5, 1.0, 2.0, 4.0)
    ]
    assert predictions == sorted(predictions)
    assert len(set(predictions)) == len(predictions)


# --- Phase 4: Censored/Contaminated Likelihoods ---

def test_upper_censored_zero_loss_when_predicted_below_observed():
    obs = TracerFitObservation(tracer="SF6", value=10.0, sigma=1.0, units="pptv", weight=1.0, likelihood="upper_censored")
    assert _standardized_observation_loss(obs, 5.0) == 0.0
    assert _standardized_observation_loss(obs, 10.0) == 0.0


def test_upper_censored_penalizes_prediction_above_observed():
    obs = TracerFitObservation(tracer="SF6", value=10.0, sigma=1.0, units="pptv", weight=1.0, likelihood="upper_censored")
    loss = _standardized_observation_loss(obs, 15.0)
    assert loss > 0.0
    assert loss == 25.0


def test_contaminated_mixture_loss_smaller_than_gaussian_for_large_residual():
    obs_gaussian = TracerFitObservation(tracer="SF6", value=10.0, sigma=1.0, units="pptv", weight=1.0, likelihood="gaussian")
    obs_robust = TracerFitObservation(tracer="SF6", value=10.0, sigma=1.0, units="pptv", weight=1.0, likelihood="contaminated_mixture")
    predicted = 50.0
    gaussian_loss = _standardized_observation_loss(obs_gaussian, predicted)
    robust_loss = _standardized_observation_loss(obs_robust, predicted)
    assert robust_loss < gaussian_loss


# --- Phase 6: Age Fraction Predictions ---

def test_young_pfm_predicts_high_anthropocene_fraction():
    fracs = age_fraction_predictions("PFM", {"mean_age_years": 5.0})
    assert fracs["anthropocene"] > 0.8
    assert fracs["pleistocene"] < 0.1


def test_old_pfm_predicts_high_pleistocene_fraction():
    fracs = age_fraction_predictions("PFM", {"mean_age_years": 20000.0})
    assert fracs["pleistocene"] > 0.8
    assert fracs["anthropocene"] < 0.1


def test_age_fractions_sum_to_approx_one():
    fracs = age_fraction_predictions("DM", {"mean_age_years": 500.0, "dispersion": 0.3})
    total = fracs["anthropocene"] + fracs["holocene"] + fracs["pleistocene"]
    assert abs(total - 1.0) < 0.05


# --- Phase 7: Gated AICc/LOO BMA ---

def test_aicc_infinite_when_n_le_k_plus_one():
    assert _aicc(aic=10.0, n=3, k=2) == float("inf")
    assert math.isfinite(_aicc(aic=10.0, n=10, k=2))


def _make_fit(model, params, aic, aicc, converged=True):
    from hydrosheaf.nuclear.joint_lpm import TracerFitResidual
    return JointLpmFit(
        model=model,
        parameters=params,
        objective=1.0,
        rmse_standardized=0.1,
        aic=aic,
        bic=aic + 2.0,
        aicc=aicc,
        effective_n_params=2,
        n_tracers=3,
        residuals=[TracerFitResidual("3H", 1.0, 1.0, 0.1, 0.0)],
        converged=converged,
    )


def test_divergent_model_ages_trigger_top_model_fallback():
    fits = [
        _make_fit("PFM", {"mean_age_years": 5.0}, aic=10.0, aicc=12.0),
        _make_fit("DM", {"mean_age_years": 5000.0}, aic=11.0, aicc=13.0),
    ]
    bma = compute_gated_bma_age(fits, max_delta_aicc=4.0, max_log_age_span=0.7)
    assert bma["bma_used"] is False
    assert "log_age_span_exceeded" in bma["bma_skip_reason"]


def test_close_model_ages_allow_bma():
    fits = [
        _make_fit("PFM", {"mean_age_years": 20.0}, aic=10.0, aicc=12.0),
        _make_fit("DM", {"mean_age_years": 25.0}, aic=11.0, aicc=13.0),
    ]
    bma = compute_gated_bma_age(fits, max_delta_aicc=4.0, max_log_age_span=0.7)
    assert bma["bma_used"] is True
    assert bma["bma_n_models"] == 2


def test_sparse_suite_uses_aic_when_aicc_is_undefined():
    fits = [
        _make_fit("PFM", {"mean_age_years": 20.0}, aic=10.0, aicc=float("inf")),
        _make_fit("DM", {"mean_age_years": 25.0}, aic=11.0, aicc=float("inf")),
    ]
    bma = compute_gated_bma_age(fits)
    assert math.isfinite(bma["age_years"])
    assert bma["bma_used"] is True
    assert bma["bma_information_criterion"] == "aic_fallback"


# --- Phase 8: BMM Refinement ---

def test_bmm_dm_dm_refinement_preserves_mean_age_2_greater_than_1():
    sample_year = 2020.0
    predicted = predict_lpm_tracers(
        "BMM-DM-DM",
        {"mean_age_1_years": 5.0, "mean_age_2_years": 50.0, "binary_fraction": 0.4, "dispersion": 0.3, "dispersion_2": 0.3},
        sample_year,
        ["3H", "SF6", "CFC12"],
        max_age_years=100.0,
    )
    observations = {
        "tritium_TU": predicted["3H"],
        "tritium_sigma_TU": 0.3,
        "sf6_pptv": predicted["SF6"],
        "sf6_sigma_pptv": 0.2,
        "cfc12_pptv": predicted["CFC12"],
        "cfc12_sigma_pptv": 5.0,
    }
    fit = fit_lpm_model(
        observations,
        sample_year=sample_year,
        model="BMM-DM-DM",
        max_age_years=100.0,
        age_steps=45,
        refine=True,
    )
    assert fit.converged
    params = fit.parameters
    assert params["mean_age_2_years"] > params["mean_age_1_years"]
    assert 0.001 <= params["binary_fraction"] <= 0.999


def test_bmm_refinement_objective_not_worse_than_grid():
    sample_year = 2020.0
    predicted = predict_lpm_tracers(
        "BMM-DM-DM",
        {"mean_age_1_years": 10.0, "mean_age_2_years": 80.0, "binary_fraction": 0.3, "dispersion": 0.1, "dispersion_2": 0.1},
        sample_year,
        ["3H", "SF6"],
        max_age_years=100.0,
    )
    observations = {
        "tritium_TU": predicted["3H"],
        "tritium_sigma_TU": 0.3,
        "sf6_pptv": predicted["SF6"],
        "sf6_sigma_pptv": 0.2,
    }
    fit_grid = fit_lpm_model(
        observations,
        sample_year=sample_year,
        model="BMM-DM-DM",
        max_age_years=100.0,
        age_steps=45,
        refine=False,
    )
    fit_refined = fit_lpm_model(
        observations,
        sample_year=sample_year,
        model="BMM-DM-DM",
        max_age_years=100.0,
        age_steps=45,
        refine=True,
    )
    assert fit_refined.converged
    assert fit_refined.objective <= fit_grid.objective * 1.05


def test_fit_lpm_models_propagates_refine_flag():
    fits = fit_lpm_models(
        {"tritium_TU": 2.0, "tritium_sigma_TU": 0.2},
        sample_year=2020.0,
        models=["PFM", "EM"],
        max_age_years=100.0,
        age_steps=20,
        refine=True,
    )
    assert all(fit.refinement_attempted for fit in fits if fit.converged)


def test_final_objective_matches_unified_robust_and_fraction_loss():
    observations = {
        "sample_year": 2020.0,
        "tritium_TU": 2.0,
        "tritium_sigma_TU": 0.2,
        "sf6_pptv": 100.0,
        "sf6_sigma_pptv": 0.2,
    }
    age_fraction_obs = {
        "anthropocene": 0.8,
        "holocene": 0.2,
        "pleistocene": 0.0,
        "sigma_fraction": 0.1,
    }
    fit = fit_lpm_model(
        observations,
        sample_year=2020.0,
        model="EM",
        max_age_years=100.0,
        age_steps=20,
        refine=True,
        age_fraction_obs=age_fraction_obs,
    )
    fit_observations = build_lpm_tracer_observations(observations)
    expected, _, valid = _evaluate_lpm_candidate(
        "EM",
        fit.parameters,
        2020.0,
        [obs.tracer for obs in fit_observations],
        fit_observations,
        max_age_years=100.0,
        age_fraction_obs=age_fraction_obs,
    )
    assert valid
    assert math.isclose(fit.objective, expected, rel_tol=1.0e-10, abs_tol=1.0e-10)


def test_site_specific_he4_rate_changes_fitted_age():
    common = {
        "he4_ccpg": 1.0e-7,
        "he4_sigma_ccpg": 1.0e-9,
        "he4_background_ccpg": 0.0,
    }
    fast = fit_lpm_model(
        {**common, "he4_accumulation_rate_ccpg_per_year": 1.0e-9},
        sample_year=2020.0,
        model="PFM",
        max_age_years=20000.0,
        age_steps=30,
        use_helium4=True,
    )
    slow = fit_lpm_model(
        {**common, "he4_accumulation_rate_ccpg_per_year": 1.0e-11},
        sample_year=2020.0,
        model="PFM",
        max_age_years=20000.0,
        age_steps=30,
        use_helium4=True,
    )
    assert fast.converged and slow.converged
    assert slow.parameters["mean_age_years"] > fast.parameters["mean_age_years"] * 20.0
