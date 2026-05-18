import math

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
