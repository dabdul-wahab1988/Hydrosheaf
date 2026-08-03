from __future__ import annotations

from dataclasses import replace

import numpy as np

from hydrosheaf.validation.kinetics_specialist_benchmark import (
    COMPETENCE_MATCHED_COMPARATOR_STATUS,
    KineticsBenchmarkConfig,
    fit_development_interval_calibrator,
    fit_effective_rate_comparator,
    fit_kinetics_specialist,
    generate_kinetics_dataset,
    run_kinetics_specialist_benchmark,
    score_kinetics_predictions,
)


def _small_config(**overrides: object) -> KineticsBenchmarkConfig:
    values = {
        "seed": 17,
        "cases_per_regime": 2,
        "n_time_points": 8,
        "ka_grid_size": 80,
    }
    values.update(overrides)
    return KineticsBenchmarkConfig(**values)


def test_dataset_is_deterministic_and_truth_blind() -> None:
    config = _small_config()
    first = generate_kinetics_dataset(config)
    second = generate_kinetics_dataset(config)

    assert config.config_hash == second.config_hash
    assert first.as_dict() == second.as_dict()
    public = first.truth_blind_view()[0]
    assert "k_per_day" not in public
    assert "surface_area" not in public
    assert "true_concentration" not in public
    assert public["truth_blind"] is True


def test_specialist_reports_effective_rate_and_structural_confounding() -> None:
    config = _small_config(surface_area_measurement_probability=0.0)
    dataset = generate_kinetics_dataset(config)
    predictions = fit_kinetics_specialist(dataset, config)
    score = score_kinetics_predictions(dataset, predictions, config)

    assert len(predictions) == 6
    assert all(prediction.effective_rate_per_day is not None for prediction in predictions)
    assert all(prediction.k_a_identifiable is False for prediction in predictions)
    assert all(prediction.parameter_abstain is True for prediction in predictions)
    assert score.identified_case_rate == 0.0
    assert score.false_commitment_rate == 0.0
    assert np.isfinite(score.predictive_rmse)


def test_independent_area_measurement_enables_separate_k_and_area_diagnostics() -> None:
    config = _small_config(
        surface_area_measurement_probability=1.0,
        missingness_rate=0.0,
    )
    dataset = generate_kinetics_dataset(config)
    predictions = fit_kinetics_specialist(dataset, config)
    score = score_kinetics_predictions(dataset, predictions, config)

    assert score.identified_case_rate == 1.0
    assert score.k_log_rmse_identified is not None
    assert score.surface_area_log_rmse_identified is not None
    assert all(prediction.k_interval is not None for prediction in predictions)


def test_interval_calibration_uses_development_cases_only() -> None:
    config = _small_config()
    dataset = generate_kinetics_dataset(config)
    predictions = fit_kinetics_specialist(dataset, config)
    calibrator = fit_development_interval_calibrator(dataset, predictions, config)

    assert calibrator.fit_scope == "development_only"
    assert set(calibrator.development_case_ids).isdisjoint(calibrator.locked_case_ids)
    assert set(calibrator.development_case_ids) | set(calibrator.locked_case_ids) == {
        observation.case_id for observation in dataset.observations
    }
    calibrated = calibrator.apply(predictions)
    assert all(
        prediction.diagnostics["interval_calibration_fit_scope"] == "development_only"
        for prediction in calibrated
    )
    locked_score = score_kinetics_predictions(
        dataset, calibrated, config, case_ids=dataset.locked_case_ids
    )
    assert locked_score.n_cases == len(dataset.locked_case_ids)
    assert np.isfinite(locked_score.predictive_rmse)

    altered_truth = tuple(
        replace(
            truth,
            true_front_m=truth.true_front_m + 1000.0,
            true_concentration=np.clip(truth.true_concentration + 0.5, 0.0, 1.0),
        )
        if truth.case_id in dataset.locked_case_ids
        else truth
        for truth in dataset._truth
    )
    altered_dataset = replace(dataset, _truth=altered_truth)
    altered_calibrator = fit_development_interval_calibrator(altered_dataset, predictions, config)
    assert altered_calibrator.as_dict() == calibrator.as_dict()


def test_comparator_matches_effective_rate_problem_but_never_claims_k_or_area() -> None:
    config = _small_config(surface_area_measurement_probability=1.0)
    dataset = generate_kinetics_dataset(config)
    specialist = fit_kinetics_specialist(dataset, config)
    comparator = fit_effective_rate_comparator(dataset, config)

    for specialist_prediction, comparator_prediction in zip(specialist, comparator):
        assert specialist_prediction.effective_rate_per_day == comparator_prediction.effective_rate_per_day
        assert comparator_prediction.k_per_day is None
        assert comparator_prediction.surface_area is None
        assert comparator_prediction.k_a_identifiable is False
        assert comparator_prediction.diagnostics["competence_matched"] is True


def test_full_report_keeps_comparison_diagnostic_and_claim_bounded() -> None:
    report = run_kinetics_specialist_benchmark(_small_config())

    assert report.comparator_status == COMPETENCE_MATCHED_COMPARATOR_STATUS
    assert report.claim_status == "ABSTAIN_COMPONENT_ONLY"
    assert "field validity" in report.claim_boundary
    assert report.as_dict()["truth_released_for_scoring"] is False
    assert report.as_dict()["dataset"]["truth_released_for_scoring"] is False
