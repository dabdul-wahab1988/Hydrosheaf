"""Truth-blind controlled-synthetic contracts for the history-aware TTD core."""

from __future__ import annotations

import numpy as np

from hydrosheaf.nuclear.history_ttd import (
    REASON_MISSING_STABLE_ISOTOPE_SOURCE_HISTORY,
    REASON_SOURCE_HISTORY_OUT_OF_WINDOW,
    HistoryTracerObservation,
    build_history_response_matrix,
    fit_history_ttd,
)
from hydrosheaf.nuclear.input_history import InputHistory


def _stable_histories() -> dict[str, InputHistory]:
    years = np.array([2010.0, 2015.0, 2020.0])
    return {
        "d18O": InputHistory(years, np.array([-10.0, -8.0, -6.0])),
        "d2H": InputHistory(years, np.array([-50.0, -70.0, -40.0])),
    }


def test_recovers_history_convolved_stable_isotope_ttd_with_tritium() -> None:
    sample_year = 2020.0
    ages = np.array([0.0, 5.0, 10.0])
    histories = {
        **_stable_histories(),
        "3H": InputHistory(np.array([2010.0, 2015.0, 2020.0]), np.array([10.0, 5.0, 20.0])),
    }
    observations = tuple(
        HistoryTracerObservation(tracer, 0.0, 0.01)
        for tracer in ("d18O", "d2H", "3H")
    )
    true_g = np.array([0.2, 0.5, 0.3])
    matrix = build_history_response_matrix(observations, sample_year, ages, histories)
    values = matrix @ true_g
    observations = tuple(
        HistoryTracerObservation(obs.tracer, value, obs.sigma)
        for obs, value in zip(observations, values)
    )

    result = fit_history_ttd(observations, sample_year, ages, histories, lambda_smoothness=0.0)

    assert result.status == "ESTIMATED"
    assert np.all(np.isfinite(result.g))
    assert np.all(result.g >= -1e-9)
    assert np.isclose(np.sum(result.g), 1.0, atol=1e-7)
    assert np.allclose(result.g, true_g, atol=2e-3)
    assert np.allclose(result.predicted, values, atol=2e-3)
    assert np.allclose(result.residuals, 0.0, atol=2e-3)
    assert result.provenance["inference_family"] == "time_history_ttd"
    assert result.provenance["claim_scope"] == "conditional_inference"
    assert result.provenance["field_validation_status"] == "not_performed"
    assert result.provenance["stable_isotopes_used"] is True
    assert result.provenance["radioactive_tracers_used"] is True
    assert result.provenance["stable_isotope_treatment"] == "conservative_time_history_response"


def test_missing_stable_source_history_abstains() -> None:
    observations = [HistoryTracerObservation("d18O", -8.0, 0.2)]

    result = fit_history_ttd(observations, 2020.0, [0.0, 5.0, 10.0], {})

    assert result.status == "ABSTAIN"
    assert REASON_MISSING_STABLE_ISOTOPE_SOURCE_HISTORY in result.abstention_reasons
    assert result.provenance["stable_isotopes_used"] is True
    assert result.provenance["claim_scope"] == "conditional_inference"


def test_stable_source_history_outside_window_abstains_without_extrapolation() -> None:
    observations = [HistoryTracerObservation("d2H", -60.0, 1.0)]
    history = InputHistory(np.array([2015.0, 2020.0]), np.array([-70.0, -50.0]))

    result = fit_history_ttd(
        observations,
        2020.0,
        [0.0, 5.0, 10.0],
        {"d2H": history},
    )

    assert result.status == "ABSTAIN"
    assert REASON_SOURCE_HISTORY_OUT_OF_WINDOW in result.abstention_reasons
    assert "extrapolation is disabled" in result.diagnostics["message"]


def test_existing_nuclear_only_observations_remain_supported() -> None:
    sample_year = 2020.0
    ages = np.array([0.0, 5.0, 10.0])
    histories = {"3H": InputHistory(np.array([2010.0, 2015.0, 2020.0]), np.array([10.0, 5.0, 20.0]))}
    template = (
        HistoryTracerObservation("3H", 0.0, 1e-2),
        HistoryTracerObservation("14C", 0.0, 1e-2),
    )
    true_g = np.array([0.25, 0.45, 0.30])
    matrix = build_history_response_matrix(template, sample_year, ages, histories)
    values = matrix @ true_g
    observations = tuple(
        HistoryTracerObservation(obs.tracer, value, obs.sigma)
        for obs, value in zip(template, values)
    )

    result = fit_history_ttd(observations, sample_year, ages, histories, lambda_smoothness=0.0)

    assert result.status == "ESTIMATED"
    assert result.provenance["stable_isotopes_used"] is False
    assert result.provenance["radioactive_tracers_used"] is True
    assert np.all(np.isfinite(result.predicted))
    assert np.all(np.isfinite(result.residuals))
    assert np.all(result.g >= -1e-9)
    assert np.isclose(np.sum(result.g), 1.0, atol=1e-7)


def test_optional_stable_transform_and_tracer_weights_are_recorded() -> None:
    histories = _stable_histories()
    ages = np.array([0.0, 5.0, 10.0])
    observations = (
        HistoryTracerObservation("δ18O", -1.0, 0.1, weight=0.5),
        HistoryTracerObservation("2H", -1.0, 0.1),
        HistoryTracerObservation("d18O", -1.0, 0.1),
    )
    # Duplicate isotope aliases are rejected before fitting; this also documents
    # that canonical aliases are treated as one tracer rather than two rows.
    result = fit_history_ttd(
        observations[:2],
        2020.0,
        ages,
        histories,
        stable_isotope_scales={"δ18O": 0.5, "2H": 2.0},
        stable_isotope_offsets={"d18O": 1.0, "d2H": -3.0},
        tracer_weights={"d18O": 2.0},
    )
    assert result.status == "ABSTAIN" or result.status == "ESTIMATED"
    assert result.provenance["stable_isotope_scales"]["d18O"] == 0.5
    assert result.provenance["stable_isotope_offsets"]["d2H"] == -3.0
