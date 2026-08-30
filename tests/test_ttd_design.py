import numpy as np

from hydrosheaf.calibration.bayesian_active_learning import AcquisitionConfig
from hydrosheaf.nuclear.ttd_design import (
    CandidateTracerMeasurement,
    TtdHypothesisEnsemble,
    rank_ttd_measurements,
    select_ttd_measurement_batch,
)


def _ensemble(*, with_probabilities=True):
    return TtdHypothesisEnsemble(
        hypothesis_ids=["modern", "old"],
        age_grid_years=[0.0, 10000.0],
        masses=np.eye(2),
        probabilities=[0.4, 0.6] if with_probabilities else None,
        probability_semantics="posterior" if with_probabilities else None,
    )


def test_design_abstains_without_probability_model():
    result = rank_ttd_measurements(
        _ensemble(with_probabilities=False),
        [CandidateTracerMeasurement("c14", "14C", "well-1", 2020.0, 1.0)],
    )

    assert result["status"] == "ABSTAIN"
    assert result["reason"] == "no_probability_model"
    assert result["rankings"] == []


def test_design_uses_candidate_specific_tracer_sensitivity_and_noise():
    candidates = [
        CandidateTracerMeasurement(
            "precise-c14", "14C", "well-1", 2020.0, standard_deviation=1.0
        ),
        CandidateTracerMeasurement(
            "noisy-c14", "14C", "well-1", 2020.0, standard_deviation=100.0
        ),
    ]
    result = rank_ttd_measurements(_ensemble(), candidates)

    assert result["status"] == "ACTIONABLE"
    assert result["selected_option_id"] == "precise-c14"
    assert result["probability_semantics"] == "posterior"
    assert result["groundwater_forward_operator"] == "tracer_response_kernel"


def test_design_abstains_for_undeclared_probability_semantics():
    ensemble = TtdHypothesisEnsemble(
        hypothesis_ids=["a", "b"],
        age_grid_years=[0.0, 10.0],
        masses=np.eye(2),
        probabilities=[0.5, 0.5],
        probability_semantics="optimizer_witnesses",
    )
    result = rank_ttd_measurements(
        ensemble,
        [CandidateTracerMeasurement("c14", "14C", "well-1", 2020.0, 1.0)],
    )

    assert result["status"] == "ABSTAIN"
    assert result["reason"] == "unsupported_probability_semantics"


def test_batch_selection_is_connected_to_groundwater_forward_operator():
    candidates = [
        CandidateTracerMeasurement("c14", "14C", "well-1", 2020.0, 1.0),
        CandidateTracerMeasurement("ar39", "39Ar", "well-1", 2020.0, 1.0),
    ]
    result = select_ttd_measurement_batch(
        _ensemble(),
        candidates,
        batch_size=1,
        config=AcquisitionConfig(batch_qmc_samples=256, random_seed=17),
    )

    assert len(result["selected_option_ids"]) == 1
    assert result["groundwater_forward_operator"] == "tracer_response_kernel"
