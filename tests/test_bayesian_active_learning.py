"""Scientific and numerical contracts for Bayesian active learning."""

from __future__ import annotations

import math

import numpy as np
import pytest

from hydrosheaf.calibration.bayesian_active_learning import (
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    expected_brier_risk_reduction,
    expected_information_gain,
    rank_measurement_options,
    select_measurement_batch,
    update_hypothesis_posterior,
)


def _option(
    option_id: str,
    means: list[float],
    *,
    standard_deviation: float = 1.0,
    cost: float = 1.0,
    second_means: list[float] | None = None,
    feasible: bool = True,
) -> MeasurementOption:
    scenarios = [
        PredictiveScenario("nominal", means, standard_deviation, weight=0.6)
    ]
    if second_means is not None:
        scenarios.append(
            PredictiveScenario(
                "stress", second_means, standard_deviation, weight=0.4
            )
        )
    return MeasurementOption(
        option_id=option_id,
        measurement_type="test",
        target_id=option_id,
        cost=cost,
        scenarios=tuple(scenarios),
        feasible=feasible,
    )


def test_eig_is_zero_for_uninformative_measurement():
    eig = expected_information_gain([0.5, 0.5], [2.0, 2.0], 0.3)
    assert eig == pytest.approx(0.0, abs=1.0e-12)


def test_eig_increases_with_signal_to_noise():
    weak = expected_information_gain([0.5, 0.5], [0.0, 0.5], 1.0)
    strong = expected_information_gain([0.5, 0.5], [0.0, 3.0], 1.0)
    assert 0.0 < weak < strong < math.log(2.0)


def test_eig_is_invariant_to_affine_unit_conversion():
    original = expected_information_gain([0.3, 0.7], [1.0, 4.0], [0.5, 0.8])
    converted = expected_information_gain(
        [0.3, 0.7],
        [12.0, 18.0],
        [1.0, 1.6],
    )
    assert converted == pytest.approx(original, rel=1.0e-12, abs=1.0e-12)


def test_eig_is_invariant_to_hypothesis_permutation():
    original = expected_information_gain(
        [0.2, 0.3, 0.5], [0.0, 1.0, 3.0], [0.5, 0.8, 1.0]
    )
    permuted = expected_information_gain(
        [0.5, 0.2, 0.3], [3.0, 0.0, 1.0], [1.0, 0.5, 0.8]
    )
    assert permuted == pytest.approx(original, rel=1.0e-12, abs=1.0e-12)


def test_eig_exactly_collapses_predictively_equivalent_hypotheses():
    expanded = expected_information_gain(
        [0.1, 0.2, 0.3, 0.4],
        [0.0, 0.0, 2.0, 2.0],
        [0.5, 0.5, 0.8, 0.8],
    )
    collapsed = expected_information_gain([0.3, 0.7], [0.0, 2.0], [0.5, 0.8])
    assert expanded == pytest.approx(collapsed, rel=1.0e-12, abs=1.0e-12)


def test_expected_brier_risk_reduction_tracks_decision_relevance():
    values = np.asarray([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    first_edge = expected_brier_risk_reduction(
        [0.25] * 4,
        [0.0, 0.0, 3.0, 3.0],
        0.4,
        values,
    )
    uninformative = expected_brier_risk_reduction(
        [0.25] * 4,
        [1.0] * 4,
        0.4,
        values,
    )
    assert first_edge > 0.10
    assert uninformative == pytest.approx(0.0, abs=1.0e-12)


def test_decision_focused_ranking_requires_and_uses_decision_values():
    options = [
        _option("edge_a", [0.0, 0.0, 3.0, 3.0], standard_deviation=0.4),
        _option("irrelevant_split", [0.0, 3.0, 0.0, 3.0], standard_deviation=0.4),
    ]
    config = AcquisitionConfig(decision_weight=1.0)
    with pytest.raises(ValueError, match="decision_values"):
        rank_measurement_options(["h0", "h1", "h2", "h3"], [0.25] * 4, options, config=config)
    result = rank_measurement_options(
        ["h0", "h1", "h2", "h3"],
        [0.25] * 4,
        options,
        decision_values=[[0.0], [0.0], [1.0], [1.0]],
        config=config,
    )
    assert result["selected_option_id"] == "edge_a"


def test_posterior_update_moves_towards_supported_hypothesis():
    option = _option("age", [0.0, 4.0], standard_deviation=0.5)
    result = update_hypothesis_posterior(
        ["absent", "present"], [0.5, 0.5], option, 3.8
    )
    assert result["posterior_probabilities"][1] > 0.999
    assert result["posterior_entropy"] < result["prior_entropy"]
    assert result["realised_information_gain"] > 0.0


def test_robust_score_penalises_scenario_fragility():
    fragile = _option(
        "fragile",
        [0.0, 5.0],
        second_means=[0.0, 0.0],
        standard_deviation=0.8,
    )
    stable = _option(
        "stable",
        [0.0, 2.5],
        second_means=[0.0, 2.0],
        standard_deviation=0.8,
    )
    result = rank_measurement_options(
        ["h0", "h1"],
        [0.5, 0.5],
        [fragile, stable],
        config=AcquisitionConfig(robustness_weight=0.9),
    )
    assert result["selected_option_id"] == "stable"
    assert result["rankings"][0]["worst_case_information_gain"] > 0.0


def test_cost_adjustment_prefers_information_efficiency():
    costly = _option("costly", [0.0, 4.0], cost=10.0)
    efficient = _option("efficient", [0.0, 2.0], cost=1.0)
    result = rank_measurement_options(
        ["h0", "h1"], [0.5, 0.5], [costly, efficient]
    )
    assert result["selected_option_id"] == "efficient"


def test_infeasible_options_are_audited_not_ranked():
    result = rank_measurement_options(
        ["h0", "h1"],
        [0.5, 0.5],
        [_option("blocked", [0.0, 5.0], feasible=False)],
    )
    assert result["status"] == "ABSTAIN"
    assert result["rankings"] == []
    assert result["excluded_options"] == [
        {"option_id": "blocked", "reason": "infeasible"}
    ]


def test_abstains_when_gain_is_below_scientific_threshold():
    result = rank_measurement_options(
        ["h0", "h1"],
        [0.5, 0.5],
        [_option("weak", [0.0, 0.01], standard_deviation=10.0)],
        config=AcquisitionConfig(minimum_expected_information_gain=0.01),
    )
    assert result["status"] == "ABSTAIN"
    assert result["selected_option_id"] is None


def test_ties_are_deterministic_by_option_id():
    result = rank_measurement_options(
        ["h0", "h1"],
        [0.5, 0.5],
        [_option("z", [0.0, 2.0]), _option("a", [0.0, 2.0])],
    )
    assert [row["option_id"] for row in result["rankings"]] == ["a", "z"]


def test_batch_selection_avoids_redundant_duplicate():
    # A and B measure the same binary split; C distinguishes the remaining
    # within-pair uncertainty. Joint EIG should choose A/C, not A/B.
    options = [
        _option("a", [0.0, 0.0, 4.0, 4.0], standard_deviation=0.4),
        _option("b_duplicate", [0.0, 0.0, 4.0, 4.0], standard_deviation=0.4),
        _option("c_complement", [0.0, 4.0, 0.0, 4.0], standard_deviation=0.4),
    ]
    result = select_measurement_batch(
        ["h0", "h1", "h2", "h3"],
        [0.25] * 4,
        options,
        batch_size=2,
        config=AcquisitionConfig(batch_qmc_samples=4096, random_seed=17),
    )
    assert result["selected_option_ids"] == ["a", "c_complement"]
    assert result["joint_robust_information_gain"] > 1.2


def test_batch_selection_is_reproducible():
    options = [
        _option("a", [0.0, 0.0, 2.0, 2.0]),
        _option("b", [0.0, 2.0, 0.0, 2.0]),
        _option("c", [0.0, 1.0, 2.0, 3.0]),
    ]
    kwargs = {
        "hypothesis_ids": ["h0", "h1", "h2", "h3"],
        "prior_probabilities": [0.25] * 4,
        "options": options,
        "batch_size": 2,
        "config": AcquisitionConfig(batch_qmc_samples=1024, random_seed=99),
    }
    assert select_measurement_batch(**kwargs) == select_measurement_batch(**kwargs)


@pytest.mark.parametrize(
    "probabilities",
    ([0.0, 0.0], [-0.1, 1.1], [np.nan, 1.0]),
)
def test_invalid_priors_fail_loudly(probabilities):
    with pytest.raises(ValueError):
        rank_measurement_options(
            ["h0", "h1"], probabilities, [_option("x", [0.0, 1.0])]
        )


def test_public_lazy_exports_are_available():
    import hydrosheaf

    assert hydrosheaf.MeasurementOption is MeasurementOption
    assert hydrosheaf.rank_measurement_options is rank_measurement_options
