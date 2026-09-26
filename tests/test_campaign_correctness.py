"""Regression checks for campaign decisions, identity, and input contracts."""

import json

import numpy as np
import pytest

from hydrosheaf.calibration.bayesian_active_learning import (
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    expected_brier_risk_reduction,
)
from hydrosheaf.calibration.well_active_learning import (
    CampaignConfig,
    HypothesisScenario,
    WellAction,
    build_well_actions,
    build_well_measurement_options,
    evaluate_abstention_gate,
    load_predictive_scenarios_file,
    select_well_campaign_batch,
)
from hydrosheaf.graph.types import Edge


def _option(option_id="A", means=(0.0, 2.0), sd=1.0, cost=2.0):
    return MeasurementOption(
        option_id=option_id,
        measurement_type="head_monitoring",
        target_id=option_id,
        cost=cost,
        scenarios=(PredictiveScenario("nominal", means, sd),),
    )


@pytest.mark.parametrize("key", ["n_edges_r_hat", "n_edges_ess"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_diagnostics_cannot_authorize_campaign(key, value):
    posterior = {"n_edges_r_hat": 1.01, "n_edges_ess": 100.0}
    posterior[key] = value
    abstain, reasons = evaluate_abstention_gate(
        posterior, prior_entropy=0.69, max_robust_eig=0.2
    )
    assert abstain, reasons


@pytest.mark.parametrize("key", ["joint_graph_entropy", "acceptance_rate"])
def test_nonfinite_posterior_quality_cannot_authorize_campaign(key):
    posterior = {
        "n_edges_r_hat": 1.01,
        "n_edges_ess": 100.0,
        "joint_graph_entropy": 0.69,
        key: float("nan"),
    }
    abstain, reasons = evaluate_abstention_gate(posterior, max_robust_eig=0.2)
    assert abstain, reasons


def test_fractional_sampling_times_keep_distinct_predictions(tmp_path):
    actions = [
        WellAction("W1", "head_monitoring", sampling_time=t) for t in (0.01, 0.04)
    ]
    assert len({action.action_id for action in actions}) == 2
    scenario_path = tmp_path / "predictions.json"
    scenario_path.write_text(
        json.dumps(
            {
                action.action_id: {"H1": {"mean": i, "sd": 0.5}}
                for i, action in enumerate(actions)
            }
        ),
        encoding="utf-8",
    )
    predictor = load_predictive_scenarios_file(str(scenario_path))
    for i, action in enumerate(actions):
        assert predictor(
            {
                "well_id": "W1",
                "measurement_type": "head_monitoring",
                "sampling_time": action.sampling_time,
                "hypothesis_id": "H1",
            }
        ) == (i, 0.5)


def test_batch_information_is_invariant_to_prior_weight_scale():
    options = [_option()]
    cfg = CampaignConfig(batch_size=1)
    normalized = select_well_campaign_batch(options, [0.9, 0.1], cfg)
    weights = select_well_campaign_batch(options, [9.0, 1.0], cfg)
    assert weights["joint_robust_information_gain"] == pytest.approx(
        normalized["joint_robust_information_gain"], abs=1e-12
    )


def test_abstention_threshold_is_invariant_to_prior_weight_scale():
    kwargs = dict(
        posterior_result={"n_edges_r_hat": 1.01, "n_edges_ess": 100.0},
        candidate_options=[_option()],
        prior_entropy=0.325,
        config=CampaignConfig(minimum_robust_eig=0.15),
    )
    normalized = evaluate_abstention_gate(prior_probabilities=[0.9, 0.1], **kwargs)
    weights = evaluate_abstention_gate(prior_probabilities=[9.0, 1.0], **kwargs)
    assert normalized[0] is True
    assert weights == normalized


def test_batch_abstains_when_predictions_do_not_distinguish_hypotheses():
    result = select_well_campaign_batch([_option(means=(2.0, 2.0))], [0.5, 0.5])
    assert result["status"] == "ABSTAIN"
    assert result["selected_options"] == []
    assert result["total_cost"] == 0.0


def test_batch_uses_declared_cost_without_inventing_travel_charges():
    result = select_well_campaign_batch([_option(cost=2.0)], [0.5, 0.5], budget=2.0)
    assert result["selected_options"] == ["A"]
    assert result["total_cost"] == 2.0


def test_batch_rejects_duplicate_action_ids():
    with pytest.raises(ValueError, match="Duplicate option_id"):
        select_well_campaign_batch([_option(), _option(means=(0.0, 3.0))], [0.5, 0.5])


def test_direct_and_converted_actions_use_the_same_predictive_model():
    action = WellAction(
        "W1",
        "head_monitoring",
        scenarios=(
            HypothesisScenario("H0", 0.0, 1.0),
            HypothesisScenario("H1", 2.0, 1.0),
        ),
    )
    options = build_well_measurement_options([action], ["H0", "H1"], [(), ("E",)], [])
    direct = select_well_campaign_batch([action], [0.5, 0.5])
    converted = select_well_campaign_batch(options, [0.5, 0.5])
    assert direct["joint_robust_information_gain"] == pytest.approx(
        converted["joint_robust_information_gain"], abs=1e-12
    )


def test_action_builder_accepts_samples_indexed_by_well():
    kwargs = dict(
        candidate_edges=[Edge("E", "W1", "W2")],
        hypothesis_ids=["H0", "H1"],
        edge_tuples=[(), ("E",)],
        config=CampaignConfig(measurement_types=("head_monitoring",)),
    )
    rows = [{"well_id": "W1"}, {"well_id": "W2"}]
    indexed = build_well_actions({row["well_id"]: row for row in rows}, **kwargs)
    sequential = build_well_actions(rows, **kwargs)
    assert indexed == sequential


def test_batch_decision_gain_is_conditional_on_previously_selected_measurements():
    # Four equally likely hypotheses encode two independent binary targets.
    # A resolves the first target. B repeats A; C resolves the second target.
    options = [
        _option("A", (0.0, 0.0, 20.0, 20.0), sd=0.1),
        _option("B", (0.0, 0.0, 20.0, 20.0), sd=0.1),
        _option("C", (0.0, 20.0, 0.0, 20.0), sd=0.1),
    ]
    cfg = CampaignConfig(
        batch_size=2,
        acquisition_config=AcquisitionConfig(decision_weight=1.0),
    )
    decisions = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)
    result = select_well_campaign_batch(
        options, [0.25] * 4, cfg, decision_values=decisions
    )
    assert result["selected_options"] == ["A", "C"]
    gains = [
        r["marginal_decision_risk_reduction"] for r in result["selected_option_details"]
    ]
    assert sum(gains) == pytest.approx(0.25, abs=1e-8)


@pytest.mark.parametrize("prior", [[0.5, 0.5], [0.9, 0.1]])
def test_joint_decision_risk_matches_gaussian_sufficient_statistic(prior):
    # Two independent equal-variance observations have a sample mean with
    # SD sigma/sqrt(2). Compare QMC batch risk against one-dimensional quadrature.
    options = [_option("A"), _option("B")]
    cfg = CampaignConfig(
        batch_size=2,
        acquisition_config=AcquisitionConfig(
            decision_weight=1.0, batch_qmc_samples=8192
        ),
    )
    result = select_well_campaign_batch(
        options, prior, cfg, decision_values=[[0.0], [1.0]]
    )
    rows = result["selected_option_details"]
    assert len(rows) == 2
    for n, row in enumerate(rows, start=1):
        reference = expected_brier_risk_reduction(
            prior, [0.0, 2.0], 1.0 / np.sqrt(n), [[0.0], [1.0]], quadrature_order=81
        )
        assert row["joint_decision_risk_reduction"] == pytest.approx(
            reference, abs=5e-4
        )
    assert (
        rows[1]["marginal_decision_risk_reduction"]
        < rows[0]["marginal_decision_risk_reduction"]
    )
