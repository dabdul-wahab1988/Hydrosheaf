import json
import math

import pytest

from hydrosheaf.validation.decision_utility import (
    CandidateMeasurementAction,
    ProspectiveMeasurementCase,
    ProspectivePolicy,
    ScenarioBelief,
    evaluate_prospective_policies,
    entropy_bits,
    select_random_measurement,
    select_declared_utility_measurement,
    select_next_measurement,
    select_specialist_measurement,
)


PRIOR = {"connected": 0.5, "not_connected": 0.5}


def _posterior_action(action_id, *, cost, posterior, feasible=True):
    return CandidateMeasurementAction(
        action_id=action_id,
        cost=cost,
        feasible=feasible,
        prior=PRIOR,
        posterior=posterior,
    )


def test_action_ranking_prefers_higher_cost_adjusted_information_gain():
    weak = _posterior_action(
        "weak_screen",
        cost=1.0,
        posterior={"connected": 0.75, "not_connected": 0.25},
    )
    strong = CandidateMeasurementAction(
        action_id="tracer_pair",
        cost=1.0,
        prior=PRIOR,
        outcome_likelihoods={
            "young": {"connected": 0.9, "not_connected": 0.1},
            "old": {"connected": 0.1, "not_connected": 0.9},
        },
    )

    report = select_next_measurement([weak, strong], min_utility_per_cost=0.0)

    assert report.decision == "MEASURE"
    assert report.selected_action_id == "tracer_pair"
    ranks = {record.action_id: record.rank for record in report.audit_records}
    assert ranks["tracer_pair"] == 1
    assert ranks["weak_screen"] == 2


def test_cost_penalty_changes_ranking_for_same_information_gain():
    cheap = _posterior_action(
        "cheap_field_ec",
        cost=2.0,
        posterior={"connected": 0.9, "not_connected": 0.1},
    )
    expensive = _posterior_action(
        "expensive_isotope_panel",
        cost=8.0,
        posterior={"connected": 0.9, "not_connected": 0.1},
    )

    report = select_next_measurement([expensive, cheap], min_utility_per_cost=0.0)

    assert report.selected_action_id == "cheap_field_ec"
    utilities = {
        record.action_id: record.robust_utility_per_cost
        for record in report.audit_records
    }
    assert utilities["cheap_field_ec"] == pytest.approx(
        4.0 * utilities["expensive_isotope_panel"]
    )


def test_declared_utility_selector_uses_only_supplied_truth_blind_scores():
    cheap = _posterior_action(
        "cheap_head",
        cost=1.0,
        posterior={"connected": 0.6, "not_connected": 0.4},
    )
    expensive = _posterior_action(
        "expensive_tracer",
        cost=2.0,
        posterior={"connected": 0.9, "not_connected": 0.1},
    )

    report = select_declared_utility_measurement(
        [cheap, expensive],
        {"cheap_head": 0.70, "expensive_tracer": 0.20},
        min_utility_per_cost=-1.0e308,
    )

    assert report.truth_blind is True
    assert report.selected_action_id == "cheap_head"
    audited = {record.action_id: record for record in report.audit_records}
    assert audited["cheap_head"].identifiability == "DECLARED_EXPECTED_UTILITY"
    assert audited["cheap_head"].robust_utility_per_cost == pytest.approx(0.70)


def test_scenario_disagreement_uses_worst_case_robust_utility():
    scenario_sensitive = CandidateMeasurementAction(
        action_id="mechanism_specific_tracer",
        cost=1.0,
        scenarios={
            "flow_prior_correct": ScenarioBelief(
                prior=PRIOR,
                posterior={"connected": 0.99, "not_connected": 0.01},
            ),
            "flow_prior_wrong": ScenarioBelief(prior=PRIOR, posterior=PRIOR),
        },
    )
    stable = CandidateMeasurementAction(
        action_id="boring_but_stable_measurement",
        cost=1.0,
        scenarios={
            "flow_prior_correct": ScenarioBelief(
                prior=PRIOR,
                posterior={"connected": 0.8, "not_connected": 0.2},
            ),
            "flow_prior_wrong": ScenarioBelief(
                prior=PRIOR,
                posterior={"connected": 0.8, "not_connected": 0.2},
            ),
        },
    )

    report = select_next_measurement(
        [scenario_sensitive, stable],
        min_utility_per_cost=0.0,
    )

    assert report.selected_action_id == "boring_but_stable_measurement"
    records = {record.action_id: record for record in report.audit_records}
    assert records["mechanism_specific_tracer"].mean_utility_per_cost > records[
        "boring_but_stable_measurement"
    ].mean_utility_per_cost
    assert records["mechanism_specific_tracer"].robust_utility_per_cost == 0.0
    assert records["mechanism_specific_tracer"].scenario_disagreement > 0.0


def test_infeasible_actions_are_excluded_but_audited():
    infeasible = _posterior_action(
        "perfect_but_blocked_well",
        cost=1.0,
        posterior={"connected": 1.0, "not_connected": 0.0},
        feasible=False,
    )
    feasible = _posterior_action(
        "available_lab_panel",
        cost=1.0,
        posterior={"connected": 0.8, "not_connected": 0.2},
    )

    report = select_next_measurement([infeasible, feasible], min_utility_per_cost=0.0)

    assert report.selected_action_id == "available_lab_panel"
    blocked = next(
        record for record in report.audit_records
        if record.action_id == "perfect_but_blocked_well"
    )
    assert blocked.status == "INFEASIBLE"
    assert blocked.rank is None
    assert blocked.robust_utility_per_cost is None


def test_no_action_abstention_below_predeclared_threshold_and_serializes():
    action = _posterior_action(
        "tiny_update",
        cost=1.0,
        posterior={"connected": 0.51, "not_connected": 0.49},
    )

    report = select_next_measurement([action], min_utility_per_cost=0.5)

    assert report.decision == "ABSTAIN"
    assert report.selected_action_id is None
    assert report.audit_records[0].selected is False
    assert report.audit_records[0].rank == 1
    payload = json.loads(report.to_json())
    assert payload["decision"] == "ABSTAIN"
    assert payload["audit_records"][0]["scenario_utilities"][0][
        "prior_entropy_bits"
    ] == entropy_bits(PRIOR)


def test_posterior_by_outcome_requires_declared_outcome_probabilities():
    action = CandidateMeasurementAction(
        action_id="bad_declaration",
        cost=1.0,
        prior=PRIOR,
        posterior_by_outcome={
            "positive": {"connected": 0.9, "not_connected": 0.1},
            "negative": {"connected": 0.1, "not_connected": 0.9},
        },
    )

    with pytest.raises(ValueError, match="outcome_probabilities"):
        select_next_measurement([action], min_utility_per_cost=0.0)


def test_entropy_ignores_scale_and_uses_bits():
    assert entropy_bits({"a": 2, "b": 2}) == pytest.approx(1.0)
    assert math.isclose(entropy_bits({"a": 1.0, "b": 0.0}), 0.0)


def test_random_and_specialist_policies_are_stable_and_truth_blind():
    actions = [
        _posterior_action("b", cost=2.0, posterior={"connected": 0.9, "not_connected": 0.1}),
        _posterior_action("a", cost=1.0, posterior={"connected": 0.6, "not_connected": 0.4}),
    ]
    first = select_random_measurement(actions, seed=11)
    second = select_random_measurement(reversed(actions), seed=11)
    assert first.to_dict() == second.to_dict()
    assert first.truth_blind is True
    specialist = select_specialist_measurement(actions, {"a": 0.1, "b": 0.9})
    assert specialist.selected_action_id == "b"


def test_prospective_scores_include_cost_adjusted_utility_regret_and_risk():
    actions = (
        _posterior_action("cheap", cost=1.0, posterior={"connected": 0.6, "not_connected": 0.4}),
        _posterior_action("strong", cost=2.0, posterior={"connected": 0.9, "not_connected": 0.1}),
    )
    case = ProspectiveMeasurementCase(
        case_id="case_1",
        actions=actions,
        true_state="connected",
        benefit_by_action_and_state={
            "cheap": {"connected": 1.0, "not_connected": 0.2},
            "strong": {"connected": 3.0, "not_connected": 0.1},
        },
    )
    benchmark = evaluate_prospective_policies(
        [case],
        [
            ProspectivePolicy("random", lambda candidates: select_random_measurement(candidates, seed=3)),
            ProspectivePolicy(
                "specialist",
                lambda candidates: select_specialist_measurement(
                    candidates, {"cheap": 0.1, "strong": 0.9}
                ),
            ),
        ],
        cost_penalty=0.5,
    )
    assert benchmark.status == "SCORED"
    assert benchmark.claim_status == "ABSTAIN_CALIBRATION_INSUFFICIENT"
    assert benchmark.calibration_sufficient is False
    assert benchmark.policies["specialist"].mean_regret == pytest.approx(0.0)
    assert benchmark.policies["random"].mean_cost_adjusted_utility is not None
    assert "random_vs_specialist" in benchmark.pairwise
    payload = benchmark.to_dict()
    assert payload["cost_penalty"] == pytest.approx(0.5)


def test_uniform_random_expectation_is_seed_stable_and_uses_all_feasible_actions():
    actions = (
        _posterior_action("cheap", cost=1.0, posterior={"connected": 0.6, "not_connected": 0.4}),
        _posterior_action("strong", cost=2.0, posterior={"connected": 0.9, "not_connected": 0.1}),
    )
    case = ProspectiveMeasurementCase(
        case_id="case_1",
        actions=actions,
        true_state="connected",
        benefit_by_action_and_state={
            "cheap": {"connected": 1.0, "not_connected": 0.2},
            "strong": {"connected": 3.0, "not_connected": 0.1},
        },
    )
    benchmark = evaluate_prospective_policies(
        [case],
        [
            ProspectivePolicy(
                "random",
                lambda candidates: select_random_measurement(candidates, seed=999),
                scoring_mode="uniform_action_expectation",
            )
        ],
        cost_penalty=0.5,
        required_policy_ids=(),
        calibration_sufficient=True,
    )
    score = benchmark.policies["random"]
    assert score.evaluation_mode == "uniform_action_expectation"
    assert score.selection_coverage == pytest.approx(1.0)
    assert score.mean_cost_adjusted_utility == pytest.approx(1.25)


def test_prospective_incomplete_outcomes_and_truth_leakage_abstain():
    action = _posterior_action(
        "measure", cost=1.0, posterior={"connected": 0.8, "not_connected": 0.2}
    )
    leaked = CandidateMeasurementAction(
        action_id="leaked",
        cost=1.0,
        prior=PRIOR,
        posterior=PRIOR,
        metadata={"truth": "connected"},
    )
    with pytest.raises(ValueError, match="truth field"):
        select_next_measurement([leaked], min_utility_per_cost=0.0)
    case = ProspectiveMeasurementCase(
        case_id="case_1",
        actions=(action,),
        true_state="connected",
        benefit_by_action_and_state={"measure": {"connected": 1.0}},
    )
    benchmark = evaluate_prospective_policies(
        [case],
        [
            ProspectivePolicy("random", lambda candidates: select_random_measurement(candidates)),
            ProspectivePolicy("specialist", lambda candidates: select_next_measurement(candidates, min_utility_per_cost=10.0)),
        ],
    )
    assert benchmark.status == "ABSTAIN_INCOMPLETE_OUTCOMES"
    assert benchmark.claim_status == "ABSTAIN"
