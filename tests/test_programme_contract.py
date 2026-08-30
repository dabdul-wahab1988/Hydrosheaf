from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    DecisionKind,
    IdentifiabilityStatus,
    ProgrammeDecision,
    ProgrammeRun,
    ProgrammeStage,
    StageStatus,
    assert_truth_blind,
)


def test_action_decision_serialises_required_fields() -> None:
    decision = ProgrammeDecision(
        decision_kind=DecisionKind.ACTION,
        identifiability=IdentifiabilityStatus.PARTIALLY_IDENTIFIED,
        reason="This action has robust positive expected utility.",
        measurement="age_tracer",
        target="well-A",
        cost=5.0,
        expected_utility=0.12,
        scenario_count=3,
    )

    payload = decision.to_dict()
    assert payload["decision_kind"] == "ACTION"
    assert payload["identifiability"] == "PARTIALLY_IDENTIFIED"
    assert payload["cost"] == 5.0
    assert payload["scenario_count"] == 3


def test_set_report_and_abstain_enforce_safe_outputs() -> None:
    set_report = ProgrammeDecision(
        decision_kind="SET_REPORT",
        reason="Two reaction families remain compatible.",
        compatible_hypotheses=("carbonate", "silicate"),
    )
    abstain = ProgrammeDecision(
        decision_kind="ABSTAIN",
        identifiability="MODEL_DISAGREEMENT",
        reason="Scenario predictions disagree materially.",
    )

    assert set_report.compatible_hypotheses == ("carbonate", "silicate")
    assert abstain.decision_kind is DecisionKind.ABSTAIN

    with pytest.raises(ValueError, match="compatible_hypotheses"):
        ProgrammeDecision(decision_kind="SET_REPORT", reason="No set supplied.")
    with pytest.raises(ValueError, match="cannot contain"):
        ProgrammeDecision(
            decision_kind="ABSTAIN",
            reason="No robust action.",
            target="well-A",
        )


def test_truth_blind_boundary_rejects_explicit_and_prefixed_truth() -> None:
    assert_truth_blind([{"site_id": "A", "hydraulic_head": 10.0}])

    with pytest.raises(ValueError, match="true_edges"):
        assert_truth_blind([{"site_id": "A", "true_edges": []}])
    with pytest.raises(ValueError, match="truth_age_years"):
        assert_truth_blind([{"site_id": "A", "truth_age_years": 20.0}])
    with pytest.raises(ValueError, match="age_years"):
        assert_truth_blind([{"site_id": "A", "age_years": 20.0}])
    with pytest.raises(ValueError, match="true_edges"):
        assert_truth_blind([{"site_id": "A", "nested": {"true_edges": []}}])


def test_run_requires_independent_generator_and_sealed_truth() -> None:
    stage = ProgrammeStage(
        name="network_fit",
        status=StageStatus.COMPLETED,
        detail="completed without generator truth",
    )
    run = ProgrammeRun(
        run_id="RUN-PROGRAMME-SMOKE-0001",
        generator="independent_synthetic_generator",
        generator_independent=True,
        stages=(stage,),
    )

    payload = run.to_dict()
    assert payload["truth_released_for_scoring"] is False
    assert payload["stages"][0]["status"] == "completed"

    with pytest.raises(ValueError, match="independent generator"):
        ProgrammeRun(
            run_id="RUN-DEPENDENT",
            generator="hydrosheaf_generator",
            generator_independent=False,
        )
    with pytest.raises(ValueError, match="sealed"):
        ProgrammeRun(
            run_id="RUN-OPEN-TRUTH",
            generator="independent_synthetic_generator",
            generator_independent=True,
            truth_released_for_scoring=True,
        )
