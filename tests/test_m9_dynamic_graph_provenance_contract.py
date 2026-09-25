"""Focused contracts for the multi-tracer dynamic-graph result metadata."""

from __future__ import annotations

from hydrosheaf.benchmarks.ttd_graph_multitracer import (
    MultiTracerBenchmarkConfig,
    generate_multitracer_case,
)
from hydrosheaf.nuclear.multi_tracer_graph_inversion import (
    MultiTracerGraphConfig,
    TRACER_CONFLICT_REASON_CODE,
    solve_joint_multitracer_node_inversion,
)


def _solve_case(panel: str):
    _, observations = generate_multitracer_case(
        MultiTracerBenchmarkConfig(panel=panel, seed=42)
    )
    train = observations.training_end_step
    calibration = {
        tracer: observations.node_tracer_observations["B"][tracer][:train]
        for tracer in observations.available_tracers
    }
    return solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=observations.time_steps[:train],
        tracer_observations=calibration,
        candidate_parents={"A": observations.node_tracer_observations["A"]},
        local_tracer_inputs=observations.local_tracer_inputs.get("B"),
        config=MultiTracerGraphConfig(),
    )


def test_recovered_result_has_controlled_synthetic_provenance_and_active_tracers():
    result = _solve_case("all")

    assert result.status == "RECOVERED"
    assert result.reason_code is None
    assert result.provenance["inference_family"] == "multi_tracer_dynamic_graph"
    assert result.provenance["claim_scope"] == "controlled_synthetic"
    assert result.provenance["field_validation_status"] == "not_performed"
    assert result.provenance["active_tracers"] == list(result.active_tracers)
    assert result.provenance["active_tracer_count"] == len(result.active_tracers)
    assert result.provenance["stable_isotopes_used"] is True
    assert result.provenance["radioactive_tracers_used"] is True
    assert result.provenance["holdout_evaluated"] is False


def test_conflict_abstention_has_stable_reason_and_legacy_human_reason():
    result = _solve_case("conflicting_tracers")

    assert result.status == "ABSTAIN"
    assert result.conflict_detected is True
    assert result.reason_code == TRACER_CONFLICT_REASON_CODE
    assert TRACER_CONFLICT_REASON_CODE in result.reason_codes
    assert "tracer_conflict_detected" in result.reason_codes
    assert any(reason.startswith("conflict:") for reason in result.reason_codes)
    assert result.provenance["claim_scope"] == "controlled_synthetic"
    assert result.provenance["field_validation_status"] == "not_performed"
    assert result.provenance["active_tracers"] == list(result.active_tracers)
    assert result.provenance["active_tracer_count"] == len(result.active_tracers)
