from __future__ import annotations

import json
from pathlib import Path

import pytest

from hydrosheaf.benchmarks.ttd_graph_evidence import (
    BenchmarkRecord,
    assess_claim_readiness,
    assert_observation_view_is_truth_blind,
    score_intervals,
    score_prediction_rmse,
    score_topology,
    write_virtual_benchmark_artifacts,
)


def test_interval_score_penalises_blanket_abstention_and_reports_agreement() -> None:
    score = score_intervals(
        truth=[10.0, 20.0, 30.0],
        lower=[8.0, None, 29.0],
        upper=[12.0, None, 31.0],
        abstained=[False, True, False],
        point_estimates=[10.5, None, 32.0],
        expected_abstention=[False, True, False],
    )

    assert score.n_scored == 2
    assert score.n_abstained == 1
    assert score.coverage_nonabstained == 1.0
    assert score.coverage_including_abstention == pytest.approx(2.0 / 3.0)
    assert score.mean_interval_width == pytest.approx(3.0)
    assert score.selective_mae == pytest.approx(1.25)
    assert score.appropriate_abstention_rate == 1.0


def test_interval_score_rejects_an_interval_attached_to_an_abstention() -> None:
    with pytest.raises(ValueError, match="Abstained targets"):
        score_intervals(
            truth=[1.0],
            lower=[0.0],
            upper=[2.0],
            abstained=[True],
        )


def test_prediction_and_directed_topology_scores_are_strict() -> None:
    assert score_prediction_rmse([1.0, 3.0], [2.0, 3.0]) == pytest.approx(2**-0.5)
    score = score_topology(
        [("A", "B"), ("B", "C")],
        [("A", "B"), ("C", "B")],
    )
    assert score["true_positive"] == 1
    assert score["false_positive"] == 1
    assert score["false_negative"] == 1
    assert score["topology_precision"] == 0.5
    assert score["topology_recall"] == 0.5


def _record(
    *,
    family: str = "analytic_particle_network",
    method: str = "local_baseline",
    scenario: str = "nominal",
    held_out: bool = True,
) -> BenchmarkRecord:
    return BenchmarkRecord(
        case_id=f"{family}:{method}:{scenario}",
        generator_family=family,
        regime="static",
        scenario=scenario,
        method=method,
        held_out=held_out,
        truth_blind=True,
        metrics={"held_out_rmse": 0.2},
    )


def test_claim_readiness_requires_every_declared_generator_control_and_scenario() -> None:
    readiness = assess_claim_readiness(
        [_record()],
        required_generator_families=["analytic_particle_network", "process_oriented_modflow_modpath"],
        required_comparators=["local_baseline", "candidate_graph"],
        required_scenarios=["nominal", "reversed_graph"],
    )
    assert readiness["execution_status"] == "FAIL"
    assert readiness["controlled_synthetic_claim_status"] == "ABSTAIN"
    assert readiness["field_validation_status"] == "DEFERRED"
    assert "generator:process_oriented_modflow_modpath" in readiness["missing_requirements"]
    assert "comparator:candidate_graph" in readiness["missing_requirements"]
    assert "scenario:reversed_graph" in readiness["missing_requirements"]


def test_claim_readiness_can_only_be_ready_for_adjudication_not_a_pass() -> None:
    records = [
        _record(family=family, method=method, scenario=scenario)
        for family in ("analytic_particle_network", "process_oriented_modflow_modpath")
        for method in ("local_baseline", "candidate_graph")
        for scenario in ("nominal", "reversed_graph")
    ]
    readiness = assess_claim_readiness(
        records,
        required_generator_families=["analytic_particle_network", "process_oriented_modflow_modpath"],
        required_comparators=["local_baseline", "candidate_graph"],
        required_scenarios=["nominal", "reversed_graph"],
    )
    assert readiness["execution_status"] == "PASS"
    assert readiness["controlled_synthetic_claim_status"] == "READY_FOR_PREREGISTERED_ADJUDICATION"


def test_truth_blind_view_rejects_truth_fields() -> None:
    with pytest.raises(ValueError, match="truth"):
        assert_observation_view_is_truth_blind([{"node_id": "A", "truth_age": 12.0}])


def test_artifact_writer_separates_truth_and_refuses_unrequested_overwrite(
    tmp_path: Path,
) -> None:
    protocol = tmp_path / "protocol.md"
    protocol.write_text("# Test protocol\n", encoding="utf-8")
    readiness = assess_claim_readiness(
        [_record()],
        required_generator_families=["analytic_particle_network"],
        required_comparators=["local_baseline"],
        required_scenarios=["nominal"],
    )
    paths = write_virtual_benchmark_artifacts(
        tmp_path / "run",
        run_id="TTD-GRAPH-TEST-001",
        protocol_path=protocol,
        config={"seed": 1},
        observations=[{"node_id": "A", "time_days": 0.0, "delta18O": -4.1}],
        records=[_record()],
        truth_for_scoring={"truth_age_days": {"A": 12.0}},
        generator_provenance={"imports_hydrosheaf": False},
        readiness=readiness,
    )
    observation_payload = json.loads(Path(paths["observations"]).read_text(encoding="utf-8"))
    truth_payload = json.loads(Path(paths["truth_scoring_only"]).read_text(encoding="utf-8"))
    manifest = json.loads(Path(paths["manifest"]).read_text(encoding="utf-8"))

    assert observation_payload == [{"delta18O": -4.1, "node_id": "A", "time_days": 0.0}]
    assert truth_payload["truth_age_days"]["A"] == 12.0
    assert manifest["generator_independent"] is True
    assert manifest["field_validation_status"] == "DEFERRED"
    with pytest.raises(FileExistsError):
        write_virtual_benchmark_artifacts(
            tmp_path / "run",
            run_id="TTD-GRAPH-TEST-002",
            protocol_path=protocol,
            config={"seed": 2},
            observations=[{"node_id": "A"}],
            records=[_record()],
            truth_for_scoring={"truth_age_days": {"A": 12.0}},
            generator_provenance={"imports_hydrosheaf": False},
            readiness=readiness,
        )
