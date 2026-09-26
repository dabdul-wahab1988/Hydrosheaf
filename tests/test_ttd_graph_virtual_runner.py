from __future__ import annotations

import json
from pathlib import Path

import pytest

from hydrosheaf.benchmarks.ttd_graph_virtual_runner import (
    load_virtual_benchmark_config,
    run_ttd_graph_virtual_benchmark,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs" / "ttd_graph_virtual_benchmark_v1.json"


def _smoke_config(tmp_path: Path) -> Path:
    config = load_virtual_benchmark_config(DEFAULT_CONFIG)
    config["static"] = dict(config["static"])
    config["dynamic"] = dict(config["dynamic"])
    config["particle"] = dict(config.get("particle", {}))
    config["static"]["seeds"] = [17]
    config["static"]["n_time_steps"] = 120
    config["dynamic"]["seeds"] = [23]
    config["dynamic"]["n_time_steps"] = 120
    config["particle"]["seeds"] = [41]
    config["particle"]["n_particles"] = 300
    config["particle"]["n_steps"] = 130
    config["particle"]["max_lag_steps"] = 24
    config["protocol"] = str(ROOT / "docs" / "TTD_GRAPH_VIRTUAL_BENCHMARK_PROTOCOL.md")
    path = tmp_path / "smoke-config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    return path


def test_runner_persists_truth_blind_artifacts_and_evaluates_multi_family_programme(
    tmp_path: Path,
) -> None:
    config_path = _smoke_config(tmp_path)
    result = run_ttd_graph_virtual_benchmark(
        tmp_path / "run", config_path=config_path
    )

    assert result["static_record_count"] == 5
    assert result["dynamic_record_count"] == 8
    assert result["particle_record_count"] == 5
    readiness = result["readiness"]
    assert readiness["static_component"]["execution_status"] == "PASS"
    assert readiness["dynamic_component"]["execution_status"] == "PASS"
    assert readiness["particle_component"]["execution_status"] == "PASS"
    assert readiness["programme"]["execution_status"] == "PASS"
    assert (
        readiness["programme"]["controlled_synthetic_claim_status"]
        == "READY_FOR_PREREGISTERED_ADJUDICATION"
    )
    assert readiness["programme"]["missing_requirements"] == []
    assert readiness["field_transfer"]["status"] == "DEFERRED"
    assert readiness["external_watres"]["status"] == "NOT_RUN"

    artifact_paths = result["artifacts"]
    observations = json.loads(Path(artifact_paths["observations"]).read_text(encoding="utf-8"))
    records = json.loads(Path(artifact_paths["records"]).read_text(encoding="utf-8"))
    manifest = json.loads(Path(artifact_paths["manifest"]).read_text(encoding="utf-8"))
    assert len(observations) == 14
    assert len(records) == 18
    assert manifest["generator_independent"] is True
    assert manifest["field_validation_status"] == "DEFERRED"
    assert all("true_edges" not in row for row in observations)

    with pytest.raises(FileExistsError):
        run_ttd_graph_virtual_benchmark(tmp_path / "run", config_path=config_path)
