"""Contract tests for the isolated two-tier age-adjacency QA audit."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "scripts"
    / "audit_age_adjacency_tiers.py"
)
SPEC = importlib.util.spec_from_file_location("audit_age_adjacency_tiers", SCRIPT)
assert SPEC and SPEC.loader
qa = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qa)


HASH = "0" * 64


def _provenance(source_ids: list[str]) -> dict[str, object]:
    return {
        "created_utc": "2026-09-08T00:00:00+00:00",
        "protocol_hash": HASH,
        "generator": {"name": "independent-generator", "version": "1.0", "revision": "gen-001"},
        "inference": {"name": "hydrosheaf", "version": "dev", "revision": "inf-001"},
        "environment": {"runtime": "python-3.14", "dependency_lock_hash": HASH},
        "inputs": [
            {"source_id": source_id, "role": "observations", "path": f"{source_id}.json", "sha256": HASH, "size_bytes": 1}
            for source_id in source_ids
        ],
    }


def _metadata(source_kind: str, claim_tier: str, *, aiken: bool = False, integrated: bool = False) -> dict[str, object]:
    return {
        "source_kind": source_kind,
        "truth_sealed": True,
        "truth_access_mode": "evaluation_only",
        "candidate_set_frozen": True,
        "integrated_scoring_allowed": integrated,
        "aiken_emulation": aiken,
        "prediction_input_columns": ["age_years", "age_sigma_years"],
        "evaluation_only_columns": ["relation_label", "is_true_edge"],
    }


def _ages() -> list[dict[str, object]]:
    return [
        {"case_id": "case-1", "node_id": "A", "age_years": 100.0, "age_sigma_years": 2.0, "age_status": "observed", "source_id": "ages"},
        {"case_id": "case-1", "node_id": "B", "age_years": 110.0, "age_sigma_years": 2.0, "age_status": "observed", "source_id": "ages"},
        {"case_id": "case-1", "node_id": "M", "age_years": 105.0, "age_sigma_years": 2.0, "age_status": "observed", "source_id": "ages"},
    ]


def _base_package(*, tier: str, source_kind: str = "independent_synthetic_truth", claim_tier: str = "controlled_synthetic_component", aiken: bool = False, integrated: bool = True) -> dict[str, object]:
    return {
        "schema": qa.PROTOCOL_ID,
        "protocol_id": qa.PROTOCOL_ID,
        "run_id": "RUN-AGE-QA-001",
        "tier": tier,
        "claim_tier": claim_tier,
        "metadata": _metadata(source_kind, claim_tier, aiken=aiken, integrated=integrated),
        "units": {
            "age_years": "years",
            "age_sigma_years": "years",
            "age_covariance_years2": "years^2",
            "travel_time_years": "years",
            "travel_sigma_years": "years",
            "probability": "1",
        },
        "provenance": _provenance(["ages", "transport"]),
        "ages": _ages(),
        "transport_hypotheses": [],
        "predictions": [],
        "truth_artifact": {"available": False, "sealed": True, "role": "unavailable", "reason": "not required for this fixture"},
    }


def _write_package(tmp_path: Path, package: dict[str, object], *, truth: dict[str, object] | None = None) -> Path:
    if truth is not None:
        truth_path = tmp_path / "truth.json"
        truth_path.write_text(json.dumps(truth, sort_keys=True) + "\n", encoding="utf-8")
        package["truth_artifact"] = {
            "available": True,
            "sealed": True,
            "role": "evaluation_only",
            "path": "truth.json",
            "sha256": hashlib.sha256(truth_path.read_bytes()).hexdigest(),
        }
    package_path = tmp_path / "package.json"
    package_path.write_text(json.dumps(package, sort_keys=True) + "\n", encoding="utf-8")
    return package_path


def test_t1_temporal_order_requires_directness_abstention(tmp_path: Path) -> None:
    package = _base_package(
        tier=qa.TIER_TEMPORAL,
        source_kind="observed_field",
        claim_tier="direction_diagnostic",
        integrated=False,
    )
    package["provenance"] = _provenance(["ages"])
    package["metadata"]["prediction_input_columns"] = ["age_years", "age_sigma_years"]
    package["transport_hypotheses"] = []
    package["truth_artifact"] = {"available": False, "sealed": True, "role": "unavailable", "reason": "no independent edge truth"}
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.95,
            "prediction_status": "ABSTAIN",
            "adjacency_status": "ABSTAIN",
            "identifiability_stratum": "endpoint_age_order_only",
        }
    ]

    report = qa.audit_package(_write_package(tmp_path, package))

    assert report["valid"] is True
    assert report["decision"] == "PASS_T1_TEMPORAL_ONLY"
    assert report["summary"]["n_scored_directness"] == 0


def test_t2_requires_independent_direct_and_indirect_hypotheses(tmp_path: Path) -> None:
    package = _base_package(tier=qa.TIER_SEGMENT_TRANSPORT)
    package["transport_hypotheses"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "hypothesis": "direct",
            "travel_time_years": 10.0,
            "travel_sigma_years": 1.0,
            "evidence_source_id": "transport",
            "independent_of_endpoint_age": True,
            "path_basis": "segment_rtd",
            "intermediate_nodes": [],
        },
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "hypothesis": "indirect",
            "travel_time_years": 30.0,
            "travel_sigma_years": 4.0,
            "evidence_source_id": "transport",
            "independent_of_endpoint_age": True,
            "path_basis": "two_segment_convolution",
            "intermediate_nodes": ["M"],
        },
    ]
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.95,
            "prediction_status": "scored",
            "adjacency_status": "scored",
            "identifiability_stratum": "edge_transport_comparison",
            "direct_probability": 0.91,
            "log_bayes_factor_direct_vs_indirect": 2.4,
        }
    ]
    truth = {
        "schema": qa.TRUTH_SCHEMA,
        "run_id": "RUN-AGE-QA-001",
        "sealed": True,
        "direct_edge_count": 1,
        "reachable_ordered_pair_count": 2,
    }

    report = qa.audit_package(_write_package(tmp_path, package, truth=truth))

    assert report["valid"] is True
    assert report["decision"] == "PASS_T2_DIRECTNESS_SCORING"
    assert report["summary"]["truth_direct_edge_count"] == 1


def test_t2_rejects_scored_row_without_indirect_hypothesis(tmp_path: Path) -> None:
    package = _base_package(tier=qa.TIER_SEGMENT_TRANSPORT)
    package["transport_hypotheses"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "hypothesis": "direct",
            "travel_time_years": 10.0,
            "travel_sigma_years": 1.0,
            "evidence_source_id": "transport",
            "independent_of_endpoint_age": True,
            "path_basis": "segment_rtd",
            "intermediate_nodes": [],
        }
    ]
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.95,
            "prediction_status": "scored",
            "adjacency_status": "scored",
            "identifiability_stratum": "edge_transport_comparison",
            "direct_probability": 0.91,
            "log_bayes_factor_direct_vs_indirect": 2.4,
        }
    ]

    report = qa.audit_package(_write_package(tmp_path, package))

    assert report["valid"] is False
    assert any("requires direct and indirect hypotheses" in item for item in report["errors"])


def test_truth_bearing_prediction_field_is_rejected(tmp_path: Path) -> None:
    package = _base_package(tier=qa.TIER_TEMPORAL, source_kind="observed_field", claim_tier="direction_diagnostic", integrated=False)
    package["provenance"] = _provenance(["ages"])
    package["truth_artifact"] = {"available": False, "sealed": True, "role": "unavailable", "reason": "not independent"}
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.9,
            "prediction_status": "ABSTAIN",
            "adjacency_status": "ABSTAIN",
            "identifiability_stratum": "endpoint_age_order_only",
            "relation_label": "direct_adjacent",
        }
    ]

    report = qa.audit_package(_write_package(tmp_path, package))

    assert report["valid"] is False
    assert any("truth-bearing field" in item for item in report["errors"])


def test_aiken_is_emulation_only_and_cannot_score_directness(tmp_path: Path) -> None:
    package = _base_package(
        tier=qa.TIER_SEGMENT_TRANSPORT,
        source_kind="calibrated_model_reference",
        claim_tier="calibrated_model_reference",
        aiken=True,
        integrated=False,
    )
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.8,
            "prediction_status": "ABSTAIN",
            "adjacency_status": "insufficient_information",
            "identifiability_stratum": "transport_censored",
            "flags": ["no_independent_direct_truth"],
        }
    ]
    package["transport_hypotheses"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "hypothesis": "direct",
            "travel_time_years": 10.0,
            "travel_sigma_years": 2.0,
            "evidence_source_id": "transport",
            "independent_of_endpoint_age": False,
            "path_basis": "calibrated_modpath_reference",
            "intermediate_nodes": [],
        },
    ]

    report = qa.audit_package(_write_package(tmp_path, package))

    assert report["valid"] is True
    assert report["decision"] == "PASS_AIKEN_EMULATION_ONLY"
    assert report["summary"]["directness_estimand"] == "ABSTAIN"

    package["predictions"][0]["prediction_status"] = "scored"
    package["predictions"][0]["adjacency_status"] = "scored"
    package["predictions"][0]["identifiability_stratum"] = "edge_transport_comparison"
    package["predictions"][0]["direct_probability"] = 0.9
    package["predictions"][0]["log_bayes_factor_direct_vs_indirect"] = 1.2
    report = qa.audit_package(_write_package(tmp_path, package))
    assert report["valid"] is False
    assert any("Aiken emulation cannot score direct adjacency" in item for item in report["errors"])


def test_truth_sidecar_hash_is_verified(tmp_path: Path) -> None:
    package = _base_package(tier=qa.TIER_SEGMENT_TRANSPORT)
    package["predictions"] = [
        {
            "case_id": "case-1",
            "edge_id": "A->B",
            "u": "A",
            "v": "B",
            "direction_probability": 0.8,
            "prediction_status": "ABSTAIN",
            "adjacency_status": "insufficient_information",
            "identifiability_stratum": "unidentifiable",
        }
    ]
    truth = {
        "schema": qa.TRUTH_SCHEMA,
        "run_id": "RUN-AGE-QA-001",
        "sealed": True,
        "direct_edge_count": 1,
        "reachable_ordered_pair_count": 2,
    }
    package_path = _write_package(tmp_path, package, truth=truth)
    loaded = json.loads(package_path.read_text(encoding="utf-8"))
    loaded["truth_artifact"]["sha256"] = HASH
    package_path.write_text(json.dumps(loaded, sort_keys=True) + "\n", encoding="utf-8")

    report = qa.audit_package(package_path)

    assert report["valid"] is False
    assert any("does not match" in item for item in report["errors"])
