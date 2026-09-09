from __future__ import annotations

import json
from pathlib import Path

from M6.m6_field_transfer_benchmark.scripts.audit_prospective_campaign_readiness import (
    run,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "M6" / "m6_field_transfer_benchmark" / "configs" / "ghana_prospective_campaign.json"


def test_protocol_declares_all_four_packages_and_abstain_gates() -> None:
    protocol = json.loads(CONFIG.read_text(encoding="utf-8"))
    assert protocol["status"] == "PROTOCOL_ONLY"
    assert {item["package_id"] for item in protocol["current_field_packages"]} == {
        "lower_anayari",
        "northen_ghana",
        "northern_ghana_new",
        "talensi_mining_area",
    }
    assert "independent_direct_adjacency_truth" not in protocol["estimands"]
    assert "direct_edge_gate" in protocol["analysis_gates"]
    assert "user_supplied_geology_dictionary" in protocol["required_measurements"]["reaction_and_geology"]


def test_readiness_audit_keeps_field_transfer_and_truth_gates_separate(tmp_path: Path) -> None:
    report = run(output=tmp_path / "RUN-GHANA-PROSPECTIVE-TEST")
    assert report["module_status"]["four_package_transfer"] == "RUN"
    assert report["module_status"]["geochemical_ratios"] == "RUN"
    assert all(row["reaction_screening"] == "RUN" for row in report["dataset_statuses"])
    assert report["module_status"]["independent_age_accuracy"] == "ABSTAIN"
    assert report["module_status"]["direct_adjacency_accuracy"] == "ABSTAIN"
    assert (tmp_path / "RUN-GHANA-PROSPECTIVE-TEST" / "readiness.json").exists()
