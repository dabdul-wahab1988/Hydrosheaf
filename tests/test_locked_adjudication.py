from __future__ import annotations

from pathlib import Path

from scripts.adjudicate_locked_ensemble import build_adjudication


def _row(
    *,
    condition: str,
    scenario: str = "complete",
    f1: float,
    mae: float,
) -> dict[str, object]:
    return {
        "phase": "locked_test",
        "family": "fixture_family",
        "seed": 1,
        "condition": condition,
        "observation_scenario": scenario,
        "topology": {"selected_f1": f1},
        "baselines": {
            "fixture_specialist": {"topology": {"selected_f1": 0.6}},
        },
        "age": {"point": {"mae": mae}},
        "reaction": {"metrics": {"accuracy": 0.25, "n": 4}},
        "decision_policy": {"decision": "ABSTAIN"},
    }


def test_adjudication_keeps_execution_separate_from_scientific_claim() -> None:
    rows = [
        _row(condition="full_sheaf", f1=0.5, mae=1.0),
        _row(condition="age_permuted", f1=0.5, mae=2.0),
        _row(
            condition="hydraulic_only",
            scenario="no_flow_head_permutation",
            f1=0.2,
            mae=1.0,
        ),
    ]
    report = build_adjudication(
        lock_dir=Path("fixture-lock"),
        manifest={
            "run_id": "fixture-run",
            "status": "PASS",
            "integrated_performance_gate": True,
            "calibration": {"target_coverage": 0.95},
            "artifacts": {},
            "source_hashes": {},
        },
        execution_gate={"status": "PASS"},
        rows=rows,
        calibration_scores={
            "fixture_family:1": {
                "calibrated_age": {"coverage": 1.0, "mean_width": 10.0},
                "locked_selective_risk_curve": [],
            },
            "locked_selective_risk_curve": [
                {"requested_acceptance": 1.0, "interval_coverage": 1.0, "accepted": 1, "total": 1}
            ],
        },
        errors=[],
    )

    case = report["case_adjudications"][0]
    assert report["execution_verdict"] == "PASS"
    assert report["system_level_claim_status"] == "ABSTAIN"
    assert case["topology"]["conditional_status"] == "FAIL"
    assert case["age"]["diagnostic_status"] == "PASS"
    assert case["negative_controls"]["age_permutation_status"] == "PASS"
    assert case["negative_controls"]["head_permutation_status"] == "PASS"
    assert case["case_claim_status"] == "ABSTAIN"


def test_common_universe_topology_component_claim_is_bounded_and_fixed() -> None:
    rows = []
    for family, offset in (
        ("analytic_lattice_v1", 0),
        ("independent_mixing_v1", 10),
        ("modflow_modpath_v3", 20),
    ):
        for replicate in range(2):
            row = _row(
                condition="full_sheaf",
                f1=0.70,
                mae=1.0,
            )
            row.update(
                {
                    "family": family,
                    "seed": offset + replicate + 1,
                    "common_candidate_universe": {
                        "algorithm": "common_all_pairs_v1",
                        "candidate_hash": "a" * 64,
                        "truth_blind": True,
                        "truth_fields_seen": [],
                    },
                    "common_topology_comparison": {
                        "hydrosheaf": {"topology": {"selected_f1": 0.70}},
                        "baselines": {
                            "hydraulic_only_directional_gradient": {
                                "selected_f1": 0.60
                            }
                        },
                    },
                }
            )
            rows.append(row)

    report = build_adjudication(
        lock_dir=Path("fixture-lock"),
        manifest={
            "run_id": "fixture-run",
            "status": "PASS_SCOPED_EXECUTION",
            "integrated_performance_gate": False,
            "calibration": {"target_coverage": 0.95},
            "artifacts": {},
            "source_hashes": {},
        },
        execution_gate={"status": "PASS"},
        rows=rows,
        calibration_scores={
            "locked_selective_risk_curve": [
                {
                    "requested_acceptance": 1.0,
                    "interval_coverage": 1.0,
                    "accepted": 1,
                    "total": 1,
                }
            ]
        },
        errors=[],
    )

    assert report["execution_verdict"] == "PASS"
    assert report["synthetic_component_claim"]["status"] == "PASS"
    assert report["synthetic_integrated_claim"]["status"] == "ABSTAIN"
    assert report["field_validation_claim"]["status"] == "DEFERRED"
