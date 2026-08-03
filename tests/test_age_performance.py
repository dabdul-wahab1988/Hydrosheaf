from __future__ import annotations

from hydrosheaf.validation.age_competent_baseline import AgeCompetentBaseline
from hydrosheaf.validation.age_performance import (
    AgePerformanceThresholds,
    evaluate_age_performance,
)


def _raw_gas_observations() -> dict[str, object]:
    return {
        "tracer_age_history": {
            "nodes": {
                "site-a": {
                    "tritium_TU": 2.6,
                    "argon39_pmc": 95.0,
                    "sf6_dissolved": 20.0,
                    "sf6_solubility_per_pptv": 2.0,
                    "sf6_excess_air_fraction": 0.25,
                    "sf6_dissolved_sigma": 1.0,
                    "sample_date": 2025.5,
                }
            }
        }
    }


def test_age_baseline_profiles_declared_excess_air_schema_without_truth() -> None:
    prediction = AgeCompetentBaseline().predict(_raw_gas_observations())["site-a"]

    assert prediction["correction_model"] == "profiled_excess_air_and_degradation"
    assert prediction["correction_scenario_count"] == 2
    assert prediction["correction_support"]["best_excess_air_fraction"] == 0.25
    assert prediction["truth_blind"] is True
    assert "true_age" not in prediction


def test_age_baseline_profiles_declared_degradation_alternative() -> None:
    observations = _raw_gas_observations()
    observations["tracer_age_history"]["nodes"]["site-a"]["sf6_degradation_fraction"] = 0.5

    prediction = AgeCompetentBaseline().predict(observations)["site-a"]

    assert prediction["correction_scenario_count"] == 4
    assert prediction["correction_support"]["best_degradation_fraction_by_tracer"]["sf6_pptv"] in {0.0, 0.5}


def test_age_performance_reports_conditional_and_unconditional_coverage() -> None:
    truth = {f"t{index}": float(10 + index * 5) for index in range(6)}
    predictions = {
        target: {
            "calibrated_mean_age_years": age + (0.5 if index % 2 else 0.0),
            "calibrated_age_low": age - 2.0,
            "calibrated_age_high": age + 2.0,
            "uncertainty_years": 2.0,
            "decision": "abstain" if index == 0 else "select",
            "calibration_status": "abstained_not_calibrated" if index == 0 else "development_fitted",
            "mean_age_years": age,
            "age_95_low": age - 2.0,
            "age_95_high": age + 2.0,
        }
        for index, (target, age) in enumerate(truth.items())
    }
    metadata = {
        target: {
            "case_id": "test-a" if index < 3 else "test-b",
            "generator_family": "analytic" if index < 3 else "mixing",
            "mechanism": "single_lpm" if index < 3 else "mixed_lpm",
            "missingness": "complete" if index != 5 else "one_tracer_missing",
        }
        for index, target in enumerate(truth)
    }
    thresholds = AgePerformanceThresholds(
        min_conditional_coverage=0.90,
        min_coverage_including_abstention=0.60,
        max_mae_years=5.0,
        max_mean_interval_width_years=10.0,
        max_selective_risk_years=5.0,
        max_false_commitment_rate=0.10,
        min_selected_rate=0.50,
        min_rows=2,
        min_rows_per_stratum=1,
    )

    report = evaluate_age_performance(
        truth,
        predictions,
        metadata=metadata,
        thresholds=thresholds,
        development_case_ids=("development-a",),
        held_out_case_ids=("test-a", "test-b"),
    )

    assert report["claim_status"] == "PASS"
    metrics = report["metrics"]
    assert metrics["n_truth"] == 6
    assert metrics["n_selected"] == 5
    assert metrics["n_abstain"] == 1
    assert metrics["conditional_interval_coverage"] == 1.0
    assert metrics["coverage_including_abstention"] == 5 / 6
    assert metrics["false_commitment_rate"] == 0.0
    assert set(report["stratified"]) == {"family", "mechanism", "missingness"}
    assert report["selective_risk_curve"]


def test_age_performance_does_not_reward_all_abstention() -> None:
    truth = {"a": 10.0, "b": 20.0, "c": 30.0, "d": 40.0}
    predictions = {
        target: {
            "mean_age_years": 60.0,
            "age_95_low": 0.0,
            "age_95_high": 120.0,
            "uncertainty_years": 60.0,
            "decision": "abstain",
        }
        for target in truth
    }
    metadata = {
        target: {
            "case_id": "test",
            "family": "analytic" if index < 2 else "mixing",
            "mechanism": "single" if index < 2 else "mixed",
            "missingness": "complete",
        }
        for index, target in enumerate(truth)
    }
    report = evaluate_age_performance(
        truth,
        predictions,
        metadata=metadata,
        use_calibrated=False,
        thresholds=AgePerformanceThresholds(min_rows=2, min_rows_per_stratum=1),
        development_case_ids=("development",),
        held_out_case_ids=("test",),
    )

    assert report["metrics"]["n_selected"] == 0
    assert report["metrics"]["coverage_including_abstention"] == 0.0
    assert any(gate["name"] == "nontrivial_selection" and gate["status"] == "FAIL" for gate in report["gates"])
    assert report["claim_status"] == "FAIL"


def test_age_performance_rejects_overlapping_development_and_test_cases() -> None:
    truth = {"a": 10.0, "b": 20.0}
    predictions = {
        target: {
            "mean_age_years": age,
            "age_95_low": age - 1.0,
            "age_95_high": age + 1.0,
            "uncertainty_years": 1.0,
            "decision": "select",
        }
        for target, age in truth.items()
    }
    metadata = {
        target: {
            "case_id": "case-overlap",
            "family": "analytic",
            "mechanism": "single",
            "missingness": "complete",
        }
        for target in truth
    }
    report = evaluate_age_performance(
        truth,
        predictions,
        metadata=metadata,
        use_calibrated=False,
        thresholds=AgePerformanceThresholds(min_rows=1, min_rows_per_stratum=1),
        development_case_ids=("case-overlap",),
        held_out_case_ids=("case-overlap",),
    )

    gate = next(gate for gate in report["gates"] if gate["name"] == "held_out_case_disjoint")
    assert gate["status"] == "FAIL"
    assert report["claim_status"] == "FAIL"
