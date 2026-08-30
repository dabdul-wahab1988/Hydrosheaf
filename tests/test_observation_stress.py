from __future__ import annotations

from hydrosheaf.validation import (
    OBSERVATION_STRESS_SCENARIOS,
    apply_observation_stress,
    consolidate_head_evidence,
)


def test_head_gls_consolidation_does_not_duplicate_the_secondary_channel() -> None:
    rows, audit = consolidate_head_evidence(
        [{"site_id": "A", "head_meas": 10.0, "hydraulic_head": 12.0}],
        model={
            "primary_sigma_m": 1.0,
            "secondary_sigma_m": 1.0,
            "measurement_error_correlation": 0.0,
            "combination": "gls",
        },
    )

    assert rows[0]["head_meas"] == 11.0
    assert "hydraulic_head" not in rows[0]
    assert rows[0]["head_evidence_source_count"] == 2
    assert audit["secondary_channel_consumed_as_independent_evidence"] is False


def test_structured_bias_is_recorded_as_discrepancy_not_independent_evidence() -> None:
    rows, audit = consolidate_head_evidence(
        [{"site_id": "A", "head_meas": 10.0, "hydraulic_head": 12.0}],
        model={
            "primary_sigma_m": 0.1,
            "secondary_sigma_m": 0.2,
            "combination": "primary_only_with_discrepancy",
        },
    )

    assert rows[0]["head_meas"] == 10.0
    assert rows[0]["head_evidence_secondary_minus_primary_m"] == 2.0
    assert "hydraulic_head" not in rows[0]
    assert audit["n_discrepancy_rows"] == 1


def test_observation_stress_is_reproducible_and_contains_real_censoring() -> None:
    rows = [
        {
            "site_id": f"S{i}",
            "tritium_TU": float(i + 1),
            "sf6_pptv": float(i + 1) / 10.0,
            "cfc11_pptv": float(i + 1),
            "cfc113_pptv": float(i + 1),
            "NO3": float(i + 1) / 10.0,
            "Fe": float(i + 1) / 100.0,
        }
        for i in range(8)
    ]

    first = apply_observation_stress(rows, "combined_stress", seed=44)
    second = apply_observation_stress(rows, "combined_stress", seed=44)

    assert first.observations == second.observations
    assert first.missing_count > 0
    assert first.censored_count > 0
    assert set(first.field_counts)
    assert set(OBSERVATION_STRESS_SCENARIOS) == {
        "complete",
        "structured_missingness",
        "left_censoring",
        "combined_stress",
    }


def test_observation_stress_rows_do_not_add_truth_fields() -> None:
    rows = [{"site_id": "S1", "tritium_TU": 1.0}]
    stressed = apply_observation_stress(rows, "structured_missingness", seed=2)

    assert all(
        not str(key).lower().startswith(("true_", "truth_"))
        for row in stressed.observations
        for key in row
    )
