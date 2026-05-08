from hydrosheaf.validation import (
    ClaimRecord,
    EvidenceLevel,
    assess_claim_records,
    interval_coverage,
    regression_metrics,
    topology_metrics,
)


def test_regression_metrics_handles_basic_benchmark_values():
    metrics = regression_metrics([1.0, 2.0, 3.0], [1.1, 1.8, 3.2])

    assert metrics["n"] == 3.0
    assert round(metrics["mae"], 6) == round((0.1 + 0.2 + 0.2) / 3.0, 6)
    assert metrics["rmse"] > 0.0
    assert metrics["nse"] < 1.0


def test_interval_coverage_reports_coverage_and_width():
    metrics = interval_coverage([10, 20, 30], [9, 21, 25], [11, 22, 35])

    assert metrics["n"] == 3.0
    assert metrics["coverage"] == 2 / 3
    assert metrics["mean_width"] > 0.0


def test_topology_metrics_use_directed_edges():
    metrics = topology_metrics(
        [("A", "B"), ("B", "C"), ("C", "D")],
        [("A", "B"), ("C", "B"), ("C", "D")],
    )

    assert metrics["tp"] == 2.0
    assert metrics["fp"] == 1.0
    assert metrics["fn"] == 1.0
    assert round(metrics["f1"], 6) == round(2 / 3, 6)


def test_claim_assessment_flags_overclaiming():
    claims = [
        ClaimRecord(
            claim_id="m3-age-equivalence",
            manuscript="M3",
            claim="Hydrosheaf reproduces full TracerLPM model families.",
            supported_level=EvidenceLevel.SCREENING,
            requested_level=EvidenceLevel.VALIDATED,
            guardrail="Use screening-level age-magnitude agreement unless model-family validation is completed.",
        )
    ]

    report = assess_claim_records(claims)

    assert report[0]["overclaim"] is True
    assert report[0]["supported_level"] == "screening"
