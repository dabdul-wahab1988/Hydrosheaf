"""Tests for the separated temporal-direction/direct-adjacency age evidence."""

from __future__ import annotations

import math

import pytest

from hydrosheaf.validation.age_adjacency import (
    compute_age_adjacency_evidence,
    compute_direction_evidence,
)


def test_reverse_age_order_is_directionally_incompatible() -> None:
    result = compute_direction_evidence(
        upstream_age_years=20.0,
        downstream_age_years=10.0,
        upstream_sigma_years=1.0,
        downstream_sigma_years=1.0,
    )

    assert result.status == "reverse_incompatible"
    assert result.age_increment_years == pytest.approx(-10.0)
    assert result.age_sigma_years == pytest.approx(math.sqrt(2.0))
    assert result.forward_compatibility_probability < 1.0e-10
    assert result.forward_compatibility_cost > 10.0


def test_direct_match_beats_a_cumulative_skip_when_travel_hypotheses_differ() -> None:
    direct_match = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=0.5,
        downstream_sigma_years=0.5,
        process_sigma_years=1.0,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
    )
    cumulative_skip = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=130.0,
        upstream_sigma_years=0.5,
        downstream_sigma_years=0.5,
        process_sigma_years=1.0,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
    )

    assert direct_match.adjacency_status == "scored"
    assert cumulative_skip.adjacency_status == "scored"
    assert direct_match.log_bayes_factor_direct_vs_indirect is not None
    assert cumulative_skip.log_bayes_factor_direct_vs_indirect is not None
    assert direct_match.log_bayes_factor_direct_vs_indirect > 0.0
    assert cumulative_skip.log_bayes_factor_direct_vs_indirect < 0.0
    # Direction remains forward-compatible in both cases; the comparison of
    # explicit travel-time hypotheses is what changes the adjacency evidence.
    assert direct_match.direction.forward_compatibility_probability > 0.5
    assert cumulative_skip.direction.forward_compatibility_probability > 0.5


def test_wider_age_and_process_uncertainty_reduces_discrimination() -> None:
    narrow = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=0.25,
        downstream_sigma_years=0.25,
        process_sigma_years=0.25,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
    )
    wide = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=8.0,
        downstream_sigma_years=8.0,
        process_sigma_years=8.0,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
    )

    assert narrow.propagated_sigma_years == pytest.approx(math.sqrt(0.25**2 * 3.0))
    assert wide.propagated_sigma_years == pytest.approx(math.sqrt(8.0**2 * 3.0))
    assert narrow.log_bayes_factor_direct_vs_indirect is not None
    assert wide.log_bayes_factor_direct_vs_indirect is not None
    assert abs(wide.log_bayes_factor_direct_vs_indirect) < abs(
        narrow.log_bayes_factor_direct_vs_indirect
    )


def test_missing_indirect_hypothesis_is_explicitly_ambiguous() -> None:
    direct_only = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=1.0,
        downstream_sigma_years=1.0,
        process_sigma_years=1.0,
        direct_travel_years=10.0,
    )
    no_path = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=1.0,
        downstream_sigma_years=1.0,
        process_sigma_years=1.0,
    )

    assert direct_only.adjacency_status == "ambiguous"
    assert direct_only.direct_log_likelihood is not None
    assert direct_only.indirect_log_likelihood is None
    assert direct_only.log_bayes_factor_direct_vs_indirect is None
    assert "indirect_travel_time_missing" in direct_only.flags
    assert no_path.adjacency_status == "insufficient_information"
    assert no_path.direct_log_likelihood is None
    assert no_path.log_bayes_factor_direct_vs_indirect is None
    assert "direct_travel_time_missing" in no_path.flags


def test_zero_propagated_uncertainty_cannot_produce_a_normal_likelihood() -> None:
    result = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
    )

    assert result.adjacency_status == "insufficient_information"
    assert result.direct_log_likelihood is None
    assert result.indirect_log_likelihood is None
    assert result.log_bayes_factor_direct_vs_indirect is None
    assert "zero_propagated_uncertainty" in result.flags


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"upstream_age_years": math.nan, "downstream_age_years": 2.0}, "finite"),
        ({"upstream_age_years": 1.0, "downstream_age_years": math.inf}, "finite"),
        (
            {
                "upstream_age_years": 1.0,
                "downstream_age_years": 2.0,
                "upstream_sigma_years": -1.0,
            },
            "non-negative",
        ),
        (
            {
                "upstream_age_years": 1.0,
                "downstream_age_years": 2.0,
                "direct_travel_years": -0.1,
            },
            "non-negative",
        ),
    ],
)
def test_nonfinite_and_negative_inputs_are_rejected(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        compute_age_adjacency_evidence(**kwargs)


def test_equal_hypotheses_are_scored_but_marked_nondiscriminating() -> None:
    result = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=1.0,
        downstream_sigma_years=1.0,
        direct_travel_years=10.0,
        indirect_travel_years=10.0,
    )

    assert result.adjacency_status == "insufficient_information"
    assert result.log_bayes_factor_direct_vs_indirect is None
    assert result.comparison_discriminating is False
    assert "equal_travel_time_hypotheses" in result.flags


def test_hypothesis_specific_dispersion_and_endpoint_covariance_are_propagated() -> None:
    result = compute_age_adjacency_evidence(
        upstream_age_years=100.0,
        downstream_age_years=110.0,
        upstream_sigma_years=2.0,
        downstream_sigma_years=2.0,
        age_covariance_years2=2.0,
        process_sigma_years=1.0,
        direct_travel_years=10.0,
        indirect_travel_years=30.0,
        direct_travel_sigma_years=1.0,
        indirect_travel_sigma_years=5.0,
    )

    # Base increment variance is 4 + 4 + 1 - 2*2 = 5.  The two travel-time
    # dispersions are then added to their respective likelihood variances.
    assert result.propagated_sigma_years == pytest.approx(math.sqrt(5.0))
    assert result.direct_sigma_years == pytest.approx(1.0)
    assert result.indirect_sigma_years == pytest.approx(5.0)
    assert result.direct_standardized_residual == pytest.approx(0.0)
    assert result.indirect_standardized_residual == pytest.approx(
        -20.0 / math.sqrt(30.0)
    )
    assert result.log_bayes_factor_direct_vs_indirect is not None
    assert result.log_bayes_factor_direct_vs_indirect > 0.0


def test_invalid_endpoint_covariance_is_rejected() -> None:
    with pytest.raises(ValueError, match="covariance"):
        compute_direction_evidence(
            upstream_age_years=1.0,
            downstream_age_years=2.0,
            upstream_sigma_years=1.0,
            downstream_sigma_years=1.0,
            age_covariance_years2=2.0,
        )
