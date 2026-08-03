from __future__ import annotations

from pathlib import Path

import pytest

from hydrosheaf.validation import generate_independent_candidate_universe


def _rows() -> list[dict[str, object]]:
    return [
        {"site_id": "A", "x_m": 0.0, "y_m": 0.0, "head_meas": 10.0},
        {"site_id": "B", "x_m": 100.0, "y_m": 0.0, "head_meas": 9.0},
        {"site_id": "C", "x_m": 200.0, "y_m": 0.0, "head_meas": 8.0},
        {"site_id": "D", "x_m": 300.0, "y_m": 0.0, "head_meas": 8.0},
    ]


def test_generation_is_order_invariant_and_hash_audited() -> None:
    first = generate_independent_candidate_universe(_rows(), max_neighbors=1)
    second = generate_independent_candidate_universe(list(reversed(_rows())), max_neighbors=1)

    assert first.input_hash == second.input_hash
    assert first.candidate_hash == second.candidate_hash
    assert [edge.edge for edge in first.edges] == [edge.edge for edge in second.edges]
    assert first.to_audit_record()["truth_blind"] is True
    assert first.to_audit_record()["truth_fields_seen"] == []


def test_head_orientation_is_not_taken_from_hydrosheaf_probability() -> None:
    universe = generate_independent_candidate_universe(_rows()[:2], max_neighbors=1)

    assert [(edge.u, edge.v) for edge in universe.edges] == [("A", "B")]
    attrs = universe.edges[0].attrs
    assert attrs["direction_rule"] == "higher_observed_head_to_lower"
    assert "p_uv" not in attrs
    assert "edge_confidence" not in attrs
    assert 0.0 <= float(attrs["independent_support"]) <= 1.0


def test_tied_or_missing_heads_keep_both_directions() -> None:
    tied = generate_independent_candidate_universe(
        [
            {"site_id": "A", "x_m": 0.0, "y_m": 0.0, "head_meas": 10.0},
            {"site_id": "B", "x_m": 20.0, "y_m": 0.0, "head_meas": 10.05},
        ],
        max_neighbors=1,
    )
    missing = generate_independent_candidate_universe(
        [
            {"site_id": "A", "x_m": 0.0, "y_m": 0.0},
            {"site_id": "B", "x_m": 20.0, "y_m": 0.0, "head_meas": 10.0},
        ],
        max_neighbors=1,
    )

    assert {edge.edge for edge in tied.edges} == {("A", "B"), ("B", "A")}
    assert {edge.edge for edge in missing.edges} == {("A", "B"), ("B", "A")}
    assert {edge.attrs["direction_rule"] for edge in tied.edges} == {"head_tie_both_directions"}
    assert {
        edge.attrs["direction_rule"] for edge in missing.edges
    } == {"missing_observed_head_both_directions"}


def test_geographic_coordinates_are_projected_without_inference_imports() -> None:
    universe = generate_independent_candidate_universe(
        [
            {"site_id": "A", "lon": -1.0, "lat": 7.0, "head_meas": 10.0},
            {"site_id": "B", "lon": -0.999, "lat": 7.0, "head_meas": 9.0},
        ],
        max_neighbors=1,
    )

    assert universe.coordinate_system == "lon_lat"
    assert universe.edges[0].attrs["distance_km"] > 0.0

    source = Path("hydrosheaf/validation/specialist_candidate_generation.py").read_text(
        encoding="utf-8"
    )
    assert "hydrosheaf.inference" not in source
    assert "from hydrosheaf.inference" not in source


def test_truth_fields_are_rejected_before_candidate_generation() -> None:
    with pytest.raises(ValueError, match="Truth/reference field"):
        generate_independent_candidate_universe(
            [
                {
                    "site_id": "A",
                    "x_m": 0.0,
                    "y_m": 0.0,
                    "head_meas": 10.0,
                    "true_age_years": 12.0,
                }
            ]
        )
