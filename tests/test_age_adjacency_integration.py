from __future__ import annotations

import pytest

from hydrosheaf import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf


def _sample(site_id: str, age_years: float) -> dict[str, object]:
    return {
        "site_id": site_id,
        "head_meas": 100.0 - age_years,
        "mean_age_years": age_years,
        "mean_age_std_years": 0.5,
        "tracer_identifiable": True,
        "lat": 0.0,
        "lon": 0.0,
    }


def _config() -> Config:
    return Config(
        edge_max_neighbors=2,
        edge_radius_km=10.0,
        sheaf_weight_head_prior=0.0,
        sheaf_weight_isotope=0.0,
        sheaf_weight_cl=0.0,
        sheaf_weight_age=1.0,
        sheaf_age_process_sigma_years=0.25,
        sheaf_age_travel_time_cv=0.0,
        sheaf_age_travel_cost_weight=0.0,
        sheaf_age_adjacency_enabled=True,
        sheaf_age_adjacency_weight=1.0,
        sheaf_max_iter=1,
    )


def test_optional_adjacency_term_ranks_direct_match_above_cumulative_skip() -> None:
    samples = [
        _sample("U", 100.0),
        _sample("V", 110.0),
        _sample("W", 130.0),
    ]
    candidates = [
        Edge(
            edge_id="U->V",
            u="U",
            v="V",
            attrs={"direct_travel_years": 10.0, "indirect_travel_years": 30.0},
        ),
        Edge(
            edge_id="U->W",
            u="U",
            v="W",
            attrs={"direct_travel_years": 10.0, "indirect_travel_years": 30.0},
        ),
    ]

    refine_edges_with_sheaf(samples, candidates, _config())

    direct = candidates[0].attrs
    skip = candidates[1].attrs
    assert direct["sheaf_age_adjacency_status"] == "scored"
    assert skip["sheaf_age_adjacency_status"] == "scored"
    assert direct["sheaf_age_direct_probability_equal_prior"] > 0.99
    assert skip["sheaf_age_direct_probability_equal_prior"] < 0.01
    assert direct["sheaf_cost_age"] < skip["sheaf_cost_age"]
    assert direct["sheaf_age_adjacency_evidence_available"] is True


def test_shared_legacy_travel_uncertainty_is_not_used_as_rtd_dispersion() -> None:
    samples = [_sample("U", 100.0), _sample("V", 110.0)]
    candidate = Edge(
        edge_id="U->V",
        u="U",
        v="V",
        attrs={"direct_travel_years": 10.0, "indirect_travel_years": 30.0},
    )

    refine_edges_with_sheaf(samples, [candidate], _config())

    evidence = candidate.attrs["sheaf_age_adjacency"]
    assert evidence["direct_sigma_years"] == pytest.approx(0.0)
    assert evidence["indirect_sigma_years"] == pytest.approx(0.0)
    assert evidence["propagated_sigma_years"] == pytest.approx(0.75)


def test_missing_indirect_hypothesis_is_annotated_but_does_not_change_age_cost() -> None:
    samples = [_sample("U", 100.0), _sample("V", 110.0)]
    candidate = Edge(
        edge_id="U->V",
        u="U",
        v="V",
        attrs={"direct_travel_years": 10.0},
    )
    config = _config()
    selected = refine_edges_with_sheaf(samples, [candidate], config)

    assert selected
    assert candidate.attrs["sheaf_age_adjacency_status"] == "ambiguous"
    assert candidate.attrs["sheaf_age_adjacency_evidence_available"] is False
    assert "indirect_travel_time_missing" in candidate.attrs["sheaf_age_adjacency_flags"]
    assert "sheaf_age_direct_probability_equal_prior" not in candidate.attrs


def test_adjacency_configuration_is_validated() -> None:
    with pytest.raises(ValueError, match="sheaf_age_adjacency_weight"):
        Config(sheaf_age_adjacency_weight=-1.0).validate()
    with pytest.raises(ValueError, match="sheaf_age_adjacency_enabled"):
        Config(sheaf_age_adjacency_enabled=1).validate()
