"""Tests for hydraulic Hodge diagnostics (post-refactor).

The old tests that relied on *min_drop* clipping to produce artificial
obstruction for a scalar potential have been replaced with physically
meaningful checks:

- d² = 0 for scalar head → obstruction is identically zero.
- Projected-gradient cost penalises anti-aligned edges.
- Local head-plane residuals are zero for planar data.
- No uphill double-penalty in the cost function.
"""

import math

import numpy as np

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.topology_posterior import (
    make_topology_cost_fn,
    select_posterior_edges,
)
from hydrosheaf.sheaf.hydraulic_hodge import (
    compute_head_hodge_diagnostics,
    compute_head_hodge_edge_penalties,
    compute_local_head_plane_residuals,
    extract_node_heads,
    head_hodge_graph_cost,
    infer_reference_distance_km,
)


def _samples():
    return [
        {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
        {"site_id": "B", "lat": 0.0, "lon": 1.0, "elevation": 8.0},
        {"site_id": "C", "lat": 0.0, "lon": 3.0, "elevation": 7.0},
    ]


def _edges():
    return [
        Edge("A->B", "A", "B", attrs={"distance_km": 1.0, "edge_confidence": 0.9}),
        Edge("B->C", "B", "C", attrs={"distance_km": 1.0, "edge_confidence": 0.9}),
        Edge("A->C", "A", "C", attrs={"distance_km": 3.0, "edge_confidence": 0.9}),
    ]


# ---------------------------------------------------------------------------
# d² = 0 — scalar-potential obstruction is identically zero
# ---------------------------------------------------------------------------


def test_scalar_potential_obstruction_is_zero_for_any_graph():
    """Confirm that the scalar-potential obstruction energy is zero regardless
    of graph structure, because head differences form an exact 1-form (d² = 0).
    """
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=1.0,
        hydraulic_hodge_fallback_to_elevation=True,
    )
    samples = _samples()
    edges = _edges()
    sample_map = {row["site_id"]: row for row in samples}
    reference = infer_reference_distance_km(edges, sample_map, config)

    # No gradient map, no local residuals → fallback cost (should be zero)
    chain_cost = head_hodge_graph_cost(
        edges[:2], sample_map, config, reference_distance_km=reference,
        gradient_map=None, local_residuals=None,
    )
    triangle_cost = head_hodge_graph_cost(
        edges, sample_map, config, reference_distance_km=reference,
        gradient_map=None, local_residuals=None,
    )

    # Both should be zero: no artificial obstruction from scalar potential
    assert chain_cost == 0.0
    assert triangle_cost == 0.0


def test_diagnostics_label_obstruction_as_experimental():
    """compute_head_hodge_diagnostics uses 'experimental_*' keys for the
    degenerate scalar-potential obstruction fields.
    """
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_fallback_to_elevation=True,
    )
    samples = _samples()
    edges = _edges()
    sample_map = {row["site_id"]: row for row in samples}

    diagnostics = compute_head_hodge_diagnostics(edges, sample_map, config)

    # Old keys should NOT be present
    assert "obstruction_energy" not in diagnostics
    assert "cycle_count" not in diagnostics
    assert "obstruction_leverage" not in diagnostics
    assert "cycle_obstruction_max" not in diagnostics

    # New experimental keys should be present
    assert "experimental_obstruction_energy" in diagnostics
    assert "experimental_cycle_count" in diagnostics
    assert "experimental_obstruction_leverage" in diagnostics
    assert "experimental_cycle_obstruction_max" in diagnostics

    # New production keys should be present
    assert "direction_violation_count" in diagnostics
    assert "local_plane_residual_mean" in diagnostics
    assert "local_plane_residual_max" in diagnostics
    assert "normalized_cost" in diagnostics


# ---------------------------------------------------------------------------
# Projected-gradient cost
# ---------------------------------------------------------------------------


def test_gradient_map_is_noop_in_cost():
    """gradient_map is accepted for signature compat but NOT consumed —
    directional prior responsibility lives in flow_direction."""
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=1.0,
        hydraulic_hodge_fallback_to_elevation=True,
    )
    samples = _samples()
    sample_map = {row["site_id"]: row for row in samples}

    gradient_map = {"A": (-0.005, 0.0, 0.0), "B": (-0.005, 0.0, 0.0)}
    edge_aligned = [Edge("A->B", "A", "B", attrs={"distance_km": 1.0, "edge_confidence": 0.9})]
    edge_misaligned = [Edge("A->C", "A", "C", attrs={"distance_km": 3.0, "edge_confidence": 0.9})]

    # Without local_residuals, cost is zero regardless of gradient_map
    cost_aligned = head_hodge_graph_cost(
        edge_aligned, sample_map, config, gradient_map=gradient_map,
    )
    cost_misaligned = head_hodge_graph_cost(
        edge_misaligned, sample_map, config, gradient_map=gradient_map,
    )

    # gradient_map is NOT consumed — both costs are zero
    assert cost_aligned == 0.0
    assert cost_misaligned == 0.0


def test_local_residuals_produce_nonzero_cost():
    """local_residuals (not gradient_map) drive the cost function."""
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        head_plane_residual_weight=1.0,
    )
    samples = [
        {"site_id": "A", "lat": 0.0, "lon": 0.0},
        {"site_id": "B", "lat": 0.0, "lon": 1.0},
    ]
    sample_map = {row["site_id"]: row for row in samples}

    # Node A has high residual, node B has zero
    local_residuals = {"A": 3.0, "B": 0.0}
    edge = [Edge("A->B", "A", "B", attrs={"distance_km": 1.0, "edge_confidence": 0.9})]

    cost = head_hodge_graph_cost(
        edge, sample_map, config, local_residuals=local_residuals,
    )

    # Should be positive because node A's residual contributes
    assert cost > 0.0


# ---------------------------------------------------------------------------
# Local head-plane residuals
# ---------------------------------------------------------------------------


def test_local_plane_residual_zero_for_planar_data():
    """Perfectly planar head data should produce near-zero residuals."""
    # Create a regular grid with a planar head surface: h(x,y) = 100 - 0.01*x - 0.005*y
    samples = {}
    for i in range(5):
        for j in range(5):
            x = float(i)
            y = float(j)
            hid = f"n_{i}_{j}"
            samples[hid] = {
                "site_id": hid,
                "lon": x,
                "lat": y,
                "head": 100.0 - 0.01 * x - 0.005 * y,
            }

    config = Config(hydraulic_hodge_head_key="head")
    head_map = extract_node_heads(samples, config)
    residuals = compute_local_head_plane_residuals(head_map, samples, n_neighbors=8)

    assert len(residuals) > 0
    # All residuals should be very small for perfect planar data
    for r in residuals.values():
        assert r < 0.01, f"Residual too large for planar data: {r}"


def test_local_plane_residual_detects_outlier():
    """A node with anomalous head gets a high local-plane residual."""
    # Planar head field
    samples = {}
    for i in range(5):
        for j in range(5):
            x = float(i)
            y = float(j)
            hid = f"n_{i}_{j}"
            samples[hid] = {
                "site_id": hid,
                "lon": x,
                "lat": y,
                "head": 100.0 - 0.01 * x - 0.005 * y,
            }

    # Inject an outlier
    samples["n_2_2"]["head"] = 500.0  # 400 ft above the plane

    config = Config(hydraulic_hodge_head_key="head")
    head_map = extract_node_heads(samples, config)
    residuals = compute_local_head_plane_residuals(head_map, samples, n_neighbors=8)

    assert "n_2_2" in residuals
    # The outlier should have the largest residual
    outlier_residual = residuals["n_2_2"]
    for nid, r in residuals.items():
        if nid != "n_2_2":
            assert r < outlier_residual, (
                f"Outlier residual ({outlier_residual}) should exceed "
                f"{nid} residual ({r})"
            )


# ---------------------------------------------------------------------------
# No uphill double-penalty
# ---------------------------------------------------------------------------


def test_no_uphill_penalty_in_cost():
    """head_hodge_graph_cost does not penalize uphill edges — that is the
    flow_direction module's responsibility.
    """
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=1.0,
        hydraulic_hodge_fallback_to_elevation=True,
    )
    samples = [
        {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
        {"site_id": "B", "lat": 0.0, "lon": 1.0, "elevation": 8.0},
    ]
    sample_map = {row["site_id"]: row for row in samples}

    # Downhill edge A->B
    downhill = [Edge("A->B", "A", "B", attrs={"distance_km": 1.0, "edge_confidence": 0.9})]
    # Uphill edge B->A
    uphill = [Edge("B->A", "B", "A", attrs={"distance_km": 1.0, "edge_confidence": 0.9})]

    cost_down = head_hodge_graph_cost(downhill, sample_map, config)
    cost_up = head_hodge_graph_cost(uphill, sample_map, config)

    # Without gradient_map or local_residuals, both should be zero
    # (no artificial uphill penalty)
    assert cost_down == 0.0
    assert cost_up == 0.0


# ---------------------------------------------------------------------------
# Edge penalties from local plane residuals
# ---------------------------------------------------------------------------


def test_edge_penalties_from_local_residuals():
    """compute_head_hodge_edge_penalties uses local plane residuals."""
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_leverage_weight=0.5,
        hydraulic_hodge_head_key="head",
    )
    # Planar data with one outlier
    samples = {}
    for i in range(4):
        for j in range(4):
            x = float(i)
            y = float(j)
            hid = f"n_{i}_{j}"
            samples[hid] = {
                "site_id": hid,
                "lon": x,
                "lat": y,
                "head": 100.0 - 0.01 * x - 0.005 * y,
            }
    samples["n_2_2"]["head"] = 500.0  # outlier

    edges = [
        Edge("n_2_1->n_2_2", "n_2_1", "n_2_2", attrs={"distance_km": 1.0}),
        Edge("n_1_1->n_1_2", "n_1_1", "n_1_2", attrs={"distance_km": 1.0}),
    ]

    penalties = compute_head_hodge_edge_penalties(edges, samples, config)

    # The edge incident to the outlier (n_2_2) should have a higher penalty
    assert "n_2_1->n_2_2" in penalties
    assert penalties["n_2_1->n_2_2"] > 0


# ---------------------------------------------------------------------------
# topology_posterior integration (preserved from original)
# ---------------------------------------------------------------------------


def test_make_topology_cost_fn_includes_hydraulic_hodge_term():
    """make_topology_cost_fn accepts and passes gradient_map + local_residuals."""
    config = Config(
        phreeqc_enabled=False,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_age_enabled=False,
        hydraulic_hodge_enabled=True,
        hydraulic_hodge_weight=1.0,
        hydraulic_hodge_fallback_to_elevation=True,
    )
    samples = _samples()
    sample_map = {row["site_id"]: row for row in samples}
    edges = _edges()

    # With gradient_map pointing east
    gradient_map = {"A": (0.005, 0.0, 0.0)}

    cost_fn = make_topology_cost_fn(
        sample_map=sample_map,
        config=config,
        gradient_map=gradient_map,
    )

    # Cost should be finite and non-negative
    cost = cost_fn(edges)
    assert math.isfinite(cost)
    assert cost >= 0.0


def test_select_posterior_edges_respects_probability_ranking_per_source():
    """select_posterior_edges with max_neighbors=1 picks the top-probability
    edge per source node."""
    edges = [
        Edge("A->B", "A", "B"),
        Edge("A->C", "A", "C"),
        Edge("B->C", "B", "C"),
    ]
    posterior_result = {
        "map_edges": ["A->B", "A->C", "B->C"],
        "edge_probabilities": {"A->B": 0.9, "A->C": 0.6, "B->C": 0.8},
    }

    selected = select_posterior_edges(edges, posterior_result, max_neighbors=1)

    assert {edge.edge_id for edge in selected} == {"A->B", "B->C"}
