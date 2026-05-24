"""Unit tests for reference graph analysis and geometry calculations in M4 pipeline."""

from __future__ import annotations

import math
import numpy as np
import pandas as pd
import pytest

from M4.m4_topology_benchmark.scripts.run_public_archive_pipeline import (
    _polygon_area,
    _convex_hull,
    _convex_iou,
    _bootstrap_ci,
    compute_topology_agreement_and_particles,
    compute_graph_priors,
    compute_topology_summary,
)


def test_polygon_area():
    """Verify polygon area calculation for simple shapes."""
    # Unit square
    square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    assert math.isclose(_polygon_area(square), 1.0)

    # Empty/invalid polygons
    assert _polygon_area([]) == 0.0
    assert _polygon_area([[0.0, 0.0], [1.0, 1.0]]) == 0.0


def test_convex_hull():
    """Verify convex hull extraction."""
    points = [[0.0, 0.0], [1.0, 0.0], [0.5, 0.5], [1.0, 1.0], [0.0, 1.0]]
    hull = _convex_hull(points)
    # The point [0.5, 0.5] is strictly inside the square and should not be on the hull
    assert len(hull) == 4
    for pt in [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]:
        assert pt in hull
    assert [0.5, 0.5] not in hull


def test_convex_iou():
    """Verify Intersection-over-Union for overlapping convex hulls."""
    square_a = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]
    square_b = [[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]]
    
    # Area A = 4, Area B = 4, Intersection = 2, Union = 6. IoU = 2 / 6 = 1/3
    iou, area_a, area_b, inter = _convex_iou(square_a, square_b)
    assert math.isclose(iou, 1.0 / 3.0)
    assert math.isclose(area_a, 4.0)
    assert math.isclose(area_b, 4.0)
    assert math.isclose(inter, 2.0)


def test_bootstrap_ci():
    """Verify bootstrap confidence interval estimation."""
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    ci = _bootstrap_ci(values, np.mean, iterations=100)
    assert ci["n"] == 5.0
    assert math.isclose(ci["estimate"], 3.0)
    assert ci["ci_low"] < ci["estimate"]
    assert ci["ci_high"] > ci["estimate"]


def test_compute_topology_agreement_and_particles():
    """Verify topology agreement and particle summary calculations."""
    endpoints = pd.DataFrame([
        {
            "particle_id": 1, "initial_cell": 10, "final_cell": 20, 
            "x0": 0.0, "y0": 0.0, "z0": 0.0, "x": 10.0, "y": 10.0, "z": 10.0,
            "time": 100.0, "status": 1
        }
    ])
    pathlines = pd.DataFrame([
        {"particle_id": 1, "x": 0.0, "y": 0.0, "z": 0.0, "time": 0.0, "cell": 10, "sequence": 1},
        {"particle_id": 1, "x": 5.0, "y": 5.0, "z": 5.0, "time": 50.0, "cell": 15, "sequence": 2},
        {"particle_id": 1, "x": 10.0, "y": 10.0, "z": 10.0, "time": 100.0, "cell": 20, "sequence": 3},
    ])

    agreement, particles = compute_topology_agreement_and_particles(endpoints, pathlines)
    
    assert not agreement.empty
    assert len(agreement) == 1
    assert agreement.loc[0, "edge_id"] == "cell_10->cell_20"
    assert agreement.loc[0, "endpoint_particle_count"] == 1
    
    assert not particles.empty
    assert len(particles) == 1
    assert particles.loc[0, "particle_id"] == 1
    assert particles.loc[0, "pathline_elapsed_time"] == 100.0
    assert particles.loc[0, "endpoint_projection_preserved"]


def test_compute_topology_agreement_detects_endpoint_pathline_mismatch():
    """Endpoint/pathline agreement must compare parsed pathline cells, not endpoint labels twice."""
    endpoints = pd.DataFrame([
        {
            "particle_id": 1, "initial_cell": 10, "final_cell": 20,
            "x0": 0.0, "y0": 0.0, "z0": 0.0, "x": 10.0, "y": 10.0, "z": 10.0,
            "time": 100.0, "status": 1
        }
    ])
    pathlines = pd.DataFrame([
        {"particle_id": 1, "x": 0.0, "y": 0.0, "z": 0.0, "time": 0.0, "cell": 10, "sequence": 1},
        {"particle_id": 1, "x": 10.0, "y": 10.0, "z": 10.0, "time": 100.0, "cell": 30, "sequence": 2},
    ])

    agreement, particles = compute_topology_agreement_and_particles(endpoints, pathlines)

    assert not particles.loc[0, "endpoint_projection_preserved"]
    assert set(agreement["classification"]) == {"FN", "FP"}


def test_compute_graph_priors():
    """Verify extraction of edge travel-time priors from endpoints."""
    endpoints = pd.DataFrame([
        {
            "particle_id": 1, "initial_cell": 10, "final_cell": 20, 
            "x0": 0.0, "y0": 0.0, "z0": 0.0, "x": 10.0, "y": 10.0, "z": 10.0,
            "time": 100.0, "status": 1
        },
        {
            "particle_id": 2, "initial_cell": 10, "final_cell": 20, 
            "x0": 0.0, "y0": 0.0, "z0": 0.0, "x": 10.0, "y": 10.0, "z": 10.0,
            "time": 200.0, "status": 1
        }
    ])
    config = {"modpath_version": 7, "endpoint_file_in_zip": "test.endpoint"}
    priors = compute_graph_priors(endpoints, config)
    
    assert not priors.empty
    assert len(priors) == 1
    assert priors.loc[0, "edge_id"] == "cell_10->cell_20"
    assert priors.loc[0, "particle_count"] == 2
    assert priors.loc[0, "travel_time_mean"] == 150.0
    assert priors.loc[0, "travel_time_p50"] == 100.0
    assert priors.loc[0, "p_uv"] == 1.0
