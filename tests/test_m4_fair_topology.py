from __future__ import annotations

import math

import pandas as pd
import pytest

import hydrosheaf.validation.m4_path_aware as m4_path_aware
from hydrosheaf.validation.m4_fair import (
    FEET_TO_METRES,
    M4_A,
    M4_B,
    SavageProjectedFrame,
    build_savage_observations,
)
from hydrosheaf.validation.topology_v2 import (
    MODEL_FEATURE_SETS,
    TopologyV2Config,
    build_topology_v2_feature_rows,
    generate_topology_v2_candidate_universe,
)


def _nodes() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"node_id": "cell_1", "x": 100.0, "y": 200.0, "z": 10.0},
            {"node_id": "cell_2", "x": 200.0, "y": 200.0, "z": 9.0},
            {"node_id": "cell_3", "x": 100.0, "y": 300.0, "z": None},
        ]
    )


def test_savage_frame_uses_projected_euclidean_distance_not_haversine():
    frame = SavageProjectedFrame(rotation_deg=-12.0)
    first = frame.project(100.0, 200.0)
    second = frame.project(200.0, 200.0)

    assert math.hypot(second[0] - first[0], second[1] - first[1]) == pytest.approx(
        100.0 * FEET_TO_METRES
    )
    assert frame.to_dict()["latitude_longitude_interpretation"] is False


def test_m4_a_excludes_heads_and_preserves_sparse_elevation_proxy():
    observations = build_savage_observations(_nodes(), mode=M4_A)

    assert all("hydraulic_head" not in row for row in observations)
    assert "elevation" in observations[0]
    assert "elevation" not in observations[2]


def test_m4_b_requires_complete_heads_and_keeps_public_budget_channels():
    context = {
        "cell_1": {
            "well_rate": 25.0,
            "river_leakage": 0.0,
            "recharge": 4.0,
            "head_boundary_flux": 0.0,
        },
        "cell_2": {
            "well_rate": -100.0,
            "river_leakage": -2.0,
            "recharge": 0.0,
            "head_boundary_flux": -3.0,
        },
        "cell_3": {
            "well_rate": 0.0,
            "river_leakage": 0.0,
            "recharge": 0.0,
            "head_boundary_flux": 0.0,
        },
    }
    heads = {"cell_1": 300.0, "cell_2": 299.0, "cell_3": 298.0}
    observations = build_savage_observations(
        _nodes(), mode=M4_B, heads=heads, budget_context=context
    )

    assert observations[0]["hydraulic_head"] == pytest.approx(300.0 * FEET_TO_METRES)
    assert observations[1]["well_rate"] == -100.0
    assert observations[1]["river_leakage"] == -2.0
    with pytest.raises(ValueError, match="actual FHD head"):
        build_savage_observations(
            _nodes(), mode=M4_B, heads={"cell_1": 300.0}, budget_context=context
        )


def test_m4_candidate_universe_is_all_pairs_and_head_gradient_is_soft_evidence():
    context = {f"cell_{index}": {"well_rate": 0.0} for index in (1, 2, 3)}
    observations = build_savage_observations(
        _nodes(),
        mode=M4_B,
        heads={"cell_1": 300.0, "cell_2": 299.0, "cell_3": 298.0},
        budget_context=context,
    )
    universe = generate_topology_v2_candidate_universe(
        observations,
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    rows = build_topology_v2_feature_rows(universe, observations)
    edge = next(row for row in rows if row.edge_id == "cell_1->cell_2")

    assert len(universe.edges) == 6
    assert universe.truth_blind is True
    assert edge.features["head_delta_m"] == pytest.approx(FEET_TO_METRES)
    assert edge.features["gradient_m_per_km"] > 0.0
    assert edge.features["target_pumping_sink_strength"] == 0.0


def test_m4_c_path_features_cover_all_pairs_and_use_public_sink_context(monkeypatch):
    monkeypatch.setattr(
        m4_path_aware,
        "aggregate_savage_cbc_face_activity",
        lambda *args, **kwargs: (
            {1: 1.0, 2: 1.0, 3: 1.0, 4: 4.0},
            {"truth_blind": True, "direction_from_cbc_sign": False},
        ),
    )
    observations = [
        {"site_id": "cell_1", "well_rate": 0.0},
        {"site_id": "cell_2", "well_rate": 0.0},
        {"site_id": "cell_4", "well_rate": -100.0},
    ]
    from hydrosheaf.physics.modflow_head import build_grid_geometry_from_params

    result = m4_path_aware.build_path_aware_features(
        observations,
        head_map={1: 10.0, 2: 9.0, 3: 8.0, 4: 7.0},
        cbc_path="unused-for-mocked-activity",
        grid=build_grid_geometry_from_params(2, 2, 1, 1.0, 1.0),
    )

    assert len(result.edge_features) == 6
    assert result.face_activity_metadata["truth_blind"] is True
    assert result.edge_features["cell_1->cell_4"]["path_endpoint_sink_support"] == 1.0
    assert "path_endpoint_sink_support" in MODEL_FEATURE_SETS["M4_C_path_aware"]
