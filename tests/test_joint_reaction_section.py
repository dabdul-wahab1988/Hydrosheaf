from types import SimpleNamespace
from unittest import mock

import pytest

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.topology_posterior import make_topology_cost_fn
from hydrosheaf.sheaf import solve_joint_reaction_section
from hydrosheaf.sheaf.directed_section import DirectedEdgeMap
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf


def _reaction_edge(
    *,
    reaction_matrix,
    signed_mask,
    transport_offset=(0.0,),
    offset=(0.0,),
):
    return DirectedEdgeMap(
        edge=Edge("u_to_v", "u", "v"),
        alpha=1.0,
        offset=list(offset),
        weight=1.0,
        objective=0.0,
        transport_model="fixed",
        endmember_id=None,
        residual_norm=0.0,
        transport_offset=list(transport_offset),
        reaction_matrix=[list(row) for row in reaction_matrix],
        reaction_labels=[f"reaction_{index}" for index in range(len(reaction_matrix))],
        reaction_penalty_scales=[1.0] * len(reaction_matrix),
        signed_reaction_mask=list(signed_mask),
        reaction_extents=[0.0] * len(reaction_matrix),
    )


def test_joint_solver_fits_node_states_and_reaction_extents_together():
    # The legacy offset is deliberately incompatible with the observations.
    # A joint solve must use the separate transport term and estimate z = 2.
    edge_map = _reaction_edge(
        reaction_matrix=[[1.0]],
        signed_mask=[False],
        offset=[999.0],
    )

    result = solve_joint_reaction_section(
        ["u", "v"],
        [edge_map],
        {"u": [1.0], "v": [3.0]},
        ["Ca"],
        obs_weight=1.0e6,
        diag_eps=0.0,
        max_iter=2000,
        tol=1.0e-9,
    )

    assert result.converged
    assert result.node_states["u"] == pytest.approx([1.0], abs=1.0e-6)
    assert result.node_states["v"] == pytest.approx([3.0], abs=1.0e-6)
    assert result.reaction_extents["u_to_v"] == pytest.approx([2.0], abs=1.0e-6)
    assert result.edge_residuals["u_to_v"] < 1.0e-6


def test_joint_solver_applies_declared_reaction_sign_constraints():
    observations = {"u": [2.0], "v": [1.0]}
    signed = solve_joint_reaction_section(
        ["u", "v"],
        [_reaction_edge(reaction_matrix=[[1.0]], signed_mask=[True])],
        observations,
        ["Ca"],
        obs_weight=1.0e6,
        diag_eps=0.0,
        max_iter=2000,
        tol=1.0e-9,
    )
    unsigned = solve_joint_reaction_section(
        ["u", "v"],
        [_reaction_edge(reaction_matrix=[[1.0]], signed_mask=[False])],
        observations,
        ["Ca"],
        obs_weight=1.0e6,
        diag_eps=0.0,
        max_iter=2000,
        tol=1.0e-9,
    )

    assert signed.converged
    assert signed.reaction_extents["u_to_v"] == pytest.approx([-1.0], abs=1.0e-6)
    assert unsigned.reaction_extents["u_to_v"][0] == pytest.approx(0.0, abs=1.0e-9)


def test_joint_solver_rejects_legacy_maps_without_separate_metadata():
    legacy_map = DirectedEdgeMap(
        edge=Edge("u_to_v", "u", "v"),
        alpha=1.0,
        offset=[0.0],
        weight=1.0,
        objective=0.0,
        transport_model="fixed",
        endmember_id=None,
        residual_norm=0.0,
    )

    with pytest.raises(ValueError, match="transport_offset"):
        solve_joint_reaction_section(
            ["u", "v"],
            [legacy_map],
            {"u": [1.0], "v": [1.0]},
            ["Ca"],
        )


def test_config_rejects_nonfinite_joint_solver_tolerance():
    config = Config()
    config.sheaf_joint_reaction_tol = float("nan")

    with pytest.raises(ValueError, match="finite and non-negative"):
        config.validate()


def test_topology_refinement_records_joint_reaction_solution():
    config = Config(
        ion_order=["Ca"],
        weights=[1.0],
        conservative_weights=[1.0],
        active_minerals=[],
        transport_models_enabled=["evap"],
        edge_max_neighbors=1,
        sheaf_max_iter=1,
        sheaf_joint_reaction_max_iter=2000,
    )
    selected = refine_edges_with_sheaf(
        [
            {"node_id": "u", "Ca": 1.0, "lat": 0.0, "lon": 0.0},
            {"node_id": "v", "Ca": 3.0, "lat": 0.0, "lon": 0.1},
        ],
        [Edge("u_to_v", "u", "v")],
        config,
    )

    assert len(selected) == 1
    attrs = selected[0].attrs
    assert attrs["sheaf_joint_reaction_scope"] == "fixed_transport_maps"
    assert attrs["sheaf_joint_reaction_status"] == "converged"
    assert attrs["sheaf_joint_reaction_iterations"] > 0
    assert isinstance(attrs["sheaf_joint_reaction_labels"], list)
    assert isinstance(attrs["sheaf_joint_reaction_extents"], list)


def test_legacy_switch_preserves_fixed_offset_section_workflow():
    config = Config(
        ion_order=["Ca"],
        weights=[1.0],
        conservative_weights=[1.0],
        active_minerals=[],
        transport_models_enabled=["evap"],
        edge_max_neighbors=1,
        sheaf_max_iter=1,
        sheaf_joint_reaction_enabled=False,
    )
    selected = refine_edges_with_sheaf(
        [
            {"node_id": "u", "Ca": 1.0, "lat": 0.0, "lon": 0.0},
            {"node_id": "v", "Ca": 3.0, "lat": 0.0, "lon": 0.1},
        ],
        [Edge("u_to_v", "u", "v")],
        config,
    )

    assert selected[0].attrs["sheaf_joint_reaction_status"] == "disabled"


def test_topology_posterior_likelihood_uses_joint_reaction_solver():
    config = Config(
        ion_order=["Ca"],
        weights=[1.0],
        conservative_weights=[1.0],
        active_minerals=[],
        transport_models_enabled=["evap"],
    )
    samples = {
        "u": {"site_id": "u", "Ca": 1.0},
        "v": {"site_id": "v", "Ca": 3.0},
    }
    edge = Edge("u_to_v", "u", "v")

    with mock.patch(
        "hydrosheaf.sheaf.joint_reaction.solve_joint_reaction_section",
        return_value=SimpleNamespace(edge_residuals={"u_to_v": 1.25}),
    ) as joint_solver:
        cost = make_topology_cost_fn(samples, config)([edge])

    joint_solver.assert_called_once()
    assert cost >= 1.25
