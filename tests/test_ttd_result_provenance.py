"""Provenance and abstention contracts for network and dynamic TTD results."""

from __future__ import annotations

import networkx as nx
import numpy as np

from hydrosheaf.nuclear.multi_tracer_graph_inversion import (
    DUPLICATE_TRACER_ALIAS_REASON_CODE,
    MISSING_TRACER_HISTORY_REASON_CODE,
    MultiTracerGraphConfig,
    solve_joint_multitracer_node_inversion,
)
from hydrosheaf.nuclear.ttd_diagnostics import CODE_MISSING_FORWARD_SYSTEM
from hydrosheaf.nuclear.ttd_grid import build_uniform_ttd_grid
from hydrosheaf.nuclear.ttd_kernel_builder import (
    MultiTracerForwardSystem,
    NodeTracerPanel,
    TracerObservation,
    build_forward_system,
)
from hydrosheaf.nuclear.ttd_network_solver import (
    SOURCE_HISTORY_UNAVAILABLE_REASON_CODE,
    solve_network_ttd,
    solve_single_node_ttd,
)


def _radioactive_system(grid):
    return build_forward_system(
        NodeTracerPanel(
            "A",
            2024.0,
            (
                TracerObservation("3H", 2.0, 0.2),
                TracerObservation("14C", 50.0, 1.0),
            ),
        ),
        grid,
    )


def test_single_node_success_and_abstention_have_direct_provenance_fields():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    system = _radioactive_system(grid)

    success = solve_single_node_ttd(system, grid)
    assert success.status == "ESTIMATED"
    assert success.inference_family == "network_ttd"
    assert success.claim_scope == "conditional_inference"
    assert success.field_validation_status == "not_performed"
    assert success.stable_isotopes_used is False
    assert success.radioactive_tracers_used is True
    assert success.source_history_required is True
    assert success.source_history_available is True
    assert success.reason_code is None

    abstained = solve_single_node_ttd(
        system,
        grid,
        lambda_smoothness=-1.0,
    )
    assert abstained.status == "ABSTAIN"
    assert abstained.reason_code == abstained.abstention_reasons[0]
    assert abstained.inference_family == "network_ttd"
    assert abstained.source_history_required is True
    assert abstained.source_history_available is True


def test_network_success_and_missing_system_abstention_preserve_provenance():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    graph = nx.DiGraph()
    graph.add_node("A")
    system = _radioactive_system(grid)

    success = solve_network_ttd(graph, {"A": system}, {}, grid)
    assert success.status == "ESTIMATED"
    assert success.inference_family == "network_ttd"
    assert success.claim_scope == "conditional_inference"
    assert success.field_validation_status == "not_performed"
    assert success.stable_isotopes_used is False
    assert success.radioactive_tracers_used is True
    assert success.source_history_required is True
    assert success.source_history_available is True
    assert success.reason_code is None
    assert success.node_results["A"].reason_code is None
    assert success.node_results["A"].source_history_available is True

    missing = nx.DiGraph([("A", "B")])
    abstained = solve_network_ttd(missing, {"A": system}, {}, grid)
    assert abstained.status == "ABSTAIN"
    assert CODE_MISSING_FORWARD_SYSTEM in abstained.abstention_reasons
    assert abstained.reason_code == CODE_MISSING_FORWARD_SYSTEM
    assert abstained.inference_family == "network_ttd"
    assert abstained.source_history_required is True
    assert abstained.source_history_available is None
    assert abstained.node_results["A"].reason_code == CODE_MISSING_FORWARD_SYSTEM
    assert abstained.node_results["B"].reason_code == CODE_MISSING_FORWARD_SYSTEM


def test_stable_isotope_forward_system_is_conditional_and_not_radioactive():
    grid = build_uniform_ttd_grid(max_age_years=2.0, dt_years=1.0)
    system = MultiTracerForwardSystem(
        node_id="iso",
        sample_year=2024.0,
        grid=grid,
        tracers=("d18O", "d2H"),
        matrix=np.asarray(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        ),
        observations=np.asarray([0.7, 0.3]),
        sigmas=np.asarray([0.1, 0.8]),
        weights=np.ones(2),
    )

    result = solve_single_node_ttd(system, grid)
    assert result.status == "ESTIMATED"
    assert result.stable_isotopes_used is True
    assert result.radioactive_tracers_used is False
    assert result.source_history_required is True
    assert result.source_history_available is None
    assert "STABLE_ISOTOPE_SOURCE_HISTORY_UNVERIFIED" in result.diagnostics["provenance_warnings"]
    assert result.diagnostics["provenance"]["claim_scope"] == "conditional_inference"


def test_dynamic_graph_missing_history_abstains_without_d_excess_conflict():
    times = np.arange(12)
    result = solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=times,
        tracer_observations={"d18O": np.linspace(-5.0, -4.0, len(times))},
        candidate_parents={"A": {}},
        config=MultiTracerGraphConfig(min_required_tracers=1),
    )

    assert result.status == "ABSTAIN"
    assert MISSING_TRACER_HISTORY_REASON_CODE in result.reason_codes
    assert result.reason_code == MISSING_TRACER_HISTORY_REASON_CODE
    assert result.provenance["stable_isotopes_used"] is True
    assert result.provenance["radioactive_tracers_used"] is False
    assert result.provenance["source_history_required"] is True
    assert result.provenance["source_history_available"] is False
    assert result.provenance["warnings"]
    assert not any("d_excess" in reason for reason in result.reason_codes)


def test_dynamic_graph_canonicalizes_tracer_aliases_before_activation():
    times = np.arange(12)
    result = solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=times,
        tracer_observations={"δ18O": np.linspace(-5.0, -4.0, len(times))},
        candidate_parents={"A": {}},
        config=MultiTracerGraphConfig(min_required_tracers=1),
    )

    assert result.status == "ABSTAIN"
    assert result.active_tracers == ("d18O",)
    assert MISSING_TRACER_HISTORY_REASON_CODE in result.reason_codes


def test_dynamic_graph_abstains_on_colliding_tracer_aliases():
    times = np.arange(12)
    result = solve_joint_multitracer_node_inversion(
        target_node="B",
        target_times=times,
        tracer_observations={
            "d18O": np.linspace(-5.0, -4.0, len(times)),
            "δ18O": np.linspace(-4.0, -3.0, len(times)),
        },
        candidate_parents={"A": {}},
        config=MultiTracerGraphConfig(min_required_tracers=1),
    )

    assert result.status == "ABSTAIN"
    assert result.reason_code == DUPLICATE_TRACER_ALIAS_REASON_CODE
    assert DUPLICATE_TRACER_ALIAS_REASON_CODE in result.reason_codes
    assert result.diagnostics["duplicate_tracer_aliases"] == ("d18O",)


def test_explicitly_unavailable_network_history_has_machine_reason_code():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    system = _radioactive_system(grid)
    # Forward systems are immutable dataclasses, so this checks the helper's
    # optional metadata contract with a small proxy carrying the same fields.
    system_with_metadata = type(
        "ForwardSystemWithHistoryMetadata",
        (),
        {
            **system.__dict__,
            "source_history_available": False,
        },
    )()

    graph = nx.DiGraph()
    graph.add_node("A")
    result = solve_network_ttd(graph, {"A": system_with_metadata}, {}, grid)

    assert result.status == "ABSTAIN"
    assert result.reason_code == SOURCE_HISTORY_UNAVAILABLE_REASON_CODE
    assert SOURCE_HISTORY_UNAVAILABLE_REASON_CODE in result.abstention_reasons
