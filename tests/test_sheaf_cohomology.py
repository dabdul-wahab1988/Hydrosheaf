"""Tests for sheaf cohomology diagnostics."""

import unittest

import numpy as np

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.sheaf.directed_section import DirectedEdgeMap, build_edge_maps
from hydrosheaf.sheaf.cohomology import (
    build_coboundary_matrix,
    build_rhs_vector,
    compute_cohomology,
    compute_cycle_obstruction,
    compute_edge_leverage,
    attach_cohomology_attrs,
)

ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]


def _make_config():
    cfg = Config()
    cfg.mixing_endmembers = {}
    cfg.transport_models_enabled = ["evap"]
    return cfg


def _make_edge(edge_id: str, u: str, v: str) -> Edge:
    return Edge(edge_id=edge_id, u=u, v=v)


def _make_consistent_edge_maps() -> list:
    """Build edge maps for a consistent 3-node chain: A -> B -> C.

    alpha=1.0 means no evaporation (conservative), offset=[0,...].
    This produces a perfectly consistent system (zero obstruction).
    """
    cfg = _make_config()
    dim = len(ION_ORDER)
    node_values = {
        "A": [1.0] * dim,
        "B": [1.0] * dim,
        "C": [1.0] * dim,
    }
    e1 = _make_edge("A->B", "A", "B")
    e1.attrs = {"edge_confidence": 1.0}
    e2 = _make_edge("B->C", "B", "C")
    e2.attrs = {"edge_confidence": 1.0}
    return list(build_edge_maps([e1, e2], node_values, cfg).values())


def _make_inconsistent_cycle() -> list:
    """Build edge maps for a 3-node cycle A -> B -> C -> A with an
    inconsistent offset on the last edge, producing positive obstruction.
    """
    cfg = _make_config()
    dim = len(ION_ORDER)
    # All nodes identical except edge C->A has a nonzero offset
    node_values = {
        "A": [1.0] * dim,
        "B": [1.0] * dim,
        "C": [1.0] * dim,
    }
    e1 = _make_edge("A->B", "A", "B")
    e1.attrs = {"edge_confidence": 1.0}
    e2 = _make_edge("B->C", "B", "C")
    e2.attrs = {"edge_confidence": 1.0}

    # Manually construct the inconsistent edge C->A
    e3 = _make_edge("C->A", "C", "A")
    e3.attrs = {"edge_confidence": 1.0}
    em = DirectedEdgeMap(
        edge=e3,
        alpha=1.0,
        offset=[2.0] + [0.0] * (dim - 1),  # Inconsistency: Ca offset
        weight=1.0,
        objective=1.0,
        transport_model="evap",
        endmember_id=None,
        residual_norm=1.0,
    )

    ems = list(build_edge_maps([e1, e2], node_values, cfg).values())
    ems.append(em)
    return ems


class TestCoboundaryMatrix(unittest.TestCase):
    def test_consistent_chain_near_zero_obstruction(self):
        ems = _make_consistent_edge_maps()
        coh = compute_cohomology(ems, dim=len(ION_ORDER))
        self.assertAlmostEqual(coh["obstruction_energy"], 0.0, delta=1e-6)
        self.assertGreaterEqual(coh["h0_dim"], 0)

    def test_inconsistent_cycle_positive_obstruction(self):
        ems = _make_inconsistent_cycle()
        coh = compute_cohomology(ems, dim=len(ION_ORDER))
        self.assertGreater(coh["obstruction_energy"], 0.1)
        self.assertGreaterEqual(coh["h1_dim"], 0)

    def test_removing_bad_edge_reduces_obstruction(self):
        ems = _make_inconsistent_cycle()
        base_energy = compute_cohomology(ems, dim=len(ION_ORDER))["obstruction_energy"]
        self.assertGreater(base_energy, 0.0)

        # Removing any edge from the cycle should reduce obstruction to zero
        for i in range(len(ems)):
            reduced = ems[:i] + ems[i + 1 :]
            reduced_energy = compute_cohomology(reduced, dim=len(ION_ORDER))["obstruction_energy"]
            self.assertLess(reduced_energy, base_energy)


class TestCycleObstruction(unittest.TestCase):
    def test_consistent_chain_no_cycles(self):
        ems = _make_consistent_edge_maps()
        info = compute_cycle_obstruction(ems, dim=len(ION_ORDER))
        self.assertEqual(info["cycle_count"], 0)
        self.assertEqual(info["cycle_obstruction_max"], 0.0)

    def test_inconsistent_cycle_detected(self):
        ems = _make_inconsistent_cycle()
        info = compute_cycle_obstruction(ems, dim=len(ION_ORDER))
        # The 3-node cycle A->B->C->A should be found
        self.assertGreaterEqual(info["cycle_count"], 1)
        self.assertGreater(info["cycle_obstruction_max"], 0.0)


class TestAttachAttrs(unittest.TestCase):
    def test_attrs_attached_to_selected_edges(self):
        ems = _make_consistent_edge_maps()
        edges = [em.edge for em in ems]
        ems_dict = {em.edge.edge_id: em for em in ems}
        attach_cohomology_attrs(edges, ems_dict, dim=len(ION_ORDER))

        for edge in edges:
            attrs = edge.attrs or {}
            self.assertIn("sheaf_h0_dim", attrs)
            self.assertIn("sheaf_h1_dim", attrs)
            self.assertIn("sheaf_obstruction_energy", attrs)
            self.assertIn("sheaf_cycle_count", attrs)


class TestCoboundaryMatrixStructure(unittest.TestCase):
    def test_chain_rank_dim1_correct_h0_h1(self):
        """Chain A->B->C with dim=1 should give rank=2, h0=1, h1=0."""
        from hydrosheaf.sheaf.cohomology import compute_cohomology
        from hydrosheaf.config import Config

        cfg = Config()
        cfg.mixing_endmembers = {}
        cfg.transport_models_enabled = ["evap"]
        dim = 1

        e1 = _make_edge("A->B", "A", "B")
        e1.attrs = {"edge_confidence": 1.0}
        e2 = _make_edge("B->C", "B", "C")
        e2.attrs = {"edge_confidence": 1.0}

        node_values = {"A": [1.0], "B": [1.0], "C": [1.0]}
        ems = list(build_edge_maps([e1, e2], node_values, cfg).values())
        coh = compute_cohomology(ems, dim=dim)

        self.assertEqual(coh["rank_D"], 2)
        self.assertEqual(coh["h0_dim"], 1)  # n_vertex_dofs(3) - rank(2) = 1
        self.assertEqual(coh["h1_dim"], 0)  # n_edge_dofs(2) - rank(2) = 0

    def test_matrix_dimensions(self):
        ems = _make_consistent_edge_maps()
        D = build_coboundary_matrix(ems, dim=len(ION_ORDER))
        n_edges = len(ems)
        n_nodes = len({em.edge.u for em in ems} | {em.edge.v for em in ems})
        self.assertEqual(D.shape[0], n_edges * len(ION_ORDER))
        self.assertEqual(D.shape[1], n_nodes * len(ION_ORDER))

    def test_rhs_shape(self):
        ems = _make_consistent_edge_maps()
        b = build_rhs_vector(ems, dim=len(ION_ORDER))
        self.assertEqual(len(b), len(ems) * len(ION_ORDER))

    def test_empty_edge_maps(self):
        coh = compute_cohomology([], dim=11)
        self.assertEqual(coh["h0_dim"], 0)
        self.assertEqual(coh["h1_dim"], 0)
        self.assertEqual(coh["obstruction_energy"], 0.0)


if __name__ == "__main__":
    unittest.main()
