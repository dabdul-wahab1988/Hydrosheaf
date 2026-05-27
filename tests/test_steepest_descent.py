"""Tests for steepest-descent flow-direction priors and MODFLOW head parser."""

import io
import math
import textwrap
from pathlib import Path

import pytest

import numpy as np

from hydrosheaf.config import Config
from hydrosheaf.graph.flow_direction import (
    apply_flow_direction_priors,
    apply_steepest_descent_priors,
    compute_projected_gradient_priors,
    compute_steepest_descent_priors,
)
from hydrosheaf.graph.types import Edge


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_edge(u: str, v: str) -> Edge:
    return Edge(edge_id=f"{u}->{v}", u=u, v=v)


def _chain_samples():
    """Three nodes in a straight downhill chain: A(head=10) -> B(head=5) -> C(head=0).

    The steepest direction from A is directly toward B (and then C).
    Edge A->C is a "shortcut" that also goes downhill but skips B.
    """
    return {
        "A": {"site_id": "A", "x": 0.0, "y": 0.0, "elevation": 10.0, "hydraulic_head": 10.0},
        "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "elevation": 5.0, "hydraulic_head": 5.0},
        "C": {"site_id": "C", "x": 2000.0, "y": 0.0, "elevation": 0.0, "hydraulic_head": 0.0},
    }


# ---------------------------------------------------------------------------
# compute_steepest_descent_priors
# ---------------------------------------------------------------------------


class TestSteepestDescentPriors:
    def test_downhill_aligned_gets_high_score(self):
        """A->B aligned with steepest descent should outscore a lateral shortcut."""
        # A at origin, B downhill in x-direction (main flow), L lateral (y-direction, small drop)
        samples = {
            "A": {"site_id": "A", "x": 0.0,    "y": 0.0,    "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0,    "hydraulic_head": 5.0},
            "L": {"site_id": "L", "x": 0.0,    "y": 1000.0, "hydraulic_head": 9.5},
        }
        edges = [_make_edge("A", "B"), _make_edge("A", "L")]
        config = Config()
        scores = compute_steepest_descent_priors(edges, samples, config)

        # A->B is along steepest descent (cos_theta=1, highest gradient)
        # A->L is perpendicular (cos_theta=0) — should get floor score
        assert scores["A->B"] > scores["A->L"]
        assert scores["A->B"] > 0.7

    def test_uphill_edge_gets_floor_score(self):
        """A reverse-direction edge (uphill) should receive the floor prior."""
        samples = _chain_samples()
        edges = [_make_edge("B", "A")]  # uphill
        config = Config()
        scores = compute_steepest_descent_priors(edges, samples, config)
        assert scores["B->A"] == pytest.approx(0.05)

    def test_missing_head_gives_uninformative_prior(self):
        """Nodes without head data should receive the uninformative 0.5 prior."""
        samples = {
            "X": {"site_id": "X", "x": 0.0, "y": 0.0},
            "Y": {"site_id": "Y", "x": 1.0, "y": 0.0},
        }
        edges = [_make_edge("X", "Y")]
        config = Config()
        scores = compute_steepest_descent_priors(edges, samples, config)
        assert scores["X->Y"] == pytest.approx(0.5)

    def test_scores_in_valid_range(self):
        """All scores must be in [0.05, 0.95]."""
        samples = _chain_samples()
        edges = [_make_edge("A", "B"), _make_edge("A", "C"), _make_edge("B", "A")]
        config = Config()
        scores = compute_steepest_descent_priors(edges, samples, config)
        for eid, score in scores.items():
            assert 0.05 <= score <= 0.95, f"{eid} has out-of-range score {score}"

    def test_lateral_edge_penalised_vs_aligned(self):
        """A lateral edge (perpendicular to flow) should score lower than aligned."""
        samples = {
            "O": {"site_id": "O", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "D": {"site_id": "D", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},   # downhill
            "L": {"site_id": "L", "x": 0.0, "y": 1000.0, "hydraulic_head": 9.0},   # lateral, slight drop
        }
        edges = [_make_edge("O", "D"), _make_edge("O", "L")]
        config = Config()
        scores = compute_steepest_descent_priors(edges, samples, config)
        assert scores["O->D"] > scores["O->L"]


# ---------------------------------------------------------------------------
# apply_steepest_descent_priors — mutation test
# ---------------------------------------------------------------------------


class TestApplySteepestDescentPriors:
    def test_mutates_edge_confidence(self):
        """apply_steepest_descent_priors must write edge_confidence to edge.attrs."""
        samples = _chain_samples()
        edge = _make_edge("A", "B")
        config = Config()
        apply_steepest_descent_priors([edge], samples, config)
        assert "edge_confidence" in edge.attrs
        assert 0.05 <= edge.attrs["edge_confidence"] <= 0.95

    def test_steepest_descent_score_attr_written(self):
        """The raw steepest_descent_score must also be written for diagnostics."""
        samples = _chain_samples()
        edge = _make_edge("A", "B")
        config = Config()
        apply_steepest_descent_priors([edge], samples, config)
        assert "steepest_descent_score" in edge.attrs


# ---------------------------------------------------------------------------
# parse_fhd — lightweight parser test with synthetic ASCII content
# ---------------------------------------------------------------------------


class TestParseFhd:
    def test_parse_single_layer(self, tmp_path: Path):
        """parse_fhd correctly reads a single-layer synthetic .fhd file."""
        from hydrosheaf.physics.modflow_head import parse_fhd

        # Minimal 2×2 layer (NCOL=2, NROW=2, ILAY=1)
        content = textwrap.dedent("""\
              1         1   3.650000E+02   3.650000E+02         2         2         1
             HEAD IN LAYER 1 AT END OF TIME STEP 1 IN STRESS PERIOD 1
              10.5  9.0  8.3  7.1
        """)
        fhd = tmp_path / "test.fhd"
        fhd.write_text(content)
        heads = parse_fhd(fhd)

        assert len(heads) == 4
        assert heads[1] == pytest.approx(10.5)
        assert heads[2] == pytest.approx(9.0)
        assert heads[3] == pytest.approx(8.3)
        assert heads[4] == pytest.approx(7.1)

    def test_parse_two_layers(self, tmp_path: Path):
        """parse_fhd handles multiple layers and assigns sequential cell IDs."""
        from hydrosheaf.physics.modflow_head import parse_fhd

        content = textwrap.dedent("""\
              1         1   365.0         365.0         2         2         1
             HEAD IN LAYER 1 AT END OF TIME STEP 1 IN STRESS PERIOD 1
              10.0  9.0  8.0  7.0
              1         1   365.0         365.0         2         2         2
             HEAD IN LAYER 2 AT END OF TIME STEP 1 IN STRESS PERIOD 1
              6.0  5.0  4.0  3.0
        """)
        fhd = tmp_path / "two_layer.fhd"
        fhd.write_text(content)
        heads = parse_fhd(fhd)

        # Layer 1: cells 1-4; layer 2: cells 5-8
        assert heads[1] == pytest.approx(10.0)
        assert heads[4] == pytest.approx(7.0)
        assert heads[5] == pytest.approx(6.0)
        assert heads[8] == pytest.approx(3.0)

    def test_dry_cells_excluded(self, tmp_path: Path):
        """Dry / inactive cells (head >= 1e28) must not appear in the result."""
        from hydrosheaf.physics.modflow_head import parse_fhd

        content = textwrap.dedent("""\
              1         1   365.0         365.0         2         2         1
             HEAD IN LAYER 1 AT END OF TIME STEP 1 IN STRESS PERIOD 1
              10.0  1e30  8.0  7.0
        """)
        fhd = tmp_path / "dry.fhd"
        fhd.write_text(content)
        heads = parse_fhd(fhd)

        assert 2 not in heads          # dry cell excluded
        assert heads[1] == pytest.approx(10.0)
        assert heads[3] == pytest.approx(8.0)


# ---------------------------------------------------------------------------
# compute_projected_gradient_priors
# ---------------------------------------------------------------------------


class TestProjectedGradientPriors:
    """Tests for the continuous-field projected-gradient priors."""

    def test_aligned_edge_gets_high_score(self):
        """Edge pointing directly along the flow direction gets a near-ceiling score.
        Requires ≥4 edges with diverse projections so the IQR-based scaling activates."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
            "L": {"site_id": "L", "x": 0.0, "y": 1000.0, "hydraulic_head": 9.0},
            "D": {"site_id": "D", "x": -1000.0, "y": 0.0, "hydraulic_head": 8.0},
            "F": {"site_id": "F", "x": 0.0, "y": -1000.0, "hydraulic_head": 7.0},
        }
        edges = [_make_edge("A", "B"), _make_edge("A", "L"), _make_edge("A", "D"), _make_edge("A", "F")]
        # Head drops eastward: dh/dx = -0.005
        gradient_map = {"A": (-0.005, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        assert "A->B" in scores
        # A->B is aligned east with strong gradient → should score above others
        assert scores["A->B"] > scores["A->L"], f"Aligned {scores['A->B']} not > lateral {scores['A->L']}"
        assert scores["A->B"] > scores["A->D"], f"Aligned {scores['A->B']} not > anti-aligned {scores['A->D']}"

    def test_perpendicular_edge_gets_low_score(self):
        """Edge perpendicular to flow direction gets a low score."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "L": {"site_id": "L", "x": 0.0, "y": 1000.0, "hydraulic_head": 9.0},
        }
        edges = [_make_edge("A", "L")]
        # Flow is east → gradient must be negative east
        # Edge A->L points north → perpendicular → proj ≈ 0
        gradient_map = {"A": (-0.005, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        assert "A->L" in scores
        # Perpendicular should score near the median (typically ~0.5) or lower
        assert scores["A->L"] < 0.6, f"Perpendicular edge got {scores['A->L']}, expected <0.6"

    def test_anti_aligned_edge_gets_floor(self):
        """Edge pointing upstream (against gradient) gets the floor score.
        Requires ≥4 edges with diverse projections for IQR-based scaling."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
            "C": {"site_id": "C", "x": 500.0, "y": 0.0, "hydraulic_head": 7.0},
            "L": {"site_id": "L", "x": 0.0, "y": 1000.0, "hydraulic_head": 6.0},
            "S": {"site_id": "S", "x": 0.0, "y": -1000.0, "hydraulic_head": 6.0},
        }
        # B->A is uphill: flow goes east, edge points west → anti-aligned
        gradient_map = {"B": (-0.005, 0.0, 0.0), "A": (-0.005, 0.0, 0.0)}
        edges = [_make_edge("B", "A"), _make_edge("B", "C"), _make_edge("B", "L"), _make_edge("B", "S")]

        # B->C: downhill east, slightly aligned
        # B->L: perpendicular north, ~neutral
        # B->A: uphill west, anti-aligned
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)
        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        assert "B->A" in scores
        assert scores["B->A"] == pytest.approx(0.05)  # floor for anti-aligned

    def test_single_neighbor_uninformative_in_low_diversity(self):
        """A node with only one neighbor and no edge diversity gets uninformative
        prior (0.5). This correctly encodes that the data don't support
        directional discrimination with a single candidate."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 100.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 50.0},
        }
        edges = [_make_edge("A", "B")]
        gradient_map = {"A": (-0.05, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        # Single edge → low diversity → uninformative for all edges
        assert scores["A->B"] == pytest.approx(0.5)

    def test_weak_gradient_produces_moderate_scores(self):
        """When head gradient is nearly flat, aligned edges get near-uninformative scores."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 9.9999},
            "D": {"site_id": "D", "x": 1000.0, "y": 100.0, "hydraulic_head": 9.9998},
        }
        edges = [_make_edge("A", "B"), _make_edge("A", "D")]
        # Very weak eastward gradient: ~1e-7 ft/ft
        gradient_map = {"A": (-1e-7, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        # Weak gradient → scores near 0.5 for the aligned edge
        assert 0.4 < scores["A->B"] < 0.6, (
            f"Weak-aligned edge got {scores['A->B']}, expected near-uniform"
        )
        # The slightly-off-diagonal edge should not be far from 0.5 either
        assert 0.4 < scores["A->D"] < 0.6, (
            f"Off-aligned edge got {scores['A->D']}, expected near-uniform"
        )

    def test_scores_in_valid_range(self):
        """All projected-gradient scores must be in [0.05, 0.95]."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
            "C": {"site_id": "C", "x": 0.0, "y": 1000.0, "hydraulic_head": 6.0},
        }
        edges = [_make_edge("A", "B"), _make_edge("A", "C"), _make_edge("B", "A")]
        gradient_map = {"A": (-0.005, 0.0, 0.0), "B": (-0.005, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        for eid, score in scores.items():
            assert 0.05 <= score <= 0.95, f"{eid} has out-of-range score {score}"

    def test_missing_gradient_gives_uninformative(self):
        """Nodes without gradient data get 0.5 (uninformative)."""
        samples = {
            "X": {"site_id": "X", "x": 0.0, "y": 0.0},
            "Y": {"site_id": "Y", "x": 1.0, "y": 0.0},
        }
        edges = [_make_edge("X", "Y")]
        gradient_map = {}  # empty — no gradient data
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        scores = compute_projected_gradient_priors(edges, samples, config, gradient_map)

        assert scores["X->Y"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# apply_flow_direction_priors — unified dispatch
# ---------------------------------------------------------------------------


class TestApplyFlowDirectionPriors:
    def test_dispatches_to_projected_gradient_when_enabled(self):
        """apply_flow_direction_priors uses projected gradient when
        projected_gradient_enabled=True and gradient_map is provided.
        Requires ≥4 edges for IQR-based scaling to activate."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
            "L": {"site_id": "L", "x": 0.0, "y": 1000.0, "hydraulic_head": 9.0},
            "D": {"site_id": "D", "x": -1000.0, "y": 0.0, "hydraulic_head": 8.0},
            "F": {"site_id": "F", "x": 0.0, "y": -1000.0, "hydraulic_head": 7.0},
        }
        edges = [_make_edge("A", "B"), _make_edge("A", "L"), _make_edge("A", "D"), _make_edge("A", "F")]
        gradient_map = {"A": (-0.005, 0.0, 0.0)}
        config = Config(projected_gradient_enabled=True, projected_gradient_sharpness=10.0)

        apply_flow_direction_priors(edges, samples, config, gradient_map)

        # All edges should have edge_confidence + projected_gradient_score
        for edge in edges:
            assert "edge_confidence" in edge.attrs
            assert "projected_gradient_score" in edge.attrs
        # Aligned edge (A->B) should have highest confidence
        scores = {e.edge_id: e.attrs["edge_confidence"] for e in edges}
        assert scores["A->B"] > scores["A->L"]
        assert scores["A->B"] > scores["A->D"]

    def test_dispatches_to_steepest_descent_when_no_gradient_map(self):
        """apply_flow_direction_priors falls back to steepest descent when
        projected_gradient_enabled=True but no gradient_map is provided."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
        }
        edge = _make_edge("A", "B")
        config = Config(projected_gradient_enabled=True)

        # No gradient_map → should fall through (no-op for projected_gradient,
        # no steepest_descent_enabled → no-op entirely)
        apply_flow_direction_priors([edge], samples, config, gradient_map=None)

        # edge_confidence should NOT have been set (neither path active)
        # This is correct — projected_gradient without gradient_map is a no-op

    def test_backward_compat_steepest_descent(self):
        """Original steepest_descent_enabled path still works."""
        samples = {
            "A": {"site_id": "A", "x": 0.0, "y": 0.0, "hydraulic_head": 10.0},
            "B": {"site_id": "B", "x": 1000.0, "y": 0.0, "hydraulic_head": 5.0},
        }
        edge = _make_edge("A", "B")
        config = Config(steepest_descent_enabled=True)

        apply_steepest_descent_priors([edge], samples, config)

        assert "edge_confidence" in edge.attrs
        assert "steepest_descent_score" in edge.attrs


# ---------------------------------------------------------------------------
# compute_head_gradient — finite-difference gradient on structured grid
# ---------------------------------------------------------------------------


class TestGridGradient:
    def test_uniform_gradient_on_regular_grid(self):
        """compute_head_gradient recovers a known uniform gradient."""
        from hydrosheaf.physics.modflow_head import (
            GridGeometry,
            compute_head_gradient,
        )

        # 5×5 grid, Δx=100, Δy=100, head(x,y) = 100 - 0.01*x - 0.005*y
        grid = GridGeometry(ncol=5, nrow=5, nlay=1, dx=100.0, dy=100.0)
        head_map = {}
        for col in range(5):
            for row in range(5):
                x = col * 100.0
                y = row * 100.0
                cell = row * 5 + col + 1
                head_map[cell] = 100.0 - 0.01 * x - 0.005 * y

        gradient = compute_head_gradient(head_map, grid, sigma=0.0)

        # All interior cells should recover the uniform gradient
        # Interior: col in [1,3], row in [1,3]
        for col in range(1, 4):
            for row in range(1, 4):
                cell = row * 5 + col + 1
                gx, gy, gz = gradient[cell]
                # h(x,y) = 100 - 0.01*x - 0.005*y, so dh/dx = -0.01, dh/dy = -0.005
                assert gx == pytest.approx(-0.01, rel=0.01), (
                    f"Cell {cell}: gx={gx}, expected -0.01"
                )
                assert gy == pytest.approx(-0.005, rel=0.01), (
                    f"Cell {cell}: gy={gy}, expected -0.005"
                )
                assert gz == 0.0

    def test_rotated_grid_gradient(self):
        """Gradient is correctly rotated from grid to world coordinates."""
        from hydrosheaf.physics.modflow_head import (
            GridGeometry,
            compute_head_gradient,
        )

        # 3×3 grid with rotation = 90° (grid x → world -y, grid y → world +x)
        grid = GridGeometry(ncol=3, nrow=3, nlay=1, dx=1.0, dy=1.0, rotation_deg=90.0)
        head_map = {}
        for col in range(3):
            for row in range(3):
                x_grid = float(col)
                y_grid = float(row)
                cell = row * 3 + col + 1
                # Head drops along grid-x: h = 100 - 5*x_grid
                head_map[cell] = 100.0 - 5.0 * x_grid

        gradient = compute_head_gradient(head_map, grid, sigma=0.0)

        # Center cell (col=1, row=1)
        gx, gy, gz = gradient[1 * 3 + 1 + 1]  # cell at col=1, row=1

        # Grid gradient: dh/dx_grid ≈ -5.0, dh/dy_grid ≈ 0
        # Rotation +90°: gx_world = 0*cos90 - (-5)*sin90 = 5
        #                  gy_world = 0*sin90 + (-5)*cos90 = 0
        # Wait, let me think again. dh/dx_grid ≈ (h[2,1] - h[0,1]) / 2
        # h[0,1] = 100, h[1,1] = 95, h[2,1] = 90
        # dh/dx_grid = (90 - 100) / 2 = -5, dh/dy_grid ≈ 0
        # Rotation by θ (=90°): gx_world = gx_grid*cosθ - gy_grid*sinθ = -5*0 - 0*1 = 0
        #                       gy_world = gx_grid*sinθ + gy_grid*cosθ = -5*1 + 0*0 = -5
        # So gradient in world coords points south (-y), which is correct because
        # the "grid-x" axis is rotated 90° → it points north in world coords,
        # and head drops along grid-x (northward in world coords).
        assert abs(gx) < 1e-6, f"gx={gx}, expected near 0"
        assert gy == pytest.approx(-5.0, rel=0.01), f"gy={gy}, expected -5.0"

    def test_dry_cells_handled_gracefully(self):
        """Cells with NaN head (dry) produce no gradient entry, and neighbors
        with at least one valid side still get one-sided gradients."""
        from hydrosheaf.physics.modflow_head import (
            GridGeometry,
            compute_head_gradient,
        )

        # 4x4 grid, dry cell at corner (3,3) so its neighbors still have
        # valid data on their other sides for one-sided differences.
        grid = GridGeometry(ncol=4, nrow=4, nlay=1, dx=1.0, dy=1.0)
        head_map = {}
        for col in range(4):
            for row in range(4):
                cell = row * 4 + col + 1
                if col == 3 and row == 3:
                    continue  # cell 16 (corner) is dry
                head_map[cell] = 100.0 - float(col) - float(row)

        gradient = compute_head_gradient(head_map, grid, sigma=0.0)

        # Cell 16 should be absent (dry)
        assert 16 not in gradient
        # Neighbors of dry cell should still get gradients (one-sided OK)
        assert 12 in gradient  # cell at (3,2), above dry corner
        assert 15 in gradient  # cell at (2,3), left of dry corner
        # Interior cells should all be present
        assert 6 in gradient   # cell at (1,1), interior
