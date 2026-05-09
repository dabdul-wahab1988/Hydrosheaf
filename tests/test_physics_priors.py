import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.network_fit import estimate_edge_residence_time_days
from hydrosheaf.physics.modpath import (
    NodeCoord,
    load_modpath_endpoint_records,
    load_modpath_pathline_points,
    priors_from_modpath_endpoints,
)
from hydrosheaf.physics.priors import PhysicsPrior, apply_physics_priors
from hydrosheaf.config import Config


class TestPhysicsPriors(unittest.TestCase):
    def test_apply_physics_priors_override(self):
        edges = [Edge(edge_id="A->B", u="A", v="B", attrs={"p_uv": 0.2})]
        priors = [
            PhysicsPrior(u="A", v="B", p_uv=0.9, tt_mean_days=12.0, tt_std_days=3.0)
        ]
        out = apply_physics_priors(edges, priors, mode="override")
        self.assertEqual(len(out), 1)
        attrs = out[0].attrs
        self.assertAlmostEqual(float(attrs["p_uv"]), 0.9)
        self.assertAlmostEqual(float(attrs["edge_confidence"]), 0.9)
        self.assertAlmostEqual(float(attrs["physics_tau_mean_days"]), 12.0)
        self.assertAlmostEqual(float(attrs["edge_residence_time_days"]), 12.0)
        self.assertAlmostEqual(float(attrs["physics_tau_std_days"]), 3.0)

    def test_apply_physics_priors_merge_and_only(self):
        edges = [Edge(edge_id="A->B", u="A", v="B", attrs={})]
        priors = [PhysicsPrior(u="B", v="C", p_uv=0.7, tt_mean_days=5.0)]
        merged = apply_physics_priors(edges, priors, mode="merge")
        self.assertEqual({e.edge_id for e in merged}, {"A->B", "B->C"})
        only = apply_physics_priors(edges, priors, mode="only")
        self.assertEqual({e.edge_id for e in only}, {"B->C"})

    def test_explicit_tau_override_used(self):
        cfg = Config()
        tau = estimate_edge_residence_time_days({"edge_residence_time_days": 42.0}, cfg)
        self.assertAlmostEqual(float(tau or 0.0), 42.0)

    def test_compact_modpath5_endpoint_fallback(self):
        content = "\n".join(
            [
                "@ [ MODPATH 5.0(COMPACT)(TREF=   0.000000E+00 ) ]",
                "1 15098 10.0 0.0 0.0 4.0 0.0 0.0 0.1 130749 1 1 11 0.0",
                "1 16745 10.0 0.0 0.0 6.0 0.0 0.0 0.1 130750 1 1 11 0.0",
            ]
        )
        with TemporaryDirectory() as tmpdir:
            endpoints = Path(tmpdir) / "compact.end"
            endpoints.write_text(content, encoding="utf-8")

            records = load_modpath_endpoint_records(str(endpoints))
            self.assertEqual(len(records), 2)
            self.assertAlmostEqual(records[0].x0, 0.0)
            self.assertAlmostEqual(records[0].x, 10.0)
            self.assertAlmostEqual(records[0].time, 4.0)

            priors = priors_from_modpath_endpoints(
                str(endpoints),
                [
                    NodeCoord("source", 0.0, 0.0, 0.1),
                    NodeCoord("sink", 10.0, 0.0, 0.0),
                ],
                max_match_distance=0.5,
                source="test_compact_modpath5",
            )
            self.assertEqual(len(priors), 1)
            self.assertEqual(priors[0].edge_id(), "source->sink")
            self.assertAlmostEqual(float(priors[0].p_uv or 0.0), 1.0)
            self.assertAlmostEqual(float(priors[0].tt_mean_days or 0.0), 5.0)

    def test_compact_modpath5_pathline_fallback(self):
        content = "\n".join(
            [
                "@ [ MODPATH 5.0(COMPACT)(TREF=   0.000000E+00 ) ]",
                "1 0.0 0.0 0.1 4.0 0.0 130749 1",
                "1 10.0 0.0 0.0 6.0 0.0 15098 1",
                "2 1.0 0.0 0.1 5.0 0.0 130749 1",
            ]
        )
        with TemporaryDirectory() as tmpdir:
            pathlines = Path(tmpdir) / "compact.path"
            pathlines.write_text(content, encoding="utf-8")

            points = load_modpath_pathline_points(str(pathlines))
            self.assertEqual(len(points), 3)
            self.assertEqual(points[0].particle_id, 1)
            self.assertEqual(points[0].cell, 130749)
            self.assertEqual(points[1].sequence, 1)
