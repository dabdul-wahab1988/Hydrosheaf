import unittest

from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.network_fit import estimate_edge_residence_time_days
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
