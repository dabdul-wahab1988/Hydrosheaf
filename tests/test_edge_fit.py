import unittest

from hydrosheaf.config import Config
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.models.reactions import build_reaction_dictionary
from hydrosheaf.models.mixing import fit_evaporation


class EdgeFitTests(unittest.TestCase):
    def test_edge_fit_evap_and_halite(self):
        config = Config(lambda_sparse=0.0)
        # Fix: Halite releases Cl, so Cl is not conservative in this synthetic case.
        # We must adjust conservative_weights to ignore Cl (index 5) and rely on others.
        # Default is Cl=1.0, others=0.01. We invert this or just set Cl=0.
        weights = [1.0] * 11
        weights[5] = 0.0  # Ignore Cl
        config.conservative_weights = weights
        
        reaction_matrix, labels, _ = build_reaction_dictionary(config)
        halite_idx = labels.index("halite")
        halite = reaction_matrix[halite_idx]

        x_u = [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0]
        gamma = 1.2
        z_true = 0.5
        x_v = [gamma * x + z_true * h for x, h in zip(x_u, halite)]

        result = fit_edge(x_u, x_v, config, edge_id="e1", u="A", v="B")

        self.assertEqual(result.transport_model, "evap")
        self.assertAlmostEqual(result.gamma, gamma, places=6)
        self.assertAlmostEqual(result.z_extents[halite_idx], z_true, places=2)
        self.assertAlmostEqual(
            sum(result.transport_probabilities.values()), 1.0, places=6
        )

    def test_edge_fit_evap_composite_reactions(self):
        config = Config(lambda_sparse=0.0)
        # Restore uniform conservative weights for this synthetic test
        # because we are simulating Halite dissolution (Cl source),
        # so Cl is NOT conservative here.
        config.conservative_weights = [1.0] * 10
        
        reaction_matrix, labels, _ = build_reaction_dictionary(config)

        halite_idx = labels.index("halite")
        gypsum_idx = labels.index("gypsum")
        halite = reaction_matrix[halite_idx]
        gypsum = reaction_matrix[gypsum_idx]

        x_u = [1.0, 0.5, 0.2, 1.0, 0.3, 0.4, 0.2, 0.1]
        gamma = 1.1
        z_halite = 0.4
        z_gypsum = 0.3
        x_v = [
            gamma * x + z_halite * h + z_gypsum * g
            for x, h, g in zip(x_u, halite, gypsum)
        ]

        result = fit_edge(x_u, x_v, config, edge_id="e2", u="A", v="B")

        self.assertEqual(result.transport_model, "evap")
        # Updated to use conservative weights logic
        transport_weights = getattr(config, "conservative_weights", config.weights)
        gamma_expected, _, _ = fit_evaporation(x_u, x_v, transport_weights)
        self.assertAlmostEqual(result.gamma, gamma_expected, places=6)
        self.assertGreater(result.z_extents[halite_idx], 0.0)

        self.assertGreater(result.z_extents[gypsum_idx], 0.0)
        self.assertLess(result.anomaly_norm, 0.2)
        self.assertAlmostEqual(
            sum(result.transport_probabilities.values()), 1.0, places=6
        )


if __name__ == "__main__":
    unittest.main()
