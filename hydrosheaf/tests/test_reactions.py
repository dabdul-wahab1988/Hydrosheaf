import unittest

from hydrosheaf.models.reactions import fit_reactions


class ReactionTests(unittest.TestCase):
    def test_fit_single_reaction(self):
        reaction_matrix = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        residual = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        fit = fit_reactions(
            residual,
            reaction_matrix,
            weights=[1.0] * 8,
            lambda_l1=0.0,
        )
        self.assertAlmostEqual(fit.extents[0], 2.0, places=4)
        self.assertTrue(all(abs(v) < 1e-4 for v in fit.residual))
        self.assertAlmostEqual(fit.l1_norm, 2.0, places=4)


if __name__ == "__main__":
    unittest.main()
