import unittest

from hydrosheaf.models.sheaf import edge_residual
from hydrosheaf.models.transport import evaporation_affine


class SheafTests(unittest.TestCase):
    def test_edge_residual_zero(self):
        x_u = [1.0, 2.0]
        gamma = 2.0
        matrix, offset = evaporation_affine(gamma, len(x_u))
        reaction_matrix = [[1.0, 0.0]]
        reaction_extents = [1.0]
        x_v = [gamma * x_u[0] + 1.0, gamma * x_u[1]]

        residual = edge_residual(x_u, x_v, matrix, offset, reaction_matrix, reaction_extents)
        self.assertAlmostEqual(residual[0], 0.0, places=6)
        self.assertAlmostEqual(residual[1], 0.0, places=6)


if __name__ == "__main__":
    unittest.main()
