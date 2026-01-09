import unittest

from hydrosheaf.models.transport import fit_evaporation, fit_mixing


class TransportTests(unittest.TestCase):
    def test_fit_evaporation(self):
        x_u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        gamma = 1.5
        x_v = [gamma * x for x in x_u]
        gamma_hat, residual, norm = fit_evaporation(x_u, x_v, [1.0] * 8)
        self.assertAlmostEqual(gamma_hat, gamma, places=6)
        self.assertAlmostEqual(norm, 0.0, places=6)
        self.assertTrue(all(abs(r) < 1e-6 for r in residual))

    def test_fit_mixing(self):
        x_u = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        x_end = [2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0]
        f = 0.25
        x_v = [u + f * (e - u) for u, e in zip(x_u, x_end)]
        f_hat, residual, norm = fit_mixing(x_u, x_v, x_end, [1.0] * 8)
        self.assertAlmostEqual(f_hat, f, places=6)
        self.assertAlmostEqual(norm, 0.0, places=6)
        self.assertTrue(all(abs(r) < 1e-6 for r in residual))


if __name__ == "__main__":
    unittest.main()
