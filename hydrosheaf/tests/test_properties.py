import unittest

from hydrosheaf.models.transport import fit_evaporation, fit_mixing


class TransportPropertyTests(unittest.TestCase):
    def test_evaporation_scaling_invariance(self):
        x_u = [1.0, 2.0, 3.0, 4.0, 0.5, 1.5, 2.5, 3.5]
        gamma = 1.8
        x_v = [gamma * x for x in x_u]
        gamma_hat, _, _ = fit_evaporation(x_u, x_v, [1.0] * 8)

        scale = 3.5
        x_u_scaled = [scale * x for x in x_u]
        x_v_scaled = [scale * x for x in x_v]
        gamma_scaled, _, _ = fit_evaporation(x_u_scaled, x_v_scaled, [1.0] * 8)

        self.assertAlmostEqual(gamma_hat, gamma_scaled, places=6)

    def test_mixing_scaling_invariance(self):
        x_u = [1.0, 1.5, 2.0, 2.5, 0.5, 0.75, 1.25, 1.0]
        x_end = [2.0, 3.0, 0.0, 4.0, 1.0, 1.5, 0.5, 2.0]
        f = 0.4
        x_v = [u + f * (e - u) for u, e in zip(x_u, x_end)]
        f_hat, _, _ = fit_mixing(x_u, x_v, x_end, [1.0] * 8)

        scale = 2.25
        x_u_scaled = [scale * x for x in x_u]
        x_end_scaled = [scale * x for x in x_end]
        x_v_scaled = [scale * x for x in x_v]
        f_scaled, _, _ = fit_mixing(x_u_scaled, x_v_scaled, x_end_scaled, [1.0] * 8)

        self.assertAlmostEqual(f_hat, f_scaled, places=6)


if __name__ == "__main__":
    unittest.main()
