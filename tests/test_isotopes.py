import unittest

from hydrosheaf.isotopes import compute_d_excess, evaporation_index, fit_lmwl, isotope_penalty


class IsotopeTests(unittest.TestCase):
    def test_compute_d_excess(self):
        self.assertAlmostEqual(compute_d_excess(-2.0, -10.0), 6.0, places=6)

    def test_fit_lmwl(self):
        samples = [
            {"18O": -2.0, "2H": -16.0},
            {"18O": -1.0, "2H": -8.0},
            {"18O": 0.0, "2H": 0.0},
        ]
        a, b = fit_lmwl(samples)
        self.assertAlmostEqual(a, 0.0, places=6)
        self.assertAlmostEqual(b, 8.0, places=6)

    def test_isotope_penalty(self):
        penalty, metrics = isotope_penalty(-2.0, -10.0, -1.0, -7.0, 0.0, 8.0, "evap", d_excess_weight=1.0)
        self.assertGreaterEqual(penalty, 0.0)
        self.assertIn("e_u", metrics)
        self.assertIn("e_v", metrics)
        self.assertIn("d_excess_penalty", metrics)


if __name__ == "__main__":
    unittest.main()
