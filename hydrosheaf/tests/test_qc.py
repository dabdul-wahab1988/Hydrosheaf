import unittest

from hydrosheaf.data.qc import qc_flags


class QcTests(unittest.TestCase):
    def test_negative_concentration_flag(self):
        values = [1.0, -0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        flags = qc_flags(values, ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F"], 0.1)
        self.assertIn("negative_concentration", flags)


if __name__ == "__main__":
    unittest.main()
