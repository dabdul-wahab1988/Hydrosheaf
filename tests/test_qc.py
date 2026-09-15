import unittest

from hydrosheaf.data.qc import charge_balance_ratio, qc_flags


class QcTests(unittest.TestCase):
    def test_negative_concentration_flag(self):
        values = [1.0, -0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        flags = qc_flags(
            values, ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F"], 0.1
        )
        self.assertIn("negative_concentration", flags)

    def test_charge_balance_includes_potassium(self):
        values = [1.0, 1.0]
        ions = ["K", "Cl"]
        self.assertEqual(charge_balance_ratio(values, ions), 0.0)
        self.assertEqual(qc_flags(values, ions, limit=0.05), [])


if __name__ == "__main__":
    unittest.main()
