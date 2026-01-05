import unittest

from hydrosheaf.config import Config
from hydrosheaf.phreeqc.runner import run_phreeqc


class PhreeqcRunnerTests(unittest.TestCase):
    def test_missing_ph_skips(self):
        samples = [
            {"sample_id": "A", "pH": None},
            {"sample_id": "B", "pH": 7.0},
        ]
        config = Config(phreeqc_enabled=True)
        results = run_phreeqc(samples, config)
        self.assertEqual(results["A"]["skipped_reason"], "missing_pH")
        self.assertFalse(results["A"]["phreeqc_ok"])

    def test_missing_database(self):
        samples = [{"sample_id": "A", "pH": 7.0}]
        config = Config(phreeqc_enabled=True, phreeqc_database="missing.dat")
        results = run_phreeqc(samples, config)
        self.assertEqual(results["A"]["skipped_reason"], "phreeqc_database_missing")


if __name__ == "__main__":
    unittest.main()
