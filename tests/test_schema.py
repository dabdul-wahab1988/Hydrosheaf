import unittest

from hydrosheaf.data.schema import parse_numeric


class SchemaTests(unittest.TestCase):
    def test_parse_detection_limit_half(self):
        self.assertAlmostEqual(parse_numeric("<0.2", "half"), 0.1, places=6)

    def test_parse_detection_limit_zero(self):
        self.assertAlmostEqual(parse_numeric("<0.2", "zero"), 0.0, places=6)

    def test_parse_detection_limit_value(self):
        self.assertAlmostEqual(parse_numeric("<0.2", "value"), 0.2, places=6)

    def test_parse_detection_limit_drop(self):
        self.assertIsNone(parse_numeric("<0.2", "drop"))


if __name__ == "__main__":
    unittest.main()
