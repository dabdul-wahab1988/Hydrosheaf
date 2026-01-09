import json
import pathlib
import unittest

from hydrosheaf.data.endmembers import load_endmembers_json


class EndmembersJsonTests(unittest.TestCase):
    def test_load_endmembers_json(self):
        payload = {
            "meta": {
                "units": "mg/L",
                "ion_order": ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F"],
            },
            "endmembers": [
                {
                    "id": "test",
                    "composition": {
                        "Ca": 40.0,
                        "Mg": 24.0,
                        "Na": 23.0,
                        "HCO3": 61.0,
                        "Cl": 35.0,
                        "SO4": 96.0,
                        "NO3": 62.0,
                        "F": 19.0,
                    },
                }
            ],
        }
        path = pathlib.Path(__file__).parent / "fixtures" / "endmembers.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")

        endmembers, meta = load_endmembers_json(str(path))
        self.assertIn("test", endmembers)
        self.assertEqual(meta["units"], "mg/L")
        self.assertEqual(len(endmembers["test"]), 8)


if __name__ == "__main__":
    unittest.main()
