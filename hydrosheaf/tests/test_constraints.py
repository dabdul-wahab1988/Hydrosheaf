import unittest

from hydrosheaf.config import Config
from hydrosheaf.models.reactions import build_reaction_dictionary
from hydrosheaf.phreeqc.constraints import build_edge_bounds


class ConstraintTests(unittest.TestCase):
    def test_si_bounds_mapping(self):
        config = Config(si_threshold_tau=0.2)
        _, labels, mineral_mask = build_reaction_dictionary(config)
        phreeqc = {
            "A": {
                "sample_id": "A",
                "phreeqc_ok": True,
                "si_calcite": -1.0,
                "si_dolomite": -0.5,
                "si_gypsum": -1.0,
                "si_halite": -1.0,
                "si_fluorite": -1.0,
            },
            "B": {
                "sample_id": "B",
                "phreeqc_ok": True,
                "si_calcite": -1.0,
                "si_dolomite": 0.5,
                "si_gypsum": -1.0,
                "si_halite": -1.0,
                "si_fluorite": -1.0,
            },
        }
        bounds = build_edge_bounds(phreeqc, [("A", "B")], labels, mineral_mask, config)
        entry = bounds["A->B"]
        calcite_idx = labels.index("calcite")
        dolomite_idx = labels.index("dolomite")
        nitrate_idx = labels.index("NO3src")

        self.assertEqual(entry["lb"][calcite_idx], 0.0)
        self.assertEqual(entry["constraints_active"]["calcite"], "dissolution_only")
        self.assertEqual(entry["ub"][dolomite_idx], 0.0)
        self.assertEqual(entry["constraints_active"]["dolomite"], "precipitation_only")
        self.assertEqual(entry["lb"][nitrate_idx], 0.0)


if __name__ == "__main__":
    unittest.main()
