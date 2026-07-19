import unittest
from unittest import mock

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import fit_network
from hydrosheaf.models.reactions import build_reaction_dictionary
from hydrosheaf.phreeqc.constraints import build_edge_bounds
from hydrosheaf.inference.edge_fit import fit_edge


class ConstraintTests(unittest.TestCase):
    def test_si_bounds_mapping(self):
        config = Config(si_threshold_tau=0.2)
        _, labels, mineral_mask, *_ = build_reaction_dictionary(config)
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
        self.assertEqual(entry["lb"][dolomite_idx], -float("inf"))
        self.assertEqual(entry["ub"][dolomite_idx], float("inf"))
        self.assertEqual(entry["constraints_active"]["dolomite"], "free")
        self.assertEqual(entry["lb"][nitrate_idx], 0.0)

    def test_precipitation_bound_enables_negative_mineral_extent(self):
        config = Config(
            ion_order=["Ca", "HCO3", "Cl"],
            weights=[1.0, 1.0, 1.0],
            conservative_weights=[0.01, 0.01, 1.0],
            active_minerals=["calcite"],
            measured_ions=["Ca", "HCO3", "Cl"],
            exchange_enabled=False,
            gibbs_enabled=False,
            phreeqc_enabled=True,
        )
        _, labels, _, _ = build_reaction_dictionary(
            config, sample={"Ca": 2.0, "HCO3": 4.0, "Cl": 1.0}
        )
        calcite_index = labels.index("calcite")
        lower = [0.0] * len(labels)
        upper = [float("inf")] * len(labels)
        lower[calcite_index] = -float("inf")
        upper[calcite_index] = 0.0
        result = fit_edge(
            [2.0, 4.0, 1.0],
            [1.0, 2.0, 1.0],
            config,
            edge_id="A->B",
            u="A",
            v="B",
            obs_u={"Ca": 2.0, "HCO3": 4.0, "Cl": 1.0},
            obs_v={"Ca": 1.0, "HCO3": 2.0, "Cl": 1.0},
            bounds={
                "lb": lower,
                "ub": upper,
                "constraints_active": {
                    label: (
                        "precipitation_only"
                        if label == "calcite"
                        else "dissolution_only"
                    )
                    for label in labels
                },
                "phreeqc_ok": True,
            },
        )
        fitted_calcite = result.z_extents[result.z_labels.index("calcite")]
        self.assertLess(fitted_calcite, -0.1)
        self.assertGreater(result.thermodynamic_constraints_active_count, 0)

    def test_network_bounds_follow_each_edge_reaction_label_order(self):
        config = Config(
            ion_order=["Ca", "HCO3", "Cl", "SO4", "NO3", "Fe"],
            weights=[1.0] * 6,
            conservative_weights=[0.01, 0.01, 1.0, 0.01, 0.01, 0.01],
            measured_ions=["Ca", "HCO3", "Cl", "SO4", "NO3", "Fe"],
            active_minerals=[],
            reaction_processes_enabled=["denitrification"],
            exchange_enabled=False,
            phreeqc_enabled=True,
            missing_policy="impute_zero",
            uncertainty_method="none",
        )
        samples = [
            {
                "site_id": "A",
                "Ca": 1.0,
                "HCO3": 2.0,
                "Cl": 0.5,
                "SO4": 0.2,
                "NO3": 0.05,
                "Fe": 0.01,
            },
            {
                "site_id": "B",
                "Ca": 1.1,
                "HCO3": 2.1,
                "Cl": 0.5,
                "SO4": 0.2,
                "NO3": 0.30,
                "Fe": 0.01,
            },
            {
                "site_id": "C",
                "Ca": 1.2,
                "HCO3": 2.2,
                "Cl": 0.5,
                "SO4": 0.2,
                "NO3": 0.10,
                "Fe": 0.01,
            },
        ]
        labels_seen = []

        def neutral_bounds(results, edges, labels, mineral_mask, cfg):
            labels_seen.append(tuple(labels))
            edge = list(edges)[0]
            return {
                edge.edge_id: {
                    "lb": [None] * len(labels),
                    "ub": [None] * len(labels),
                    "constraints_active": {},
                    "phreeqc_ok": True,
                }
            }

        phreeqc = {
            sample["site_id"]: {
                "sample_id": sample["site_id"],
                "phreeqc_ok": True,
            }
            for sample in samples
        }
        with mock.patch(
            "hydrosheaf.inference.network_fit.build_edge_bounds",
            side_effect=neutral_bounds,
        ):
            results = fit_network(
                samples,
                [("A", "B"), ("B", "C")],
                config,
                phreeqc_results=phreeqc,
            )

        self.assertEqual(len(results), 2)
        self.assertNotIn("denit", labels_seen[0])
        self.assertIn("denit", labels_seen[1])
        self.assertEqual(tuple(results[0].z_labels), labels_seen[0])
        self.assertEqual(tuple(results[1].z_labels), labels_seen[1])


if __name__ == "__main__":
    unittest.main()
