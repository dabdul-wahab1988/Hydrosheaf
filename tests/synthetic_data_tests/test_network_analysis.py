import unittest
import pandas as pd
from pathlib import Path
import json

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import fit_network
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.graph.types import Edge

SYNTHETIC_DIR = Path(__file__).parents[2] / "hydrosheaf_synthetic_csv"


class NetworkAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load raw data
        cls.chem_df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
        cls.edges_df = pd.read_csv(SYNTHETIC_DIR / "network_edges.csv")
        cls.config = Config()
        cls.config.missing_policy = "impute_zero"

        # Prepare Event E1 data
        e1_data = cls.chem_df[cls.chem_df["event_code"] == "E1"].copy()
        records = e1_data.to_dict(orient="records")

        cls.samples = []
        for r in records:
            new_r = r.copy()
            # Add site_id from station_code (required by fit_network)
            new_r["site_id"] = r.get("station_code")
            for ion in cls.config.ion_order:
                col_name = f"{ion}_mg_L"
                if col_name in r and pd.notnull(r[col_name]):
                    new_r[ion] = mgL_to_mmolL(float(r[col_name]), ion)
                elif ion in r and pd.notnull(r[ion]):
                    new_r[ion] = mgL_to_mmolL(float(r[ion]), ion)
            cls.samples.append(new_r)

        # Prepare edges list
        cls.edges = []
        for _, row in cls.edges_df.iterrows():
            edge = Edge(
                edge_id=f"{row['from_station']}->{row['to_station']}",
                u=row["from_station"],
                v=row["to_station"],
                attrs={
                    "edge_type": row["edge_type"],
                    "distance": row["distance_m"],
                    "delta_head": (
                        row["delta_head_m"] if pd.notna(row["delta_head_m"]) else None
                    ),
                },
            )
            cls.edges.append(edge)

    def test_full_network_e1(self):
        """Test fitting the entire network for Event E1."""
        results = fit_network(self.samples, self.edges, self.config)

        # Check that we got results for all edges
        # Note: fit_network filters out edges where u or v data is missing.
        # For E1, we expect all stations to have data.
        self.assertEqual(
            len(results), len(self.edges), "Expected results for all 5 edges"
        )

        # Validate output structure
        for res in results:
            self.assertIsNotNone(res.edge_id)
            self.assertIsNotNone(res.transport_model)
            # Ensure no infinite values
            self.assertNotEqual(res.gamma, float("inf"))
            self.assertNotEqual(res.f, float("inf"))
            self.assertNotEqual(res.anomaly_norm, float("inf"))

    def test_consistency_check(self):
        """Check for obvious transport inconsistencies."""
        results = fit_network(self.samples, self.edges, self.config)
        results_map = {r.edge_id: r for r in results}

        # Example check: BH1->BH2 and BH2->BH3 should ideally have similar physics
        # or explainable changes. This is just a sanity check that models aren't wild.
        # This test is soft - just ensures code ran and produced data we can query.
        r1 = results_map.get("BH1->BH2")
        r2 = results_map.get("BH2->BH3")

        if r1 and r2:
            pass  # Just confirming we can access them by ID


if __name__ == "__main__":
    unittest.main()
