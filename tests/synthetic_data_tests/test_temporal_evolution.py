import unittest
import pandas as pd
from pathlib import Path
from datetime import datetime

from hydrosheaf.config import Config
from hydrosheaf.api import fit_temporal_edges
from hydrosheaf.temporal import TemporalNode, TemporalEdgeResult, TimeSeriesSample
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.graph.types import Edge

SYNTHETIC_DIR = Path(__file__).parents[2] / "hydrosheaf_synthetic_csv"


class TemporalEvolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chem_df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
        cls.config = Config()
        cls.config.missing_policy = "impute_zero"

        # Build TemporalNodes: Map[station, TemporalNode]
        # TemporalNode contains time-series data for a station
        cls.temporal_nodes = {}

        station_groups = cls.chem_df.groupby("station_code")
        for station, group in station_groups:
            # Sort by date
            group = group.sort_values("collection_date")

            # Convert to TimeSeriesSample objects
            time_series_samples = []
            for idx, (_, row) in enumerate(group.iterrows()):
                timestamp = datetime.strptime(row["collection_date"], "%Y-%m-%d")

                # Build concentrations vector in ion_order
                concentrations = []
                for ion in cls.config.ion_order:
                    col_name = f"{ion}_mg_L"
                    value = 0.0
                    if col_name in row and pd.notnull(row[col_name]):
                        value = mgL_to_mmolL(float(row[col_name]), ion)
                    elif ion in row and pd.notnull(row[ion]):
                        value = mgL_to_mmolL(float(row[ion]), ion)
                    concentrations.append(value)

                sample = TimeSeriesSample(
                    sample_id=f"{station}_{row['event_code']}",
                    node_id=station,
                    timestamp=timestamp,
                    concentrations=concentrations,
                    temperature_c=(
                        row.get("temp_C") if pd.notnull(row.get("temp_C")) else None
                    ),
                    ph=row.get("pH") if pd.notnull(row.get("pH")) else None,
                )
                time_series_samples.append(sample)

            node = TemporalNode(node_id=station, samples=time_series_samples)
            cls.temporal_nodes[station] = node

    def test_temporal_fitting(self):
        """Test residence time estimation for edges over 7 events."""
        # Edge: BH1 -> BH2
        edges = [
            Edge(
                edge_id="BH1->BH2",
                u="BH1",
                v="BH2",
                attrs={"edge_type": "groundwater_flow", "distance": 109.5},
            )
        ]

        # Fit temporal - returns tuple (List[TemporalEdgeResult], Dict[str, float])
        results, residence_overrides = fit_temporal_edges(
            self.temporal_nodes, edges, self.config
        )

        self.assertIsInstance(results, list)
        if len(results) > 0:
            res = results[0]
            self.assertIsInstance(res, TemporalEdgeResult)
            self.assertEqual(res.edge_id, "BH1->BH2")

            # Check for residence time
            self.assertIsNotNone(res.residence_time_days)
            self.assertGreaterEqual(res.residence_time_days, 0)

    def test_seasonal_pattern(self):
        """Check if we have enough events to see seasonality."""
        l1_node = self.temporal_nodes["L1"]  # Lysimeter
        self.assertEqual(len(l1_node.samples), 7, "Expected 7 events for L1")

        # Check date range
        start = l1_node.samples[0].timestamp
        end = l1_node.samples[-1].timestamp
        delta = end - start

        # Feb 2024 to Sep 2025 is ~1.5 years (~550 days)
        self.assertGreater(delta.days, 500)


if __name__ == "__main__":
    unittest.main()
