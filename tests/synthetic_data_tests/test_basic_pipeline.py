import unittest
import pandas as pd
from pathlib import Path
import json

from hydrosheaf.config import Config
from hydrosheaf.inference.edge_fit import fit_edge, EdgeResult
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.data.schema import vector_from_sample
from hydrosheaf.graph.types import Edge

SYNTHETIC_DIR = Path(__file__).parents[2] / "hydrosheaf_synthetic_csv"


class BasicPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load raw data
        cls.chem_df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
        cls.config = Config()

        # Prepare Event E1 data (Dry season)
        e1_data = cls.chem_df[cls.chem_df["event_code"] == "E1"].copy()

        # Convert to list of dicts and transform units
        records = e1_data.to_dict(orient="records")
        # Fix conversion: Map X_mg_L to X and convert mg/L to mmol/L
        cls.samples = []
        for r in records:
            new_r = r.copy()
            for ion in cls.config.ion_order:
                # Check for ion_mg_L column
                col_name = f"{ion}_mg_L"
                if col_name in r and pd.notnull(r[col_name]):
                    # Convert found column
                    new_r[ion] = mgL_to_mmolL(float(r[col_name]), ion)
                elif ion in r and pd.notnull(r[ion]):
                    # If already has short name (unlikely but safe)
                    new_r[ion] = mgL_to_mmolL(float(r[ion]), ion)

            cls.samples.append(new_r)

        # Index by station code
        cls.samples_map = {s["station_code"]: s for s in cls.samples}

    def test_single_edge_vadose_recharge(self):
        """Test fitting a single edge: L1 -> BH2 (Vadose recharge path)."""
        u = self.samples_map.get("L1")
        v = self.samples_map.get("BH2")

        if not u or not v:
            self.skipTest("Missing E1 data for L1 or BH2")

        edge = Edge(
            edge_id="L1->BH2",
            u="L1",
            v="BH2",
            attrs={"edge_type": "vadose_to_groundwater", "distance": 104.6},
        )

        # Run fit
        x_u, _ = vector_from_sample(
            u,
            self.config.ion_order,
            missing_policy="impute_zero",
            detection_policy="half",
        )
        x_v, _ = vector_from_sample(
            v,
            self.config.ion_order,
            missing_policy="impute_zero",
            detection_policy="half",
        )
        result = fit_edge(
            x_u, x_v, self.config, edge_id=edge.edge_id, u=edge.u, v=edge.v
        )

        # Checks
        self.assertIsInstance(result, EdgeResult)
        self.assertIn(result.transport_model, ["evap", "mix"])

        # For vadose to groundwater, mixing is often the model if concentration decreases,
        # or evap if it increases. L1->BH2 might be mixing with groundwater.
        # Just ensure it calculated something valid.
        self.assertGreater(result.objective_score, 0.0)

        # Ensure we have reaction results
        self.assertTrue(len(result.z_labels) > 0)
        self.assertTrue(len(result.z_extents) > 0)

        # Check parameters based on model
        if result.transport_model == "evap":
            self.assertFalse(pd.isna(result.gamma))
        elif result.transport_model == "mix":
            self.assertFalse(pd.isna(result.f))

    def test_single_edge_groundwater_flow(self):
        """Test fitting a single edge: BH1 -> BH2 (Groundwater flow)."""
        u = self.samples_map.get("BH1")
        v = self.samples_map.get("BH2")

        if not u or not v:
            self.skipTest("Missing E1 data for BH1 or BH2")

        edge = Edge(
            edge_id="BH1->BH2",
            u="BH1",
            v="BH2",
            attrs={
                "edge_type": "groundwater_flow",
                "distance": 109.5,
                "delta_head": 0.6,
            },
        )

        # Run fit
        x_u, _ = vector_from_sample(
            u,
            self.config.ion_order,
            missing_policy="impute_zero",
            detection_policy="half",
        )
        x_v, _ = vector_from_sample(
            v,
            self.config.ion_order,
            missing_policy="impute_zero",
            detection_policy="half",
        )
        result = fit_edge(
            x_u, x_v, self.config, edge_id=edge.edge_id, u=edge.u, v=edge.v
        )

        self.assertIsInstance(result, EdgeResult)

        # Verify residuals are reasonable (not huge)
        # Normalized residual should ideally be small
        self.assertLess(
            result.anomaly_norm,
            10.0,
            "Residual norm surprisingly high for synthetic data",
        )


if __name__ == "__main__":
    unittest.main()
