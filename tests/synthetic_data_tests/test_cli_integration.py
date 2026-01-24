import unittest
import subprocess
import sys
import json
from pathlib import Path
import tempfile
import os
import pandas as pd

SYNTHETIC_DIR = Path(__file__).parents[2] / "data/synthetic"


class CLIIntegrationTests(unittest.TestCase):
    def test_cli_basic_run(self):
        """Run the full CLI pipeline on synthetic data and check output."""

        # Temporary output file
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
            output_path = tmp.name

        # Create temporary samples CSV with proper column names (site_id, standard ion names)
        with tempfile.NamedTemporaryFile(
            suffix="_samples.csv", delete=False, mode="w", newline=""
        ) as samples_tmp:
            samples_path = samples_tmp.name

        # Create temporary edges CSV with proper column names (u, v instead of from_station, to_station)
        with tempfile.NamedTemporaryFile(
            suffix="_edges.csv", delete=False, mode="w", newline=""
        ) as edges_tmp:
            edges_path = edges_tmp.name

        try:
            # Prepare samples CSV with required columns
            chem_df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
            # Add site_id from station_code
            chem_df["site_id"] = chem_df["station_code"]
            chem_df["sample_id"] = chem_df["station_code"] + "_" + chem_df["event_code"]
            # Map column names: Ca_mg_L -> Ca, etc.
            ion_mapping = {
                "Ca_mg_L": "Ca",
                "Mg_mg_L": "Mg",
                "Na_mg_L": "Na",
                "K_mg_L": "K",
                "Cl_mg_L": "Cl",
                "SO4_mg_L": "SO4",
                "HCO3_mg_L": "HCO3",
                "NO3_mg_L": "NO3",
                "EC_uS_cm": "EC",
                "TDS_mg_L": "TDS",
            }
            for old_col, new_col in ion_mapping.items():
                if old_col in chem_df.columns:
                    chem_df[new_col] = chem_df[old_col]
            # Add missing ions with zeros
            for ion in ["F", "Fe", "PO4"]:
                if ion not in chem_df.columns:
                    chem_df[ion] = 0.0
            chem_df.to_csv(samples_path, index=False)

            # Prepare edges CSV with correct column names
            edges_df = pd.read_csv(SYNTHETIC_DIR / "network_edges.csv")
            edges_df["u"] = edges_df["from_station"]
            edges_df["v"] = edges_df["to_station"]
            edges_df.to_csv(edges_path, index=False)

            # Construct command
            cmd = [
                sys.executable,
                "-m",
                "hydrosheaf.cli",
                "--samples",
                samples_path,
                "--output",
                output_path,
                "--edges",
                edges_path,
                "--missing-policy",
                "impute_zero",
                "--unit",
                "mg/L",
            ]

            # Run
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print("CLI stdout:", result.stdout)
                print("CLI stderr:", result.stderr)

            self.assertEqual(
                result.returncode, 0, f"CLI execution failed: {result.stderr}"
            )

            # Verify Output
            self.assertTrue(os.path.exists(output_path), "Output JSON not created")

            if os.path.exists(output_path):
                with open(output_path, "r") as f:
                    data = json.load(f)

                # Check basic structure of results
                # Should be a list of EdgeResult dicts
                self.assertIsInstance(data, list)
                self.assertGreater(len(data), 0)
                first = data[0]
                self.assertIn("edge_id", first)
                self.assertIn("transport_model", first)

        finally:
            for path in [output_path, samples_path, edges_path]:
                if os.path.exists(path):
                    try:
                        os.remove(path)
                    except OSError:
                        pass


if __name__ == "__main__":
    unittest.main()
