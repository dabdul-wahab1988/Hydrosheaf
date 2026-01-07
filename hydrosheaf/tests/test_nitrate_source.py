import math
import unittest
import pandas as pd
from unittest.mock import MagicMock, patch

from hydrosheaf.nitrate_source_v2 import (
    fit_robust_stats, infer_node_posteriors, NitrateStats, compute_evidence
)
from hydrosheaf.inference.edge_fit import EdgeResult

class NitrateSourceTests(unittest.TestCase):
    def test_stats_computation(self):
        data = [
            {"NO3": 10.0, "Cl": 1.0, "site_id": "A"},
            {"NO3": 20.0, "Cl": 2.0, "site_id": "B"},
            {"NO3": 100.0, "Cl": 1.0, "site_id": "C"},  # High ratio outlier
        ]
        df = pd.DataFrame(data)
        stats = fit_robust_stats(df)
        self.assertAlmostEqual(stats.ln_no3_cl_median, math.log(10.0), places=1)
        
    def test_evidence_high_no3_cl(self):
        stats = NitrateStats()
        stats.ln_no3_cl_median = math.log(1.0) 
        stats.ln_no3_cl_mad = 0.5
        sample = {"NO3": 100.0, "Cl": 1.0} 
        weights = {"w1_no3_cl": 1.0}
        logit, evidence = compute_evidence(sample, stats, weights)
        self.assertLess(logit, -1.0)
        self.assertIn("NO3/Cl_high_fert", evidence)
        
    @patch("hydrosheaf.nitrate_source_v2.fit_robust_stats")
    def test_denitrification_boost(self, mock_fit):
        # Setup specific stats to ensure non-zero Z
        fixed_stats = NitrateStats()
        fixed_stats.denit_median = 0.0
        fixed_stats.denit_mad = 1.0
        # Also set other medians to avoid noise
        fixed_stats.ln_no3_cl_median = math.log(5.0) # Baseline
        fixed_stats.ln_no3_cl_mad = 1.0
        mock_fit.return_value = fixed_stats

        e = MagicMock(spec=EdgeResult)
        e.u = "u"
        e.v = "v"
        e.z_labels = ["denit"]
        e.z_extents = [2.0] # Strong removal (Z=2)
        e.transport_model = "mix"
        
        nodes_data = [
            {"site_id": "u", "NO3": 50.0, "HCO3": 100.0, "Cl": 10.0, "K": 5.0, "Ca": 20.0, "Mg": 10.0, "Na": 10.0, "SO4": 10.0},
            {"site_id": "v", "NO3": 10.0, "HCO3": 140.0, "Cl": 10.0, "K": 5.0, "Ca": 20.0, "Mg": 10.0, "Na": 10.0, "SO4": 10.0},
        ]
        df = pd.DataFrame(nodes_data).set_index("site_id", drop=False)
        
        overrides = {"weights": {"w5_denitrif": 2.0}, "prior_prob": 0.5}
        
        results = infer_node_posteriors(df, [e], overrides)
        res_v = results["v"]
        
        self.assertGreater(res_v.p_manure, 0.7)
        self.assertIn("denitrif_strong", res_v.top_evidence)
        
    @patch("hydrosheaf.nitrate_source_v2.fit_robust_stats")
    def test_evap_gate(self, mock_fit):
        # Setup stats ensuring NO3/Cl gives High Z
        fixed_stats = NitrateStats()
        fixed_stats.ln_no3_cl_median = math.log(1.0)
        fixed_stats.ln_no3_cl_mad = 0.5
        # Set Default P25
        fixed_stats.d_excess_p25 = 10.0
        mock_fit.return_value = fixed_stats
        
        sample = {
            "NO3": 100.0, "Cl": 1.0, # Ratio 100, ln=4.6, Z high
            "d_excess": 5.0, # < 10.0 -> Activates flag
            "site_id": "A",
            "Ca":1,"Mg":1,"Na":1,"K":1,"HCO3":1,"SO4":1 
        } 
        df = pd.DataFrame([sample]).set_index("site_id", drop=False)
        
        # Override threshold just to be sure (it overwrites stats.d_excess_p25)
        overrides = {
            "weights": {"w1_no3_cl": 10.0}, 
            "evap_gate_factor": 0.1,
            "nitrate_source_d_excess_p25": 10.0
        }
        
        results = infer_node_posteriors(df, [], overrides)
        res = results["A"]
        
        self.assertIn("low_d_excess", res.gating_flags)
        
        # Compare with ungated
        overrides_ungated = {
            "weights": {"w1_no3_cl": 10.0}, 
            "evap_gate_factor": 1.0,
            "nitrate_source_d_excess_p25": 10.0
        }
        results_ungated = infer_node_posteriors(df, [], overrides_ungated)
        
        self.assertLess(abs(res.logit_score), abs(results_ungated["A"].logit_score))

if __name__ == "__main__":
    unittest.main()
