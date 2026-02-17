import unittest
import numpy as np
from hydrosheaf.sheaf.directed_section import solve_coupled_section, DirectedEdgeMap
from hydrosheaf.graph.types import Edge
from hydrosheaf.isotopes import rayleigh_fractionation
from hydrosheaf.calibration.adapters_iso import WaterIsotopeMixingAdapter, WaterEndmember
from hydrosheaf.reactive_transport.kinetic_phreeqc import build_kinetic_block
from hydrosheaf.reactive_transport.rate_laws import KineticParameters

class TestCritiqueImprovements(unittest.TestCase):

    def test_rayleigh_fractionation(self):
        """Test accuracy of Rayleigh fractionation function."""
        # R = R0 * f^(alpha - 1)
        # alpha = 0.5 (extreme fractionation), f = 0.5
        # R should be R0 * 0.5^(-0.5) = R0 * sqrt(2) = R0 * 1.414...
        
        # Using deltas
        # delta = 0 => R/Rstd = 1.0
        # R_new = 1.0 * (0.5)^(-0.5) = 1.41421356
        # delta_new = (1.41421356 - 1) * 1000 = 414.21356
        
        res = rayleigh_fractionation(0.0, 0.5, 0.5)
        self.assertAlmostEqual(res, 414.21356, places=4)
        
        # Test alpha=1 (no fractionation)
        res = rayleigh_fractionation(10.0, 0.5, 1.0)
        self.assertAlmostEqual(res, 10.0)

    def test_coupled_solver_charge_balance(self):
        """Test that coupled solver reduces charge imbalance."""
        nodes = ["A", "B"]
        species = ["Ca", "Cl"]
        
        # A: Ca=1.0, Cl=2.0 (Balanced)
        # B: Observed Ca=1.0, Cl=1.8 (Unbalanced)
        
        edge = Edge("e1", "A", "B")
        edge_map = DirectedEdgeMap(
            edge=edge,
            alpha=1.0,
            offset=[0.0, 0.0],
            weight=1.0,
            objective=0.0,
            transport_model="conservative",
            endmember_id=None,
            residual_norm=0.0
        )
        
        node_obs = {
            "A": [1.0, 2.0],
            "B": [1.0, 1.8]
        }
        
        # 1. Uncoupled (low/no weight)
        res_uncoupled = solve_coupled_section(
            nodes, [edge_map], node_obs, species_names=species, 
            obs_weight=100.0, 
            charge_balance_weight=0.0
        )
        # Should stay close to observations
        self.assertAlmostEqual(res_uncoupled["B"][0], 1.0, places=2)
        self.assertAlmostEqual(res_uncoupled["B"][1], 1.8, places=2)
        
        # 2. Coupled (high weight)
        res_coupled = solve_coupled_section(
            nodes, [edge_map], node_obs, species_names=species, 
            obs_weight=1.0,
            charge_balance_weight=10000.0
        )
        # Should adjust to balance charge: 2*Ca - 1*Cl = 0  => Cl = 2*Ca
        ca = res_coupled["B"][0]
        cl = res_coupled["B"][1]
        imbalance = 2*ca - cl
        self.assertLess(abs(imbalance), 0.01, "Charge imbalance should be minimized")

    def test_isotope_adapter_evaporation(self):
        """Test WaterIsotopeMixingAdapter with evaporation enabled."""
        endmembers = [
            WaterEndmember("Rain", -10.0, -70.0),
            WaterEndmember("Gw", -5.0, -30.0)
        ]
        
        obs = {
            "Sample1": {"18O": -5.0, "2H": -50.0}
        }
        group_map = {"Sample1": "G1"}
        
        adapter = WaterIsotopeMixingAdapter(endmembers, obs, group_map, allow_evaporation=True)
        params = adapter.get_parameters()
        
        # Check if f_evap parameter exists
        param_names = [p.name for p in params]
        self.assertIn("f_evap_G1", param_names)
        
        # Case: 100% Rain, 20% Evap
        vals = {
            "theta_G1_Rain": 10.0, # High -> 100% Rain
            "f_evap_G1": 0.2       # 20% Evap => f=0.8
        }
        
        res = adapter.run_model(vals)
        pred_18 = res["Sample1_18O"]
        self.assertGreater(pred_18, -10.0)

    def test_kinetic_initial_moles(self):
        """Test that build_kinetic_block respects initial_moles."""
        params = KineticParameters(reaction_name="Calcite", rate_constant=1e-9, surface_area=1.0)
        # Inject initial_moles
        params.initial_moles = 123.456
        
        kinetic_params = {"Calcite": params}
        
        block = build_kinetic_block(
            reaction_labels=["Calcite"],
            extents=[1.0],
            residence_time_days=1.0,
            kinetic_params=kinetic_params
        )
        
        self.assertIn("-m0 1.234560e+02", block)

if __name__ == "__main__":
    unittest.main()
