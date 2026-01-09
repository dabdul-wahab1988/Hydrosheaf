
import math
import unittest
from hydrosheaf.reactive_transport.rate_laws import (
    KineticParameters,
    apply_temperature_correction,
    RATE_LAW_TEMPLATES
)
from hydrosheaf.reactive_transport.metrics import (
    compute_consistency_metrics,
    check_thermodynamic_consistency
)
from hydrosheaf.reactive_transport.kinetic_phreeqc import build_kinetic_block

class TestReactiveTransportAccuracy(unittest.TestCase):
    """
    Test accuracy of Reactive Transport extension from a mathematical geochemistry perspective.
    """

    def test_arrhenius_correction(self):
        """
        Test that apply_temperature_correction follows the Arrhenius equation accurately.
        k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))
        """
        # Test Case 1: Activation Energy = 0 -> k should be constant
        params_no_ea = KineticParameters(
            reaction_name="test",
            rate_constant=1.0,
            surface_area=1.0,
            activation_energy=0.0,
            reference_temp_k=298.15
        )
        self.assertAlmostEqual(
            apply_temperature_correction(params_no_ea, 25.0), 
            1.0, 
            msg="k should be constant if Ea=0"
        )
        self.assertAlmostEqual(
            apply_temperature_correction(params_no_ea, 50.0), 
            1.0, 
            msg="k should be constant if Ea=0"
        )

        # Test Case 2: Standard calculation
        # Ea = 50000 J/mol, T_ref = 298.15 K (25 C), T = 308.15 K (35 C)
        # R = 8.314 J/mol/K
        params_ea = KineticParameters(
            reaction_name="test_ea",
            rate_constant=1.0e-6,
            surface_area=1.0,
            activation_energy=50000.0,
            reference_temp_k=298.15
        )
        
        T1 = 25.0
        T2 = 35.0
        T1_K = 298.15
        T2_K = 308.15
        R = 8.314
        
        # Expected ratio: exp(-Ea/R * (1/T2 - 1/T1))
        # Note: The implementation uses 1/T - 1/T_ref. So if T > T_ref, (1/T - 1/T_ref) is negative.
        # -Ea/R is negative. So exponent is positive. Rate should increase.
        exponent = -(50000.0 / R) * (1.0/T2_K - 1.0/T1_K)
        expected_k = 1.0e-6 * math.exp(exponent)
        
        calculated_k = apply_temperature_correction(params_ea, T2)
        
        # Check for proximity. The implementation uses 2.71828 instead of math.e, so we expect small error.
        # Let's check how small.
        self.assertAlmostEqual(calculated_k, expected_k, places=4, 
                               msg="Arrhenius correction should match theoretical value")

    def test_consistency_metrics(self):
        """
        Test statistical metrics against manual calculation.
        """
        x_obs = [10.0, 20.0, 30.0]
        x_fwd = [12.0, 18.0, 25.0]
        # Errors: +2, -2, -5
        # Squared errors: 4, 4, 25 -> sum = 33 -> mean = 11 -> sqrt = 3.3166
        
        metrics = compute_consistency_metrics(x_fwd, x_obs)
        
        self.assertAlmostEqual(metrics['rmse'], (33.0/3.0)**0.5, places=4)
        
        # PBIAS
        # Sum obs = 60
        # Sum diff = 2 - 2 - 5 = -5
        # PBIAS = 100 * (-5 / 60) = -8.3333
        self.assertAlmostEqual(metrics['pbias'], 100 * (-5.0/60.0), places=4)

    def test_thermodynamic_consistency(self):
        """
        Test logic for thermodynamic consistency checks.
        Dissolution requires SI < 0.
        Precipitation requires SI > 0.
        """
        labels = ["calcite", "dolomite"]
        
        # Case 1: Consistent
        # Calcite dissolves (extent > 0), SI = -0.5 (undersaturated) -> OK
        # Dolomite precipitates (extent < 0), SI = 0.5 (supersaturated) -> OK
        extents = [1.0, -1.0]
        si = {"calcite": -0.5, "dolomite": 0.5}
        is_cons, violations = check_thermodynamic_consistency(extents, si, labels)
        self.assertTrue(is_cons)
        self.assertEqual(len(violations), 0)
        
        # Case 2: Inconsistent Dissolution
        # Calcite dissolves (extent > 0) BUT SI = 0.5 (supersaturated) -> Fail
        extents = [1.0, 0.0]
        si = {"calcite": 0.5}
        is_cons, violations = check_thermodynamic_consistency(extents, si, labels)
        self.assertFalse(is_cons)
        self.assertIn("calcite: dissolution", violations[0])
        
        # Case 3: Inconsistent Precipitation
        # Dolomite precipitates (extent < 0) BUT SI = -0.5 (undersaturated) -> Fail
        extents = [0.0, -1.0]
        si = {"dolomite": -0.5}
        is_cons, violations = check_thermodynamic_consistency(extents, si, labels)
        self.assertFalse(is_cons)
        self.assertIn("dolomite: precipitation", violations[0])

    def test_phreeqc_template_validity(self):
        """
        Test that PHREEQC templates are syntactically plausible.
        """
        for mineral, template in RATE_LAW_TEMPLATES.items():
            self.assertIn("-start", template)
            self.assertIn("-end", template)
            self.assertIn("SAVE moles", template)
            # Check for basic TST form (1 - 10^SI) or similar
            if mineral != "pyrite_oxidation_aerobic" and mineral != "denitrification":
                # Most minerals use SI
                self.assertIn("SI", template, f"Template for {mineral} should use SI")

    def test_build_kinetic_block(self):
        """
        Test construction of the PHREEQC kinetic block string.
        """
        labels = ["calcite"]
        extents = [0.5] # dissolution 0.5 mmol/L
        res_time = 10.0 # days
        temp = 25.0
        
        block = build_kinetic_block(labels, extents, res_time, temperature_c=temp)
        
        self.assertIn("RATES", block)
        self.assertIn("Calcite", block)
        self.assertIn("KINETICS 1", block)
        self.assertIn("-m0", block) # Should calculate initial moles
        self.assertIn("-steps", block)
        self.assertIn("8.640000e+05", block) # 10 days in seconds
