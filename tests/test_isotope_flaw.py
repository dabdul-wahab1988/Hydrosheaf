import pytest
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.config import Config
from hydrosheaf.isotopes import isotope_penalty

class TestIsotopeMixing:
    def test_isotope_penalty_geometry(self):
        """
        Test that the isotope_penalty function correctly calculates the geometric
        deviation from a mixing line.
        """
        # U = (0, 0)
        # End = (10, 80) (Slope 8, GMWL-like)
        # V = (5, 40) -> On the line. Penalty should be 0.
        
        pen, metrics = isotope_penalty(0, 0, 5, 40, 0, 0, "mix", endmember_iso=(10, 80))
        assert pen == pytest.approx(0.0, abs=1e-6)
        
        # V = (5, 50) -> Off the line (above).
        # Vector U->End is (10, 80). Length = sqrt(100 + 6400) = 80.62
        # Vector U->V is (5, 50).
        # Cross product: 10*50 - 80*5 = 500 - 400 = 100.
        # Deviation = 100 / 80.62 = 1.24
        
        pen, metrics = isotope_penalty(0, 0, 5, 50, 0, 0, "mix", endmember_iso=(10, 80))
        assert pen == pytest.approx(1.2403, abs=1e-3)
        assert metrics["mix_deviation"] == pytest.approx(1.2403, abs=1e-3)

    def test_missing_endmember_bypass(self):
        """
        Reproduction Test:
        If endmember_iso is missing (None), the code currently returns 0 penalty
        even if the point is wildly off any reasonable mixing trajectory.
        """
        # U = (-5, -30)
        # V = (0, 0) (Wildly different, chemically plausible if mixing with seawater,
        # but if we don't check the isotope line, we accept it blindly).
        
        # This should theoretically fail or warn, but currently passes with 0 penalty.
        pen, metrics = isotope_penalty(-5, -30, 0, 0, 8, 10, "mix", endmember_iso=None)
        
        # The flaw: It returns 0.0 because the check is skipped.
        assert pen == 0.0
        assert metrics["mix_deviation"] == 0.0

    def test_edge_fit_integration_flaw(self):
        """
        Reproduction Test in fit_edge:
        Show that fit_edge fails to pass the endmember isotopes even if available
        in the configuration, leading to the bypass demonstrated above.
        """
        config = Config(
            transport_models_enabled=["mix"],
            isotope_enabled=True,
            lmwl_defined=True,
            isotope_weight=1000.0 # Huge weight to ensure rejection if checked
        )
        
        # Fresh Water U
        x_u = [0.0] * 11
        x_u[5] = 10.0 # Cl
        iso_u = {"18O": -10.0, "2H": -70.0} # Classic meteoric
        
        # Seawater Endmember
        x_end = [0.0] * 11
        x_end[5] = 100.0 # Cl
        iso_end = {"18O": 0.0, "2H": 0.0}
        
        # Config setup
        config.mixing_endmembers = {"Seawater": x_end}
        config.mixing_endmembers_isotopes = {"Seawater": (iso_end["18O"], iso_end["2H"])}
        
        # OBS V: Chemically a perfect 50/50 mix
        x_v = [0.0] * 11
        x_v[5] = 55.0 
        
        # OBS V Isotopes: IMPOSSIBLE
        # Instead of (-5, -35), let's say (-5, -10) -> Huge deuterium excess, evaporated?
        # But certainly not on the mixing line between U(-10, -70) and End(0, 0).
        iso_v = {"18O": -5.0, "2H": -10.0} 
        
        # Create full dictionaries
        obs_u = {config.ion_order[i]: x for i, x in enumerate(x_u)}
        obs_u.update(iso_u)
        
        obs_v = {config.ion_order[i]: x for i, x in enumerate(x_v)}
        obs_v.update(iso_v)
        
        result = fit_edge(
            x_u, x_v, config, 
            edge_id="flaw_test", u="U", v="V",
            obs_u=obs_u, obs_v=obs_v
        )
        
        # The Isotope Penalty should now be HUGE because V is off the mixing line.
        # Before fix, it was 0.0. Now we assert it's > 10.0 (config.weight is 1000).
        assert result.isotope_penalty > 10.0

