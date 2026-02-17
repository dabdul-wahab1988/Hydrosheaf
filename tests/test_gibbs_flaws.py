import pytest
from hydrosheaf.models.gibbs import classify_gibbs_dominance, gibbs_evaporation_penalty

class TestGibbsFuzzy:
    def test_high_tds_rock_correct_classification(self):
        """
        Verify that high TDS (e.g. from Gypsum dissolution) is correctly 
        classified as 'rock' if Na/Cl ratios are low.
        """
        # TDS = 1500 (High), but low Na and Cl ratios (0.1)
        sample = {
            "Na": 10,
            "Ca": 90,
            "Cl": 10,
            "HCO3": 90,
            "TDS": 1500
        }
        
        dom = classify_gibbs_dominance(sample)
        # CORRECT BEHAVIOR: Returns "rock" because ratios are low despite high TDS
        assert dom == "rock" 
        
        # Penalty should be high (close to 1.0) because P(Evap) is near 0
        pen = gibbs_evaporation_penalty(sample)
        assert pen > 0.9

    def test_smooth_penalty_transition(self):
        """
        Verify that the penalty transition is smooth, not brittle.
        """
        # Ratios high (0.9) to make evaporation possible
        s1 = {"Na": 90, "Ca": 10, "Cl": 90, "HCO3": 10, "TDS": 990}
        s2 = {"Na": 90, "Ca": 10, "Cl": 90, "HCO3": 10, "TDS": 1010}
        
        p1 = gibbs_evaporation_penalty(s1)
        p2 = gibbs_evaporation_penalty(s2)
        
        # Difference should be small (continuous), not a 0.3 jump
        assert abs(p1 - p2) < 0.05
        # And both should be low penalties because ratios are high and TDS is near 1000
        assert p1 < 0.6
        assert p2 < 0.6

