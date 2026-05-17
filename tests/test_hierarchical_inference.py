import pytest
from hydrosheaf.config import Config
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.models.reactions import build_reaction_dictionary

class TestHierarchicalInference:
    def test_pure_mixing_no_reactions(self):
        """
        Verify that pure physical mixing is identified without invoking chemical reactions,
        even if a chemical reaction could theoretically explain the data (equifinality check).
        """
        config = Config(lambda_sparse=0.1) # Apply penalty to reactions
        # Standard conservative weights (Cl is high)
        config.conservative_weights = [0.01] * 11
        config.conservative_weights[5] = 1.0 # Cl is conservative anchor
        
        # Setup:
        # Endmember: Cl = 100, Na = 100 (Saline)
        # Upstream: Cl = 10, Na = 10 (Fresh)
        # Downstream: Cl = 55, Na = 55 (Exact 50-50 Mix)
        
        # Ion order: ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]
        # Na is index 2, Cl is index 5.
        
        x_u = [0.0] * 11
        x_u[2] = 10.0 # Na
        x_u[5] = 10.0 # Cl
        
        x_end = [0.0] * 11
        x_end[2] = 100.0 # Na
        x_end[5] = 100.0 # Cl
        
        x_v = [0.0] * 11
        x_v[2] = 55.0 # Na
        x_v[5] = 55.0 # Cl
        
        config.mixing_endmembers = {"Seawater": x_end}
        config.transport_models_enabled = ["mix"]
        
        result = fit_edge(x_u, x_v, config, edge_id="test_mix", u="Fresh", v="Mix")
        
        # Assertions
        assert result.transport_model == "mix"
        assert result.f == pytest.approx(0.5, abs=0.01)
        assert result.endmember_id == "Seawater"
        
        # Crucially: No reactions should be invoked because mixing explains everything perfectly.
        # The L1 penalty on reactions should suppress any noise.
        # Check reaction extents (z)
        total_reaction = sum(abs(z) for z in result.z_extents)
        assert total_reaction < 1e-6, f"Model invoked reactions {total_reaction} when pure mixing was sufficient."

    def test_mixing_plus_reaction(self):
        """
        Verify that the model can identify a superimposed chemical reaction
        after correctly establishing the mixing baseline.
        """
        config = Config(lambda_sparse=0.001)
        config.conservative_weights = [0.01] * 11
        config.conservative_weights[5] = 1.0 # Cl is conservative anchor
        
        # Setup:
        # Endmember: Cl = 100, Na = 100, Ca = 0
        # Upstream: Cl = 10, Na = 10, Ca = 10
        # Transport: 50-50 Mix -> Cl = 55, Na = 55, Ca = 5.
        # Reaction: Dissolve Calcite (CaCO3) -> Adds Ca + HCO3.
        # Let's add 5 mmol/L Calcite.
        # Final Ca = 5 + 5 = 10.
        # Final HCO3 = 0 + 5 = 5.
        
        x_u = [0.0] * 11
        x_u[0] = 10.0 # Ca
        x_u[2] = 10.0 # Na
        x_u[5] = 10.0 # Cl
        
        x_end = [0.0] * 11
        x_end[2] = 100.0 # Na
        x_end[5] = 100.0 # Cl
        
        # Expected after mixing (f=0.5):
        # Cl = 55, Na = 55, Ca = 5.
        
        # Add Calcite dissolution (+5 Ca, +5 HCO3)
        x_v = [0.0] * 11
        x_v[0] = 10.0 # Ca (5 from mix + 5 reaction)
        x_v[2] = 55.0 # Na
        x_v[4] = 5.0  # HCO3
        x_v[5] = 55.0 # Cl
        
        config.mixing_endmembers = {"Brine": x_end}
        config.transport_models_enabled = ["mix"]
        config.active_minerals = ["calcite_closed"]
        config.exchange_enabled = False # Disable exchange to isolate Calcite
        
        result = fit_edge(x_u, x_v, config, edge_id="test_mix_react", u="Fresh", v="Evolved")
        
        # Assertions
        # 1. Transport should be anchored by Cl (which didn't react)
        assert result.transport_model == "mix"
        assert result.f == pytest.approx(0.5, abs=0.01)
    
        # 2. Reaction should be identified
        # Check Calcite extent
        _, labels, *_ = build_reaction_dictionary(config)
        calcite_idx = labels.index("calcite_closed")
        calcite_extent = result.z_extents[calcite_idx]        
        assert calcite_extent == pytest.approx(5.0, abs=0.1)
