import pytest
import numpy as np
from hydrosheaf.models.redox import classify_redox, get_redox_constraints
from hydrosheaf.data.minerals import MINERAL_LIBRARY, get_mineral_stoich

def test_redox_classification():
    # Oxic (High NO3)
    assert classify_redox({"NO3": 0.1, "Fe": 0.0}) == "oxic"
    # Reducing (Low NO3, High Fe)
    assert classify_redox({"NO3": 0.0, "Fe": 0.05}) == "reducing"
    # Ambiguous
    assert classify_redox({"NO3": 0.01, "Fe": 0.001}) == "ambiguous"

def test_redox_constraints():
    labels = ["calcite", "pyrite_oxidation_aerobic", "halite"]
    # Reducing sample
    sample_v = {"NO3": 0.0, "Fe": 0.1}
    overrides = get_redox_constraints(sample_v, labels)
    assert "pyrite_oxidation_aerobic" in overrides
    assert overrides["pyrite_oxidation_aerobic"] == (0.0, 0.0)
    
    # Oxic sample - should not have overrides
    sample_v_oxic = {"NO3": 0.5, "Fe": 0.0}
    overrides_oxic = get_redox_constraints(sample_v_oxic, labels)
    assert "pyrite_oxidation_aerobic" not in overrides_oxic

def test_biotite_stoichiometry():
    stoich = get_mineral_stoich("biotite")
    assert stoich["K"] == 1.0
    assert stoich["Mg"] == 1.5
    assert stoich["F"] == 0.2
    assert stoich["Fe"] == 1.5

def test_mass_balance_logic():
    # Simple check that minerals in library have required ions
    for name, stoich in MINERAL_LIBRARY.items():
        assert len(stoich) > 0
        for ion in stoich:
            assert isinstance(ion, str)
            assert isinstance(stoich[ion], (int, float))
