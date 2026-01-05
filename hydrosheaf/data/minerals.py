"""Standardized mineral library for geochemical reactions."""

from typing import Dict, Mapping

# Type alias for reaction stoichiometry: ion -> coefficient
Stoich = Mapping[str, float]

# Complete Mineral Library (IUPAC/Phreeqc Compliant Standards)
MINERAL_LIBRARY: Dict[str, Stoich] = {
    # Carbonates
    "calcite": {"Ca": 1, "HCO3": 2},  # CaCO3 (Alkalinity balanced: Ca + CO3 + H = Ca + HCO3)
    "aragonite": {"Ca": 1, "HCO3": 2},  # CaCO3 polymorph
    "dolomite": {"Ca": 1, "Mg": 1, "HCO3": 4},  # CaMg(CO3)2
    "magnesite": {"Mg": 1, "HCO3": 2},  # MgCO3
    "siderite": {"Fe": 1, "HCO3": 2},  # FeCO3

    # Evaporites
    "gypsum": {"Ca": 1, "SO4": 1},  # CaSO4:2H2O
    "anhydrite": {"Ca": 1, "SO4": 1},  # CaSO4
    "halite": {"Na": 1, "Cl": 1},  # NaCl
    "sylvite": {"K": 1, "Cl": 1},  # KCl
    "fluorite": {"Ca": 1, "F": 2},  # CaF2
    "carnallite": {"K": 1, "Mg": 1, "Cl": 3},  # KMgCl3:6H2O
    "apatite": {"Ca": 5, "PO4": 3, "F": 1},    # Fluorapatite Ca5(PO4)3F
    "goethite": {"Fe": 1},                     # FeOOH (Iron sink for pyrite oxidation)

    # Silicates (Weathering Proxies)
    # Modeled as pure water weathering reactions usually transforming to Kaolinite/Clay
    # Albite: 2 NaAlSi3O8 + 2 CO2 + 11 H2O -> Al2Si2O5(OH)4 (Kaolinite) + 2 Na+ + 2 HCO3- + 4 H4SiO4
    # Simplified Ion Proxy: Na:1, HCO3:1
    "albite": {"Na": 1, "HCO3": 1},
    
    # Anorthite: CaAl2Si2O8 + 2 CO2 + 3 H2O -> Al2Si2O5(OH)4 + Ca++ + 2 HCO3-
    # Simplified Ion Proxy: Ca:1, HCO3:2
    "anorthite": {"Ca": 1, "HCO3": 2},
    
    # K-Feldspar (Microcline/Orthoclase): 2 KAlSi3O8 + ...
    # Simplified Ion Proxy: K:1, HCO3:1
    "k_feldspar": {"K": 1, "HCO3": 1},
    "microcline": {"K": 1, "HCO3": 1},
    
    # Olivine (Forsterite): Mg2SiO4 + 4CO2 + 4H2O -> 2Mg++ + 4HCO3- + H4SiO4
    "forsterite": {"Mg": 2, "HCO3": 4},

    # Crystalline Basement Weathering (Granite/Gneiss)
    # Biotite: K(Mg,Fe)3AlSi3O10(OH)2 -> Releases K, Mg, Fe, and often Fluoride
    # Simplified Stoichiometry: 1 K, 1.5 Mg, 1.5 Fe, 0.2 F, 7.2 HCO3
    "biotite": {"K": 1, "Mg": 1.5, "Fe": 1.5, "F": 0.2, "HCO3": 7.2},
    # Chlorite (Weathering product/Metamorphic): (Mg,Fe)5Al2Si3O10(OH)8
    # Simplified Stoichiometry: 3 Mg, 2 Fe, 10 HCO3
    "chlorite": {"Mg": 3, "Fe": 2, "HCO3": 10},

    # Sulfides & Redox
    # Pyrite (Aerobic Oxidation): FeS2 + ... -> Fe(OH)3 + 2 SO4-- + 4 H+
    # Acid attacks carbonate buffer: 4 H+ + 4 CaCO3 -> 4 Ca++ + 8 HCO3- (Net: 4Ca, 8HCO3, 2SO4?)
    # Simplified Model without coupling: SO4: 2 (assuming Fe precipitates as Hydroxide)
    # Or coupled with alkalinity consumption if explicit H+ not tracked.
    # Current Hydrosheaf model often separates acidity. Here we assume Fe precipitates.
    "pyrite_oxidation_aerobic": {"SO4": 2, "Fe": 1},
    
    # Pyrite (Denitrification Coupled):
    # 5 FeS2 + 14 NO3- + 4 H+ -> 7 N2 + 10 SO4-- + 5 Fe++ + 2 H2O
    # Stoich per mole Pyrite: SO4: 2, NO3: -2.8, Fe: 1
    "pyrite_oxidation_denit": {"SO4": 2, "NO3": -2.8, "Fe": 1},

    # Generic/Legacy Proxies (for backward compatibility)
    "NaSil": {"Na": 1, "HCO3": 1},
    "CaMgSil": {"Ca": 1, "Mg": 1, "HCO3": 4},  # Updated to match typical Anorthite/Pyroxene mix
}

def get_mineral_stoich(name: str) -> Stoich:
    """Retrieve stoichiometry for a named mineral."""
    normalized = name.lower().replace(" ", "_")
    if normalized in MINERAL_LIBRARY:
        return MINERAL_LIBRARY[normalized]
    raise ValueError(f"Mineral '{name}' not found in library.")
