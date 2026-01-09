"""
Kinetic rate law templates and default parameters for reactive transport.
"""

from typing import Dict
import math

from . import KineticParameters


# PHREEQC RATES block templates
RATE_LAW_TEMPLATES = {
    "calcite": """
Calcite
    -start
    10 SI_cal = SI("Calcite")
    20 k = PARM(1)  # rate constant mol/m²/s
    30 A = PARM(2)  # surface area m²/L
    40 rate = k * A * (1 - 10^SI_cal)
    50 moles = rate * TIME
    60 SAVE moles
    -end
""",
    "dolomite": """
Dolomite
    -start
    10 SI_dol = SI("Dolomite")
    20 k = PARM(1)
    30 A = PARM(2)
    40 rate = k * A * (1 - 10^(0.5*SI_dol))
    50 moles = rate * TIME
    60 SAVE moles
    -end
""",
    "gypsum": """
Gypsum
    -start
    10 SI_gyp = SI("Gypsum")
    20 k = PARM(1)
    30 A = PARM(2)
    40 rate = k * A * (1 - 10^SI_gyp)
    50 moles = rate * TIME
    60 SAVE moles
    -end
""",
    "halite": """
Halite
    -start
    10 SI_hal = SI("Halite")
    20 k = PARM(1)
    30 A = PARM(2)
    40 rate = k * A * (1 - 10^SI_hal)
    50 moles = rate * TIME
    60 SAVE moles
    -end
""",
    "fluorite": """
Fluorite
    -start
    10 SI_flu = SI("Fluorite")
    20 k = PARM(1)
    30 A = PARM(2)
    40 rate = k * A * (1 - 10^SI_flu)
    50 moles = rate * TIME
    60 SAVE moles
    -end
""",
    "albite": """
Albite
    -start
    10 SI_alb = SI("Albite")
    20 k = PARM(1)
    30 A = PARM(2)
    40 IF SI_alb > 0 THEN GOTO 100
    50 rate = k * A * (1 - 10^SI_alb)
    60 moles = rate * TIME
    70 SAVE moles
    100 END
    -end
""",
    "anorthite": """
Anorthite
    -start
    10 SI_ano = SI("Anorthite")
    20 k = PARM(1)
    30 A = PARM(2)
    40 IF SI_ano > 0 THEN GOTO 100
    50 rate = k * A * (1 - 10^SI_ano)
    60 moles = rate * TIME
    70 SAVE moles
    100 END
    -end
""",
    "pyrite_oxidation_aerobic": """
Pyrite_oxidation
    -start
    # FeS2 + 3.75 O2 + 3.5 H2O -> Fe(OH)3 + 2 SO4 + 4 H+
    10 O2 = TOT("O(0)")
    20 IF O2 < 1e-6 THEN GOTO 100
    30 k = PARM(1)  # rate constant 1/s
    40 A = PARM(2)  # surface area
    50 rate = k * A * O2^0.5
    60 moles = rate * TIME
    70 SAVE moles
    100 END
    -end
""",
    "denitrification": """
Denitrification
    -start
    # 5 CH2O + 4 NO3- -> 2 N2 + 4 HCO3- + CO2 + 3 H2O
    10 NO3 = TOT("N(5)")
    20 IF NO3 < 1e-6 THEN GOTO 100
    30 k = PARM(1)  # first-order rate constant 1/s
    40 rate = -k * NO3
    50 moles = rate * TIME
    60 SAVE moles
    100 END
    -end
""",
}


# Default kinetic parameters from literature
DEFAULT_KINETIC_PARAMS: Dict[str, KineticParameters] = {
    "calcite": KineticParameters(
        reaction_name="calcite",
        rate_constant=1e-6,  # mol/m²/s at 25°C (Plummer et al., 1978)
        surface_area=0.1,  # m²/L (typical for aquifer material)
        activation_energy=41840,  # J/mol
        reference_temp_k=298.15,
    ),
    "dolomite": KineticParameters(
        reaction_name="dolomite",
        rate_constant=1e-8,  # mol/m²/s (slower than calcite)
        surface_area=0.05,  # m²/L
        activation_energy=52000,  # J/mol
        reference_temp_k=298.15,
    ),
    "gypsum": KineticParameters(
        reaction_name="gypsum",
        rate_constant=1e-4,  # mol/m²/s (fast dissolution)
        surface_area=0.1,  # m²/L
        activation_energy=15000,  # J/mol
        reference_temp_k=298.15,
    ),
    "halite": KineticParameters(
        reaction_name="halite",
        rate_constant=1e-3,  # mol/m²/s (very fast)
        surface_area=0.1,  # m²/L
        activation_energy=10000,  # J/mol
        reference_temp_k=298.15,
    ),
    "fluorite": KineticParameters(
        reaction_name="fluorite",
        rate_constant=1e-7,  # mol/m²/s
        surface_area=0.05,  # m²/L
        activation_energy=45000,  # J/mol
        reference_temp_k=298.15,
    ),
    "albite": KineticParameters(
        reaction_name="albite",
        rate_constant=1e-12,  # mol/m²/s (very slow silicate dissolution)
        surface_area=0.1,  # m²/L
        activation_energy=70000,  # J/mol
        reference_temp_k=298.15,
    ),
    "anorthite": KineticParameters(
        reaction_name="anorthite",
        rate_constant=1e-11,  # mol/m²/s
        surface_area=0.1,  # m²/L
        activation_energy=65000,  # J/mol
        reference_temp_k=298.15,
    ),
    "pyrite_oxidation_aerobic": KineticParameters(
        reaction_name="pyrite_oxidation_aerobic",
        rate_constant=1e-9,  # 1/s (oxygen-dependent)
        surface_area=0.05,  # m²/L
        activation_energy=56000,  # J/mol
        reference_temp_k=298.15,
    ),
    "denitrification": KineticParameters(
        reaction_name="denitrification",
        rate_constant=1e-7,  # 1/s (microbially mediated)
        surface_area=1.0,  # conceptual (not surface-limited)
        activation_energy=69000,  # J/mol (Arrhenius for microbial)
        reference_temp_k=298.15,
    ),
}


def get_rate_law_template(reaction_name: str) -> str:
    """
    Get PHREEQC RATES block template for a reaction.

    Parameters
    ----------
    reaction_name : str
        Name of reaction (e.g., "calcite", "dolomite")

    Returns
    -------
    str
        PHREEQC RATES block text

    Raises
    ------
    KeyError
        If reaction_name not in RATE_LAW_TEMPLATES
    """
    if reaction_name not in RATE_LAW_TEMPLATES:
        raise KeyError(f"No rate law template for '{reaction_name}'")
    return RATE_LAW_TEMPLATES[reaction_name]


def get_default_kinetic_params(reaction_name: str) -> KineticParameters:
    """
    Get default kinetic parameters for a reaction.

    Parameters
    ----------
    reaction_name : str
        Name of reaction

    Returns
    -------
    KineticParameters
        Default parameters from literature

    Raises
    ------
    KeyError
        If reaction_name not in DEFAULT_KINETIC_PARAMS
    """
    if reaction_name not in DEFAULT_KINETIC_PARAMS:
        raise KeyError(f"No default kinetic parameters for '{reaction_name}'")
    return DEFAULT_KINETIC_PARAMS[reaction_name]


def apply_temperature_correction(
    params: KineticParameters, temperature_c: float
) -> float:
    """
    Apply Arrhenius temperature correction to rate constant.

    k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))

    Parameters
    ----------
    params : KineticParameters
        Kinetic parameters with reference temperature
    temperature_c : float
        Temperature in Celsius

    Returns
    -------
    float
        Temperature-corrected rate constant
    """
    if params.activation_energy == 0:
        return params.rate_constant

    R = 8.314  # J/mol/K (gas constant)
    T = temperature_c + 273.15  # Convert to Kelvin
    T_ref = params.reference_temp_k

    exponent = -(params.activation_energy / R) * (1.0 / T - 1.0 / T_ref)

    k_T = params.rate_constant * math.exp(exponent)

    return k_T
