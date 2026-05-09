"""
Nuclear data and definitions for groundwater tracers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

@dataclass(frozen=True)
class Nuclide:
    """
    Radioactive nuclide properties.
    """
    name: str
    symbol: str
    half_life_years: float
    default_units: str
    decay_constant_per_day: float
    description: str

    @classmethod
    def from_half_life(cls, name: str, symbol: str, half_life_years: float, units: str, desc: str) -> Nuclide:
        # lambda = ln(2) / t_1/2
        # t_1/2 in days = half_life_years * 365.25
        # lambda_day = ln(2) / (half_life_years * 365.25)
        import math
        half_life_days = half_life_years * 365.25
        decay = math.log(2) / half_life_days
        return cls(name, symbol, half_life_years, units, decay, desc)


# Standard groundwater dating nuclides
TRITIUM = Nuclide.from_half_life(
    name="Tritium",
    symbol="3H",
    half_life_years=12.32,
    units="TU",
    desc="Hydrogen-3, produced by cosmic rays and 1950s/60s bomb tests."
)

CARBON14 = Nuclide.from_half_life(
    name="Radiocarbon",
    symbol="14C",
    half_life_years=5730.0,
    units="pmc",
    desc="Carbon-14, used for dating older groundwater (requires correction)."
)

ARGON39 = Nuclide.from_half_life(
    name="Argon-39",
    symbol="39Ar",
    half_life_years=269.0,
    units="pmc",
    desc="Noble gas isotope for intermediate ages (50-1000 years)."
)

KRYPTON85 = Nuclide.from_half_life(
    name="Krypton-85",
    symbol="85Kr",
    half_life_years=10.76,
    units="dpm/cc",
    desc="Noble gas tracer for young groundwater."
)

# Registry of known nuclides
REGISTRY: Dict[str, Nuclide] = {
    "3H": TRITIUM,
    "H3": TRITIUM,
    "tritium": TRITIUM,
    "14C": CARBON14,
    "C14": CARBON14,
    "carbon14": CARBON14,
    "39Ar": ARGON39,
    "Ar39": ARGON39,
    "argon39": ARGON39,
    "85Kr": KRYPTON85,
    "Kr85": KRYPTON85,
}

def get_nuclide(name: str) -> Optional[Nuclide]:
    """Retrieve nuclide data by name or symbol (case-insensitive)."""
    key = name.strip()
    # Direct lookup first
    if key in REGISTRY:
        return REGISTRY[key]
    
    # Normalize
    norm = key.lower().replace("_", "")
    for k, v in REGISTRY.items():
        if k.lower().replace("_", "") == norm:
            return v
    return None
