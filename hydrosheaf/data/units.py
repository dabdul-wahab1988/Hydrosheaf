"""Unit conversion helpers and dynamic chemical element registry."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class ChemicalSpecies:
    """Represents a chemical element or species with speciation rules."""

    symbol: str
    name: str
    molar_mass_g_mol: float
    base_valence: float
    charge_equiv: float
    who_guideline_mg_l: Optional[float] = None
    speciation_rule: Optional[Callable[[Optional[float], Optional[float]], float]] = None

    def get_effective_valence(
        self,
        ph: Optional[float] = None,
        do_mg_l: Optional[float] = None,
    ) -> float:
        """Calculate effective charge valence based on pH and dissolved oxygen."""
        if self.speciation_rule is not None:
            return self.speciation_rule(ph, do_mg_l)
        return self.base_valence


def _boron_speciation(ph: Optional[float], _do: Optional[float]) -> float:
    # Continuous Henderson-Hasselbalch speciation:
    # B(OH)3 + H2O <=> B(OH)4^- + H+ (pKa = 9.24)
    if ph is not None:
        alpha = 1.0 / (1.0 + 10.0 ** (9.24 - float(ph)))
        return -alpha
    return 0.0


def _silica_speciation(_ph: Optional[float], _do: Optional[float]) -> float:
    # H4SiO4^0 is neutral (z = 0) across typical groundwater pH (5.5 - 8.5).
    return 0.0


def _arsenic_speciation(ph: Optional[float], do_mg_l: Optional[float]) -> float:
    # Under reducing conditions (DO < 1.0 mg/L), neutral arsenite H3AsO3^0 dominates (valence = 0)
    # Under oxic conditions (DO >= 1.0 mg/L), arsenate H2AsO4^- / HAsO4^2- dominates (pKa2 = 6.96)
    is_reducing = do_mg_l is not None and float(do_mg_l) < 1.0
    if is_reducing:
        return 0.0
    if ph is not None:
        alpha = 1.0 / (1.0 + 10.0 ** (6.96 - float(ph)))
        return -(1.0 + alpha)
    return -1.0


class ChemicalRegistry:
    """Dynamic Element Registry for hydrogeochemical species."""

    def __init__(self) -> None:
        self._species: Dict[str, ChemicalSpecies] = {}
        self._init_defaults()

    def _init_defaults(self) -> None:
        defaults = [
            # Major ions
            ChemicalSpecies("Ca", "Calcium", 40.078, base_valence=2, charge_equiv=2),
            ChemicalSpecies("Mg", "Magnesium", 24.305, base_valence=2, charge_equiv=2),
            ChemicalSpecies("Na", "Sodium", 22.990, base_valence=1, charge_equiv=1),
            ChemicalSpecies("K", "Potassium", 39.098, base_valence=1, charge_equiv=1),
            ChemicalSpecies("HCO3", "Bicarbonate", 61.017, base_valence=-1, charge_equiv=1),
            ChemicalSpecies("Cl", "Chloride", 35.453, base_valence=-1, charge_equiv=1),
            ChemicalSpecies("SO4", "Sulfate", 96.064, base_valence=-2, charge_equiv=2),
            ChemicalSpecies("NO3", "Nitrate", 62.005, base_valence=-1, charge_equiv=1, who_guideline_mg_l=50.0),
            ChemicalSpecies("F", "Fluoride", 18.998, base_valence=-1, charge_equiv=1, who_guideline_mg_l=1.5),
            ChemicalSpecies("Fe", "Iron", 55.845, base_valence=2, charge_equiv=2),
            ChemicalSpecies("PO4", "Phosphate", 94.971, base_valence=-2, charge_equiv=2),
            # Diagnostic Tracers & Trace Elements
            ChemicalSpecies(
                "B",
                "Boron",
                10.811,
                base_valence=0,
                charge_equiv=0,
                who_guideline_mg_l=2.4,
                speciation_rule=_boron_speciation,
            ),
            ChemicalSpecies("Br", "Bromide", 79.904, base_valence=-1, charge_equiv=1),
            ChemicalSpecies(
                "SiO2",
                "Silica",
                60.0843,
                base_valence=0,
                charge_equiv=0,
                speciation_rule=_silica_speciation,
            ),
            ChemicalSpecies("Sr", "Strontium", 87.62, base_valence=2, charge_equiv=2),
            ChemicalSpecies("Mn", "Manganese", 54.938, base_valence=2, charge_equiv=2, who_guideline_mg_l=0.4),
            ChemicalSpecies(
                "As",
                "Arsenic",
                74.922,
                base_valence=0,
                charge_equiv=0,
                who_guideline_mg_l=0.01,
                speciation_rule=_arsenic_speciation,
            ),
            ChemicalSpecies("Li", "Lithium", 6.941, base_valence=1, charge_equiv=1),
            ChemicalSpecies("Ba", "Barium", 137.327, base_valence=2, charge_equiv=2, who_guideline_mg_l=1.3),
        ]
        for sp in defaults:
            self.register(sp)

    def register(self, species: ChemicalSpecies) -> None:
        self._species[species.symbol] = species

    def get(self, symbol: str) -> Optional[ChemicalSpecies]:
        return self._species.get(symbol)

    def __contains__(self, symbol: str) -> bool:
        return symbol in self._species

    def __iter__(self) -> Iterable[str]:
        return iter(self._species)

    @property
    def molar_mass_dict(self) -> Dict[str, float]:
        return {k: sp.molar_mass_g_mol for k, sp in self._species.items()}

    @property
    def charge_equiv_dict(self) -> Dict[str, float]:
        return {k: sp.charge_equiv for k, sp in self._species.items()}

    def get_effective_charge_equiv(
        self,
        symbol: str,
        ph: Optional[float] = None,
        do_mg_l: Optional[float] = None,
    ) -> float:
        sp = self._species.get(symbol)
        if sp is None:
            raise KeyError(f"Unknown chemical species: {symbol}")
        return abs(sp.get_effective_valence(ph=ph, do_mg_l=do_mg_l))

    def is_who_exceeded(self, symbol: str, value_mg_l: float) -> Optional[bool]:
        sp = self._species.get(symbol)
        if sp is None or sp.who_guideline_mg_l is None:
            return None
        return float(value_mg_l) > sp.who_guideline_mg_l


# Default global registry singleton
REGISTRY = ChemicalRegistry()

# Synchronized mapping dictionaries for full backwards compatibility
MOLAR_MASS_G_MOL: Dict[str, float] = REGISTRY.molar_mass_dict
CHARGE_EQUIV: Dict[str, float] = REGISTRY.charge_equiv_dict


def _live_species(symbol: str) -> ChemicalSpecies:
    species = REGISTRY.get(symbol)
    if species is None:
        raise KeyError(f"Unknown chemical species: {symbol}")
    return species


def mgL_to_mmolL(value: float, ion: str) -> float:
    molar_mass = _live_species(ion).molar_mass_g_mol
    # g/mol == mg/mmol. value (mg/L) / M (mg/mmol) = mmol/L
    return value / molar_mass


def mmolL_to_mgL(value: float, ion: str) -> float:
    """Convert concentration from mmol/L to mg/L."""
    molar_mass = _live_species(ion).molar_mass_g_mol
    # mmol/L * mg/mmol = mg/L
    return value * molar_mass


def mmolL_to_meqL(
    value: float,
    ion: str,
    ph: Optional[float] = None,
    do_mg_l: Optional[float] = None,
) -> float:
    """Convert concentration from mmol/L to meq/L, with optional pH/DO speciation."""
    if ion in REGISTRY:
        equiv = REGISTRY.get_effective_charge_equiv(ion, ph=ph, do_mg_l=do_mg_l)
        return value * equiv
    if ion not in CHARGE_EQUIV:
        raise KeyError(f"Unknown ion: {ion}")
    return value * abs(CHARGE_EQUIV[ion])


def meqL_to_mmolL(
    value: float,
    ion: str,
    ph: Optional[float] = None,
    do_mg_l: Optional[float] = None,
) -> float:
    """Convert concentration from meq/L to mmol/L, with optional pH/DO speciation."""
    if ion in REGISTRY:
        charge = REGISTRY.get_effective_charge_equiv(ion, ph=ph, do_mg_l=do_mg_l)
    elif ion in CHARGE_EQUIV:
        charge = abs(CHARGE_EQUIV[ion])
    else:
        raise KeyError(f"Unknown ion: {ion}")
    if charge == 0:
        return 0.0
    return value / charge


def row_mgL_to_mmolL(row: Dict[str, float], ion_order: Iterable[str]) -> List[float]:
    return [mgL_to_mmolL(float(row[ion]), ion) for ion in ion_order]


def row_meqL_to_mmolL(values: Iterable[float], ion_order: Iterable[str]) -> List[float]:
    return [meqL_to_mmolL(float(v), ion) for v, ion in zip(values, ion_order)]


def row_mmolL_to_meqL(values: Iterable[float], ion_order: Iterable[str]) -> List[float]:
    return [mmolL_to_meqL(float(v), ion) for v, ion in zip(values, ion_order)]


SPECIES_REGISTRY: ChemicalRegistry = REGISTRY


def get_species_molar_mass(symbol: str) -> float:
    if symbol in REGISTRY:
        sp = REGISTRY.get(symbol)
        if sp is not None:
            return sp.molar_mass_g_mol
    if symbol in MOLAR_MASS_G_MOL:
        return MOLAR_MASS_G_MOL[symbol]
    raise KeyError(f"Unknown chemical species: {symbol}")


def get_species_charge_equivalent(
    symbol: str,
    ph: Optional[float] = None,
    do_mg_l: Optional[float] = None,
) -> float:
    if symbol in REGISTRY:
        sp = REGISTRY.get(symbol)
        if sp is not None:
            return sp.get_effective_valence(ph=ph, do_mg_l=do_mg_l)
    if symbol in CHARGE_EQUIV:
        return float(CHARGE_EQUIV[symbol])
    raise KeyError(f"Unknown chemical species: {symbol}")


def get_who_drinking_water_limit(symbol: str) -> Optional[float]:
    sp = REGISTRY.get(symbol)
    if sp is not None:
        return sp.who_guideline_mg_l
    return None
