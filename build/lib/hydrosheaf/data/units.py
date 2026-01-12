"""Unit conversion helpers."""

from typing import Dict, Iterable, List

MOLAR_MASS_G_MOL: Dict[str, float] = {
    "Ca": 40.078,
    "Mg": 24.305,
    "Na": 22.990,
    "K": 39.0983,
    "HCO3": 61.0168,
    "Cl": 35.45,
    "SO4": 96.06,
    "NO3": 62.0049,
    "F": 19.0,
    "Fe": 55.845,
    "PO4": 94.9714,
}

CHARGE_EQUIV: Dict[str, int] = {
    "Ca": 2,
    "Mg": 2,
    "Na": 1,
    "K": 1,
    "HCO3": 1,
    "Cl": 1,
    "SO4": 2,
    "NO3": 1,
    "F": 1,
    "Fe": 2,    # Assuming Fe(II) for charge balance estimation
    "PO4": 3,   # Assuming PO4(3-)
}


def mgL_to_mmolL(value: float, ion: str) -> float:
    if ion not in MOLAR_MASS_G_MOL:
        raise KeyError(f"Unknown ion: {ion}")
    return value / (MOLAR_MASS_G_MOL[ion] * 1000.0)


def mmolL_to_meqL(value: float, ion: str) -> float:
    if ion not in CHARGE_EQUIV:
        raise KeyError(f"Unknown ion: {ion}")
    return value * abs(CHARGE_EQUIV[ion])


def meqL_to_mmolL(value: float, ion: str) -> float:
    if ion not in CHARGE_EQUIV:
        raise KeyError(f"Unknown ion: {ion}")
    charge = abs(CHARGE_EQUIV[ion])
    if charge == 0:
        return 0.0
    return value / charge


def row_mgL_to_mmolL(row: Dict[str, float], ion_order: Iterable[str]) -> List[float]:
    return [mgL_to_mmolL(float(row[ion]), ion) for ion in ion_order]


def row_meqL_to_mmolL(values: Iterable[float], ion_order: Iterable[str]) -> List[float]:
    return [meqL_to_mmolL(float(v), ion) for v, ion in zip(values, ion_order)]


def row_mmolL_to_meqL(values: Iterable[float], ion_order: Iterable[str]) -> List[float]:
    return [mmolL_to_meqL(float(v), ion) for v, ion in zip(values, ion_order)]
