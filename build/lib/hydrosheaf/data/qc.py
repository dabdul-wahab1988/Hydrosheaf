"""Quality control utilities."""

from typing import Dict, Iterable, List

CHARGE_EQUIV: Dict[str, int] = {
    "Ca": 2,
    "Mg": 2,
    "Na": 1,
    "HCO3": -1,
    "Cl": -1,
    "SO4": -2,
    "NO3": -1,
    "F": -1,
}


def charge_balance_ratio(values: Iterable[float], ion_order: Iterable[str]) -> float:
    cations = 0.0
    anions = 0.0
    for value, ion in zip(values, ion_order):
        charge = CHARGE_EQUIV.get(ion)
        if charge is None:
            continue
        if charge > 0:
            cations += value * charge
        else:
            anions += value * abs(charge)
    denom = cations + anions
    if denom == 0:
        return 0.0
    return (cations - anions) / denom


def nonnegative(values: Iterable[float]) -> bool:
    return all(v >= 0 for v in values)


def qc_flags(values: Iterable[float], ion_order: Iterable[str], limit: float = 0.1) -> List[str]:
    flags = []
    if not nonnegative(values):
        flags.append("negative_concentration")
    if abs(charge_balance_ratio(values, ion_order)) > limit:
        flags.append("charge_balance")
    return flags
