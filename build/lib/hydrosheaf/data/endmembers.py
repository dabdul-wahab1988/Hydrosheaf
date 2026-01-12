"""Endmember JSON parsing."""

import json
from typing import Dict, List, Mapping, Tuple

from .units import mgL_to_mmolL, meqL_to_mmolL


def _convert_value(value: float, ion: str, units: str) -> float:
    if units == "mg/L":
        return mgL_to_mmolL(value, ion)
    if units == "meq/L":
        return meqL_to_mmolL(value, ion)
    return float(value)


def load_endmembers_json(path: str) -> Tuple[Dict[str, List[float]], Dict[str, object]]:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    meta = payload.get("meta") or {}
    units = meta.get("units", "mg/L")
    ion_order = meta.get("ion_order") or []
    if units not in {"mg/L", "mmol/L", "meq/L"}:
        raise ValueError("meta.units must be mg/L, mmol/L, or meq/L")
    if len(ion_order) != 8:
        raise ValueError("meta.ion_order must list 8 ions")

    endmembers_list = payload.get("endmembers") or []
    if not endmembers_list:
        raise ValueError("endmembers list is empty")

    endmembers: Dict[str, List[float]] = {}
    for entry in endmembers_list:
        end_id = entry.get("id")
        if not end_id:
            raise ValueError("endmember id is required")
        composition = entry.get("composition") or {}
        if any(ion not in composition for ion in ion_order):
            missing = [ion for ion in ion_order if ion not in composition]
            raise ValueError(f"endmember '{end_id}' missing ions: {missing}")
        vector = [
            _convert_value(float(composition[ion]), ion, units)
            for ion in ion_order
        ]
        endmembers[str(end_id)] = vector

    return endmembers, meta
