"""Parsing helpers for PHREEQC outputs."""

import csv
from typing import Dict, Iterable, List, Mapping, Optional


def _to_float(value: str) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_selected_output(csv_path: str, id_map: Mapping[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = id_map.get(str(row.get("sample_id")), row.get("sample_id"))
            rows.append(
                {
                    "sample_id": sample_id,
                    "pH": _to_float(row.get("pH")),
                    "ionic_strength": _to_float(row.get("I")),
                    "charge_balance": _to_float(row.get("charge_balance")),
                    "si_calcite": _to_float(row.get("SI_Calcite")),
                    "si_dolomite": _to_float(row.get("SI_Dolomite")),
                    "si_gypsum": _to_float(row.get("SI_Gypsum")),
                    "si_anhydrite": _to_float(row.get("SI_Anhydrite")),
                    "si_halite": _to_float(row.get("SI_Halite")),
                    "si_fluorite": _to_float(row.get("SI_Fluorite")),
                    "phreeqc_ok": True,
                }
            )
    return rows
