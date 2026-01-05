"""PHREEQC execution helpers."""

import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping

from ..config import Config
from .outputs import parse_selected_output
from .templates import build_selected_output_block, build_solution_block


def _default_result(sample_id: str, reason: str) -> Dict[str, object]:
    return {
        "sample_id": sample_id,
        "phreeqc_ok": False,
        "skipped_reason": reason,
        "ionic_strength": None,
        "charge_balance": None,
        "si_calcite": None,
        "si_dolomite": None,
        "si_gypsum": None,
        "si_anhydrite": None,
        "si_halite": None,
        "si_fluorite": None,
        "si_siderite": None,
        "si_apatite": None,
        "si_goethite": None,
        "si_sylvite": None,
    }


def _ensure_sample_id(sample: Mapping[str, object]) -> str:
    sample_id = sample.get("sample_id") or sample.get("site_id")
    if sample_id is None:
        raise ValueError("Each sample must include sample_id or site_id.")
    return str(sample_id)


def _resolve_database_path(config: Config) -> str:
    if not config.phreeqc_database:
        return ""
    path = Path(config.phreeqc_database)
    if path.exists():
        return str(path)
    candidate = Path(__file__).resolve().parent.parent / "databases" / path.name
    if candidate.exists():
        return str(candidate)
    return str(path)


def _parse_float(value: object) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _add_value(
    composition: Dict[str, object],
    key: str,
    value: object,
    suffix: str | None = None,
) -> None:
    numeric = _parse_float(value)
    if numeric is None:
        return
    if suffix:
        composition[key] = f"{numeric} {suffix}"
    else:
        composition[key] = numeric


def _build_solution_composition(
    sample: Mapping[str, object],
    temp_default_c: float,
) -> Dict[str, object]:
    temp_c = sample.get("temp_c")
    if temp_c in (None, ""):
        temp_c = temp_default_c
    composition: Dict[str, object] = {"units": "mmol/L"}  # Hydrosheaf internal units are mmol/L
    _add_value(composition, "temp", temp_c)
    _add_value(composition, "pH", sample.get("pH"))

    for key in ("Ca", "Mg", "Na", "Cl", "F", "Fe"):
        _add_value(composition, key, sample.get(key))

    _add_value(composition, "S(6)", sample.get("SO4"))
    _add_value(composition, "N(5)", sample.get("NO3"))
    _add_value(composition, "P", sample.get("PO4"))
    _add_value(composition, "Alkalinity", sample.get("HCO3"))

    return composition


def run_phreeqc(samples: Iterable[Mapping[str, object]], config: Config) -> Dict[str, Dict[str, object]]:
    sample_list = list(samples)
    results: Dict[str, Dict[str, object]] = {}
    for sample in sample_list:
        sample_id = _ensure_sample_id(sample)
        results[sample_id] = _default_result(sample_id, "phreeqc_disabled")

    if not config.phreeqc_enabled:
        return results

    missing_ph = [sample for sample in sample_list if sample.get("pH") in (None, "")]
    for sample in missing_ph:
        sample_id = _ensure_sample_id(sample)
        results[sample_id] = _default_result(sample_id, "missing_pH")

    valid_samples = [sample for sample in sample_list if sample.get("pH") not in (None, "")]
    if not valid_samples:
        return results

    db_path = _resolve_database_path(config)
    if config.phreeqc_database and not Path(db_path).exists():
        for sample in valid_samples:
            sample_id = _ensure_sample_id(sample)
            results[sample_id] = _default_result(sample_id, "phreeqc_database_missing")
        return results

    if config.phreeqc_mode == "phreeqpython":
        try:
            import phreeqpython  # type: ignore
        except ImportError:
            for sample in valid_samples:
                sample_id = _ensure_sample_id(sample)
                results[sample_id] = _default_result(sample_id, "phreeqpython_unavailable")
            return results
        try:
            if db_path:
                db_file = Path(db_path)
                engine = phreeqpython.PhreeqPython(
                    database=db_file.name,
                    database_directory=db_file.parent,
                )
            else:
                engine = phreeqpython.PhreeqPython()
        except Exception:
            for sample in valid_samples:
                sample_id = _ensure_sample_id(sample)
                results[sample_id] = _default_result(sample_id, "phreeqpython_failed")
            return results

        for sample in valid_samples:
            sample_id = _ensure_sample_id(sample)
            try:
                composition = _build_solution_composition(sample, config.temp_default_c)
                solution = engine.add_solution(composition)
                entry = _default_result(sample_id, "phreeqpython_ok")
                entry["phreeqc_ok"] = True
                entry["skipped_reason"] = None
                entry["ionic_strength"] = solution.I
                entry["si_calcite"] = solution.si("Calcite")
                entry["si_dolomite"] = solution.si("Dolomite")
                entry["si_gypsum"] = solution.si("Gypsum")
                entry["si_anhydrite"] = solution.si("Anhydrite")
                entry["si_halite"] = solution.si("Halite")
                entry["si_fluorite"] = solution.si("Fluorite")
                entry["si_siderite"] = solution.si("Siderite")
                entry["si_apatite"] = solution.si("Fluorapatite")
                entry["si_goethite"] = solution.si("Goethite")
                entry["si_sylvite"] = solution.si("Sylvite")
                results[sample_id] = entry
            except Exception as exc:
                reason = "phreeqpython_failed"
                if "No database is loaded" in str(exc):
                    reason = "phreeqpython_database_unavailable"
                results[sample_id] = _default_result(sample_id, reason)
        return results

    if not config.phreeqc_executable:
        for sample in valid_samples:
            sample_id = _ensure_sample_id(sample)
            results[sample_id] = _default_result(sample_id, "phreeqc_executable_missing")
        return results

    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = f"{temp_dir}/input.pqi"
        output_path = f"{temp_dir}/selected_output.csv"
        id_map: Dict[str, str] = {}
        blocks: List[str] = []
        if db_path:
            blocks.append(f"DATABASE {db_path}")
        for idx, sample in enumerate(valid_samples, start=1):
            sample_id = _ensure_sample_id(sample)
            id_map[str(idx)] = sample_id
            blocks.append(build_solution_block(sample, idx, config.temp_default_c))
        blocks.append(build_selected_output_block())

        with open(input_path, "w", encoding="utf-8") as handle:
            handle.write("\n\n".join(blocks))

        try:
            import subprocess

            subprocess.run(
                [config.phreeqc_executable, input_path, output_path],
                check=True,
                capture_output=True,
            )
            parsed = parse_selected_output(output_path, id_map)
        except Exception:
            for sample in valid_samples:
                sample_id = _ensure_sample_id(sample)
                results[sample_id] = _default_result(sample_id, "phreeqc_failed")
            return results

    for row in parsed:
        sample_id = str(row["sample_id"])
        results[sample_id] = row

    return results
