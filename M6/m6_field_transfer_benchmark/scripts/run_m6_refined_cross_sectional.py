"""Run a provenance-checked descriptive audit of the refined field cohorts.

This replacement for the legacy seasonal field branch uses only the completed
Central Region and Upper East Region workbooks. It does not infer a flow graph,
temporal transfer, or field reaction truth from cross-sectional samples.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
FIELD_ROOT = REPO_ROOT / "data" / "FieldData"
DEFAULT_OUTPUT = (
    REPO_ROOT
    / ".codex_work"
    / "runs"
    / "REFINED-CR-UER-20260922-01"
    / "M6_field_input_QA_v2"
)

SOURCES = {
    "central_region": {
        "workbook": Path("CRdata/CentralRegion_completed.xlsx"),
        "mirror": Path("derived/cr_field_integration_dataset.csv"),
        "n_expected": 252,
        "sha256_expected": "6f8c65bbd4691727762f9824977a3e7ba2dc9b4626cbe399a7f10746e17be4e4",
        "fluoride_in_stored_cbe": False,
    },
    "upper_east_region": {
        "workbook": Path("UERdata/compiled UER data_new_completed.xlsx"),
        "mirror": Path("derived/uer_field_integration_dataset.csv"),
        "n_expected": 237,
        "sha256_expected": "a802efa607100269f42b52e428ed43f20cd19cfd46cd0e2d2502c1c10420dff7",
        "fluoride_in_stored_cbe": True,
    },
}

MAJOR_IONS = ("Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3")
ION_SPEC = {
    "Ca": (40.078, 2.0),
    "Mg": (24.305, 2.0),
    "Na": (22.98977, 1.0),
    "K": (39.0983, 1.0),
    "HCO3": (61.0168, -1.0),
    "Cl": (35.453, -1.0),
    "SO4": (96.06, -2.0),
    "NO3": (62.0049, -1.0),
    "F": (18.9984, -1.0),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )


def _read_source(dataset: str, spec: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    workbook = FIELD_ROOT / spec["workbook"]
    mirror = FIELD_ROOT / spec["mirror"]
    if not workbook.is_file() or not mirror.is_file():
        raise FileNotFoundError(
            f"Required approved source or its validation mirror is missing: {workbook} / {mirror}"
        )
    source_sha256 = _sha256(workbook)
    if source_sha256 != spec["sha256_expected"]:
        raise ValueError(
            f"Approved source hash changed for {dataset}: {source_sha256}; "
            f"expected {spec['sha256_expected']}. Review and version the new source before analysis."
        )
    frame = pd.read_excel(workbook, sheet_name="GW_Field_Integration")
    mirror_frame = pd.read_csv(mirror)
    try:
        pd.testing.assert_frame_equal(
            frame.reset_index(drop=True),
            mirror_frame.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=1e-9,
            atol=1e-9,
        )
    except AssertionError as exc:
        raise ValueError(
            f"Completed workbook and derived validation mirror disagree for {dataset}: {exc}"
        ) from exc

    if len(frame) != int(spec["n_expected"]):
        raise ValueError(
            f"Unexpected record count for {dataset}: {len(frame)}; "
            f"expected {spec['n_expected']}"
        )
    for column in ("node_id", "site_id", "sample_no"):
        if column not in frame:
            raise ValueError(f"{dataset} source is missing identifier field {column!r}")
        if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"{dataset} has missing {column} identifiers")
        if frame[column].duplicated().any():
            raise ValueError(f"{dataset} has duplicate {column} identifiers")

    source_info = {
        "source_of_truth": str(workbook.relative_to(REPO_ROOT)),
        "source_sheet": "GW_Field_Integration",
        "source_sha256": source_sha256,
        "validation_mirror": str(mirror.relative_to(REPO_ROOT)),
        "validation_mirror_sha256": _sha256(mirror),
        "workbook_mirror_rowwise_match": True,
        "n_rows": int(len(frame)),
        "n_unique_node_ids": int(frame["node_id"].nunique()),
        "n_unique_site_ids": int(frame["site_id"].nunique()),
        "source_dataset_label_ignored_for_cohort_identity": "dataset" in frame.columns,
    }
    return frame, source_info


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _cbe_percent(frame: pd.DataFrame, *, include_fluoride: bool) -> pd.Series:
    """Independently recalculate CBE from mg/L and declared molar masses."""

    cations = pd.Series(0.0, index=frame.index)
    anions = pd.Series(0.0, index=frame.index)
    required = list(MAJOR_IONS)
    if include_fluoride and "f_mg_L" in frame:
        required.append("F")
    for ion in required:
        column = f"{ion.lower()}_mg_L"
        values = _numeric(frame, column)
        if ion == "F":
            # Reproduce the workbook's documented fluoride contribution
            # convention for reconciliation only; strict 8-ion QC is separate.
            values = values.fillna(0.0)
        mass, charge = ION_SPEC[ion]
        equivalents = values / mass * abs(charge)
        if charge > 0:
            cations = cations.add(equivalents, fill_value=np.nan)
        else:
            anions = anions.add(equivalents, fill_value=np.nan)
    # The Upper East workbook's documented CBE recipe inserts zero for
    # unmeasured fluoride. Preserve that workbook recipe only for reconciliation;
    # the standard eight-major-ion CBE below remains separately reported.
    totals = cations + anions
    result = 100.0 * (cations - anions) / totals.where(totals > 0.0)
    return result


def _recalculate_cbe_without_imputation(frame: pd.DataFrame) -> pd.Series:
    complete = pd.Series(True, index=frame.index)
    cations = pd.Series(0.0, index=frame.index)
    anions = pd.Series(0.0, index=frame.index)
    for ion in MAJOR_IONS:
        values = _numeric(frame, f"{ion.lower()}_mg_L")
        complete &= values.notna() & values.ge(0.0)
        mass, charge = ION_SPEC[ion]
        equivalents = values / mass * abs(charge)
        if charge > 0:
            cations += equivalents.fillna(0.0)
        else:
            anions += equivalents.fillna(0.0)
    totals = cations + anions
    cbe = 100.0 * (cations - anions) / totals.where(totals > 0.0)
    return cbe.where(complete)


def _cohort_outputs(dataset: str, frame: pd.DataFrame, spec: dict[str, Any]) -> tuple[dict[str, Any], pd.DataFrame, list[dict[str, Any]], list[dict[str, Any]]]:
    measured = {ion: _numeric(frame, f"{ion.lower()}_mg_L") for ion in MAJOR_IONS}
    complete_ions = pd.concat(measured, axis=1).notna().all(axis=1)
    valid_ions = pd.concat(measured, axis=1).ge(0.0).all(axis=1)
    complete_valid = complete_ions & valid_ions
    cbe_major8 = _recalculate_cbe_without_imputation(frame)
    cbe_reproduced = _cbe_percent(
        frame,
        include_fluoride=bool(spec["fluoride_in_stored_cbe"]),
    )
    stored_cbe = _numeric(frame, "cbe_percent")
    cbe_frame = pd.DataFrame(
        {
            "dataset": dataset,
            "node_id": frame["node_id"].astype(str),
            "site_id": frame["site_id"].astype(str),
            "complete_major8_panel": complete_ions,
            "nonnegative_major8_panel": valid_ions,
            "cbe_major8_independent_percent": cbe_major8,
            "cbe_source_recipe_reproduced_percent": cbe_reproduced,
            "cbe_workbook_percent": stored_cbe,
            "cbe_source_recipe_abs_difference_points": (cbe_reproduced - stored_cbe).abs(),
            "fluoride_observed": _numeric(frame, "f_mg_L").notna(),
            "cbe_source_recipe_includes_fluoride": bool(spec["fluoride_in_stored_cbe"]),
        }
    )
    cbe_frame["cbe_source_recipe_matches_0_02_points"] = (
        cbe_frame["cbe_source_recipe_abs_difference_points"] <= 0.02
    )

    summary: dict[str, Any] = {
        "dataset": dataset,
        "n_samples": int(len(frame)),
        "n_unique_sites": int(frame["site_id"].nunique()),
        "n_complete_major8": int(complete_ions.sum()),
        "n_valid_nonnegative_major8": int(complete_valid.sum()),
        "n_incomplete_major8": int((~complete_ions).sum()),
        "n_workbook_cbe": int(stored_cbe.notna().sum()),
        "n_independent_major8_cbe": int(cbe_major8.notna().sum()),
        "n_major8_cbe_abs_le_5_percent": int(cbe_major8.abs().le(5.0).sum()),
        "n_major8_cbe_abs_le_10_percent": int(cbe_major8.abs().le(10.0).sum()),
        "median_abs_major8_cbe_percent": (
            float(cbe_major8.abs().median()) if cbe_major8.notna().any() else np.nan
        ),
        "n_source_recipe_cbe_reconciled_within_0_02_points": int(
            cbe_frame["cbe_source_recipe_matches_0_02_points"].sum()
        ),
        "n_fluoride_measured": int(_numeric(frame, "f_mg_L").notna().sum()),
        "n_stable_water_isotope_pairs": int(
            (_numeric(frame, "d18O_permil").notna() & _numeric(frame, "d2H_permil").notna()).sum()
        ),
        "n_tritium_measurements": int(_numeric(frame, "tritium_TU").notna().sum()),
        "n_nitrate_isotope_pairs": int(
            (_numeric(frame, "d15N_NO3_permil_air").notna() & _numeric(frame, "d18O_NO3_permil_VSMOW").notna()).sum()
        ),
        "n_measured_hydraulic_head_field_values": int(_numeric(frame, "hydraulic_head_m").notna().sum()),
        "n_samples_with_sample_date": int(
            frame["sample_date"].notna().sum() if "sample_date" in frame else 0
        ),
        "n_utm_coordinate_pairs": int(
            (_numeric(frame, "utm_easting_m").notna() & _numeric(frame, "utm_northing_m").notna()).sum()
        ),
        "stored_cbe_formula_note": (
            "Eight major ions; no fluoride term."
            if not spec["fluoride_in_stored_cbe"]
            else "Workbook recipe includes fluoride and treats an unmeasured fluoride contribution as zero; see both source-recipe and independently complete eight-major-ion CBE outputs."
        ),
        "tds_source_counts": (
            frame["tds_source"].fillna("missing").astype(str).value_counts().to_dict()
            if "tds_source" in frame
            else {}
        ),
    }

    ion_rows: list[dict[str, Any]] = []
    for ion in MAJOR_IONS + (("F",) if "f_mg_L" in frame else ()):
        values = _numeric(frame, f"{ion.lower()}_mg_L").dropna()
        ion_rows.append(
            {
                "dataset": dataset,
                "ion": ion,
                "unit": "mg/L",
                "n_observed": int(values.size),
                "mean": float(values.mean()) if not values.empty else np.nan,
                "median": float(values.median()) if not values.empty else np.nan,
                "q25": float(values.quantile(0.25)) if not values.empty else np.nan,
                "q75": float(values.quantile(0.75)) if not values.empty else np.nan,
                "maximum": float(values.max()) if not values.empty else np.nan,
            }
        )

    facies_rows: list[dict[str, Any]] = []
    if "facies_label" in frame:
        counts = frame["facies_label"].fillna("unclassified").astype(str).value_counts()
        facies_rows.extend(
            {"dataset": dataset, "facies_label": str(label), "n": int(count)}
            for label, count in counts.items()
        )
    return summary, cbe_frame, ion_rows, facies_rows


def run(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty run directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    source_manifest: dict[str, Any] = {}
    cohort_summaries: list[dict[str, Any]] = []
    cbe_tables: list[pd.DataFrame] = []
    ion_summaries: list[dict[str, Any]] = []
    facies_summaries: list[dict[str, Any]] = []
    for dataset, spec in SOURCES.items():
        frame, source_info = _read_source(dataset, spec)
        source_manifest[dataset] = source_info
        summary, cbe_frame, ion_rows, facies_rows = _cohort_outputs(dataset, frame, spec)
        cohort_summaries.append(summary)
        cbe_tables.append(cbe_frame)
        ion_summaries.extend(ion_rows)
        facies_summaries.extend(facies_rows)

    cbe_all = pd.concat(cbe_tables, ignore_index=True)
    cohort_frame = pd.DataFrame(cohort_summaries)
    cohort_frame.to_csv(output / "cohort_coverage_and_qc.csv", index=False)
    cbe_all.to_csv(output / "sample_charge_balance_reconciliation.csv", index=False)
    pd.DataFrame(ion_summaries).to_csv(output / "major_ion_descriptive_summary.csv", index=False)
    pd.DataFrame(facies_summaries).to_csv(output / "workbook_facies_counts.csv", index=False)

    manifest = {
        "schema": "hydrosheaf.refined-field-input-qa.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_script_sha256": _sha256(Path(__file__)),
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "field_sources_used": list(SOURCES),
        "source_manifest": source_manifest,
        "output_directory": str(output),
        "validation": {
            "derived_csv_used_as_analysis_input": False,
            "derived_csv_used_only_for_rowwise_source_consistency_check": True,
            "chemistry_zero_imputation": False,
            "independent_cbe_uses_complete_eight_major_ion_panel": True,
            "workbook_cbe_recipe_reconciliation_reported_separately": True,
        },
        "unsupported_inferences_not_run": [
            "Northern Ghana workbook or any derivative of it",
            "Lower Anayari or Talensi packages",
            "wet/dry seasonal transfer",
            "cross-sectional temporal hold-forward",
            "flow direction from elevation or chemistry",
            "field reaction-family accuracy without independent reaction labels",
            "pooled Central-versus-Upper-East inferential tests",
        ],
        "claim_boundary": (
            "This is a provenance, completeness, charge-balance reconciliation, "
            "and descriptive-screening run. It is not independent validation of "
            "flow paths, temporal prediction, groundwater age, or reaction truth. "
            "Cohorts are summarized separately."
        ),
    }
    _write_json(output / "run_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(args.output)
    print(json.dumps(report["validation"], indent=2))
    print(f"Wrote refined-cohort QA outputs to {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
