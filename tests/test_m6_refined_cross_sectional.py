from __future__ import annotations

import json
from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M6" / "m6_field_transfer_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_m6_refined_cross_sectional as refined  # noqa: E402


def test_completed_workbooks_match_validation_mirrors_and_expected_cohorts() -> None:
    assert set(refined.SOURCES) == {"central_region", "upper_east_region"}
    sizes = {}
    for dataset, spec in refined.SOURCES.items():
        frame, manifest = refined._read_source(dataset, spec)
        sizes[dataset] = len(frame)
        assert manifest["workbook_mirror_rowwise_match"] is True
        assert manifest["n_unique_node_ids"] == len(frame)
    assert sizes == {"central_region": 252, "upper_east_region": 237}


def test_strict_major_ion_cbe_never_zero_fills_missing_rows() -> None:
    frame = pd.DataFrame(
        {
            "ca_mg_L": [40.078, 40.078],
            "mg_mg_L": [24.305, 24.305],
            "na_mg_L": [22.98977, 22.98977],
            "k_mg_L": [39.0983, 39.0983],
            "hco3_mg_L": [122.0336, 122.0336],
            "cl_mg_L": [35.453, 35.453],
            "so4_mg_L": [96.06, 96.06],
            "no3_mg_L": [62.0049, None],
        }
    )
    result = refined._recalculate_cbe_without_imputation(frame)
    assert result.iloc[0] == pytest.approx(0.0, abs=1e-10)
    assert pd.isna(result.iloc[1])


def test_run_outputs_only_separate_descriptive_qa_and_refuses_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "refined-field-QA"
    manifest = refined.run(output)
    assert manifest["field_sources_used"] == ["central_region", "upper_east_region"]
    assert manifest["validation"]["chemistry_zero_imputation"] is False
    assert "flow direction from elevation or chemistry" in manifest[
        "unsupported_inferences_not_run"
    ]
    assert (output / "sample_charge_balance_reconciliation.csv").exists()
    assert (output / "major_ion_descriptive_summary.csv").exists()
    assert not (output / "field_prequential_predictions.csv").exists()

    cohort = pd.read_csv(output / "cohort_coverage_and_qc.csv").set_index("dataset")
    assert cohort.loc["central_region", "n_samples"] == 252
    assert cohort.loc["central_region", "n_complete_major8"] == 168
    assert cohort.loc["upper_east_region", "n_samples"] == 237
    assert cohort.loc["upper_east_region", "n_complete_major8"] == 237

    loaded_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    assert set(loaded_manifest["source_manifest"]) == {"central_region", "upper_east_region"}
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        refined.run(output)
