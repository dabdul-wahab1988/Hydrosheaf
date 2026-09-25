from __future__ import annotations

from pathlib import Path

import pytest

from hydrosheaf.data.field import (
    FIELD_DATA_RELATIVE_PATHS,
    flatten_field_datasets,
    load_all_field_datasets,
    load_field_dataset,
)
from hydrosheaf.models.ratios import (
    compare_ratio_diagnostics,
    compute_geochemical_ratios,
)


ROOT = Path(__file__).resolve().parents[1]
APPROVED = {"central_region", "upper_east_region"}


def test_only_refined_central_and_upper_east_cohorts_are_loadable() -> None:
    datasets = load_all_field_datasets(field_root=ROOT / "data" / "FieldData")
    assert set(FIELD_DATA_RELATIVE_PATHS) == APPROVED
    assert set(datasets) == APPROVED
    assert datasets["central_region"].n_records == 252
    assert datasets["upper_east_region"].n_records == 237

    # Missing measurements remain missing; no zero-imputation is permitted.
    assert datasets["central_region"].coverage()["Ca"] == 188
    assert datasets["upper_east_region"].coverage()["Ca"] == 237
    assert datasets["upper_east_region"].coverage()["F"] == 188
    assert datasets["upper_east_region"].records[0]["geology_join_status"] == "MATCHED"

    # Upper East topography is not promoted to a measured hydraulic head.
    assert datasets["upper_east_region"].records[0]["elevation"] is not None
    assert datasets["upper_east_region"].records[0]["hydraulic_head"] is None

    with pytest.raises(ValueError, match="Unknown field dataset"):
        load_field_dataset("northern_ghana")
    with pytest.raises(ValueError, match="Unknown field dataset"):
        load_field_dataset("northern_ghana_new")


def test_ratios_use_mmol_and_preserve_missingness() -> None:
    diagnostics = compute_geochemical_ratios(
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 2.0, "HCO3": 6.0}
    )
    assert diagnostics.values["Na_Cl"] == pytest.approx(2.0)
    assert diagnostics.values["HCO3_CaMg_equiv"] == pytest.approx(1.0)
    assert "Ca_Sr" in diagnostics.missing
    assert diagnostics.log_values["Na_Cl"] == pytest.approx(0.69314718056)


def test_flattening_keeps_the_two_refined_cohorts_separate() -> None:
    datasets = load_all_field_datasets(field_root=ROOT / "data" / "FieldData")
    rows = flatten_field_datasets(datasets)
    assert len(rows) == 489
    assert len({row["node_id"] for row in rows}) == 489
    assert {row["dataset"] for row in rows} == APPROVED
    assert all(row["season"] == "unknown" for row in rows)


def test_external_sidecar_joins_are_not_implicitly_loaded(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="external sidecar joins are not accepted"):
        load_field_dataset(
            "upper_east_region",
            field_root=ROOT / "data" / "FieldData",
            geology_join_path=tmp_path / "unverified_join.csv",
        )


def test_ratio_comparison_does_not_impute_unmeasured_trace_pairs() -> None:
    comparison = compare_ratio_diagnostics(
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 2.0},
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 4.0},
    )
    assert comparison["n_pairs"] >= 2
    assert comparison["similarity"] is not None
    assert comparison["missing_upstream"]
