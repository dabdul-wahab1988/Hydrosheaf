from __future__ import annotations

from pathlib import Path

import pytest

from hydrosheaf.data.field import (
    flatten_field_datasets,
    load_all_field_datasets,
)
from hydrosheaf.models.ratios import (
    compare_ratio_diagnostics,
    compute_geochemical_ratios,
)


ROOT = Path(__file__).resolve().parents[1]


def test_all_four_field_packages_are_harmonised_without_source_mutation() -> None:
    datasets = load_all_field_datasets(field_root=ROOT / "data" / "FieldData")
    assert set(datasets) == {
        "lower_anayari",
        "northen_ghana",
        "northern_ghana_new",
        "talensi_mining_area",
    }
    assert [datasets[key].n_records for key in datasets] == [41, 320, 237, 63]
    assert datasets["northen_ghana"].coverage()["SiO2"] == 320
    assert datasets["northen_ghana"].coverage()["Sr"] == 320
    assert datasets["lower_anayari"].coverage()["Fe"] == 41
    assert datasets["talensi_mining_area"].coverage()["Fe"] == 63
    assert datasets["northern_ghana_new"].coverage()["F"] == 188
    assert datasets["northern_ghana_new"].records[0]["geology_join_status"] == "MATCHED"
    assert datasets["northen_ghana"].records[0]["well_depth"] == pytest.approx(93.5)
    assert datasets["northern_ghana_new"].auxiliary_tables["rain"]


def test_ratios_use_mmol_and_preserve_missingness() -> None:
    diagnostics = compute_geochemical_ratios(
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 2.0, "HCO3": 6.0}
    )
    assert diagnostics.values["Na_Cl"] == pytest.approx(2.0)
    assert diagnostics.values["HCO3_CaMg_equiv"] == pytest.approx(1.0)
    assert "Ca_Sr" in diagnostics.missing
    assert diagnostics.log_values["Na_Cl"] == pytest.approx(0.69314718056)


def test_flattening_preserves_all_four_packages_and_seasonal_nodes() -> None:
    datasets = load_all_field_datasets(field_root=ROOT / "data" / "FieldData")
    rows = flatten_field_datasets(datasets)
    assert len(rows) == 661
    assert len({row["node_id"] for row in rows}) == 661
    northern_nodes = [
        row["node_id"]
        for row in rows
        if row["dataset"] == "northen_ghana" and row["site_id"] == "NG_NGW-001"
    ]
    assert sorted(northern_nodes) == ["NG_NGW-001_dry", "NG_NGW-001_wet"]


def test_ratio_comparison_does_not_impute_unmeasured_trace_pairs() -> None:
    comparison = compare_ratio_diagnostics(
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 2.0},
        {"Ca": 2.0, "Mg": 1.0, "Na": 4.0, "Cl": 4.0},
    )
    assert comparison["n_pairs"] >= 2
    assert comparison["similarity"] is not None
    assert comparison["missing_upstream"]
