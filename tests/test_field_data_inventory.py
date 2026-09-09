from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from field_data_inventory import (  # noqa: E402
    inventory_field_data,
    render_markdown,
)


def _datasets(report: dict[str, object]) -> dict[str, dict[str, object]]:
    return {
        str(dataset["dataset"]): dataset
        for dataset in report["datasets"]  # type: ignore[index]
    }


def _sheet(dataset: dict[str, object], name: str) -> dict[str, object]:
    return next(
        sheet for sheet in dataset["sheets"]  # type: ignore[index]
        if sheet["sheet"] == name
    )


def test_field_inventory_covers_the_four_live_sources_without_writing() -> None:
    report = inventory_field_data(ROOT)
    datasets = _datasets(report)
    assert set(datasets) == {
        "LowerAnayari/manu",
        "NorthenGhana/NorthernGhana",
        "NorthernGhanaNew/compiled UER data_new",
        "Talensi_MiningArea/talensi",
    }
    assert all(dataset["exists"] for dataset in datasets.values())

    lower = datasets["LowerAnayari/manu"]
    assert lower["effective_rows"] == 41
    lower_fields = lower["requested_field_availability"]
    assert lower_fields["SiO2"]["available_in_any_sheet"] is False
    assert lower_fields["Sr"]["available_in_any_sheet"] is False
    assert lower_fields["Fe"]["available_in_any_sheet"] is True
    assert lower_fields["Fe"]["unit_status"] == "not_declared"
    assert any(
        proxy["name"] == "elevation_gradient_context"
        for proxy in lower["flow_proxies"]
    )

    northern = datasets["NorthenGhana/NorthernGhana"]
    assert northern["effective_rows"] == 320
    northern_fields = northern["requested_field_availability"]
    for field in ("SiO2", "Sr", "depth"):
        assert northern_fields[field]["available_in_any_sheet"] is True
    assert northern_fields["SiO2"]["declared_units"] == ["mg/L"]
    assert northern_fields["Sr"]["declared_units"] == ["mg/L"]
    assert northern_fields["depth"]["declared_units"] == ["m"]
    assert northern_fields["Fe"]["available_in_any_sheet"] is False
    overlap = northern["cross_sheet_sample_id_overlaps"]
    assert overlap[0]["sample_id_overlap_count"] == 160
    for name in ("Dry", "Wet"):
        sheet = _sheet(northern, name)
        assert sheet["raw_rows"] == 160
        assert sheet["effective_rows"] == 160
        assert sheet["duplicates"]["exact_duplicate_extra_rows"] == 0

    uer = datasets["NorthernGhanaNew/compiled UER data_new"]
    assert uer["effective_rows"] == 408
    gw = _sheet(uer, "GW")
    assert gw["raw_rows"] == 241
    assert gw["effective_rows"] == 240
    assert gw["unit_row_detected"] is True
    assert gw["unit_row_index"] == 0
    assert _sheet(uer, "rain")["canonical_fields"]["date"]["available"] is True
    monitoring = _sheet(uer, "monitoring wells")
    assert monitoring["canonical_fields"]["date"]["available"] is True
    assert any(
        proxy["name"] == "monitoring_time_series"
        for proxy in monitoring["flow_proxies"]
    )

    talensi = datasets["Talensi_MiningArea/talensi"]
    assert talensi["effective_rows"] == 63
    talensi_fields = talensi["requested_field_availability"]
    assert talensi_fields["Fe"]["available_in_any_sheet"] is True
    assert talensi_fields["Fe"]["unit_status"] == "not_declared"
    assert talensi_fields["SiO2"]["available_in_any_sheet"] is False
    assert talensi_fields["Sr"]["available_in_any_sheet"] is False


def test_inventory_does_not_invent_units_and_markdown_is_explicit() -> None:
    report = inventory_field_data(ROOT)
    markdown = render_markdown(report)
    assert "undocumented units are not inferred" in markdown
    assert "LowerAnayari/manu" in markdown
    assert "NorthenGhana/NorthernGhana" in markdown
    assert "NorthernGhanaNew/compiled UER data_new" in markdown
    assert "Talensi_MiningArea/talensi" in markdown
    assert "yes (manu: 41/41; not declared)" in markdown
    assert "yes (Dry: 160/160; Wet: 160/160; mg/L)" in markdown
    assert "unit row 0 removed" in markdown
    assert "160 (repeated identifiers across sheets" in markdown


def test_missing_source_is_reported_without_raising(tmp_path: Path) -> None:
    report = inventory_field_data(
        tmp_path,
        {"missing": Path("does-not-exist.csv")},
    )
    assert report["datasets"] == [
        {
            "dataset": "missing",
            "path": str((tmp_path / "does-not-exist.csv").resolve()),
            "exists": False,
            "error": "source file is missing",
        }
    ]
