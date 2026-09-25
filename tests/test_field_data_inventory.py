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


def test_inventory_contains_only_the_two_approved_completed_workbooks() -> None:
    report = inventory_field_data(ROOT)
    datasets = _datasets(report)
    assert set(datasets) == {"central_region", "upper_east_region"}
    assert all(dataset["exists"] for dataset in datasets.values())
    assert _sheet(datasets["central_region"], "GW_Field_Integration")["effective_rows"] == 252
    assert _sheet(datasets["upper_east_region"], "GW_Field_Integration")["effective_rows"] == 237
    assert "NorthernGhana" not in str(report)
    assert "LowerAnayari" not in str(report)
    assert "Talensi_MiningArea" not in str(report)


def test_inventory_rendering_preserves_unit_and_flow_proxy_caveats() -> None:
    report = inventory_field_data(ROOT)
    markdown = render_markdown(report)
    assert "undocumented units are not inferred" in markdown
    assert "central_region" in markdown
    assert "upper_east_region" in markdown
    assert "cannot establish well-to-well flow direction" in markdown


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
