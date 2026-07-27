from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "M6" / "m6_field_transfer_benchmark"
SCRIPTS = BENCHMARK / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import m6_common as m6  # noqa: E402
from run_m6_q1 import ANALYSIS_STEPS, R_STEPS  # noqa: E402


def test_seasonal_well_pairs_keeps_only_complete_wells() -> None:
    frame = pd.DataFrame(
        [
            {"site_id": "A", "season": "Wet", "value": 1},
            {"site_id": "A", "season": "Dry", "value": 2},
            {"site_id": "B", "season": "Dry", "value": 3},
        ]
    )
    pairs = m6.seasonal_well_pairs(frame)
    assert set(pairs) == {"A"}
    wet, dry = pairs["A"]
    assert wet["value"] == 1
    assert dry["value"] == 2


def test_q1_runner_contains_every_authoritative_stage() -> None:
    assert ANALYSIS_STEPS == (
        "run_m6_field_transfer.py",
        "run_m6_synthetic_validation.py",
        "run_m6_robustness_diagnostics.py",
        "run_m6_null_sensitivity.py",
        "export_m6_figure_data.py",
        "make_m6_tables.py",
    )
    assert R_STEPS == (
        "plot_m6_publication_figures.R",
        "plot_m6_validation_figures.R",
        "plot_m6_supplementary_figures.R",
    )


def test_q1_display_assets_are_complete() -> None:
    main = [
        "figure1_dataset_tier_design",
        "figure2_workflow",
        "figure3_ng_stability",
        "figure4_tier_ablation",
        "figure5_field_prequential",
        "figure6_limitation_map",
    ]
    extended = [
        "figureED1_external_transfer",
        "figureED2_synthetic_validation",
        "figureED3_circularity_sensitivity",
    ]
    formats = (".pdf", ".png", ".tif")
    for stem in main:
        assert all(
            (BENCHMARK / "figures" / "r_publication" / f"{stem}{suffix}").exists()
            for suffix in formats
        )
    for stem in extended:
        assert all(
            (BENCHMARK / "figures" / "extended_data" / f"{stem}{suffix}").exists()
            for suffix in formats
        )
    assert (
        BENCHMARK
        / "figures"
        / "r_publication"
        / "supplementary"
        / "figureS11_null_sensitivity.pdf"
    ).exists()
    assert (BENCHMARK / "tables" / "tableS9_null_sensitivity.csv").exists()
