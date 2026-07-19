from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_strong_integration" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from field_prequential import ION_COLUMNS, run_prequential_frames  # noqa: E402
from independent_modflow_generator import (  # noqa: E402
    ION_ORDER,
    _chemistry_step,
)
from run_m7_2_strong_integration import _predict_fusion  # noqa: E402


def test_independent_generator_does_not_import_hydrosheaf() -> None:
    source = (SCRIPT_DIR / "independent_modflow_generator.py").read_text(
        encoding="utf-8"
    )
    import_lines = [
        line.strip()
        for line in source.splitlines()
        if line.strip().startswith(("import ", "from "))
    ]
    assert not any("hydrosheaf" in line.lower() for line in import_lines)


def test_independent_sulfate_reduction_has_expected_net_signature() -> None:
    state = np.ones(len(ION_ORDER), dtype=float)
    rng = np.random.default_rng(7)
    reacted = _chemistry_step(
        state,
        travel_years=20.0,
        process="sulfate_reduction",
        rng=rng,
    )
    index = {ion: idx for idx, ion in enumerate(ION_ORDER)}
    assert reacted[index["SO4"]] < state[index["SO4"]]
    assert reacted[index["HCO3"]] > state[index["HCO3"]]
    assert reacted[index["NO3"]] < state[index["NO3"]]


def test_independent_denitrification_has_expected_net_signature() -> None:
    state = np.ones(len(ION_ORDER), dtype=float)
    rng = np.random.default_rng(9)
    reacted = _chemistry_step(
        state,
        travel_years=20.0,
        process="denitrification",
        rng=rng,
    )
    index = {ion: idx for idx, ion in enumerate(ION_ORDER)}
    assert reacted[index["NO3"]] < state[index["NO3"]]
    assert reacted[index["HCO3"]] > state[index["HCO3"]]


def test_independent_iron_reduction_has_expected_net_signature() -> None:
    state = np.ones(len(ION_ORDER), dtype=float)
    state[ION_ORDER.index("NO3")] = 0.01
    state[ION_ORDER.index("Fe")] = 0.01
    rng = np.random.default_rng(11)
    reacted = _chemistry_step(
        state,
        travel_years=20.0,
        process="iron_reduction",
        rng=rng,
    )
    index = {ion: idx for idx, ion in enumerate(ION_ORDER)}
    assert reacted[index["Fe"]] > state[index["Fe"]]
    assert reacted[index["HCO3"]] > state[index["HCO3"]]


def test_confirmatory_age_gate_only_suppresses_incompatible_edges() -> None:
    frame = pd.DataFrame(
        {
            "hydraulic_logit": [0.0, 0.0],
            "negative_chemistry_log_objective": [0.0, 0.0],
            "age_cost": [0.01, 0.10],
        }
    )
    model = {
        "feature_names": [
            "hydraulic_logit",
            "negative_chemistry_log_objective",
        ],
        "means": [0.0, 0.0],
        "scales": [1.0, 1.0],
        "coefficients": [0.0, 0.0],
        "intercept": 0.0,
        "kind": "age_compatibility_gate",
        "age_cost_max": 0.05,
        "incompatible_probability": 1.0e-6,
    }
    probability = _predict_fusion(frame, model)
    assert probability[0] == 0.5
    assert probability[1] == 1.0e-6


def test_field_predictions_cannot_access_future_wet_batches() -> None:
    workbook = (
        ROOT
        / "data"
        / "NorthenGhana"
        / "Aquifers_Dataset_Mendeley.xlsx"
    )
    wells = pd.read_excel(workbook, sheet_name="Wells_Nodes")
    hydro = pd.read_excel(workbook, sheet_name="Hydrochemistry_Seasonal")
    original = run_prequential_frames(wells, hydro).predictions

    altered_hydro = hydro.copy()
    dates = pd.to_datetime(altered_hydro["Sampling_Date"])
    future = (
        altered_hydro["Season"].astype(str).str.lower().eq("wet")
        & (dates > pd.Timestamp("2025-08-10"))
    )
    altered_hydro.loc[future, list(ION_COLUMNS)] *= 1000.0
    altered = run_prequential_frames(wells, altered_hydro).predictions

    cutoff = pd.Timestamp("2025-08-10")
    columns = [
        "issue_date",
        "well_id",
        "ion",
        "method",
        "prediction_log1p",
    ]
    original_early = (
        original.loc[original["issue_date"] <= cutoff, columns]
        .sort_values(columns[:-1])
        .reset_index(drop=True)
    )
    altered_early = (
        altered.loc[altered["issue_date"] <= cutoff, columns]
        .sort_values(columns[:-1])
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(original_early, altered_early)
