from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from field_prequential import ION_COLUMNS, run_prequential_frames  # noqa: E402
from independent_modflow_generator import (  # noqa: E402
    ION_ORDER,
    _chemistry_step,
)
from run_supporting_validation import (  # noqa: E402
    BASELINE_FEATURES,
    _age_permuted_frame,
    _fit_fusion_model,
    _fit_monotone_age_model,
    _predict_fusion,
)


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


def test_continuous_age_penalty_cannot_reverse_physical_direction() -> None:
    frame = pd.DataFrame(
        {
            "seed": [1, 1, 1, 1],
            "is_true_edge": [1, 0, 1, 0],
            "hydraulic_logit": [1.0, 0.0, 1.0, 0.0],
            "negative_chemistry_log_objective": [0.0, 0.0, 0.0, 0.0],
            "age_cost": [0.0, 0.25, 0.05, 0.75],
            "age_evidence_available": [True, True, True, True],
        }
    )
    baseline = _fit_fusion_model(frame, BASELINE_FEATURES)
    model = _fit_monotone_age_model(frame, baseline)
    assert model["kind"] == "monotone_age_penalty"
    assert float(model["age_cost_coefficient"]) >= 0.0

    same_baseline = frame.iloc[[0, 1]].copy()
    same_baseline["hydraulic_logit"] = 0.5
    same_baseline["negative_chemistry_log_objective"] = -0.1
    same_baseline["age_cost"] = [0.0, 1.0]
    probabilities = _predict_fusion(same_baseline, model)
    assert probabilities[0] >= probabilities[1]

    missing_age = same_baseline.copy()
    missing_age["age_cost"] = [100.0, 0.0]
    missing_age["age_evidence_available"] = [False, True]
    missing_probabilities = _predict_fusion(missing_age, model)
    assert missing_probabilities[0] == pytest.approx(missing_probabilities[1])


def test_age_permutation_rebuilds_the_monotone_transform() -> None:
    frame = pd.DataFrame(
        {
            "seed": [1, 1, 2, 2],
            "age_cost": [0.0, 0.5, 0.1, 0.9],
            "negative_age_cost": [0.0, -np.log1p(0.5), -np.log1p(0.1), -np.log1p(0.9)],
            "age_evidence_available": [True, False, True, True],
        }
    )
    permuted = _age_permuted_frame(
        frame,
        np.random.default_rng(13),
        "logistic",
    )
    np.testing.assert_allclose(
        permuted["negative_age_cost"].to_numpy(float),
        -np.log1p(permuted["age_cost"].to_numpy(float)),
    )


def test_prequential_helper_has_no_future_label_leakage_on_synthetic_rows() -> None:
    n_wells = 40
    wells = pd.DataFrame(
        {
            "Well_ID": [f"S{i:03d}" for i in range(n_wells)],
            "Latitude": [5.0 + i * 0.001 for i in range(n_wells)],
            "Longitude": [-1.0 - i * 0.001 for i in range(n_wells)],
            "Elevation_m": [100.0 + i for i in range(n_wells)],
            "Static_Water_Level_m": [10.0 + (i % 5) for i in range(n_wells)],
            "Borehole_Depth_m": [45.0 + i % 7 for i in range(n_wells)],
            "Distance_River_km": [1.0 + i % 3 for i in range(n_wells)],
            "Distance_Farm_km": [0.5 + i % 4 for i in range(n_wells)],
            "Distance_Settlement_km": [0.7 + i % 6 for i in range(n_wells)],
            "Region": ["A" if i % 2 else "B" for i in range(n_wells)],
        }
    )
    base_ions = {
        "Ca_mg_L": 20.0,
        "Mg_mg_L": 10.0,
        "Na_mg_L": 10.0,
        "K_mg_L": 5.0,
        "HCO3_mg_L": 90.0,
        "Cl_mg_L": 20.0,
        "SO4_mg_L": 10.0,
        "NO3_mg_L": 5.0,
        "F_mg_L": 0.2,
        "Sr_mg_L": 0.1,
        "SiO2_mg_L": 15.0,
    }
    hydro_rows = []
    for i in range(n_wells):
        for season, multiplier in (("Dry", 1.0), ("Wet", 1.02 + 0.0005 * i)):
            row = {
                "Well_ID": f"S{i:03d}",
                "Season": season,
                **{ion: value * multiplier for ion, value in base_ions.items()},
                "pH": 7.0 + (i % 4) * 0.05,
                "EC_uS_cm": 350.0 + i,
                "TDS_mg_L": 225.0 + i,
                "Temperature_C": 25.0,
                "d18O_permil": -4.0 + i * 0.001,
                "d2H_permil": -25.0 + i * 0.01,
            }
            hydro_rows.append(row)
    hydro = pd.DataFrame(hydro_rows)
    original = run_prequential_frames(wells, hydro, n_batches=8).predictions

    # Batch assignment depends only on which wells are eligible (a
    # scale-invariant charge-balance screen), not on chemistry magnitude, so
    # grossly altering wet-season values for later-batch wells cannot change
    # which wells fall in earlier batches. Cutting off at the midpoint batch
    # index isolates "not yet revealed" wells from "already revealed" ones.
    cutoff_batch = int(original["issue_batch_index"].max()) // 2
    later_wells = set(
        original.loc[original["issue_batch_index"] > cutoff_batch, "well_id"]
    )
    altered_hydro = hydro.copy()
    future = (
        altered_hydro["Season"].astype(str).str.lower().eq("wet")
        & altered_hydro["Well_ID"].astype(str).isin(later_wells)
    )
    altered_hydro.loc[future, list(ION_COLUMNS)] *= 1000.0
    altered = run_prequential_frames(wells, altered_hydro, n_batches=8).predictions

    columns = [
        "issue_batch_index",
        "well_id",
        "ion",
        "method",
        "prediction_log1p",
    ]
    original_early = (
        original.loc[original["issue_batch_index"] <= cutoff_batch, columns]
        .sort_values(columns[:-1])
        .reset_index(drop=True)
    )
    altered_early = (
        altered.loc[altered["issue_batch_index"] <= cutoff_batch, columns]
        .sort_values(columns[:-1])
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(original_early, altered_early)
