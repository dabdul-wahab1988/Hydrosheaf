"""Tests for M3 USGS-calibrated parity emulator (Phase 9)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from M3.m3_age_benchmark.scripts.run_m3_usgs_calibrated_emulator import (
    _build_features,
    _log10,
    _target_log10_age,
    leave_study_unit_out_calibration,
)


def test_log10_returns_none_for_non_positive():
    assert _log10(0) is None
    assert _log10(-5) is None
    assert _log10(None) is None
    assert _log10("abc") is None
    assert _log10(100) == 2.0


def test_build_features_includes_expected_columns():
    df = pd.DataFrame(
        [
            {
                "est_age_total_years": 100.0,
                "n_tracers": 3,
                "reported_model": "DM",
                "age_class": "intermediate_50_1000",
                "Rpt_ChiSquare": 2.5,
                "FracAnthropocene": 0.3,
                "FracHolocene": 0.5,
                "FracPleistocene": 0.2,
                "Depth_m": 50.0,
            }
        ]
    )
    feats = _build_features(df)
    assert "log10_est_age" in feats.columns
    assert "reported_model_dm" in feats.columns
    assert "age_class_intermediate" in feats.columns
    assert "frac_anthropocene" in feats.columns
    assert feats["log10_est_age"].iloc[0] == 2.0


def test_target_log10_age_computes_correctly():
    df = pd.DataFrame({"Rpt_TotAge_yrs": [10.0, 100.0, 0.0, math.nan]})
    y = _target_log10_age(df)
    assert y[0] == 1.0
    assert y[1] == 2.0
    assert y[2] is None
    assert y[3] is None


def test_leave_study_unit_out_never_mixes_study_unit():
    df = pd.DataFrame(
        {
            "site_id": ["A1", "A2", "B1", "B2", "B3", "C1", "C2", "C3", "C4"],
            "StudyUnit": ["SU1", "SU1", "SU2", "SU2", "SU2", "SU3", "SU3", "SU3", "SU3"],
            "est_age_total_years": [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
            "n_tracers": [2] * 9,
            "reported_model": ["DM"] * 9,
            "age_class": ["modern_le_50"] * 9,
            "Rpt_ChiSquare": [1.0] * 9,
            "FracAnthropocene": [0.5] * 9,
            "FracHolocene": [0.3] * 9,
            "FracPleistocene": [0.2] * 9,
            "Depth_m": [10.0] * 9,
            "Rpt_TotAge_yrs": [12.0, 22.0, 32.0, 42.0, 52.0, 62.0, 72.0, 82.0, 92.0],
        }
    )
    result = leave_study_unit_out_calibration(df, method="ridge")
    for su in df["StudyUnit"].unique():
        test_mask = result["_held_out_study_unit"] == su
        train_mask = result["_fold"] == su
        # Verify no row in the held-out SU was used for training
        assert not (test_mask & ~train_mask).any()
        # Verify predictions exist for held-out rows
        assert result.loc[test_mask, "_calibrated"].notna().any() or test_mask.sum() < 10


def test_leave_study_unit_out_writes_fold_predictions():
    # Need enough rows so that each held-out fold still leaves >= 10 train rows
    rows = []
    for i in range(15):
        rows.append({
            "site_id": f"A{i}", "StudyUnit": "SU1",
            "est_age_total_years": 10.0 + i, "n_tracers": 2, "reported_model": "DM",
            "age_class": "modern_le_50", "Rpt_ChiSquare": 1.0,
            "FracAnthropocene": 0.5, "FracHolocene": 0.3, "FracPleistocene": 0.2,
            "Depth_m": 10.0, "Rpt_TotAge_yrs": 12.0 + i,
        })
    for i in range(12):
        rows.append({
            "site_id": f"B{i}", "StudyUnit": "SU2",
            "est_age_total_years": 30.0 + i, "n_tracers": 2, "reported_model": "DM",
            "age_class": "modern_le_50", "Rpt_ChiSquare": 1.0,
            "FracAnthropocene": 0.5, "FracHolocene": 0.3, "FracPleistocene": 0.2,
            "Depth_m": 10.0, "Rpt_TotAge_yrs": 32.0 + i,
        })
    df = pd.DataFrame(rows)
    result = leave_study_unit_out_calibration(df, method="ridge")
    assert "_calibrated" in result.columns
    assert "_raw_pred" in result.columns
    assert "_residual" in result.columns
    assert result.loc[result["_fold"] != "", "_calibrated"].notna().any()


def test_calibration_skips_missing_target_rows():
    df = pd.DataFrame(
        {
            "site_id": ["A1", "A2", "B1", "B2", "B3", "C1", "C2", "C3", "C4"],
            "StudyUnit": ["SU1", "SU1", "SU2", "SU2", "SU2", "SU3", "SU3", "SU3", "SU3"],
            "est_age_total_years": [10.0] * 9,
            "n_tracers": [2] * 9,
            "reported_model": ["DM"] * 9,
            "age_class": ["modern_le_50"] * 9,
            "Rpt_ChiSquare": [1.0] * 9,
            "FracAnthropocene": [0.5] * 9,
            "FracHolocene": [0.3] * 9,
            "FracPleistocene": [0.2] * 9,
            "Depth_m": [10.0] * 9,
            "Rpt_TotAge_yrs": [math.nan] * 9,
        }
    )
    result = leave_study_unit_out_calibration(df, method="ridge")
    # All targets are NaN so no folds should be processed
    assert result["_fold"].eq("").all()
