"""Pure-function tests for the M7.6 auxiliary M3 diagnostic."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import m7_6_analysis as analysis  # noqa: E402


def _row(site: str, *, x: float, age: float, e: float) -> dict[str, object]:
    row: dict[str, object] = {
        "site_id": site,
        "hydraulic_head": 100.0 - 0.01 * x,
        "elevation": 112.0 - 0.01 * x,
        "x_m": x,
        "y_m": 100.0,
        "pH": 7.0 + 0.01 * age,
        "temp_c": 25.0,
        "d18O_permil": e,
        "d2H_permil": 8.0 * e + 10.0,
        "sr87_sr86": 0.710 + 1.0e-4 * x,
        "d13C_permil": -23.0,
        "sample_date": 2025.5,
        "c14_pmc": 100.0 * np.exp(-np.log(2.0) * age / 5730.0),
        "cfc11_pptv": 200.0,
        "cfc12_pptv": 400.0,
        "cfc113_pptv": 70.0,
        "sf6_pptv": 8.0,
        "he4_ccpg": 4.5e-8 + 1.0e-10 * age,
        "h3_he3_TU": 0.2,
        "tritium_TU": 6.2 * np.exp(-age / 17.0),
        "argon39_pmc": 100.0 * np.exp(-np.log(2.0) * age / 269.0),
    }
    for ion, value in zip(analysis.ION_FEATURES, np.linspace(0.1, 1.2, len(analysis.ION_FEATURES))):
        row[ion] = float(value + 0.01 * age)
    return row


def test_all_four_streams_have_exactly_15_nonempty_subsets():
    assert len(analysis.SUBSETS) == 15
    assert set(analysis.SUBSETS) == {
        "H", "C", "N", "E", "HC", "HN", "HE", "CN", "CE", "NE",
        "HCN", "HCE", "HNE", "CNE", "HCNE",
    }


def test_environmental_isotope_age_mask_is_empty_and_not_accidentally_reused():
    assert analysis._feature_names("T2", "E") == []
    assert analysis._feature_names("T2", "NE") == analysis._feature_names("T2", "N")


def test_edge_features_are_truth_free_when_truth_is_omitted():
    observations = [_row("MF1_P0_M0", x=0.0, age=1.5, e=-3.0), _row("MF1_P0_M1", x=100.0, age=15.0, e=-3.0)]
    frame = analysis.build_edge_features(observations)
    assert len(frame) == 2
    assert frame["is_true_edge"].sum() == 0
    assert set(("H_score", "C_score", "N_score", "E_score")) <= set(frame.columns)


def test_m3_pair_diagnostic_returns_family_rows_without_truth_labels():
    observations = [_row(f"MF1_P0_M{i}", x=100.0 * i, age=1.5 + 10.0 * i, e=-3.0) for i in range(4)]
    node, pair = analysis.ttd_feasibility_diagnostics(
        observations,
        redox_by_node={row["site_id"]: "nonreducing" for row in observations},
        age_grid=np.linspace(0.25, 100.0, 40),
    )
    assert len(node) == 4
    assert {"cfc_cfc", "tritium_pair"} <= set(pair["pair_family"])
    assert set(pair["pair_infeasible"].unique()) <= {0, 1}


def test_case_bootstrap_is_case_block_paired():
    frame = pd.DataFrame(
        [
            {"seed": 1, "target": "T1", "condition": "native", "nuisance_level": "none", "panel": "H", "pr_auc": 0.2},
            {"seed": 1, "target": "T1", "condition": "native", "nuisance_level": "none", "panel": "HC", "pr_auc": 0.3},
            {"seed": 2, "target": "T1", "condition": "native", "nuisance_level": "none", "panel": "H", "pr_auc": 0.4},
            {"seed": 2, "target": "T1", "condition": "native", "nuisance_level": "none", "panel": "HC", "pr_auc": 0.5},
        ]
    )
    result = analysis.bootstrap_case_contrasts(
        frame,
        metric_columns=("pr_auc",),
        comparisons=(("HC", "H", "chemistry_added"),),
        n_bootstrap=100,
    )
    assert len(result) == 1
    assert result.loc[0, "n_cases"] == 2
    assert result.loc[0, "mean_difference"] == pytest.approx(0.1)
