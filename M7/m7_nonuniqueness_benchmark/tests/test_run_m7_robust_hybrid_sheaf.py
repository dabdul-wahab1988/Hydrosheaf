from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from independent_sheaf_graph_generator import generate_independent_section_case  # noqa: E402
from run_m7_robust_hybrid_sheaf import (  # noqa: E402
    GROUP_FOLDS,
    _crossfit,
    _feature_frame,
    _fold_assignments,
    _matrix,
    _representation,
    _robust_loo_section,
)


def test_loo_robust_estimator_is_finite_and_weight_bounded():
    case = generate_independent_section_case(
        202, "heterogeneous_affine", n_nodes=8, false_edge_ratio=1.0
    )
    residuals, weights, runtime_ms = _robust_loo_section(case, "affine")
    assert set(residuals) == {candidate.edge_id for candidate in case.candidates}
    assert np.isfinite(list(residuals.values())).all()
    assert all(1.0e-4 <= value <= 1.0 for value in weights.values())
    assert runtime_ms >= 0.0


def test_hybrid_boundaries_and_missing_endpoint_rule():
    frame = pd.DataFrame(
        {
            "local_affine_residual": [3.0, np.nan],
            "original_affine_residual": [7.0, 11.0],
        }
    )
    local_boundary, _ = _representation(frame, "original_hybrid", 1.0)
    global_boundary, _ = _representation(frame, "original_hybrid", 0.0)
    assert local_boundary[0] == np.log1p(3.0)
    assert local_boundary[1] == np.log1p(11.0)
    np.testing.assert_allclose(global_boundary, np.log1p([7.0, 11.0]))


def test_all_arms_have_equal_feature_dimension():
    case = generate_independent_section_case(
        203, "identity_limit", n_nodes=8, false_edge_ratio=1.0
    )
    frame, _ = _feature_frame(case)
    for arm, local_weight in (
        ("edge_local", None),
        ("identity_graph", None),
        ("original_affine_global", None),
        ("original_hybrid", 0.5),
        ("robust_affine_global", None),
        ("robust_hybrid", 0.5),
        ("robust_hybrid_permuted", 0.5),
    ):
        matrix, _ = _matrix(frame, arm, local_weight)
        assert matrix.shape == (len(frame), 3)


def test_group_folds_keep_each_seed_intact():
    frame = pd.DataFrame({"seed": np.repeat(np.arange(100, 116), 4)})
    assignments = _fold_assignments(frame)
    assert len(set(assignments.values())) == GROUP_FOLDS
    assert frame["seed"].map(assignments).groupby(frame["seed"]).nunique().max() == 1


def test_crossfit_is_deterministic_and_finite():
    rows = []
    for seed in range(300, 316):
        for index in range(6):
            truth = int(index < 2)
            rows.append(
                {
                    "seed": seed,
                    "scenario": "identity_limit",
                    "prior_logit": 1.0 if truth else -0.5,
                    "local_missing": int(index == 5),
                    "local_affine_residual": np.nan if index == 5 else float(index + 1),
                    "is_true_edge": truth,
                }
            )
    frame = pd.DataFrame(rows)
    first, first_score = _crossfit(frame, "edge_local", None, 1.0)
    second, second_score = _crossfit(frame, "edge_local", None, 1.0)
    np.testing.assert_allclose(first, second, rtol=0.0, atol=0.0)
    assert np.isfinite(first).all()
    assert first_score == second_score
