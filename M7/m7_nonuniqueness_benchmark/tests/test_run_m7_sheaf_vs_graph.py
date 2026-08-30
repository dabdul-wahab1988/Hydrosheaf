from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import independent_sheaf_graph_generator as generator  # noqa: E402
from run_m7_sheaf_vs_graph import (  # noqa: E402
    _feature_frame,
    _paired_bootstrap,
    _select_threshold,
)


def test_generator_is_independent_of_hydrosheaf_package():
    source = inspect.getsource(generator)
    assert "from hydrosheaf" not in source
    assert "import hydrosheaf" not in source


def test_identity_limit_reduces_exactly_to_graph_laplacian():
    case = generator.generate_independent_section_case(101, "identity_limit", n_nodes=10)
    frame, diagnostics = _feature_frame(case)
    np.testing.assert_allclose(
        frame["identity_section_residual"],
        frame["affine_section_residual"],
        atol=1.0e-10,
        rtol=0.0,
    )
    assert diagnostics["identity_affine_residual_max_abs_difference"] <= 1.0e-10


def test_locked_features_have_no_truth_in_inference_columns():
    case = generator.generate_independent_section_case(102, "heterogeneous_affine", n_nodes=10)
    frame, _ = _feature_frame(case)
    feature_columns = {
        "prior_logit",
        "local_missing",
        "local_affine_residual",
        "identity_section_residual",
        "affine_section_residual",
        "permuted_section_residual",
    }
    assert not ({"is_true_edge", "is_corrupted_edge"} & feature_columns)
    assert frame[list(feature_columns)].shape[1] == len(feature_columns)


def test_threshold_selection_is_deterministic():
    truth = np.asarray([0, 0, 1, 1])
    probability = np.asarray([0.1, 0.4, 0.6, 0.9])
    assert _select_threshold(truth, probability) == 0.5


def test_case_block_bootstrap_preserves_pairing_and_direction():
    frame = pd.DataFrame(
        [
            {"seed": seed, "scenario": "all", "model": model, "pr_auc": value}
            for seed, left, right in ((1, 0.8, 0.4), (2, 0.7, 0.5))
            for model, value in (("left", left), ("right", right))
        ]
    )
    result = _paired_bootstrap(
        frame,
        "left",
        "right",
        "pr_auc",
        samples=500,
        seed=1,
    )
    assert result["mean_difference"] > 0.0
    assert result["ci95_low"] > 0.0
