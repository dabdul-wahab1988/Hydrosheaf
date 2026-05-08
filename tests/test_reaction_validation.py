from __future__ import annotations

from hydrosheaf.config import Config
from hydrosheaf.validation import (
    l1_penalty_sensitivity,
    missing_ion_sensitivity,
    thermodynamic_bound_violations,
    validate_sparse_inverse_reaction_model,
)


def test_l1_penalty_sensitivity_reports_sparser_high_lambda_fit():
    matrix = [
        [1.0, 0.0],
        [1.0, 1.0],
    ]
    labels = ["calcite", "gypsum"]
    residual = [1.0, 1.0]
    weights = [1.0, 1.0]

    rows = l1_penalty_sensitivity(
        residual,
        matrix,
        labels,
        weights,
        lambda_grid=[0.0, 0.5, 4.0],
    )

    assert rows[0]["lambda_l1"] == 0.0
    assert rows[-1]["l1_norm"] <= rows[0]["l1_norm"]
    assert rows[-1]["n_selected_reactions"] <= rows[0]["n_selected_reactions"]


def test_missing_ion_sensitivity_uses_zero_weight_dropouts():
    matrix = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    labels = ["calcite", "gypsum"]
    residual = [1.0, 1.0]
    rows = missing_ion_sensitivity(
        residual,
        matrix,
        labels,
        ion_order=["Ca", "SO4"],
        weights=[1.0, 1.0],
        missing_ion_sets=[["SO4"]],
        lambda_l1=0.0,
    )

    assert rows[0]["missing_ions"] == ["SO4"]
    selected = {row["reaction"] for row in rows[0]["selected_reactions"]}
    assert "gypsum" not in selected


def test_thermodynamic_bound_violations_detect_si_incompatible_extent():
    violations = thermodynamic_bound_violations(
        labels=["calcite", "gypsum"],
        extents=[1.0, -0.5],
        bounds={
            "lb": [0.0, 0.0],
            "ub": [float("inf"), float("inf")],
            "constraints_active": {"calcite": "dissolution_only", "gypsum": "dissolution_only"},
        },
    )

    assert len(violations) == 1
    assert violations[0]["reaction"] == "gypsum"


def test_sparse_inverse_reaction_validation_reports_m5_guardrails():
    config = Config(
        ion_order=["Ca", "SO4"],
        weights=[1.0, 1.0],
        conservative_weights=[1.0, 1.0],
    )
    report = validate_sparse_inverse_reaction_model(
        upstream=[1.0, 1.0],
        downstream=[2.0, 2.0],
        post_transport=[1.0, 1.0],
        reaction_matrix=[
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        reaction_labels=["calcite", "gypsum"],
        config=config,
        lambda_grid=[0.0, 0.5],
        missing_ion_sets=[["SO4"]],
        phreeqc_bounds={
            "lb": [0.0, 0.0],
            "ub": [float("inf"), float("inf")],
            "constraints_active": {
                "calcite": "dissolution_only",
                "gypsum": "dissolution_only",
            },
        },
    )

    assert report["not_a_fully_coupled_nonlinear_phreeqc_inverse_solver"] is True
    assert "sparse linear inverse reaction model" in report["model_framing"]
    assert report["transport_reaction_separation"]["transport_residual"] == [1.0, 1.0]
    assert report["best_fit"]["n_selected_reactions"] >= 1
    assert "not describe it as a fully coupled" in report["claim_guardrail"]
