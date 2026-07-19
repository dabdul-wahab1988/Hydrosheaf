from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
M7_2_SCRIPTS = ROOT / "M7" / "m7_strong_integration" / "scripts"
M7_3_SCRIPTS = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
for script_dir in (M7_3_SCRIPTS, M7_2_SCRIPTS):
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))

from m7_3_analysis import (  # noqa: E402
    CORE_IONS,
    apply_evidence_condition,
    audit_ghana_workbook,
    fit_evidence_models,
    normalized_binary_entropy,
    predict_evidence_model,
    reaction_support_summary,
    topology_age_sensitivity,
)
from strong_inference import ION_ORDER, strong_config  # noqa: E402


def _development_frame() -> pd.DataFrame:
    rows = []
    for seed in (1, 2, 3):
        for index in range(12):
            truth = int(index < 3)
            rows.append(
                {
                    "seed": seed,
                    "is_true_edge": truth,
                    "hydraulic_logit": 2.0 * truth + 0.1 * index,
                    "negative_age_cost": -0.1 * (1 - truth) - 0.01 * index,
                    "negative_chemistry_log_objective": (1.2 * truth - 0.05 * index),
                }
            )
    return pd.DataFrame(rows)


def test_evidence_models_fit_all_locked_panels() -> None:
    frame = _development_frame()
    models = fit_evidence_models(frame)
    assert set(models) == {"H", "A", "C", "HA", "HC", "AC", "HAC"}
    probability = predict_evidence_model(frame, models["HAC"])
    assert np.all((probability > 0.0) & (probability < 1.0))
    assert (
        probability[frame["is_true_edge"].eq(1)].mean()
        > probability[frame["is_true_edge"].eq(0)].mean()
    )
    assert models["HAC"]["class_weight"] is None


def test_negative_control_permutation_is_casewise_and_deterministic() -> None:
    frame = _development_frame()
    altered = apply_evidence_condition(frame, "joint_misspecified")
    repeated = apply_evidence_condition(frame, "joint_misspecified")
    pd.testing.assert_frame_equal(altered, repeated)
    for seed in frame["seed"].unique():
        original = frame[frame["seed"] == seed]
        changed = altered[altered["seed"] == seed]
        assert sorted(original["hydraulic_logit"]) == sorted(changed["hydraulic_logit"])
        assert sorted(original["negative_age_cost"]) == sorted(
            changed["negative_age_cost"]
        )
    pd.testing.assert_series_equal(
        frame["negative_chemistry_log_objective"],
        altered["negative_chemistry_log_objective"],
    )


def test_normalized_binary_entropy_has_expected_limits() -> None:
    entropy = normalized_binary_entropy(np.asarray([0.5, 1.0e-8, 1.0 - 1.0e-8]))
    assert entropy[0] == pytest.approx(1.0)
    assert entropy[1] < 1.0e-5
    assert entropy[2] < 1.0e-5


def test_correct_and_reversed_topology_move_age_order_oppositely() -> None:
    true_ages = {"A": 3.0, "B": 20.0, "C": 42.0}
    observations = []
    for node, age in true_ages.items():
        observations.append(
            {
                "site_id": node,
                "tritium_TU": 6.2 * math.exp(-math.log(2.0) * age / 12.32),
                "argon39_pmc": 97.0 * math.exp(-math.log(2.0) * age / 269.0),
            }
        )
    result = topology_age_sensitivity(
        observations=observations,
        true_ages=true_ages,
        true_edges=(("A", "B"), ("B", "C")),
        seed=17,
        regime="tritium_only",
        n_particles=20_000,
    ).set_index("graph_condition")
    assert (
        result.loc["complete_true", "true_edge_order_violation_probability"]
        < result.loc["none", "true_edge_order_violation_probability"]
    )
    assert (
        result.loc["reversed", "true_edge_order_violation_probability"]
        > result.loc["complete_true", "true_edge_order_violation_probability"]
    )
    assert (
        result.loc["reversed", "log_mean_topology_weight"]
        < result.loc["complete_true", "log_mean_topology_weight"]
    )


def test_reaction_support_reports_entropy_and_carbonate_rows() -> None:
    frame = pd.DataFrame(
        {
            "seed": [1, 1, 1, 1],
            "edge_id": ["A->B"] * 4,
            "tier": ["core"] * 4,
            "true_family": ["carbonate"] * 4,
            "true_process": ["carbonate_weathering"] * 4,
            "predicted_family": [
                "carbonate",
                "carbonate",
                "silicate_exchange",
                "carbonate",
            ],
        }
    )
    support, edges, summary = reaction_support_summary(frame)
    carbonate = support[support["predicted_family"] == "carbonate"].iloc[0]
    assert carbonate["support_probability"] == pytest.approx(0.75)
    edge = edges.iloc[0]
    assert edge["true_family_probability"] == pytest.approx(0.75)
    assert edge["family_support_entropy"] > 0.0
    assert 1.0 < edge["effective_supported_families"] < 2.0
    row = summary[
        (summary["tier"] == "core")
        & (summary["true_process"] == "carbonate_weathering")
    ].iloc[0]
    assert row["modal_family_accuracy"] == 1.0
    assert 1.0 < row["mean_effective_supported_families"] < 2.0


def test_core_chemistry_config_zero_weights_excluded_ions() -> None:
    config = strong_config(phreeqc_enabled=True, measured_ions=CORE_IONS)
    by_ion = dict(zip(ION_ORDER, config.weights))
    assert all(by_ion[ion] == 1.0 for ion in CORE_IONS)
    assert all(by_ion[ion] == 0.0 for ion in set(ION_ORDER) - set(CORE_IONS))
    with pytest.raises(ValueError, match="Unknown chemistry ions"):
        strong_config(measured_ions=["not_an_ion"])


def test_ghana_scope_audit_enforces_data_limited_claim() -> None:
    workbook = ROOT / "data" / "NorthenGhana" / "Aquifers_Dataset_Mendeley.xlsx"
    audit = audit_ghana_workbook(workbook)
    assert audit["n_wells"] == 160
    assert audit["n_hydrochemistry_rows"] == 320
    assert audit["environmental_age_tracer_panel_available"] is False
    assert audit["screen_intervals_available"] is False
    assert audit["single_occasion_head_proxy_possible"] is True
    assert audit["time_varying_head_series_available"] is False
    assert audit["coordinates_masked"] is True
    assert audit["independent_field_connectivity_truth_available"] is False


def test_locked_result_decisions_trace_to_raw_artifacts() -> None:
    result = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "results" / "m7_3_locked"
    manifest = json.loads((result / "manifest.json").read_text(encoding="utf-8"))
    decision = json.loads(
        (result / "confirmatory_decision.json").read_text(encoding="utf-8")
    )
    contrasts = pd.read_csv(result / "evidence_case_bootstrap_contrasts.csv")
    topology = pd.read_csv(result / "topology_age_sensitivity.csv")

    assert manifest["protocol_commit_before_5301_series"] == (
        decision["protocol_commit_before_test_generation"]
    )
    assert manifest["locked_test_seeds"] == list(range(5301, 5313))

    age_pr_auc = contrasts[
        (contrasts["contrast"] == "native_incremental_age")
        & (contrasts["metric"] == "pr_auc")
    ].iloc[0]
    assert decision["decisions"]["age_adds_topology_ranking_value"][
        "pr_auc_difference_hac_minus_hc"
    ] == pytest.approx(age_pr_auc["mean_difference"])

    correct = topology[topology["graph_condition"] == "complete_true"]
    assert correct["importance_stable_ess_ge_400"].all()
    reversed_tritium = topology[
        (topology["graph_condition"] == "reversed")
        & (topology["tracer_regime"] == "tritium_only")
    ]
    assert (~reversed_tritium["importance_stable_ess_ge_400"]).sum() == 8

    case_directories = [path for path in (result / "cases").iterdir() if path.is_dir()]
    assert len(case_directories) == 18
    required = {
        "blind_observations.csv",
        "heldout_truth.csv",
        "modpath_pathline_truth.csv",
        "provenance.json",
        "diagnostics.json",
    }
    assert all(
        required <= {path.name for path in case.iterdir()} for case in case_directories
    )
