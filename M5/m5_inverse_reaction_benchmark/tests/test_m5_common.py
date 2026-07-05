"""Focused tests for the M5 identifiability analysis helpers."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
REPO_ROOT = Path(__file__).resolve().parents[3]
for path in (SCRIPT_DIR, REPO_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from m5_common import (  # noqa: E402
    ION_ORDER,
    REACTION_LABELS,
    class_set,
    equivalence_classes,
    fit_inverse,
    matrix_diagnostics,
    reaction_matrix,
    support_from_extents,
)
import run_m5_inverse_reaction_benchmark as m5  # noqa: E402


def test_expected_signed_equivalence_classes() -> None:
    class_map, _ = equivalence_classes()
    assert class_map["calcite"] == class_map["anorthite"]
    assert class_map["CaNa_exch"] == class_map["NaCa_exch"]
    assert class_map["MgNa_exch"] == class_map["NaMg_exch"]


def test_noise_free_fit_recovers_equivalence_class_and_mass_balance() -> None:
    truth = np.zeros(len(REACTION_LABELS))
    truth[REACTION_LABELS.index("anorthite")] = 0.2
    residual = reaction_matrix().T @ truth
    fit = fit_inverse(
        residual,
        ION_ORDER,
        "elastic_net",
        lambda_l1=0.00001,
        lambda_l2=0.005,
    )
    class_map, _ = equivalence_classes()
    recovered = support_from_extents(fit["extents"], threshold=0.01)
    assert class_map["anorthite"] in class_set(recovered, class_map)
    assert fit["rmse"] < 1e-3


def test_missing_ions_reduce_or_preserve_rank() -> None:
    full = matrix_diagnostics(ION_ORDER)
    core = matrix_diagnostics(["Ca", "Mg", "Na", "HCO3", "Cl"])
    assert core["rank"] <= full["rank"]
    assert core["nullity"] >= full["nullity"]


def test_detection_threshold_is_separate_from_truth_threshold() -> None:
    truth = np.zeros(len(REACTION_LABELS))
    truth[REACTION_LABELS.index("calcite")] = 0.02
    residual = reaction_matrix().T @ truth
    scenario = {"true_extents": truth, "true_delta": residual}
    class_map, _ = equivalence_classes()

    loose = m5.evaluate_fit(
        scenario,
        truth.copy(),
        residual,
        ION_ORDER,
        class_map,
        support_threshold=0.015,
    )
    strict = m5.evaluate_fit(
        scenario,
        truth.copy(),
        residual,
        ION_ORDER,
        class_map,
        support_threshold=0.03,
    )

    assert loose["phase_f1"] == 1.0
    assert strict["phase_recall"] == 0.0


def test_phreeqc_inverse_baseline_uses_parseable_phase_names_and_fallback() -> None:
    for phase_name, _ in m5.PHREEQC_INVERSE_PHASES.values():
        assert len(phase_name) <= 15

    labels, exchanges = m5._phreeqc_inverse_fallback_candidates(
        ["albite", "halite"],
        ["NaX"],
    )

    assert "dolomite" in labels
    assert "gypsum" in labels
    assert {"CaX2", "NaX", "MgX2"}.issubset(set(exchanges))


def test_penalty_scales_can_prioritize_equivalent_reaction_members() -> None:
    truth = np.zeros(len(REACTION_LABELS))
    truth[REACTION_LABELS.index("calcite")] = 0.2
    residual = reaction_matrix().T @ truth
    penalty_scales = np.full(len(REACTION_LABELS), 3.0)
    penalty_scales[REACTION_LABELS.index("calcite")] = 0.05
    penalty_scales[REACTION_LABELS.index("anorthite")] = 4.0

    fit = fit_inverse(
        residual,
        ION_ORDER,
        "elastic_net",
        lambda_l1=0.01,
        lambda_l2=0.005,
        penalty_scales=penalty_scales,
    )
    extents = np.asarray(fit["extents"])

    assert abs(extents[REACTION_LABELS.index("calcite")]) > abs(
        extents[REACTION_LABELS.index("anorthite")]
    )


def test_hydrosheaf_core_evidence_uses_cai_exchange_direction() -> None:
    upstream = np.asarray(
        [1.0, 0.5, 4.0, 0.1, 2.0, 3.0, 0.3, 0.2, 0.02, 0.002, 0.01],
        dtype=float,
    )
    residual = np.zeros(len(ION_ORDER), dtype=float)
    residual[ION_ORDER.index("Ca")] = 0.20
    residual[ION_ORDER.index("Mg")] = 0.08
    residual[ION_ORDER.index("Na")] = -0.45
    residual[ION_ORDER.index("Cl")] = -0.05

    scales, rows = m5.hydrosheaf_core_evidence(
        upstream,
        residual,
        ION_ORDER,
        {},
    )
    by_reaction = {str(row["reaction"]): row for row in rows}

    assert len(scales) == len(REACTION_LABELS)
    assert (
        by_reaction["CaNa_exch"]["hydrosheaf_core_evidence_score"]
        > by_reaction["NaCa_exch"]["hydrosheaf_core_evidence_score"]
    )
    assert by_reaction["CaNa_exch"]["penalty_scale"] < by_reaction["NaCa_exch"][
        "penalty_scale"
    ]


def test_plus_lite_diagnostics_prioritize_silicate_over_carbonate() -> None:
    diagnostics = ["SiO2", "Sr", "water_isotope_evaporation"]
    observed = np.asarray([0.4, 0.02, 0.0], dtype=float)

    scales, rows = m5._optional_diagnostic_evidence(diagnostics, observed)
    by_reaction = {str(row["reaction"]): row for row in rows}

    assert len(scales) == len(REACTION_LABELS)
    assert (
        by_reaction["anorthite"]["optional_evidence_score"]
        > by_reaction["calcite"]["optional_evidence_score"]
    )
    assert by_reaction["anorthite"]["optional_penalty_scale"] < by_reaction[
        "calcite"
    ]["optional_penalty_scale"]


def test_evidence_lifted_resolution_reports_class_preference_not_uniqueness() -> None:
    class_map, _ = equivalence_classes()
    diagnostics = ["SiO2", "Sr", "water_isotope_evaporation"]
    observed = np.asarray([0.4, 0.02, 0.0], dtype=float)
    _, optional_rows = m5._optional_diagnostic_evidence(diagnostics, observed)
    frame_rows = []
    for row in optional_rows:
        reaction = str(row["reaction"])
        frame_rows.append(
            {
                "scenario_id": "unit",
                "archetype": "unit",
                "replicate": 0,
                "data_tier": "plus_lite",
                "reaction": reaction,
                "equivalence_class": class_map[reaction],
                "combined_evidence_score": m5._combine_evidence_scores(
                    0.5,
                    row["optional_evidence_score"],
                ),
                "true_active": reaction == "anorthite",
                "recovered_active": reaction == "anorthite",
            }
        )

    resolution = m5._evidence_lifted_resolution_frame(
        m5.pd.DataFrame(frame_rows),
        class_map,
        score_column="combined_evidence_score",
        group_columns=["scenario_id", "archetype", "replicate", "data_tier"],
        evidence_source="unit_test",
    )
    ca_hco3 = resolution[resolution["members"].eq("anorthite;calcite")].iloc[0]

    assert ca_hco3["top_member"] == "anorthite"
    assert ca_hco3["resolution_status"] in {
        "conditionally_preferred",
        "evidence_lifted_resolved",
    }
    assert ca_hco3["top_member_true_active"]


def test_data_tier_simulation_is_deterministic_for_same_scenario() -> None:
    scenario = {
        "scenario_index": 7,
        "true_extents": np.zeros(len(REACTION_LABELS), dtype=float),
    }
    scenario["true_extents"][REACTION_LABELS.index("halite")] = 0.1

    first = m5._simulate_optional_diagnostics(scenario, "enhanced")
    second = m5._simulate_optional_diagnostics(scenario, "enhanced")

    assert first[0] == second[0]
    np.testing.assert_allclose(first[1], second[1])
    np.testing.assert_allclose(first[2], second[2])
