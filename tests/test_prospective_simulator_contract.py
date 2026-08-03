"""Mutation tests for the prospective synthetic measurement contract."""

from __future__ import annotations

import sys
from pathlib import Path

from hydrosheaf.validation.decision_utility import CandidateMeasurementAction

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from run_programme_generator_ensemble import _prospective_realised_outcomes  # noqa: E402


def test_outcome_simulator_is_invariant_to_selector_quality_mutation():
    base_action = CandidateMeasurementAction(
        action_id="head:A->B",
        cost=1.0,
        metadata={"observation_informativeness": 0.05},
    )
    mutated_action = CandidateMeasurementAction(
        action_id="head:A->B",
        cost=1.0,
        metadata={"observation_informativeness": 0.99},
    )

    base = _prospective_realised_outcomes(
        [{"site_id": "A", "head_meas": 10.0}],
        (base_action,),
        (("A", "B"),),
        seed=1234,
    )
    mutated = _prospective_realised_outcomes(
        [{"site_id": "A", "head_meas": -999.0}],
        (mutated_action,),
        (("A", "B"),),
        seed=1234,
    )

    assert base == mutated
    assert base[2]["reuses_selection_quality"] is False
    assert base[2]["rule"] == "independent_state_conditioned_measurement_v2"
