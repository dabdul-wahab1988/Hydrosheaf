from __future__ import annotations

from collections import Counter

import pytest

from hydrosheaf.validation.controls import (
    apply_field_permutation,
    make_no_flow_control,
)
from hydrosheaf.validation.programme_contract import assert_truth_blind


def _rows() -> list[dict[str, object]]:
    return [
        {"site_id": "A", "x_m": 0.0, "head_meas": 100.0, "tritium_TU": 1.0},
        {"site_id": "B", "x_m": 1.0, "head_meas": 101.0, "tritium_TU": 2.0},
        {"site_id": "C", "x_m": 2.0, "head_meas": 102.0, "tritium_TU": None},
        {"site_id": "D", "x_m": 3.0, "head_meas": 103.0},
    ]


def test_permutation_is_deterministic_and_row_safe() -> None:
    first = apply_field_permutation(_rows(), ["head_meas"], seed=17)
    second = apply_field_permutation(_rows(), ["head_meas"], seed=17)
    assert first.observations == second.observations
    assert [row["site_id"] for row in first.observations] == ["A", "B", "C", "D"]
    assert [row["x_m"] for row in first.observations] == [0.0, 1.0, 2.0, 3.0]
    assert Counter(row["head_meas"] for row in first.observations) == Counter(
        row["head_meas"] for row in _rows()
    )
    assert first.to_dict()["changed_values"] >= 0


def test_no_flow_control_preserves_marginals_and_truth_blindness() -> None:
    result = make_no_flow_control(_rows(), seed=91, head_fields=("head_meas",))
    assert result.scenario == "no_flow_head_permutation"
    assert result.fields == ("head_meas",)
    assert_truth_blind(result.observations)
    assert all(row["observation_flags"]["control_seed"] == 91 for row in result.observations)


def test_missing_field_is_allowed_but_empty_field_is_rejected() -> None:
    result = apply_field_permutation(_rows(), ["not_present"], seed=3)
    assert result.preserved_marginal_counts == {"not_present": 0}
    with pytest.raises(ValueError, match="At least one"):
        apply_field_permutation(_rows(), [], seed=3)


def test_truth_fields_are_rejected_before_control() -> None:
    with pytest.raises(ValueError, match="truth fields"):
        apply_field_permutation([{"site_id": "A", "true_age": 4.0}], ["site_id"], seed=1)
