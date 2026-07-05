from __future__ import annotations

import numpy as np

from hydrosheaf.models.evidence_lifted import (
    evidence_lifted_resolution,
    stoichiometric_equivalence_classes,
)


def test_equal_evidence_keeps_equivalence_class_unresolved() -> None:
    result = evidence_lifted_resolution(
        ["calcite", "anorthite"],
        {"calcite": 0.5, "anorthite": 0.5},
        class_id="EC05",
    )

    assert result.status == "unresolved_equivalence_class"
    assert result.evidence_lifted_resolution_index == 0.0
    assert result.probabilities["calcite"] == result.probabilities["anorthite"]


def test_dominant_evidence_resolves_equivalence_class_conditionally() -> None:
    result = evidence_lifted_resolution(
        ["calcite", "anorthite"],
        {"calcite": 0.1, "anorthite": 0.9},
        class_id="EC05",
    )

    assert result.top_member == "anorthite"
    assert result.status == "evidence_lifted_resolved"
    assert result.evidence_lifted_resolution_index > 0.5


def test_stoichiometric_equivalence_classes_groups_signed_opposites() -> None:
    matrix = np.asarray(
        [
            [1.0, 2.0],
            [1.0, 2.0],
            [1.0, -2.0],
            [-1.0, 2.0],
        ]
    )
    labels = ["calcite", "anorthite", "CaNa_exch", "NaCa_exch"]

    class_map, rows = stoichiometric_equivalence_classes(matrix, labels)

    assert class_map["calcite"] == class_map["anorthite"]
    assert class_map["CaNa_exch"] == class_map["NaCa_exch"]
    assert sum(1 for row in rows if row["ambiguous"]) == 2
