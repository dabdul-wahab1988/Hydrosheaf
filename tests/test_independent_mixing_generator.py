from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from independent_mixing_generator import generate_independent_mixing  # noqa: E402


def test_mixing_generator_is_reproducible_and_multilayered() -> None:
    first = generate_independent_mixing(9501)
    second = generate_independent_mixing(9501)

    assert first.observations == second.observations
    assert first.true_edges == second.true_edges
    assert len(first.observations) == 11
    assert len(first.true_edges) == 14
    assert first.provenance["imports_hydrosheaf"] is False
    assert first.provenance["generator_family"] == "independent_mixing"
    assert len({row["aquifer_layer"] for row in first.observations}) == 3
    assert all(
        not any(str(key).lower().startswith(("true_", "truth_")) for key in row)
        for row in first.observations
    )


def test_mixing_generator_changes_seed_and_contains_shortcut_processes() -> None:
    first = generate_independent_mixing(9501)
    other = generate_independent_mixing(9502)

    assert first.observations != other.observations
    assert len(first.true_processes) == 14
    assert "mixing" in set(first.true_processes.values())
    assert set(first.true_processes.values()) >= {
        "carbonate_weathering",
        "silicate_weathering",
        "denitrification",
        "sulfate_reduction",
        "iron_reduction",
    }
