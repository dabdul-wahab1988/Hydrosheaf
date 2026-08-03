from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from independent_lattice_generator import generate_independent_lattice  # noqa: E402


def test_lattice_generator_is_reproducible_and_truth_blind() -> None:
    first = generate_independent_lattice(9301)
    second = generate_independent_lattice(9301)

    assert first.observations == second.observations
    assert first.true_edges == second.true_edges
    assert len(first.observations) == 12
    assert len(first.true_edges) == 13
    assert first.provenance["imports_hydrosheaf"] is False
    assert first.provenance["generator_family"] == "analytic_lattice"
    assert all(
        not any(str(key).lower().startswith(("true_", "truth_")) for key in row)
        for row in first.observations
    )


def test_lattice_generator_changes_with_seed_and_has_process_truth() -> None:
    first = generate_independent_lattice(9301)
    other = generate_independent_lattice(9302)

    assert first.observations != other.observations
    assert len(first.true_processes) == 13
    assert set(first.true_processes.values()) >= {
        "carbonate_weathering",
        "silicate_weathering",
        "denitrification",
        "sulfate_reduction",
        "iron_reduction",
    }
