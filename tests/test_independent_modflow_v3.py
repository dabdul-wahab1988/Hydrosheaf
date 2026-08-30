from __future__ import annotations

import ast
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from independent_modflow_generator_v3 import (  # noqa: E402
    generate_independent_aquifer_v3,
)


MF6 = ROOT / ".codex_work" / "modflow-bin" / "mf6.exe"
MP7 = ROOT / ".codex_work" / "modflow-bin" / "mp7.exe"


def test_v3_source_is_independent_of_inference_and_historical_generators() -> None:
    source_path = SCRIPT_DIR / "independent_modflow_generator_v3.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imported_modules = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.append(str(node.module or ""))
    assert not any(module == "hydrosheaf" or module.startswith("hydrosheaf.") for module in imported_modules)
    assert "independent_modflow_generator" not in imported_modules
    assert "independent_modflow_generator_v2" not in imported_modules


@pytest.mark.skipif(
    not MF6.exists() or not MP7.exists(),
    reason="configured MODFLOW 6/MODPATH 7 executables are unavailable",
)
def test_v3_solver_case_has_independent_layered_branch_merge_truth(tmp_path: Path) -> None:
    case = generate_independent_aquifer_v3(9101, tmp_path / "case", MF6, MP7)
    rows = {str(row["site_id"]): row for row in case.observations}

    assert len(rows) == 6
    assert len({(float(row["x_m"]), float(row["y_m"])) for row in rows.values()}) == 6
    assert {int(row["aquifer_layer"]) for row in rows.values()} == {1, 2}
    assert int(case.provenance["nlay"]) == 2
    assert bool(case.provenance["layered_flow"])
    assert bool(case.provenance["vertical_exchange"])
    assert float(case.provenance["vertical_exchange_abs_flow_sum"]) > 0.0
    assert case.provenance.get("base_generator_module") is None
    assert case.provenance["mf6_sha256"]
    assert case.provenance["mp7_sha256"]

    outgoing: dict[str, set[str]] = defaultdict(set)
    incoming: dict[str, set[str]] = defaultdict(set)
    for source, target in case.true_edges:
        outgoing[source].add(target)
        incoming[target].add(source)
        assert source in rows and target in rows
        assert case.true_ages_years[target] > case.true_ages_years[source]
        assert float(rows[source]["head_meas"]) > float(rows[target]["head_meas"])
    assert any(len(targets) >= 2 for targets in outgoing.values())
    assert any(len(sources) >= 2 for sources in incoming.values())

    path_edges: set[tuple[str, str]] = set()
    by_particle: dict[int, list[dict[str, object]]] = defaultdict(list)
    for pathline_row in case.pathline_rows:
        by_particle[int(pathline_row["particle"])].append(dict(pathline_row))
    for particle_rows in by_particle.values():
        particle_rows.sort(key=lambda row: int(row["milestone"]))
        path_edges.update(
            (str(left["node_id"]), str(right["node_id"]))
            for left, right in zip(particle_rows, particle_rows[1:])
        )
    assert set(case.true_edges).issubset(path_edges)
    assert {int(row["layer_zero_based"]) for row in case.pathline_rows} == {0, 1}
    assert any(
        int(value) > 0
        for value in case.provenance["pathline_layer_transitions"].values()
    )

    suspicious_tokens = ("age", "travel", "particle", "pathline", "process", "edge", "truth")
    assert not any(
        any(token in str(key).lower() for token in suspicious_tokens)
        for row in case.observations
        for key in row
    )

    for node, component_ages in case.true_age_components.items():
        weights = np.asarray(case.true_mixture_weights[node], dtype=float)
        assert len(component_ages) == len(weights)
        assert np.isclose(weights.sum(), 1.0)
        assert np.isclose(
            case.true_ages_years[node],
            np.average(np.asarray(component_ages, dtype=float), weights=weights),
        )


@pytest.mark.skipif(
    not MF6.exists() or not MP7.exists(),
    reason="configured MODFLOW 6/MODPATH 7 executables are unavailable",
)
def test_v3_same_seed_reproduces_and_alternate_seed_changes(tmp_path: Path) -> None:
    first = generate_independent_aquifer_v3(9101, tmp_path / "first", MF6, MP7)
    repeat = generate_independent_aquifer_v3(9101, tmp_path / "repeat", MF6, MP7)
    alternate = generate_independent_aquifer_v3(9102, tmp_path / "alternate", MF6, MP7)

    assert first.observations == repeat.observations
    assert first.true_edges == repeat.true_edges
    assert first.true_ages_years == repeat.true_ages_years
    assert first.observations != alternate.observations
