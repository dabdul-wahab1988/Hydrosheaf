"""Guard against publication-asset output-stem collisions.

Two distinct artifacts once shared the stem ``table1_benchmark_design``: the
hand-authored seven-audit design-and-claim map required by reviewer comment
M7-20260729-R04 and cited as TAB-1, and an unrelated three-column scale table
emitted by ``make_m7_3_publication_assets.py``.  Re-running the asset script
silently replaced the reviewer-requested table, and the manuscript still
assembled without error because the token mapping only resolves a path.

These tests make that class of failure loud:

1. no two asset generators may claim the same output stem;
2. no generator may write a hand-maintained stem; and
3. every stem the manuscript assembler resolves must actually exist.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1]
SCRIPTS = PACKAGE / "scripts"
TABLES = PACKAGE / "tables" / "publication"

ASSET_SCRIPTS = (
    "make_m7_3_publication_assets.py",
    "make_m7_sheaf_vs_graph_assets.py",
    "make_m7_robust_hybrid_assets.py",
)

# Authored by hand; no generator may produce them.
HAND_MAINTAINED_STEMS = frozenset(
    {
        "table1_benchmark_design",
        "table2_primary_m7_3_decision",
    }
)


def _written_stems(path: Path) -> set[str]:
    """Stems a script writes, via write_table(...) or a direct table*.csv path."""
    source = path.read_text(encoding="utf-8")
    stems: set[str] = set()

    tree = ast.parse(source)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = getattr(func, "id", None) or getattr(func, "attr", None)
        if name != "write_table":
            continue
        # write_table(frame, stem, title) - stem is the second positional arg.
        if len(node.args) >= 2 and isinstance(node.args[1], ast.Constant):
            value = node.args[1].value
            if isinstance(value, str):
                stems.add(value)

    for match in re.finditer(r"[\"'](table[A-Za-z0-9_]+)\.(?:csv|md)[\"']", source):
        stems.add(match.group(1))
    return stems


@pytest.fixture(scope="module")
def stems_by_script() -> dict[str, set[str]]:
    return {name: _written_stems(SCRIPTS / name) for name in ASSET_SCRIPTS}


def test_no_two_generators_claim_the_same_stem(stems_by_script):
    owners: dict[str, list[str]] = {}
    for script, stems in stems_by_script.items():
        for stem in stems:
            owners.setdefault(stem, []).append(script)
    collisions = {s: sorted(o) for s, o in owners.items() if len(o) > 1}
    assert not collisions, (
        "Publication asset stems are written by more than one generator; "
        f"re-running either would silently clobber the other: {collisions}"
    )


def test_no_generator_writes_a_hand_maintained_stem(stems_by_script):
    violations = {
        script: sorted(stems & HAND_MAINTAINED_STEMS)
        for script, stems in stems_by_script.items()
        if stems & HAND_MAINTAINED_STEMS
    }
    assert not violations, (
        "A generator writes a hand-authored table. TAB-1 "
        "(table1_benchmark_design) is the seven-audit claim map required by "
        f"reviewer M7-20260729-R04 and must not be generated: {violations}"
    )


def test_write_table_refuses_hand_maintained_stems():
    """The runtime guard, not just the static scan."""
    pytest.importorskip("pandas")
    import sys

    sys.path.insert(0, str(SCRIPTS))
    import pandas as pd

    from make_m7_3_publication_assets import (  # noqa: E402
        HAND_MAINTAINED_STEMS as runtime_protected,
        write_table,
    )

    assert "table1_benchmark_design" in runtime_protected
    with pytest.raises(RuntimeError, match="hand-maintained"):
        write_table(pd.DataFrame({"a": [1]}), "table1_benchmark_design", "x")


def test_hand_maintained_tables_exist_and_are_not_the_generated_variant():
    table1 = TABLES / "table1_benchmark_design.md"
    assert table1.exists(), "TAB-1 source table is missing"
    text = table1.read_text(encoding="utf-8")
    # The generated variant's title and column header; either means it was
    # clobbered by make_m7_3_publication_assets.py.
    assert "Locked benchmark design" not in text, (
        "table1_benchmark_design.md has been overwritten by the generated "
        "scale table; restore the seven-audit design-and-claim map."
    )
    assert "Design item" not in text, (
        "table1_benchmark_design.md carries the generated table's column "
        "header; it has been clobbered."
    )
    assert "Audit" in text and "Permitted claim" in text


def test_assembler_table_tokens_resolve_to_existing_files():
    assembler = (SCRIPTS / "assemble_m7_manuscript.py").read_text(encoding="utf-8")
    referenced = set(re.findall(r"[\"'](table[A-Za-z0-9_]+)\.md[\"']", assembler))
    assert referenced, "no table stems found in the assembler"
    missing = sorted(s for s in referenced if not (TABLES / f"{s}.md").exists())
    assert not missing, f"assembler references missing table files: {missing}"
