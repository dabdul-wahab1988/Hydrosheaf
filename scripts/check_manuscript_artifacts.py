"""Guard against silent artifact rot in the manuscript benchmarks.

Three failure modes were found during the 2026-07-27 audit, all of the same kind:
a committed table or figure no longer corresponded to the code that claims to
produce it, and nothing detected it.

  1. M4's topology-posterior columns were produced before a 2026-07-19 change to
     ``hydrosheaf/inference/topology_posterior.py`` and no longer reproduced.
  2. M3's Table 6 silently kept stale values because the statistics function
     returns early when its per-scenario input files are missing.
  3. The M3 table builders concatenate per-scenario CSVs *before* the
     design-matrix output with ``keep="first"``, so a leftover per-scenario file
     silently overrides freshly regenerated results.

This check re-runs the table builders into a temporary directory and diffs the
output against the tracked tables.  It does not re-run the benchmarks: it
verifies that the *tables* are the ones the current code produces from the
*committed* results.  That is the link that broke in all three cases.

Usage
-----
    python scripts/check_manuscript_artifacts.py            # check M3 and M4
    python scripts/check_manuscript_artifacts.py --module M4
    python scripts/check_manuscript_artifacts.py --update   # refresh tracked tables

Exit status is non-zero when a tracked table does not match, so this can be wired
into CI or a pre-commit hook.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import pandas as pd

# Comparison uses .clip(); silence the unrelated pandas downcasting notice.
warnings.simplefilter("ignore", FutureWarning)

REPO = Path(__file__).resolve().parents[1]

MODULES = {
    "M3": {
        "builder": REPO / "M3/m3_age_benchmark/scripts/make_publication_tables.py",
        "tables": REPO / "M3/m3_age_benchmark/tables/Manuscript_Ready",
    },
    "M4": {
        "builder": REPO / "M4/m4_topology_benchmark/scripts/make_publication_tables.py",
        "tables": REPO / "M4/m4_topology_benchmark/tables/Manuscript_Ready",
    },
}

# Columns that legitimately vary between runs (timestamps, absolute paths).
VOLATILE = {"run_utc", "generated_at", "path", "source_path"}
# Tolerance for float comparison; tighter than any reported precision.
RTOL = 1e-9


def _compare(tracked: Path, fresh: Path) -> list[str]:
    problems: list[str] = []
    try:
        a = pd.read_csv(tracked)
        b = pd.read_csv(fresh)
    except Exception as exc:  # unreadable file is itself a failure
        return [f"{tracked.name}: unreadable ({exc})"]

    a = a.drop(columns=[c for c in a.columns if c in VOLATILE], errors="ignore")
    b = b.drop(columns=[c for c in b.columns if c in VOLATILE], errors="ignore")

    if list(a.columns) != list(b.columns):
        problems.append(
            f"{tracked.name}: column set changed\n"
            f"    tracked: {list(a.columns)}\n"
            f"    fresh  : {list(b.columns)}"
        )
        return problems
    if len(a) != len(b):
        problems.append(f"{tracked.name}: row count {len(a)} tracked vs {len(b)} regenerated")
        return problems

    for col in a.columns:
        ca, cb = a[col], b[col]
        # Booleans are numeric to pandas but do not support subtraction, so
        # compare them as strings alongside the other non-float columns.
        both_float = (
            pd.api.types.is_numeric_dtype(ca)
            and pd.api.types.is_numeric_dtype(cb)
            and not pd.api.types.is_bool_dtype(ca)
            and not pd.api.types.is_bool_dtype(cb)
        )
        if both_float:
            close = ((ca - cb).abs() <= RTOL * cb.abs().clip(lower=1e-30)) | (ca.isna() & cb.isna())
            bad = (~close).to_numpy().nonzero()[0]
        else:
            bad = (ca.astype(str).ne(cb.astype(str)) & ~(ca.isna() & cb.isna())).to_numpy().nonzero()[0]
        if len(bad):
            i = int(bad[0])
            problems.append(
                f"{tracked.name}: column '{col}' differs at row {i} "
                f"(tracked={ca.iloc[i]!r}, regenerated={cb.iloc[i]!r}; {len(bad)} row(s) differ)"
            )
    return problems


def check(module: str, update: bool = False) -> int:
    cfg = MODULES[module]
    builder, tracked_dir = cfg["builder"], cfg["tables"]
    if not builder.exists():
        print(f"[{module}] SKIP - builder not found: {builder}")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        # The builders write to a fixed location, so snapshot the tracked tables,
        # run the builder, compare, then put the snapshot back.
        #
        # Restoration copies file-by-file and NEVER deletes the tracked directory.
        # An earlier version used shutil.rmtree() + copytree(); on Windows the
        # rmdir failed after the contents had already been removed, destroying
        # the tracked tables. A checking tool must not be able to lose data.
        backup = Path(tmp) / "backup"
        backup.mkdir(parents=True)
        for f in tracked_dir.glob("*.csv"):
            shutil.copy2(f, backup / f.name)

        def restore() -> None:
            for saved in backup.glob("*.csv"):
                shutil.copy2(saved, tracked_dir / saved.name)
            # Remove only files the builder newly created, leaving the dir intact.
            for current in tracked_dir.glob("*.csv"):
                if not (backup / current.name).exists():
                    current.unlink()

        proc = subprocess.run(
            [sys.executable, str(builder)], capture_output=True, text=True, cwd=str(REPO)
        )
        if proc.returncode != 0:
            restore()
            print(f"[{module}] FAIL - table builder exited {proc.returncode}")
            print(proc.stderr[-1500:])
            return 1

        problems: list[str] = []
        regenerated = sorted(tracked_dir.glob("*.csv"))
        for fresh in regenerated:
            old = backup / fresh.name
            if not old.exists():
                problems.append(f"{fresh.name}: produced by the builder but NOT tracked")
                continue
            problems.extend(_compare(old, fresh))
        for old in sorted(backup.glob("*.csv")):
            if not (tracked_dir / old.name).exists():
                problems.append(f"{old.name}: tracked but the builder no longer produces it")

        if update:
            print(f"[{module}] UPDATED {len(regenerated)} table(s) from current code")
            return 0

        restore()

    if problems:
        print(f"[{module}] FAIL - {len(problems)} artifact(s) do not reproduce:")
        for p in problems:
            print("  -", p)
        return 1
    print(f"[{module}] OK - all tracked tables reproduce from committed results")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--module", choices=sorted(MODULES), action="append",
                    help="limit to one module (repeatable); default is all")
    ap.add_argument("--update", action="store_true",
                    help="rewrite tracked tables from current code instead of checking")
    args = ap.parse_args(argv)

    status = 0
    for module in (args.module or sorted(MODULES)):
        status |= check(module, update=args.update)
    if status and not args.update:
        print("\nA tracked table no longer matches what the current code produces from the "
              "committed results. Either the results were regenerated without refreshing "
              "the tables, a stale per-scenario file is shadowing them, or the builder "
              "changed. Investigate before trusting any manuscript number.")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
