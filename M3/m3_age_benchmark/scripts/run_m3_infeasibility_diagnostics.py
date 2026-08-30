"""Locate the cause of identified-TTD infeasibility in the M3 local baseline.

Two complementary diagnostics are computed over the same rows and the same
frozen protocol used by ``run_m3_identified_ttd_benchmark.py``:

ENVELOPE (individual reproducibility).  For tracer ``i`` the forward response
``A[i, :] @ w`` over the unit simplex is bounded by
``[min_j A[i, j], max_j A[i, j]]``.  An observation interval
``[obs - k*sigma, obs + k*sigma]`` that misses that range is reproducible by no
non-negative unit-mass TTD, at any k.  Censored likelihoods are one-sided and
are handled exactly as ``_compile_constraints`` handles them.

MINIMAL INFEASIBLE SUBSET (where joint inconsistency lives).  For every site
carrying at least three constraints, all singletons and all pairs are tested
for feasibility, and for infeasible full panels the smallest infeasible subset
is found by exhaustive search over increasing subset size.  Mass at size 2
indicates pairwise tracer conflict; mass at size >= 3 would indicate that no
common TTD exists.

Every row and every fold that is not examined is counted against an explicit
skip reason and reconciled against the runner-equivalent eligible-fold count.
A diagnostic that quantifies failure rates may not discard rows silently.

This script performs no inference beyond feasibility testing and authorizes no
scientific claim; it is a development diagnostic for
``docs/m3_identified_ttd_infeasibility_audit_20260731.md``.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import linprog

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.nuclear.ttd_identified import (  # noqa: E402
    TracerConstraint,
    _compile_constraints,
    build_tracer_constraints,
)
import run_m3_identified_ttd_benchmark as runner  # noqa: E402
import run_m3_usgs_benchmark as usgs  # noqa: E402

BENCHMARK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = runner.DEFAULT_PROTOCOL
DEFAULT_OUTPUT = BENCHMARK_DIR / "results" / "m3_infeasibility_diagnostics.json"
DEFAULT_K_VALUES: tuple[float, ...] = (1.96, 6.0)

# Minimum panel size for the minimal-infeasible-subset diagnostic.  Below three
# constraints the subset search cannot distinguish pairwise from higher-order
# conflict, which is the question the diagnostic exists to answer.
MIS_MINIMUM_PANEL_SIZE = 3
# Pairs seen fewer times than this are not reported as rates.
PAIR_REPORTING_MINIMUM = 20


# --------------------------------------------------------------------------
# Pure geometry: no dataset, no protocol.  These are the unit-tested core.
# --------------------------------------------------------------------------


def response_envelope(constraint: TracerConstraint) -> tuple[float, float]:
    """Achievable range of ``response @ w`` over the unit simplex."""
    response = np.asarray(constraint.response, dtype=float)
    return float(response.min()), float(response.max())


def violates_envelope(constraint: TracerConstraint, sigma_multiplier: float) -> bool:
    """Is the declared observation interval disjoint from the achievable range?

    The interval semantics mirror ``_compile_constraints`` exactly: a Gaussian
    observation is two-sided, ``upper_censored`` constrains only the upper
    limit ``observed + margin``, and ``lower_censored`` only the lower limit
    ``observed - margin``.
    """
    lower, upper = response_envelope(constraint)
    margin = float(sigma_multiplier) * float(constraint.sigma)
    observed = float(constraint.observed)
    likelihood = str(constraint.likelihood).strip().lower()
    if likelihood == "gaussian":
        return (observed + margin < lower) or (observed - margin > upper)
    if likelihood == "upper_censored":
        return observed + margin < lower
    if likelihood == "lower_censored":
        return observed - margin > upper
    return False


def is_feasible(
    constraints: Sequence[TracerConstraint],
    n_age_bins: int,
    sigma_multiplier: float,
) -> bool:
    """Does a non-negative unit-mass ``w`` satisfy every declared interval?"""
    a_ub, b_ub = _compile_constraints(
        constraints, n_age_bins, float(sigma_multiplier)
    )
    result = linprog(
        np.zeros(n_age_bins, dtype=float),
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=np.ones((1, n_age_bins), dtype=float),
        b_eq=np.ones(1, dtype=float),
        bounds=(0.0, None),
        method="highs",
    )
    return bool(result.success)


def minimal_infeasible_subset(
    constraints: Sequence[TracerConstraint],
    n_age_bins: int,
    sigma_multiplier: float,
) -> tuple[str, ...] | None:
    """Smallest infeasible subset of ``constraints``, or ``None`` if feasible.

    Searched by increasing subset size, so the first infeasible subset found is
    minimal by cardinality.
    """
    panel = list(constraints)
    for size in range(1, len(panel) + 1):
        for combination in combinations(range(len(panel)), size):
            subset = [panel[index] for index in combination]
            if not is_feasible(subset, n_age_bins, sigma_multiplier):
                return tuple(panel[index].tracer for index in combination)
    return None


# --------------------------------------------------------------------------
# Dataset-facing diagnostics
# --------------------------------------------------------------------------


class SkipLedger:
    """Counted, named skip reasons; nothing leaves a diagnostic unaccounted."""

    def __init__(self) -> None:
        self._counts: Counter[str] = Counter()

    def record(self, reason: str, count: int = 1) -> None:
        self._counts[reason] += int(count)

    def record_exception(self, stage: str, exc: BaseException) -> str:
        reason = f"{stage}_error:{type(exc).__name__}"
        self._counts[reason] += 1
        return reason

    @property
    def total(self) -> int:
        return int(sum(self._counts.values()))

    def as_dict(self) -> dict[str, int]:
        return {key: int(value) for key, value in sorted(self._counts.items())}


def _build_row_constraints(
    row: Mapping[str, Any],
    observations: Sequence[Any],
    provenance: Mapping[str, Any],
    protocol: Mapping[str, Any],
    ages: np.ndarray,
    sample_year: float,
) -> tuple[TracerConstraint, ...]:
    """Constraint construction identical to the runner's inference path."""
    histories, _ = usgs._get_site_histories(row)
    constraints, _excluded = build_tracer_constraints(
        observations,
        ages,
        sample_year,
        histories=histories,
        use_helium4=bool(protocol["observations"].get("use_helium4", False)),
        initial_c14_pmc=provenance["initial_c14_pmc"],
        helium4_background_ccpg=(
            runner._finite_float(row.get("he4_background_ccpg")) or 0.0
        ),
        helium4_accumulation_rate_ccpg_per_year=(
            runner._finite_float(row.get("he4_accumulation_rate_ccpg_per_year"))
            or 1.0e-11
        ),
    )
    return constraints


def run_envelope_diagnostic(
    rows: Sequence[Mapping[str, Any]],
    protocol: Mapping[str, Any],
    ages: np.ndarray,
    k_values: Sequence[float],
) -> dict[str, Any]:
    """Per-fold envelope violations with a full skip reconciliation."""
    held_out_tracers = list(protocol["observations"]["held_out_tracers"])
    minimum_conditioning = int(
        protocol["observations"].get("minimum_conditioning_tracers", 1)
    )
    ledger = SkipLedger()
    row_ledger = SkipLedger()

    folds_examined = 0
    # Runner-equivalent eligibility: a fold is eligible when the row has a
    # finite sample year and the held-out observation exists.  When the row's
    # observations cannot be prepared at all the runner still reports every
    # tracer for that row as an eligible ABSTAIN, so they are counted here too.
    eligible_folds = 0
    total_constraints = Counter()  # keyed by k
    violating_constraints = Counter()
    folds_with_violation = Counter()
    violations_by_tracer: dict[float, Counter] = {k: Counter() for k in k_values}
    seen_by_tracer: dict[float, Counter] = {k: Counter() for k in k_values}

    for row in rows:
        sample_year = runner._finite_float(row.get("sample_year"))
        if sample_year is None:
            row_ledger.record("missing_sample_year")
            continue
        try:
            observations, provenance = runner.prepare_row_observations(row, protocol)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            reason = row_ledger.record_exception("prepare_row_observations", exc)
            # The runner would emit one eligible ABSTAIN fold per tracer here.
            eligible_folds += len(held_out_tracers)
            ledger.record(reason, len(held_out_tracers))
            continue

        for held_out_tracer in held_out_tracers:
            conditioning, held_out = runner.split_held_out_observation(
                observations, held_out_tracer
            )
            if held_out is None:
                # Not eligible in the runner either; not a skip.
                continue
            eligible_folds += 1
            if len(conditioning) < minimum_conditioning:
                ledger.record("insufficient_conditioning_tracers")
                continue
            try:
                constraints = _build_row_constraints(
                    row, conditioning, provenance, protocol, ages, sample_year
                )
            except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
                ledger.record_exception("build_tracer_constraints", exc)
                continue
            if not constraints:
                ledger.record("no_linear_interval_constraints")
                continue

            folds_examined += 1
            for k in k_values:
                n_violations = 0
                for constraint in constraints:
                    seen_by_tracer[k][constraint.tracer] += 1
                    if violates_envelope(constraint, k):
                        n_violations += 1
                        violations_by_tracer[k][constraint.tracer] += 1
                total_constraints[k] += len(constraints)
                violating_constraints[k] += n_violations
                if n_violations:
                    folds_with_violation[k] += 1

    by_k: dict[str, Any] = {}
    for k in k_values:
        by_k[_k_key(k)] = {
            "sigma_multiplier": float(k),
            "folds_examined": int(folds_examined),
            "folds_with_envelope_violation": int(folds_with_violation[k]),
            "violating_constraints": int(violating_constraints[k]),
            "total_constraints": int(total_constraints[k]),
            "by_tracer": {
                tracer: {
                    "violations": int(violations_by_tracer[k][tracer]),
                    "seen": int(seen),
                }
                for tracer, seen in sorted(seen_by_tracer[k].items())
            },
        }
    accounted = ledger.total
    return {
        "by_k": by_k,
        "accounting": {
            "n_rows": int(len(rows)),
            "row_skips": row_ledger.as_dict(),
            "eligible_folds": int(eligible_folds),
            "folds_examined": int(folds_examined),
            "folds_skipped": int(accounted),
            "fold_skips": ledger.as_dict(),
            "reconciled": bool(eligible_folds == folds_examined + accounted),
        },
    }


def run_subset_diagnostic(
    rows: Sequence[Mapping[str, Any]],
    protocol: Mapping[str, Any],
    ages: np.ndarray,
    k_values: Sequence[float],
) -> dict[str, Any]:
    """Singleton, pairwise and minimal-infeasible-subset structure per site."""
    n_bins = int(ages.size)
    ledger = SkipLedger()
    panels: list[tuple[str, tuple[TracerConstraint, ...]]] = []

    for row in rows:
        sample_year = runner._finite_float(row.get("sample_year"))
        if sample_year is None:
            ledger.record("missing_sample_year")
            continue
        try:
            observations, provenance = runner.prepare_row_observations(row, protocol)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            ledger.record_exception("prepare_row_observations", exc)
            continue
        try:
            constraints = _build_row_constraints(
                row, observations, provenance, protocol, ages, sample_year
            )
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            ledger.record_exception("build_tracer_constraints", exc)
            continue
        if len(constraints) < MIS_MINIMUM_PANEL_SIZE:
            ledger.record(
                f"fewer_than_{MIS_MINIMUM_PANEL_SIZE}_constraints"
            )
            continue
        panels.append((str(row.get("site_id", "")), constraints))

    by_k: dict[str, Any] = {}
    for k in k_values:
        singleton_seen: Counter[str] = Counter()
        singleton_infeasible: Counter[str] = Counter()
        pair_seen: Counter[tuple[str, str]] = Counter()
        pair_infeasible: Counter[tuple[str, str]] = Counter()
        mis_size: Counter[int] = Counter()
        mis_combination: Counter[tuple[str, ...]] = Counter()
        n_full_infeasible = 0
        n_infeasible_without_mis = 0

        for _site_id, constraints in panels:
            for constraint in constraints:
                singleton_seen[constraint.tracer] += 1
                if not is_feasible([constraint], n_bins, k):
                    singleton_infeasible[constraint.tracer] += 1
            for i, j in combinations(range(len(constraints)), 2):
                key = tuple(sorted((constraints[i].tracer, constraints[j].tracer)))
                pair_seen[key] += 1
                if not is_feasible([constraints[i], constraints[j]], n_bins, k):
                    pair_infeasible[key] += 1
            if is_feasible(list(constraints), n_bins, k):
                continue
            n_full_infeasible += 1
            subset = minimal_infeasible_subset(list(constraints), n_bins, k)
            if subset is None:
                # Only reachable if feasibility is not reproducible between the
                # full-panel solve and the subset search; counted, never hidden.
                n_infeasible_without_mis += 1
                continue
            mis_size[len(subset)] += 1
            mis_combination[tuple(sorted(subset))] += 1

        by_k[_k_key(k)] = {
            "sigma_multiplier": float(k),
            "panels_examined": int(len(panels)),
            "panels_infeasible": int(n_full_infeasible),
            "panels_infeasible_without_minimal_subset": int(n_infeasible_without_mis),
            "minimal_infeasible_subset_size_counts": {
                str(size): int(count) for size, count in sorted(mis_size.items())
            },
            "minimal_infeasible_subset_combinations": {
                "+".join(combo): int(count)
                for combo, count in sorted(
                    mis_combination.items(), key=lambda item: (-item[1], item[0])
                )
            },
            "singleton": {
                tracer: {
                    "infeasible": int(singleton_infeasible[tracer]),
                    "seen": int(seen),
                }
                for tracer, seen in sorted(singleton_seen.items())
            },
            "pairwise": {
                "+".join(key): {
                    "infeasible": int(pair_infeasible[key]),
                    "seen": int(seen),
                }
                for key, seen in sorted(pair_seen.items())
                if seen >= PAIR_REPORTING_MINIMUM
            },
            "pairwise_reporting_minimum": PAIR_REPORTING_MINIMUM,
        }

    return {
        "by_k": by_k,
        "accounting": {
            "n_rows": int(len(rows)),
            "panels_examined": int(len(panels)),
            "rows_skipped": int(ledger.total),
            "row_skips": ledger.as_dict(),
            "reconciled": bool(len(rows) == len(panels) + ledger.total),
        },
    }


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _k_key(k: float) -> str:
    return f"k={float(k):g}"


def _rate(numerator: int, denominator: int) -> float | None:
    return float(numerator) / float(denominator) if denominator else None


def flatten_to_records(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Long-format rows for the CSV rendering of the same result."""
    records: list[dict[str, Any]] = []

    def emit(
        diagnostic: str,
        k: Any,
        metric: str,
        key: str,
        count: int,
        total: int | None,
    ) -> None:
        records.append(
            {
                "diagnostic": diagnostic,
                "sigma_multiplier": k,
                "metric": metric,
                "key": key,
                "count": int(count),
                "total": ("" if total is None else int(total)),
                "rate": ("" if total in (None, 0) else _rate(count, int(total))),
            }
        )

    envelope = payload.get("envelope", {})
    for block in envelope.get("by_k", {}).values():
        k = block["sigma_multiplier"]
        emit(
            "envelope",
            k,
            "folds_with_envelope_violation",
            "all",
            block["folds_with_envelope_violation"],
            block["folds_examined"],
        )
        emit(
            "envelope",
            k,
            "violating_constraints",
            "all",
            block["violating_constraints"],
            block["total_constraints"],
        )
        for tracer, stats in block["by_tracer"].items():
            emit(
                "envelope",
                k,
                "violating_constraints_by_tracer",
                tracer,
                stats["violations"],
                stats["seen"],
            )
    for reason, count in envelope.get("accounting", {}).get("fold_skips", {}).items():
        emit("envelope", "", "fold_skip", reason, count, None)
    for reason, count in envelope.get("accounting", {}).get("row_skips", {}).items():
        emit("envelope", "", "row_skip", reason, count, None)
    if envelope:
        accounting = envelope["accounting"]
        emit("envelope", "", "fold_total", "eligible", accounting["eligible_folds"], None)
        emit("envelope", "", "fold_total", "examined", accounting["folds_examined"], None)
        emit("envelope", "", "fold_total", "skipped", accounting["folds_skipped"], None)

    subset = payload.get("subset", {})
    for block in subset.get("by_k", {}).values():
        k = block["sigma_multiplier"]
        emit(
            "subset",
            k,
            "panels_infeasible",
            "all",
            block["panels_infeasible"],
            block["panels_examined"],
        )
        n_mis = sum(block["minimal_infeasible_subset_size_counts"].values())
        for size, count in block["minimal_infeasible_subset_size_counts"].items():
            emit("subset", k, "minimal_infeasible_subset_size", size, count, n_mis)
        for combo, count in block["minimal_infeasible_subset_combinations"].items():
            emit("subset", k, "minimal_infeasible_subset_combination", combo, count, n_mis)
        for tracer, stats in block["singleton"].items():
            emit(
                "subset", k, "singleton_infeasible", tracer, stats["infeasible"], stats["seen"]
            )
        for pair, stats in block["pairwise"].items():
            emit("subset", k, "pair_infeasible", pair, stats["infeasible"], stats["seen"])
    for reason, count in subset.get("accounting", {}).get("row_skips", {}).items():
        emit("subset", "", "row_skip", reason, count, None)
    if subset:
        emit(
            "subset",
            "",
            "panel_total",
            "examined",
            subset["accounting"]["panels_examined"],
            None,
        )
    return records


def print_report(payload: Mapping[str, Any]) -> None:
    envelope = payload.get("envelope")
    if envelope:
        accounting = envelope["accounting"]
        print("\n=== fold accounting (envelope diagnostic) ===")
        print(f"  rows loaded              : {accounting['n_rows']}")
        for reason, count in accounting["row_skips"].items():
            print(f"    row skip {reason:38s} {count}")
        print(f"  eligible folds           : {accounting['eligible_folds']}")
        print(f"  folds examined           : {accounting['folds_examined']}")
        print(f"  folds skipped            : {accounting['folds_skipped']}")
        for reason, count in accounting["fold_skips"].items():
            print(f"    fold skip {reason:37s} {count}")
        print(f"  reconciled               : {accounting['reconciled']}")
        for block in envelope["by_k"].values():
            folds = block["folds_examined"]
            violating = block["folds_with_envelope_violation"]
            print(f"\n=== envelope, k = {block['sigma_multiplier']:g} ===")
            print(
                f"  folds with >=1 violation : {violating}/{folds}"
                f"  ({100 * (_rate(violating, folds) or 0.0):.1f}%)"
            )
            total = block["total_constraints"]
            print(
                f"  violating constraints    : {block['violating_constraints']}/{total}"
                f"  ({100 * (_rate(block['violating_constraints'], total) or 0.0):.1f}%)"
            )
            for tracer, stats in block["by_tracer"].items():
                if not stats["violations"]:
                    continue
                print(
                    f"    {tracer:9s} {stats['violations']:5d}/{stats['seen']:5d}"
                    f"  ({100 * (_rate(stats['violations'], stats['seen']) or 0.0):.1f}%)"
                )

    subset = payload.get("subset")
    if subset:
        print("\n=== panel accounting (subset diagnostic) ===")
        accounting = subset["accounting"]
        print(f"  rows loaded              : {accounting['n_rows']}")
        print(f"  panels examined          : {accounting['panels_examined']}")
        for reason, count in accounting["row_skips"].items():
            print(f"    row skip {reason:38s} {count}")
        print(f"  reconciled               : {accounting['reconciled']}")
        for block in subset["by_k"].values():
            print(f"\n=== minimal infeasible subsets, k = {block['sigma_multiplier']:g} ===")
            print(
                f"  panels infeasible        : {block['panels_infeasible']}"
                f"/{block['panels_examined']}"
            )
            for size, count in block["minimal_infeasible_subset_size_counts"].items():
                print(f"    MIS size {size}: {count}")
            print("  singleton infeasible (infeasible/seen):")
            for tracer, stats in block["singleton"].items():
                print(f"    {tracer:9s} {stats['infeasible']:5d}/{stats['seen']:5d}")
            print("  pairwise infeasible (infeasible/seen):")
            pairs = sorted(
                block["pairwise"].items(),
                key=lambda item: -(_rate(item[1]["infeasible"], item[1]["seen"]) or 0.0),
            )
            for pair, stats in pairs:
                print(
                    f"    {pair:22s} {stats['infeasible']:5d}/{stats['seen']:5d}"
                    f"  ({100 * (_rate(stats['infeasible'], stats['seen']) or 0.0):.1f}%)"
                )


# --------------------------------------------------------------------------
# Provenance and CLI
# --------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def write_output(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        pd.DataFrame.from_records(flatten_to_records(payload)).to_csv(
            path, index=False
        )
    else:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
        )


def build_manifest(
    args: argparse.Namespace,
    protocol: Mapping[str, Any],
    k_values: Sequence[float],
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    source_files = sorted(
        path for path in usgs.M2_USGS_DATA.glob("Table_*.txt") if path.is_file()
    )
    counts: dict[str, Any] = {}
    if "envelope" in payload:
        counts["envelope"] = payload["envelope"]["accounting"]
    if "subset" in payload:
        counts["subset"] = payload["subset"]["accounting"]
    return {
        "schema_version": "1.0",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostic_id": "m3-identified-ttd-infeasibility-diagnostics",
        "protocol_id": protocol["protocol_id"],
        "protocol_status": protocol["status"],
        "claim_authority": protocol.get("claim_authority"),
        "git_commit": _git_commit(),
        "python_version": platform.python_version(),
        "source_mode": args.sources,
        "max_rows": args.max_rows,
        "diagnostics": args.diagnostic,
        "sigma_multipliers": [float(k) for k in k_values],
        "input_hashes": {
            str(path.relative_to(REPO_ROOT)): _sha256(path) for path in source_files
        },
        "protocol_path": str(args.protocol),
        "protocol_sha256": _sha256(Path(args.protocol)),
        "script_sha256": _sha256(Path(__file__)),
        "output_path": str(args.output),
        "output_sha256": _sha256(Path(args.output)),
        "graph_mode": "disabled_local_only",
        "claim_boundary": (
            "development diagnostic; quantifies local infeasibility structure "
            "only and authorizes no scientific claim"
        ),
        "row_counts": counts,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--sources", default="national")
    parser.add_argument(
        "--k",
        type=float,
        action="append",
        dest="k_values",
        metavar="SIGMA_MULTIPLIER",
        help=(
            "Sigma multiplier for the observation intervals; repeatable. "
            f"Defaults to {' and '.join(str(k) for k in DEFAULT_K_VALUES)}."
        ),
    )
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--diagnostic",
        choices=["both", "envelope", "subset"],
        default="both",
        help="Which diagnostic to compute (default: both).",
    )
    args = parser.parse_args(argv)
    if not args.k_values:
        args.k_values = list(DEFAULT_K_VALUES)
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    protocol = runner.load_protocol(args.protocol)
    ages = runner.build_protocol_age_grid(protocol)
    k_values = [float(k) for k in dict.fromkeys(args.k_values)]

    frame = usgs.load_benchmark_dataset(sources=args.sources)
    if args.max_rows is not None:
        frame = frame.head(max(0, int(args.max_rows)))
    rows = [row for _, row in frame.iterrows()]

    payload: dict[str, Any] = {
        "schema_version": "1.0",
        "diagnostic_id": "m3-identified-ttd-infeasibility-diagnostics",
        "protocol_id": protocol["protocol_id"],
        "protocol_status": protocol["status"],
        "source_mode": args.sources,
        "sigma_multipliers": k_values,
        "n_rows": int(len(rows)),
    }
    if args.diagnostic in ("both", "envelope"):
        print(f"Running envelope diagnostic over {len(rows)} rows...", flush=True)
        payload["envelope"] = run_envelope_diagnostic(rows, protocol, ages, k_values)
    if args.diagnostic in ("both", "subset"):
        print(f"Running subset diagnostic over {len(rows)} rows...", flush=True)
        payload["subset"] = run_subset_diagnostic(rows, protocol, ages, k_values)

    write_output(Path(args.output), payload)
    manifest = build_manifest(args, protocol, k_values, payload)
    manifest_path = Path(args.output).with_name(f"{Path(args.output).stem}_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print_report(payload)
    print(f"\nWrote {args.output}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
