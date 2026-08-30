"""Leakage-guarded M3 evaluation of set-valued TTD evidence.

This runner withholds one tracer at a time, constructs the local identified set
from the remaining observations, and evaluates the withheld concentration
against sharp predictive bounds.  Reported USGS ages and age fractions are not
read by the inference path.  The development protocol authorizes an
implementation benchmark only, not a field-validation or true-TTD claim.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, MutableMapping, Sequence

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.nuclear.joint_lpm import (  # noqa: E402
    TracerFitObservation,
    build_lpm_tracer_observations,
)
from hydrosheaf.nuclear.old_groundwater import prepare_c14_observation  # noqa: E402
from hydrosheaf.nuclear.ttd_graph import (  # noqa: E402
    MassTransportMap,
    audit_ttd_graph_compatibility,
)
from hydrosheaf.nuclear.ttd_identified import (  # noqa: E402
    TtdEvidenceReport,
    assess_held_out_tracer,
    infer_ttd_evidence,
    standard_age_functionals,
)
import run_m3_usgs_benchmark as usgs  # noqa: E402

BENCHMARK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = BENCHMARK_DIR / "configs" / "identified_ttd_protocol.yaml"
DEFAULT_OUTPUT = BENCHMARK_DIR / "results" / "m3_identified_ttd_heldout.csv"

# The per-fold inference path is local-only in every mode.  Supplying declared
# edges enables a downstream compatibility audit that reads frozen local
# evidence; it never feeds back into constraint construction or bounds, so the
# per-row ``graph_mode`` column stays ``GRAPH_MODE_LOCAL_ONLY`` unconditionally.
GRAPH_MODE_LOCAL_ONLY = "disabled_local_only"
GRAPH_MODE_COMPATIBILITY_AUDIT = "compatibility_audit_only"
SUPPORTED_TRANSPORT_GENERATORS = ("identity",)

_WITHHOLD_EQUIVALENCE = {
    "3H": {"3H", "3H/3He"},
    "SF6": {"SF6"},
    "CFC11": {"CFC11"},
    "CFC12": {"CFC12"},
    "CFC113": {"CFC113"},
    "14C": {"14C"},
}


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_protocol(path: Path = DEFAULT_PROTOCOL) -> dict[str, Any]:
    protocol = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(protocol, dict):
        raise ValueError("Identified-TTD protocol must be a YAML mapping.")
    required = {"schema_version", "protocol_id", "status", "inference", "observations"}
    missing = sorted(required.difference(protocol))
    if missing:
        raise ValueError(f"Identified-TTD protocol is missing: {missing}.")
    return protocol


def build_protocol_age_grid(protocol: Mapping[str, Any]) -> np.ndarray:
    segments = protocol["inference"]["age_grid_segments"]
    values: list[np.ndarray] = []
    for segment in segments:
        start = float(segment["start"])
        stop = float(segment["stop"])
        step = float(segment["step"])
        if step <= 0.0 or stop <= start:
            raise ValueError("Every age-grid segment needs stop > start and step > 0.")
        count = int(math.floor((stop - start) / step + 1.0e-12))
        values.append(start + step * np.arange(count + 1, dtype=float))
        if values[-1][-1] < stop - 1.0e-9:
            values.append(np.asarray([stop], dtype=float))
    ages = np.unique(np.concatenate(values))
    if ages[0] != 0.0:
        raise ValueError("The identified-TTD age grid must start at zero.")
    return ages


@dataclass(frozen=True)
class DeclaredEdge:
    """One externally declared transport edge between two benchmark sites."""

    edge_id: str
    source_site_id: str
    target_site_id: str
    transport: MassTransportMap


def _build_transport_matrix(entry: Mapping[str, Any], n_age_bins: int) -> np.ndarray:
    """Resolve an explicit matrix or a named generator to a dense operator."""
    has_matrix = entry.get("transport_matrix") is not None
    has_generator = entry.get("transport_generator") is not None
    if has_matrix == has_generator:
        raise ValueError(
            "Each edge needs exactly one of 'transport_matrix' or "
            "'transport_generator'."
        )
    if has_generator:
        generator = str(entry["transport_generator"]).strip().lower()
        if generator not in SUPPORTED_TRANSPORT_GENERATORS:
            raise ValueError(
                f"Unsupported transport_generator {generator!r}; "
                f"supported: {list(SUPPORTED_TRANSPORT_GENERATORS)}."
            )
        return np.eye(n_age_bins, dtype=float)
    matrix = np.asarray(entry["transport_matrix"], dtype=float)
    if matrix.shape != (n_age_bins, n_age_bins):
        raise ValueError(
            "transport_matrix must have shape (n_target_bins, n_source_bins) = "
            f"({n_age_bins}, {n_age_bins}); got {matrix.shape}."
        )
    return matrix


def load_graph_edges(path: Path, n_age_bins: int) -> list[DeclaredEdge]:
    """Load and validate a declared edge file (YAML or JSON).

    Every edge is validated against ``MassTransportMap``'s own constructor
    invariants (finite, non-negative, columns summing to one) plus the shape
    implied by the protocol age grid.  No edge is synthesised or defaulted:
    an unusable declaration raises rather than silently degrading.
    """
    document = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(document, Mapping):
        raise ValueError("Graph edge file must be a mapping with an 'edges' key.")
    raw_edges = document.get("edges")
    if not isinstance(raw_edges, Sequence) or isinstance(raw_edges, (str, bytes)):
        raise ValueError("Graph edge file 'edges' must be a list.")
    if not raw_edges:
        raise ValueError("Graph edge file declares no edges.")

    required = (
        "edge_id",
        "source_site_id",
        "target_site_id",
        "provenance_tier",
        "source",
    )
    edges: list[DeclaredEdge] = []
    seen: set[str] = set()
    for index, entry in enumerate(raw_edges):
        if not isinstance(entry, Mapping):
            raise ValueError(f"Edge at position {index} must be a mapping.")
        missing = [key for key in required if not str(entry.get(key, "")).strip()]
        if missing:
            raise ValueError(f"Edge at position {index} is missing: {missing}.")
        edge_id = str(entry["edge_id"]).strip()
        if edge_id in seen:
            raise ValueError(f"Duplicate edge_id {edge_id!r} in graph edge file.")
        seen.add(edge_id)
        matrix = _build_transport_matrix(entry, n_age_bins)
        metadata = entry.get("metadata") or {}
        if not isinstance(metadata, Mapping):
            raise ValueError(f"Edge {edge_id!r} metadata must be a mapping.")
        transport = MassTransportMap(
            edge_id=edge_id,
            matrix=matrix,
            provenance_tier=str(entry["provenance_tier"]).strip(),
            source=str(entry["source"]).strip(),
            metadata=dict(metadata),
        )
        edges.append(
            DeclaredEdge(
                edge_id=edge_id,
                source_site_id=str(entry["source_site_id"]).strip(),
                target_site_id=str(entry["target_site_id"]).strip(),
                transport=transport,
            )
        )
    return edges


def _edge_abstention(
    edge: DeclaredEdge, held_out_tracer: str, reason: str
) -> dict[str, Any]:
    return {
        "edge_id": edge.edge_id,
        "source_site_id": edge.source_site_id,
        "target_site_id": edge.target_site_id,
        "held_out_tracer": held_out_tracer,
        "status": "ABSTAIN",
        "minimum_l1_gap": np.nan,
        "maximum_bin_gap": np.nan,
        "provenance_tier": edge.transport.provenance_tier,
        "transport_source": edge.transport.source,
        "tightening_authorized": False,
        "message": reason,
        "reason": reason,
        "graph_mode": GRAPH_MODE_COMPATIBILITY_AUDIT,
    }


def run_graph_compatibility_audits(
    edges: Sequence[DeclaredEdge],
    reports: Mapping[tuple[str, str], TtdEvidenceReport],
    fold_reasons: Mapping[tuple[str, str], str],
    held_out_tracers: Sequence[str],
) -> list[dict[str, Any]]:
    """Audit declared edges against frozen local evidence, never tightening it.

    Abstentions are emitted, never dropped: a missing or infeasible endpoint
    report yields an ABSTAIN row carrying the reason.
    """
    records: list[dict[str, Any]] = []
    for edge in edges:
        for tracer in held_out_tracers:
            source_key = (edge.source_site_id, tracer)
            target_key = (edge.target_site_id, tracer)
            source_report = reports.get(source_key)
            target_report = reports.get(target_key)
            if source_report is None or target_report is None:
                missing = []
                if source_report is None:
                    detail = fold_reasons.get(source_key, "site_not_in_run")
                    missing.append(f"source={edge.source_site_id}({detail})")
                if target_report is None:
                    detail = fold_reasons.get(target_key, "site_not_in_run")
                    missing.append(f"target={edge.target_site_id}({detail})")
                records.append(
                    _edge_abstention(
                        edge,
                        tracer,
                        "endpoint_report_unavailable: " + "; ".join(missing),
                    )
                )
                continue
            try:
                audit = audit_ttd_graph_compatibility(
                    source_report, target_report, edge.transport
                )
            except Exception as exc:  # preserved as an auditable abstention
                records.append(
                    _edge_abstention(
                        edge,
                        tracer,
                        f"audit_error:{type(exc).__name__}:{exc}",
                    )
                )
                continue
            records.append(
                {
                    "edge_id": audit.edge_id,
                    "source_site_id": edge.source_site_id,
                    "target_site_id": edge.target_site_id,
                    "held_out_tracer": tracer,
                    "status": audit.status,
                    "minimum_l1_gap": float(audit.minimum_l1_gap),
                    "maximum_bin_gap": float(audit.maximum_bin_gap),
                    "provenance_tier": audit.provenance_tier,
                    "transport_source": audit.transport_source,
                    "tightening_authorized": False,
                    "message": audit.message,
                    "reason": (audit.message if audit.status == "ABSTAIN" else ""),
                    "graph_mode": GRAPH_MODE_COMPATIBILITY_AUDIT,
                }
            )
    return records


def summarize_graph_audits(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    frame = pd.DataFrame.from_records(list(records))
    if frame.empty:
        return {"n_edge_folds": 0, "status_counts": {}, "tightening_authorized": False}
    gaps = pd.to_numeric(frame.get("minimum_l1_gap"), errors="coerce").dropna()
    return {
        "n_edge_folds": int(len(frame)),
        "n_edges": int(frame["edge_id"].nunique()),
        "status_counts": frame["status"].value_counts().to_dict(),
        "median_minimum_l1_gap": float(gaps.median()) if len(gaps) else None,
        "maximum_minimum_l1_gap": float(gaps.max()) if len(gaps) else None,
        "tightening_authorized": False,
    }


def _base_observations(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "sample_year": row.get("sample_year"),
        "tritium_TU": row.get("tritium_TU"),
        "tritium_sigma_TU": row.get("tritium_sigma_TU"),
        "he3_trit_TU": row.get("he3_trit_TU"),
        "he3_trit_sigma_TU": row.get("he3_trit_sigma_TU"),
        "sf6_pptv": row.get("sf6_pptv"),
        "sf6_sigma_pptv": row.get("sf6_sigma_pptv"),
        "cfc11_pptv": row.get("cfc11_pptv"),
        "cfc11_sigma_pptv": row.get("cfc11_sigma_pptv"),
        "cfc12_pptv": row.get("cfc12_pptv"),
        "cfc12_sigma_pptv": row.get("cfc12_sigma_pptv"),
        "cfc113_pptv": row.get("cfc113_pptv"),
        "cfc113_sigma_pptv": row.get("cfc113_sigma_pptv"),
        "c14_pmc": row.get("c14_pmc"),
        "c14_sigma_pmc": row.get("c14_sigma_pmc"),
        "corrected_c14_pmc": row.get("corrected_c14_pmc"),
        "corrected_a0_pmc": row.get("corrected_a0_pmc"),
        "he4_ccpg": row.get("he4_ccpg"),
        "he4_sigma_ccpg": row.get("he4_sigma_ccpg"),
    }


def prepare_row_observations(
    row: Mapping[str, Any], protocol: Mapping[str, Any]
) -> tuple[list[TracerFitObservation], dict[str, Any]]:
    """Apply only declared corrections and return native observation objects."""
    settings = protocol["observations"]
    factors = {
        "gas_correction_mode": settings["gas_correction_mode"],
        "he4_mode": "calibrated" if settings.get("use_helium4", False) else "disabled",
    }
    observations, factors = usgs._apply_design_factors(
        _base_observations(row), row, factors
    )
    observations, initial_c14_pmc, c14_diagnostics = prepare_c14_observation(
        observations,
        mode=str(settings["c14_correction_mode"]),
        candidate_corrected_pmcs=row.get("c14_candidate_corrected_pmc_json"),
        candidate_initial_pmcs=row.get("c14_candidate_a0_pmc_json"),
        candidate_models=row.get("c14_candidate_models_json"),
    )
    fit_observations = build_lpm_tracer_observations(
        observations, use_helium4=bool(settings.get("use_helium4", False))
    )
    provenance = {
        "gas_correction_mode": factors.get("factor_gas_correction_mode"),
        "gas_correction_input_source": factors.get("gas_correction_input_source"),
        "c14_correction_mode": settings["c14_correction_mode"],
        "c14_effective_source": c14_diagnostics.get("c14_effective_source"),
        "initial_c14_pmc": float(initial_c14_pmc),
        "use_helium4": bool(settings.get("use_helium4", False)),
    }
    return fit_observations, provenance


def split_held_out_observation(
    observations: Sequence[TracerFitObservation], held_out_tracer: str
) -> tuple[list[TracerFitObservation], TracerFitObservation | None]:
    """Remove the target and any direct algebraic equivalent before fitting."""
    target = str(held_out_tracer)
    if target not in _WITHHOLD_EQUIVALENCE:
        raise ValueError(f"Unsupported held-out tracer: {held_out_tracer!r}.")
    excluded_names = _WITHHOLD_EQUIVALENCE[target]
    conditioning: list[TracerFitObservation] = []
    held_out: TracerFitObservation | None = None
    for observation in observations:
        if observation.tracer in excluded_names:
            if observation.tracer == target:
                held_out = observation
            continue
        conditioning.append(observation)
    return conditioning, held_out


def evaluate_held_out_row(
    row: Mapping[str, Any],
    held_out_tracer: str,
    protocol: Mapping[str, Any],
    *,
    report_sink: MutableMapping[tuple[str, str], TtdEvidenceReport] | None = None,
) -> dict[str, Any]:
    """Evaluate one predeclared local-only held-out tracer fold.

    ``report_sink`` is a write-only side channel used by the opt-in graph
    compatibility audit.  It never influences inference or the returned record.
    """
    site_id = str(row.get("site_id", ""))
    sample_year = _finite_float(row.get("sample_year"))
    base = {
        "site_id": site_id,
        "StudyUnit": str(row.get("StudyUnit", "")),
        "AqGroup": str(row.get("AqGroup", "")),
        "held_out_tracer": held_out_tracer,
        "graph_mode": GRAPH_MODE_LOCAL_ONLY,
        "reference_fields_used": False,
    }
    if sample_year is None:
        return {**base, "status": "NOT_ELIGIBLE", "reason": "missing_sample_year"}

    observations, correction_provenance = prepare_row_observations(row, protocol)
    conditioning, held_out = split_held_out_observation(observations, held_out_tracer)
    if held_out is None:
        return {
            **base,
            "status": "NOT_ELIGIBLE",
            "reason": "held_out_observation_unavailable",
        }
    minimum = int(protocol["observations"].get("minimum_conditioning_tracers", 1))
    if len(conditioning) < minimum:
        return {
            **base,
            "status": "ABSTAIN",
            "reason": "insufficient_conditioning_tracers",
            "n_conditioning_tracers": len(conditioning),
        }

    ages = build_protocol_age_grid(protocol)
    thresholds = protocol["inference"].get("maximum_fraction_widths", {})
    functionals = standard_age_functionals(
        ages,
        cutoffs_years=protocol["inference"]["age_fraction_cutoffs_years"],
        maximum_fraction_widths=thresholds,
    )
    histories, history_provenance = usgs._get_site_histories(row)
    report = infer_ttd_evidence(
        conditioning,
        sample_year,
        ages,
        functionals=functionals,
        histories=histories,
        use_helium4=bool(protocol["observations"].get("use_helium4", False)),
        initial_c14_pmc=correction_provenance["initial_c14_pmc"],
        helium4_background_ccpg=(_finite_float(row.get("he4_background_ccpg")) or 0.0),
        helium4_accumulation_rate_ccpg_per_year=(
            _finite_float(row.get("he4_accumulation_rate_ccpg_per_year")) or 1.0e-11
        ),
        sigma_multiplier=float(protocol["inference"]["sigma_multiplier"]),
        feasibility_tolerance=float(protocol["inference"]["feasibility_tolerance"]),
        scenario_name=str(protocol["protocol_id"]),
        assumptions=(
            "fine-grid non-negative unit-mass TTD",
            "declared observation intervals",
            "local inference with graph disabled",
        ),
        provenance={
            **correction_provenance,
            **history_provenance,
            "protocol_id": protocol["protocol_id"],
        },
    )
    if report_sink is not None:
        report_sink[(site_id, held_out_tracer)] = report
    result: dict[str, Any] = {
        **base,
        "status": report.status,
        "reason": report.message,
        "n_conditioning_tracers": len(report.constraints),
        "conditioning_tracers": "|".join(
            constraint.tracer for constraint in report.constraints
        ),
        "response_rank": report.response_rank,
        "nullity": report.nullity,
        "excluded_observations_json": json.dumps(
            [dict(item) for item in report.excluded_observations], sort_keys=True
        ),
        "held_out_observed": float(held_out.value),
        "held_out_sigma": float(held_out.sigma),
        "held_out_units": held_out.units,
    }
    for name, bound in report.bounds.items():
        result[f"{name}_lower"] = bound.lower
        result[f"{name}_upper"] = bound.upper
        result[f"{name}_width"] = bound.width
        result[f"{name}_status"] = bound.status
    if report.feasible_witness is None:
        return result

    check = assess_held_out_tracer(
        report,
        held_out,
        sample_year,
        histories=histories,
        initial_c14_pmc=correction_provenance["initial_c14_pmc"],
        helium4_background_ccpg=(_finite_float(row.get("he4_background_ccpg")) or 0.0),
        helium4_accumulation_rate_ccpg_per_year=(
            _finite_float(row.get("he4_accumulation_rate_ccpg_per_year")) or 1.0e-11
        ),
    )
    if check["status"] == "ABSTAIN":
        result.update(
            {
                "held_out_status": "ABSTAIN",
                "held_out_compatible": np.nan,
                "held_out_reason": check.get(
                    "reason", "unsupported_held_out_likelihood"
                ),
                "conditioned_on_held_out_observation": False,
            }
        )
        return result
    prediction_width = check["prediction_upper"] - check["prediction_lower"]
    result.update(
        {
            "held_out_status": check["status"],
            "held_out_compatible": check["status"] == "COMPATIBLE",
            "prediction_lower": check["prediction_lower"],
            "prediction_upper": check["prediction_upper"],
            "prediction_width": prediction_width,
            "prediction_width_sigma": prediction_width / float(held_out.sigma),
            "observed_interval_lower": check["observed_interval_lower"],
            "observed_interval_upper": check["observed_interval_upper"],
            "conditioned_on_held_out_observation": check[
                "conditioned_on_held_out_observation"
            ],
        }
    )
    return result


def evaluate_held_out_row_safely(
    row: Mapping[str, Any],
    held_out_tracer: str,
    protocol: Mapping[str, Any],
    *,
    report_sink: MutableMapping[tuple[str, str], TtdEvidenceReport] | None = None,
) -> dict[str, Any]:
    """Preserve unexpected row-level failures as auditable abstentions."""
    try:
        return evaluate_held_out_row(
            row, held_out_tracer, protocol, report_sink=report_sink
        )
    except Exception as exc:
        return {
            "site_id": str(row.get("site_id", "")),
            "StudyUnit": str(row.get("StudyUnit", "")),
            "AqGroup": str(row.get("AqGroup", "")),
            "held_out_tracer": held_out_tracer,
            "status": "ABSTAIN",
            "reason": f"evaluation_error:{type(exc).__name__}:{exc}",
            "graph_mode": GRAPH_MODE_LOCAL_ONLY,
            "reference_fields_used": False,
            "conditioned_on_held_out_observation": False,
        }


def summarize_results(results: pd.DataFrame) -> dict[str, Any]:
    if results.empty or "status" not in results:
        return {
            "n_rows": 0,
            "n_eligible": 0,
            "n_predictions": 0,
            "abstention_rate": None,
            "held_out_coverage": None,
            "median_prediction_width": None,
            "p90_prediction_width": None,
            "status_counts": {},
            "reason_counts": {},
        }
    eligible = results[results["status"] != "NOT_ELIGIBLE"]
    predicted = (
        eligible[eligible["held_out_compatible"].notna()]
        if ("held_out_compatible" in eligible)
        else eligible.iloc[0:0]
    )
    widths = pd.to_numeric(predicted.get("prediction_width"), errors="coerce")
    reason_counts = eligible["reason"].fillna("").value_counts().to_dict()
    return {
        "n_rows": int(len(results)),
        "n_eligible": int(len(eligible)),
        "n_predictions": int(len(predicted)),
        "abstention_rate": (
            float((eligible["status"] == "ABSTAIN").mean()) if len(eligible) else None
        ),
        "held_out_coverage": (
            float(predicted["held_out_compatible"].astype(bool).mean())
            if len(predicted)
            else None
        ),
        "median_prediction_width": (
            float(widths.median()) if len(widths.dropna()) else None
        ),
        "p90_prediction_width": (
            float(widths.quantile(0.9)) if len(widths.dropna()) else None
        ),
        "status_counts": results["status"].value_counts().to_dict(),
        "reason_counts": reason_counts,
    }


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


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sources", default="national")
    parser.add_argument(
        "--withhold-tracer", choices=[*_WITHHOLD_EQUIVALENCE, "all"], default="all"
    )
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--development-overwrite", action="store_true")
    parser.add_argument(
        "--graph-edges",
        type=Path,
        default=None,
        help=(
            "Optional declared edge file (YAML or JSON) enabling the opt-in "
            "graph compatibility audit.  Local inference is unchanged: the "
            "audit only reports obstructions and never tightens a local set."
        ),
    )
    parser.add_argument(
        "--graph-output",
        type=Path,
        default=None,
        help=(
            "Edge-level audit CSV path.  Defaults to '<output stem>"
            "_graph_edges.csv'.  Used only with --graph-edges."
        ),
    )
    args = parser.parse_args(argv)

    if args.graph_output is not None and args.graph_edges is None:
        parser.error("--graph-output requires --graph-edges.")

    protocol = load_protocol(args.protocol)
    if args.output.exists() and not args.development_overwrite:
        raise FileExistsError(
            f"Refusing to overwrite completed result {args.output}; use a new path "
            "or --development-overwrite for a declared development run."
        )
    if args.development_overwrite and protocol["status"] != "development":
        raise ValueError("Overwrite is allowed only by a development protocol.")

    edges: list[DeclaredEdge] = []
    if args.graph_edges is not None:
        edges = load_graph_edges(
            args.graph_edges, build_protocol_age_grid(protocol).size
        )

    data = usgs.load_benchmark_dataset(sources=args.sources)
    if args.max_rows is not None:
        data = data.head(max(0, int(args.max_rows)))
    tracers = (
        list(protocol["observations"]["held_out_tracers"])
        if args.withhold_tracer == "all"
        else [args.withhold_tracer]
    )
    report_sink: dict[tuple[str, str], TtdEvidenceReport] | None = {} if edges else None
    records: list[dict[str, Any]] = []
    total_folds = len(data) * len(tracers)
    for _, row in data.iterrows():
        row_mapping = dict(row)
        for tracer in tracers:
            records.append(
                evaluate_held_out_row_safely(
                    row_mapping, tracer, protocol, report_sink=report_sink
                )
            )
            completed = len(records)
            if completed % 250 == 0 or completed == total_folds:
                print(
                    f"Completed {completed}/{total_folds} held-out folds",
                    flush=True,
                )
    results = pd.DataFrame.from_records(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)
    summary = summarize_results(results)

    graph_output_path: Path | None = None
    graph_summary: dict[str, Any] | None = None
    if edges:
        fold_reasons = {
            (str(record.get("site_id", "")), str(record.get("held_out_tracer", ""))): (
                f"{record.get('status', '')}:{record.get('reason', '')}"
            )
            for record in records
        }
        edge_records = run_graph_compatibility_audits(
            edges, report_sink or {}, fold_reasons, tracers
        )
        graph_summary = summarize_graph_audits(edge_records)
        graph_output_path = args.graph_output or args.output.with_name(
            f"{args.output.stem}_graph_edges.csv"
        )
        graph_output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame.from_records(edge_records).to_csv(graph_output_path, index=False)

    source_files = sorted(
        path for path in usgs.M2_USGS_DATA.glob("Table_*.txt") if path.is_file()
    )
    manifest = {
        "schema_version": "1.0",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_id": protocol["protocol_id"],
        "protocol_status": protocol["status"],
        "claim_authority": protocol.get("claim_authority"),
        "git_commit": _git_commit(),
        "python_version": platform.python_version(),
        "source_mode": args.sources,
        "withheld_tracers": tracers,
        "input_hashes": {
            str(path.relative_to(REPO_ROOT)): _sha256(path) for path in source_files
        },
        "protocol_path": str(args.protocol),
        "protocol_sha256": _sha256(args.protocol),
        "script_sha256": _sha256(Path(__file__)),
        "output_path": str(args.output),
        "output_sha256": _sha256(args.output),
        "leakage_control": (
            "held-out tracer and direct algebraic equivalent removed before "
            "constraint construction; reported ages and age fractions unused"
        ),
        "graph_mode": (
            GRAPH_MODE_COMPATIBILITY_AUDIT if edges else GRAPH_MODE_LOCAL_ONLY
        ),
        "summary": summary,
    }
    if edges:
        manifest["graph_audit"] = {
            "graph_edges_path": str(args.graph_edges),
            "graph_edges_sha256": _sha256(args.graph_edges),
            "n_declared_edges": len(edges),
            "graph_output_path": str(graph_output_path),
            "graph_output_sha256": _sha256(graph_output_path),
            "conditioning_performed": False,
            "tightening_authorized": False,
            "note": (
                "Version 1 compatibility/obstruction audit only. Local "
                "identified sets are frozen; no bound was narrowed by any edge."
            ),
            "summary": graph_summary,
        }
    manifest_path = args.output.with_name(f"{args.output.stem}_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.output}")
    print(f"Wrote {manifest_path}")
    if graph_output_path is not None:
        print(json.dumps(graph_summary, indent=2))
        print(f"Wrote {graph_output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
