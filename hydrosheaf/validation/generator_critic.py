"""Adversarial, failure-facing audits for synthetic groundwater generators.

The critic is intentionally stricter than an execution test.  It checks the
generator source and withheld truth separately from the inference rows, looking
for circular imports, hidden truth, trivial topology, duplicated evidence,
weak physical directionality, absent stress structure, and non-reproducibility.
Major findings do not disappear merely because HydroSheaf can process the
rows; they remain a generator-design review item.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from enum import IntEnum
import json
import math
from pathlib import Path
from typing import Mapping, Sequence


class CriticSeverity(IntEnum):
    """Ordered severity for adversarial generator findings."""

    INFO = 0
    MINOR = 1
    MAJOR = 2
    BLOCKER = 3


@dataclass(frozen=True)
class CriticFinding:
    """One auditable generator criticism."""

    check_id: str
    severity: CriticSeverity
    title: str
    evidence: str
    remediation: str

    def to_dict(self) -> dict[str, object]:
        return {
            "check_id": self.check_id,
            "severity": self.severity.name,
            "title": self.title,
            "evidence": self.evidence,
            "remediation": self.remediation,
        }


@dataclass(frozen=True)
class GeneratorCriticReport:
    """Aggregated adversarial review for one generator case."""

    family: str
    seed: int
    verdict: str
    findings: Sequence[CriticFinding]
    summary: Mapping[str, int]

    @property
    def blockers(self) -> tuple[CriticFinding, ...]:
        return tuple(
            finding
            for finding in self.findings
            if finding.severity is CriticSeverity.BLOCKER
        )

    @property
    def majors(self) -> tuple[CriticFinding, ...]:
        return tuple(
            finding
            for finding in self.findings
            if finding.severity is CriticSeverity.MAJOR
        )

    @property
    def critic_gate(self) -> bool:
        return not self.blockers and not self.majors

    def to_dict(self) -> dict[str, object]:
        return {
            "family": self.family,
            "seed": self.seed,
            "verdict": self.verdict,
            "critic_gate": self.critic_gate,
            "summary": dict(self.summary),
            "findings": [finding.to_dict() for finding in self.findings],
        }


def _finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return True


def _nonnegative_tracer(value: object) -> bool:
    """Accept finite values and explicit one-sided detection limits."""

    if value is None:
        return True
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        numeric = None
    if numeric is not None:
        return math.isfinite(numeric) and numeric >= 0.0
    text = str(value).strip()
    if text[:1] in {"<", ">"}:
        try:
            limit = float(text[1:].strip())
            return math.isfinite(limit) and limit >= 0.0
        except (TypeError, ValueError):
            return True
    return True


def _stress_counts(rows: Sequence[Mapping[str, object]]) -> tuple[int, int]:
    """Count missing and censored observations without using generator truth."""

    missing = 0
    censored = 0
    for row in rows:
        flags = row.get("observation_flags")
        if isinstance(flags, Mapping):
            missing += len({str(item) for item in flags.get("missing_fields", ())})
            censored += len({str(item) for item in flags.get("censored_fields", ())})
            continue
        for value in row.values():
            if value is None:
                missing += 1
            elif isinstance(value, str) and value.strip()[:1] in {"<", ">"}:
                censored += 1
    return missing, censored


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, default=str, separators=(",", ":"))


def _edge_pair(edge: object) -> tuple[str, str] | None:
    if isinstance(edge, Mapping):
        if "u" in edge and "v" in edge:
            return str(edge["u"]), str(edge["v"])
        edge_id = str(edge.get("edge_id", ""))
        if "->" in edge_id:
            return tuple(edge_id.split("->", 1))  # type: ignore[return-value]
    if isinstance(edge, Sequence) and not isinstance(edge, (str, bytes)) and len(edge) >= 2:
        return str(edge[0]), str(edge[1])
    return None


def _has_cycle(nodes: set[str], edges: Sequence[tuple[str, str]]) -> bool:
    adjacency: dict[str, list[str]] = {node: [] for node in nodes}
    for source, target in edges:
        adjacency.setdefault(source, []).append(target)
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str) -> bool:
        if node in visiting:
            return True
        if node in visited:
            return False
        visiting.add(node)
        if any(visit(child) for child in adjacency.get(node, ())):
            return True
        visiting.remove(node)
        visited.add(node)
        return False

    return any(visit(node) for node in nodes)


def _pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) != len(right) or len(left) < 3:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_scale = math.sqrt(sum((a - left_mean) ** 2 for a in left))
    right_scale = math.sqrt(sum((b - right_mean) ** 2 for b in right))
    if left_scale <= 0.0 or right_scale <= 0.0:
        return None
    return numerator / (left_scale * right_scale)


def _finding(
    findings: list[CriticFinding],
    check_id: str,
    severity: CriticSeverity,
    title: str,
    evidence: str,
    remediation: str,
) -> None:
    findings.append(
        CriticFinding(
            check_id=check_id,
            severity=severity,
            title=title,
            evidence=evidence,
            remediation=remediation,
        )
    )


def audit_source_independence(source_path: Path) -> list[CriticFinding]:
    """Audit actual imports rather than trusting provenance metadata."""

    findings: list[CriticFinding] = []
    path = Path(source_path)
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except (OSError, SyntaxError) as exc:
        _finding(
            findings,
            "SOURCE_PARSE",
            CriticSeverity.BLOCKER,
            "Generator source cannot be parsed",
            f"{path}: {type(exc).__name__}: {exc}",
            "Repair the source before using the generator in a locked benchmark.",
        )
        return findings

    imported_modules: list[str] = []
    imported_names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.append(str(node.module or ""))
            imported_names.extend(alias.name for alias in node.names)
    forbidden_modules = tuple(
        module
        for module in imported_modules
        if module == "hydrosheaf" or module.startswith("hydrosheaf.")
    )
    forbidden_names = tuple(
        name
        for name in imported_names
        if name in {"fit_network_pipeline", "run_programme_synthetic_gate"}
    )
    if forbidden_modules or forbidden_names:
        _finding(
            findings,
            "SOURCE_INDEPENDENCE",
            CriticSeverity.BLOCKER,
            "Generator imports inference code",
            f"modules={forbidden_modules}, names={forbidden_names}",
            "Move the forward model into an independent source package and expose only observations plus sealed truth.",
        )
    return findings


def audit_generator_case(
    family: str,
    case: object,
    *,
    source_path: Path,
    repeat_case: object | None = None,
    alternate_case: object | None = None,
    observation_scenarios: Mapping[str, Sequence[Mapping[str, object]]] | None = None,
    provenance_override: Mapping[str, object] | None = None,
    covariance_consumed: bool = False,
) -> GeneratorCriticReport:
    """Run the adversarial review against one generated case."""

    findings: list[CriticFinding] = audit_source_independence(source_path)
    observations_raw = getattr(case, "observations", ())
    observations = [dict(row) for row in observations_raw]
    node_ids = [str(row.get("site_id", "")) for row in observations]
    node_set = set(node_ids)
    if not observations:
        _finding(
            findings,
            "OBSERVATION_COUNT",
            CriticSeverity.BLOCKER,
            "Generator returned no observations",
            "observation_count=0",
            "Return a non-empty observed case before any inference stage.",
        )
    if len(node_ids) != len(node_set) or any(not node for node in node_ids):
        _finding(
            findings,
            "OBSERVATION_IDS",
            CriticSeverity.BLOCKER,
            "Observation identifiers are missing or duplicated",
            f"rows={len(node_ids)}, unique_ids={len(node_set)}",
            "Assign one stable, unique site identifier to every observed node.",
        )

    truth_tokens = ("true", "truth", "age", "travel", "pathline", "modpath", "particle", "process", "edge")
    leaked_keys = sorted(
        {
            str(key)
            for row in observations
            for key in row
            if any(token in str(key).lower() for token in truth_tokens)
        }
    )
    if leaked_keys:
        _finding(
            findings,
            "OBSERVATION_TRUTH_LEAK",
            CriticSeverity.BLOCKER,
            "Observation rows expose latent truth or route metadata",
            f"suspicious_keys={leaked_keys}",
            "Remove latent ages, process labels, edge IDs, particle IDs, and route fields from inference rows.",
        )

    numeric_bad = [
        str(key)
        for row in observations
        for key, value in row.items()
        if not _finite(value)
    ]
    if numeric_bad:
        _finding(
            findings,
            "FINITE_OBSERVATIONS",
            CriticSeverity.BLOCKER,
            "Observation values contain non-finite numbers",
            f"fields={sorted(set(numeric_bad))}",
            "Emit finite values or explicit censored/missing records with a declared observation model.",
        )

    nonnegative_fields = (
        "tritium_TU",
        "argon39_pmc",
        "c14_pmc",
        "cfc11_pptv",
        "cfc12_pptv",
        "cfc113_pptv",
        "sf6_pptv",
        "he4_ccpg",
        "h3_he3_TU",
    )
    negative_fields = sorted(
        {
            field
            for row in observations
            for field in nonnegative_fields
            if field in row and not _nonnegative_tracer(row[field])
        }
    )
    if negative_fields:
        _finding(
            findings,
            "NONNEGATIVE_TRACERS",
            CriticSeverity.BLOCKER,
            "Tracer concentrations are negative",
            f"fields={negative_fields}",
            "Apply a declared censoring/detection-limit model instead of silently emitting negative concentrations.",
        )

    required_fields = ("site_id", "head_meas", "hydraulic_head", "sample_date", "pH")
    missing_required = sorted(field for field in required_fields if any(field not in row for row in observations))
    if missing_required:
        _finding(
            findings,
            "OBSERVATION_SCHEMA",
            CriticSeverity.MAJOR,
            "Core observation schema is incomplete",
            f"missing_fields={missing_required}",
            "Declare a per-family schema and a controlled missingness experiment before scoring.",
        )

    true_edges = []
    for edge in getattr(case, "true_edges", ()):
        parsed = _edge_pair(edge)
        if parsed is not None:
            true_edges.append(parsed)
    invalid_edges = [edge for edge in true_edges if edge[0] not in node_set or edge[1] not in node_set]
    duplicate_edges = len(true_edges) != len(set(true_edges))
    if invalid_edges or duplicate_edges:
        _finding(
            findings,
            "TRUTH_GRAPH_INTEGRITY",
            CriticSeverity.BLOCKER,
            "Withheld truth graph is invalid",
            f"invalid_edges={invalid_edges[:3]}, duplicate_edges={duplicate_edges}",
            "Validate every truth edge against the observation registry and remove duplicate edges.",
        )
    if _has_cycle(node_set, true_edges):
        _finding(
            findings,
            "TRUTH_GRAPH_CYCLE",
            CriticSeverity.BLOCKER,
            "Directed groundwater truth contains a cycle",
            "The flow truth is not acyclic despite age-directed edges.",
            "Use a physically justified recirculation model or enforce an acyclic advective truth graph.",
        )

    true_ages = {str(key): float(value) for key, value in getattr(case, "true_ages_years", {}).items()}
    if not true_ages or not all(math.isfinite(value) and value >= 0.0 for value in true_ages.values()):
        _finding(
            findings,
            "AGE_TRUTH",
            CriticSeverity.BLOCKER,
            "Age truth is absent or invalid",
            f"age_count={len(true_ages)}",
            "Return finite non-negative ages for every scored node, outside the inference rows.",
        )
    age_order = [true_ages[v] > true_ages[u] for u, v in true_edges if u in true_ages and v in true_ages]
    if age_order and sum(age_order) / len(age_order) < 0.90:
        _finding(
            findings,
            "AGE_DIRECTION",
            CriticSeverity.MAJOR,
            "Most true edges do not increase groundwater age",
            f"age_increasing_fraction={sum(age_order) / len(age_order):.3f}",
            "Make edge direction, travel time, and age evolution consistent or declare mixing/recirculation explicitly.",
        )

    observation_by_id = {str(row["site_id"]): row for row in observations if "site_id" in row}
    head_order = [
        float(observation_by_id[u]["head_meas"]) > float(observation_by_id[v]["head_meas"])
        for u, v in true_edges
        if u in observation_by_id and v in observation_by_id and "head_meas" in observation_by_id[u] and "head_meas" in observation_by_id[v]
    ]
    if head_order and sum(head_order) / len(head_order) < 0.80:
        _finding(
            findings,
            "HEAD_FLOW_DIRECTION",
            CriticSeverity.MAJOR,
            "Head observations do not support the declared flow direction",
            f"head_decreasing_fraction={sum(head_order) / len(head_order):.3f}",
            "Reconcile the head field and directed truth edges, including observation noise and boundary conditions.",
        )

    outgoing: dict[str, set[str]] = {node: set() for node in node_set}
    incoming: dict[str, set[str]] = {node: set() for node in node_set}
    for source, target in true_edges:
        outgoing.setdefault(source, set()).add(target)
        incoming.setdefault(target, set()).add(source)
    branch_nodes = sum(len(targets) > 1 for targets in outgoing.values())
    merge_nodes = sum(len(sources) > 1 for sources in incoming.values())
    if branch_nodes == 0 and merge_nodes == 0:
        _finding(
            findings,
            "TOPOLOGY_COMPLEXITY",
            CriticSeverity.MAJOR,
            "Truth topology is only a set of disjoint chains",
            f"nodes={len(node_set)}, edges={len(true_edges)}, branch_nodes=0, merge_nodes=0",
            "Add branching, merging, heterogeneous connectivity, and at least one topology family with nontrivial alternatives.",
        )

    coordinates = [
        (float(row["x_m"]), float(row["y_m"]))
        for row in observations
        if "x_m" in row and "y_m" in row
    ]
    if len(coordinates) != len(set(coordinates)):
        _finding(
            findings,
            "GEOMETRY_DUPLICATION",
            CriticSeverity.MAJOR,
            "Multiple observed nodes share identical coordinates",
            f"coordinate_count={len(coordinates)}, unique_coordinates={len(set(coordinates))}",
            "Use distinct locations or explicitly model co-located multi-screen observations.",
        )
    layer_values = {row.get("aquifer_layer") for row in observations if "aquifer_layer" in row}
    if len(layer_values) < 2:
        _finding(
            findings,
            "LAYER_HETEROGENEITY",
            CriticSeverity.MAJOR,
            "Generator has no vertical/layer heterogeneity",
            f"unique_aquifer_layers={len(layer_values)}",
            "Add multilayer or stratigraphically heterogeneous cases before making 3-D connectivity claims.",
        )

    dynamic_fields = ("tritium_TU", "argon39_pmc", "c14_pmc")
    weak_dynamic = []
    for field in dynamic_fields:
        values = [float(row[field]) for row in observations if field in row]
        if len(set(values)) < 4 or max(values, default=0.0) - min(values, default=0.0) <= 0.0:
            weak_dynamic.append(field)
    if weak_dynamic:
        _finding(
            findings,
            "TRACER_DYNAMIC_RANGE",
            CriticSeverity.MAJOR,
            "Age-tracer panel has insufficient dynamic range",
            f"weak_fields={weak_dynamic}",
            "Vary age, recharge history, mixing, and censoring so tracer identifiability is genuinely stress-tested.",
        )

    ion_fields = ("Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "Fe", "SiO2")
    varying_ions = sum(
        len({float(row[field]) for row in observations if field in row}) >= 4
        for field in ion_fields
    )
    if varying_ions < 4:
        _finding(
            findings,
            "CHEMISTRY_DYNAMIC_RANGE",
            CriticSeverity.MAJOR,
            "Chemistry panel is nearly invariant",
            f"varying_ion_count={varying_ions}",
            "Use independent reaction paths, endmembers, dilution, and analytical noise that create measurable but physically bounded variation.",
        )

    missing_values = sum(
        value is None or (isinstance(value, float) and not math.isfinite(value))
        for row in observations
        for value in row.values()
    )
    stress_has_missingness = False
    if observation_scenarios:
        for scenario_name, scenario_rows_raw in observation_scenarios.items():
            scenario_rows = [dict(row) for row in scenario_rows_raw]
            scenario_missing, scenario_censored = _stress_counts(scenario_rows)
            if str(scenario_name) != "complete" and (
                scenario_missing > 0 or scenario_censored > 0
            ):
                stress_has_missingness = True
            if not scenario_rows:
                _finding(
                    findings,
                    "OBSERVATION_STRESS_EMPTY",
                    CriticSeverity.BLOCKER,
                    "Declared observation-stress scenario is empty",
                    f"scenario={scenario_name!r}",
                    "Generate a non-empty blind observation variant before scoring it.",
                )
            scenario_leaked = sorted(
                {
                    str(key)
                    for row in scenario_rows
                    for key in row
                    if any(token in str(key).lower() for token in truth_tokens)
                }
            )
            if scenario_leaked:
                _finding(
                    findings,
                    "OBSERVATION_STRESS_TRUTH_LEAK",
                    CriticSeverity.BLOCKER,
                    "Observation-stress rows expose latent truth",
                    f"scenario={scenario_name!r}, suspicious_keys={scenario_leaked}",
                    "Apply missingness and censoring only to blind observation rows.",
                )
    if missing_values == 0 and not stress_has_missingness:
        _finding(
            findings,
            "MISSINGNESS_STRESS",
            CriticSeverity.MAJOR,
            "Generator emits complete uncensored observations only",
            "missing_or_censored_value_count=0",
            "Add preregistered missingness, detection limits, censoring, and structured dropout scenarios.",
        )

    if all("head_meas" in row and "hydraulic_head" in row for row in observations):
        head_measured = [float(row["head_meas"]) for row in observations if "hydraulic_head" in row]
        head_hydraulic = [float(row["hydraulic_head"]) for row in observations if "head_meas" in row]
        correlation = _pearson(head_measured, head_hydraulic)
        if correlation is not None and correlation > 0.995:
            provenance = provenance_override or getattr(case, "provenance", {})
            covariance_declared = isinstance(provenance, Mapping) and bool(
                provenance.get("head_channel_covariance_model")
            )
            covariance_is_consumed = covariance_declared and covariance_consumed
            _finding(
                findings,
                "REDUNDANT_HEAD_EVIDENCE",
                (
                    CriticSeverity.INFO if covariance_is_consumed else
                    CriticSeverity.MINOR if covariance_declared else
                    CriticSeverity.MAJOR
                ),
                "Two head channels are nearly duplicate evidence",
                f"pearson_correlation={correlation:.5f}",
                (
                    "The declared covariance/discrepancy model was consumed and "
                    "the secondary channel was not passed as independent evidence."
                    if covariance_is_consumed
                    else "Ensure the inference stage consumes the declared covariance "
                    "rather than treating the channels as independent evidence."
                    if covariance_declared
                    else "Specify whether the channels are independent instruments, "
                    "repeated measurements, or one derived field; avoid double-counting them."
                ),
            )

    if repeat_case is not None:
        same_case = _canonical(getattr(case, "observations", ())) == _canonical(
            getattr(repeat_case, "observations", ())
        )
        if not same_case:
            _finding(
                findings,
                "REPRODUCIBILITY",
                CriticSeverity.BLOCKER,
                "Same seed does not reproduce the same observations",
                "Repeated generation produced different observation payloads.",
                "Freeze the random-stream design and test exact regeneration before locked runs.",
            )
    else:
        _finding(
            findings,
            "REPRODUCIBILITY_AUDIT",
            CriticSeverity.MINOR,
            "Same-seed regeneration was not supplied to the critic",
            "The case could not be independently regenerated in this audit call.",
            "Always provide a same-seed regeneration record for the generator gate.",
        )

    if alternate_case is not None:
        same_seed_different_case = _canonical(getattr(case, "observations", ())) == _canonical(
            getattr(alternate_case, "observations", ())
        )
        if same_seed_different_case:
            _finding(
                findings,
                "SEED_VARIATION",
                CriticSeverity.BLOCKER,
                "Different seeds produce identical observations",
                "Alternate seed generated the same observation payload.",
                "Repair seed handling and include seed-diversity tests in the locked generator registry.",
            )

    provenance = provenance_override or getattr(case, "provenance", {})
    if isinstance(provenance, Mapping) and provenance.get("base_generator_module"):
        _finding(
            findings,
            "FORWARD_MODEL_REUSE",
            CriticSeverity.MAJOR,
            "Generator family reuses an earlier family core",
            f"base_generator_module={provenance.get('base_generator_module')}",
            "Treat this family as an extension/baseline and add a genuinely independent held-out process/geometry family.",
        )

    if not findings:
        verdict = "PASS"
    else:
        maximum = max(finding.severity for finding in findings)
        verdict = "BLOCKED" if maximum is CriticSeverity.BLOCKER else (
            "REVIEW_REQUIRED" if maximum is CriticSeverity.MAJOR else "PASS_WITH_NOTES"
        )
    summary = {severity.name: sum(finding.severity is severity for finding in findings) for severity in CriticSeverity}
    return GeneratorCriticReport(
        family=str(family),
        seed=int(getattr(case, "seed", -1)),
        verdict=verdict,
        findings=tuple(findings),
        summary=summary,
    )


__all__ = [
    "CriticFinding",
    "CriticSeverity",
    "GeneratorCriticReport",
    "audit_generator_case",
    "audit_source_independence",
]
