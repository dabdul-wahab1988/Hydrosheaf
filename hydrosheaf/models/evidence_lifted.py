"""Evidence-lifted identifiability diagnostics for reaction dictionaries.

These helpers separate structural identifiability from conditional evidence.
If two reactions have the same measured-ion vector, mass balance alone cannot
distinguish them. Evidence-lifted resolution quantifies whether independent
or contextual evidence gives one member of that equivalence class stronger
support, without converting that preference into a false uniqueness claim.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class EvidenceLiftedResolution:
    """Resolution summary for one stoichiometric equivalence class."""

    class_id: str
    members: tuple[str, ...]
    probabilities: Mapping[str, float]
    evidence_scores: Mapping[str, float]
    entropy: float
    normalized_entropy: float
    evidence_lifted_resolution_index: float
    top_member: str
    top_probability: float
    top_margin: float
    status: str
    independent_evidence_available: bool = False
    evidence_source_count: int = 0
    missing_evidence_members: tuple[str, ...] = ()

    def as_row(self) -> dict[str, object]:
        """Return a flat, CSV-friendly representation."""
        return {
            "class_id": self.class_id,
            "members": ";".join(self.members),
            "n_members": len(self.members),
            "top_member": self.top_member,
            "top_probability": self.top_probability,
            "top_margin": self.top_margin,
            "entropy": self.entropy,
            "normalized_entropy": self.normalized_entropy,
            "evidence_lifted_resolution_index": (
                self.evidence_lifted_resolution_index
            ),
            "resolution_status": self.status,
            "independent_evidence_available": self.independent_evidence_available,
            "evidence_source_count": self.evidence_source_count,
            "missing_evidence_members": ";".join(self.missing_evidence_members),
            "member_probabilities": ";".join(
                f"{member}:{self.probabilities[member]:.6g}"
                for member in self.members
            ),
            "member_evidence_scores": ";".join(
                f"{member}:{self.evidence_scores[member]:.6g}"
                for member in self.members
            ),
        }


def _connected_components(edges: Mapping[str, set[str]]) -> list[list[str]]:
    remaining = set(edges)
    components: list[list[str]] = []
    while remaining:
        seed = remaining.pop()
        stack = [seed]
        component = {seed}
        while stack:
            node = stack.pop()
            for neighbour in edges[node]:
                if neighbour not in component:
                    component.add(neighbour)
                    remaining.discard(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))
    return sorted(components, key=lambda values: (values[0], len(values)))


def stoichiometric_equivalence_classes(
    matrix: Sequence[Sequence[float]],
    labels: Sequence[str],
    *,
    cosine_threshold: float = 0.999999,
    signed: bool = True,
) -> tuple[dict[str, str], list[dict[str, object]]]:
    """Group reactions that are indistinguishable in the supplied ion space.

    Parameters
    ----------
    matrix:
        Reaction-by-ion matrix in the measured ion space.
    labels:
        Reaction labels corresponding to matrix rows.
    cosine_threshold:
        Absolute cosine-similarity threshold used for equivalence.
    signed:
        If true, opposite vectors are grouped as the same reversible process.
    """
    arr = np.asarray(matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if arr.shape[0] != len(labels):
        raise ValueError("labels length must match the number of matrix rows")
    if not np.isfinite(arr).all():
        raise ValueError("matrix must contain only finite values")
    clean_labels = [str(label).strip() for label in labels]
    if any(not label for label in clean_labels):
        raise ValueError("reaction labels must not be empty")
    if len(set(clean_labels)) != len(clean_labels):
        raise ValueError("reaction labels must be unique")
    cosine_threshold = float(cosine_threshold)
    if not math.isfinite(cosine_threshold) or not 0.0 <= cosine_threshold <= 1.0:
        raise ValueError("cosine_threshold must be finite and in [0, 1]")

    graph = {label: set() for label in clean_labels}
    for i, left in enumerate(clean_labels):
        for j in range(i + 1, len(labels)):
            right = clean_labels[j]
            norm = float(np.linalg.norm(arr[i]) * np.linalg.norm(arr[j]))
            cosine = float(np.dot(arr[i], arr[j]) / norm) if norm else 0.0
            similarity = abs(cosine) if signed else cosine
            if similarity >= cosine_threshold:
                graph[left].add(right)
                graph[right].add(str(left))

    class_map: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    for index, members in enumerate(_connected_components(graph), start=1):
        class_id = f"EC{index:02d}"
        for member in members:
            class_map[member] = class_id
        rows.append(
            {
                "class_id": class_id,
                "members": ";".join(members),
                "n_members": len(members),
                "ambiguous": len(members) > 1,
            }
        )
    return class_map, rows


def evidence_lifted_resolution(
    members: Sequence[str],
    evidence_scores: Mapping[str, float | None],
    *,
    class_id: str = "",
    score_floor: float = 1e-9,
    weight_transform: str = "odds",
    evidence_sources: Mapping[str, Iterable[str]] | None = None,
    conditional_probability: float = 0.60,
    conditional_margin: float = 0.10,
    resolved_probability: float = 0.80,
    resolved_margin: float = 0.30,
) -> EvidenceLiftedResolution:
    """Calculate entropy-normalised resolution inside an equivalence class.

    Equal evidence scores yield a resolution index of zero. A single dominant
    member yields values approaching one. By default, bounded evidence scores
    are converted to support odds, score / (1 - score), before normalisation.
    The resulting probabilities are relative evidence weights, not calibrated
    posterior probabilities.
    """
    clean_members = tuple(str(member).strip() for member in members)
    if not clean_members:
        raise ValueError("members must contain at least one reaction")
    if any(not member for member in clean_members):
        raise ValueError("members must not contain empty reaction labels")
    if len(set(clean_members)) != len(clean_members):
        raise ValueError("members must be unique")

    score_floor = float(score_floor)
    if not math.isfinite(score_floor) or not 0.0 < score_floor < 0.5:
        raise ValueError("score_floor must be finite and in (0, 0.5)")
    for name, threshold in (
        ("conditional_probability", conditional_probability),
        ("conditional_margin", conditional_margin),
        ("resolved_probability", resolved_probability),
        ("resolved_margin", resolved_margin),
    ):
        value = float(threshold)
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be finite and in [0, 1]")

    scores: dict[str, float] = {}
    missing_score_members: list[str] = []
    for member in clean_members:
        raw_value = evidence_scores.get(member)
        try:
            value = float(raw_value) if raw_value is not None else 0.5
        except (TypeError, ValueError):
            missing_score_members.append(member)
            value = 0.5
        if not math.isfinite(value):
            missing_score_members.append(member)
            value = 0.5
        elif raw_value is None:
            missing_score_members.append(member)
        scores[member] = float(np.clip(value, score_floor, 1.0 - score_floor))

    if weight_transform not in {"odds", "linear"}:
        raise ValueError("weight_transform must be 'odds' or 'linear'")
    if weight_transform == "odds":
        weights = {
            member: scores[member] / max(score_floor, 1.0 - scores[member])
            for member in clean_members
        }
    else:
        weights = scores

    total = sum(weights.values())
    if total <= 0.0:
        probabilities = {member: 1.0 / len(clean_members) for member in clean_members}
    else:
        probabilities = {
            member: weights[member] / total for member in clean_members
        }

    entropy = -sum(
        probability * math.log(max(score_floor, probability))
        for probability in probabilities.values()
    )
    max_entropy = math.log(len(clean_members)) if len(clean_members) > 1 else 0.0
    normalized_entropy = entropy / max_entropy if max_entropy else 0.0
    resolution_index = (
        1.0 - normalized_entropy if len(clean_members) > 1 else 1.0
    )
    resolution_index = float(np.clip(resolution_index, 0.0, 1.0))

    ranked = sorted(
        probabilities.items(),
        key=lambda item: (-item[1], item[0]),
    )
    top_member, top_probability = ranked[0]
    second_probability = ranked[1][1] if len(ranked) > 1 else 0.0
    top_margin = top_probability - second_probability

    source_gate_active = evidence_sources is not None
    source_counts: dict[str, int] = {}
    if source_gate_active:
        assert evidence_sources is not None
        for member in clean_members:
            raw_sources = evidence_sources.get(member, ())
            if isinstance(raw_sources, str):
                source_values = (raw_sources,)
            else:
                try:
                    source_values = tuple(raw_sources)
                except TypeError:
                    source_values = ()
            source_counts[member] = len(
                {
                    str(source).strip()
                    for source in source_values
                    if str(source).strip()
                }
            )
    else:
        source_counts = {member: 0 for member in clean_members}
    top_source_count = source_counts.get(top_member, 0)
    independent_evidence_available = bool(top_source_count) if source_gate_active else False
    evidence_source_count = int(sum(source_counts.values()))

    if len(clean_members) == 1:
        status = "structurally_unique"
    elif source_gate_active and not independent_evidence_available:
        # Unequal scores without a declared independent source are a ranking,
        # not evidence-lifted mechanism resolution.  In particular, missing
        # or unmatched contextual metadata cannot manufacture a preference.
        status = "unresolved_equivalence_class"
    elif (
        top_probability >= resolved_probability
        and top_margin >= resolved_margin
    ):
        status = "evidence_lifted_resolved"
    elif (
        top_probability >= conditional_probability
        and top_margin >= conditional_margin
    ):
        status = "conditionally_preferred"
    else:
        status = "unresolved_equivalence_class"

    return EvidenceLiftedResolution(
        class_id=class_id,
        members=clean_members,
        probabilities=probabilities,
        evidence_scores=scores,
        entropy=float(entropy),
        normalized_entropy=float(normalized_entropy),
        evidence_lifted_resolution_index=resolution_index,
        top_member=top_member,
        top_probability=float(top_probability),
        top_margin=float(top_margin),
        status=status,
        independent_evidence_available=independent_evidence_available,
        evidence_source_count=evidence_source_count,
        missing_evidence_members=tuple(missing_score_members),
    )


def evidence_score_map(
    labels: Iterable[str],
    scores: Iterable[float],
) -> dict[str, float]:
    """Create a reaction-score mapping with float coercion."""
    return {str(label): float(score) for label, score in zip(labels, scores)}


def reaction_identifiability_diagnostics(
    matrix: Sequence[Sequence[float]],
    labels: Sequence[str],
    *,
    cosine_threshold: float = 0.999999,
    signed: bool = True,
    evidence_scores: Mapping[str, float | None] | None = None,
    evidence_sources: Mapping[str, Iterable[str]] | None = None,
) -> dict[str, Any]:
    """Report structural rank/coherence and optional evidence-lifted status.

    The matrix is interpreted as reaction-by-ion in the supplied measured-ion
    space.  The report describes what that space can identify; it never
    promotes a preferred member of an equivalence class to a unique chemical
    mechanism without explicit evidence sources.
    """

    arr = np.asarray(matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if not np.isfinite(arr).all():
        raise ValueError("matrix must contain only finite values")
    if arr.shape[0] != len(labels):
        raise ValueError("labels length must match the number of matrix rows")
    clean_labels = [str(label).strip() for label in labels]
    if any(not label for label in clean_labels):
        raise ValueError("reaction labels must not be empty")
    if len(set(clean_labels)) != len(clean_labels):
        raise ValueError("reaction labels must be unique")
    cosine_threshold = float(cosine_threshold)
    if not math.isfinite(cosine_threshold) or not 0.0 <= cosine_threshold <= 1.0:
        raise ValueError("cosine_threshold must be finite and in [0, 1]")

    n_reactions, n_ions = (int(value) for value in arr.shape)
    if n_reactions and n_ions:
        singular_values = np.linalg.svd(arr, compute_uv=False)
        rank = int(np.linalg.matrix_rank(arr))
        condition_number = float(np.linalg.cond(arr))
    else:
        singular_values = np.asarray([], dtype=float)
        rank = 0
        condition_number = float("nan")

    nonzero_rows = [row for row in arr if float(np.linalg.norm(row)) > 0.0]
    coherences: list[float] = []
    for index, left in enumerate(nonzero_rows):
        left_norm = float(np.linalg.norm(left))
        for right in nonzero_rows[index + 1 :]:
            right_norm = float(np.linalg.norm(right))
            coherences.append(
                abs(float(np.dot(left, right) / (left_norm * right_norm)))
            )
    maximum_coherence = max(coherences, default=0.0)

    class_map, equivalence_classes = stoichiometric_equivalence_classes(
        arr,
        clean_labels,
        cosine_threshold=cosine_threshold,
        signed=signed,
    )
    ambiguous_classes = [
        row for row in equivalence_classes if int(row["n_members"]) > 1
    ]

    resolution_rows: list[dict[str, Any]] = []
    if evidence_scores is not None:
        for row in equivalence_classes:
            members = tuple(str(value) for value in str(row["members"]).split(";"))
            resolution = evidence_lifted_resolution(
                members,
                evidence_scores,
                class_id=str(row["class_id"]),
                evidence_sources=evidence_sources,
            )
            resolution_rows.append(resolution.as_row())

    structural_status = (
        "unique_in_supplied_ion_space"
        if not ambiguous_classes
        else "equivalence_classes_present"
    )
    if evidence_scores is None:
        evidence_status = "not_supplied"
    elif not ambiguous_classes:
        evidence_status = "not_needed_for_structural_uniqueness"
    elif evidence_sources is None:
        evidence_status = "relative_scores_only"
    elif any(row["independent_evidence_available"] for row in resolution_rows):
        evidence_status = "evidence_lifted_partial_or_conditional"
    else:
        evidence_status = "unresolved_without_independent_evidence"

    return {
        "status": (
            "structurally_unique"
            if not ambiguous_classes
            else "structurally_non_unique"
        ),
        "structural_status": structural_status,
        "evidence_status": evidence_status,
        "n_reactions": n_reactions,
        "n_ions": n_ions,
        "rank": rank,
        "nullity": max(0, n_reactions - rank),
        "underdetermined": bool(n_reactions > n_ions),
        "rank_deficient": bool(rank < min(n_reactions, n_ions)),
        "condition_number": condition_number,
        "singular_values": [float(value) for value in singular_values],
        "maximum_absolute_coherence": float(maximum_coherence),
        "equivalence_class_count": len(equivalence_classes),
        "ambiguous_equivalence_class_count": len(ambiguous_classes),
        "equivalence_classes": equivalence_classes,
        "class_map": class_map,
        "evidence_lifted_resolution": resolution_rows,
        "claim_guardrail": (
            "A preferred reaction member remains conditional evidence; rank, "
            "coherence, and equivalence classes limit unique mechanism claims."
        ),
    }


def reaction_panel_diagnostics(
    matrix: Sequence[Sequence[float]],
    labels: Sequence[str],
    *,
    ion_order: Sequence[str] | None = None,
    observed_ions: Sequence[str] | None = None,
    cosine_threshold: float = 0.999999,
    signed: bool = True,
    evidence_scores: Mapping[str, float | None] | None = None,
    evidence_sources: Mapping[str, Iterable[str]] | None = None,
) -> dict[str, Any]:
    """Diagnose one observed reaction panel and gate unresolved carbonates.

    ``matrix`` is reaction-by-ion.  If ``ion_order`` is supplied, its length
    must match the matrix columns.  ``observed_ions`` selects an explicitly
    observed subset by name; no unobserved column is imputed or reconstructed.
    If the panel contains an ambiguous equivalence class involving a carbonate
    reaction and there is no independently sourced evidence-lifted resolution,
    the returned top-level ``status`` is ``"ABSTAIN"``.

    The status is a reporting gate, not a claim that carbonate chemistry is
    impossible.  It means this panel cannot support a unique carbonate
    mechanism under the supplied evidence contract.
    """

    arr = np.asarray(matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if not np.isfinite(arr).all():
        raise ValueError("matrix must contain only finite values")
    n_columns = int(arr.shape[1])
    if ion_order is None:
        resolved_ion_order = [f"ion_{index}" for index in range(n_columns)]
    else:
        resolved_ion_order = [str(ion).strip() for ion in ion_order]
        if len(resolved_ion_order) != n_columns:
            raise ValueError("ion_order length must match matrix columns")
        if any(not ion for ion in resolved_ion_order):
            raise ValueError("ion_order must not contain empty names")
        if len(set(resolved_ion_order)) != len(resolved_ion_order):
            raise ValueError("ion_order must contain unique names")

    if observed_ions is None:
        panel = list(resolved_ion_order)
        panel_matrix = arr
    else:
        panel = [str(ion).strip() for ion in observed_ions]
        if not panel or any(not ion for ion in panel):
            raise ValueError("observed_ions must contain non-empty names")
        if len(set(panel)) != len(panel):
            raise ValueError("observed_ions must contain unique names")
        unknown = [ion for ion in panel if ion not in resolved_ion_order]
        if unknown:
            raise ValueError(
                "observed_ions contains names absent from ion_order: "
                + ", ".join(unknown)
            )
        indices = [resolved_ion_order.index(ion) for ion in panel]
        panel_matrix = arr[:, indices]

    diagnostics = reaction_identifiability_diagnostics(
        panel_matrix,
        labels,
        cosine_threshold=cosine_threshold,
        signed=signed,
        evidence_scores=evidence_scores,
        evidence_sources=evidence_sources,
    )
    resolution_by_class = {
        str(row["class_id"]): row
        for row in diagnostics["evidence_lifted_resolution"]
    }

    unresolved_carbonate_classes: list[dict[str, Any]] = []
    for equivalence in diagnostics["equivalence_classes"]:
        members = tuple(
            value for value in str(equivalence["members"]).split(";") if value
        )
        carbonate_members = tuple(
            member
            for member in members
            if any(
                token in member.casefold()
                for token in ("calcite", "dolomite", "magnesite", "aragonite", "carbonate")
            )
        )
        if len(members) <= 1 or not carbonate_members:
            continue
        resolution = resolution_by_class.get(str(equivalence["class_id"]))
        resolved_by_independent_evidence = bool(
            resolution
            and resolution.get("resolution_status") == "evidence_lifted_resolved"
            and resolution.get("independent_evidence_available") is True
        )
        if not resolved_by_independent_evidence:
            unresolved_carbonate_classes.append(
                {
                    **equivalence,
                    "carbonate_members": list(carbonate_members),
                    "resolution_status": (
                        resolution.get("resolution_status")
                        if resolution
                        else "not_supplied"
                    ),
                }
            )

    panel_status = "ABSTAIN" if unresolved_carbonate_classes else "RUN"
    if unresolved_carbonate_classes:
        status_reason = (
            "carbonate equivalence class remains unresolved in the observed "
            "ion panel; unique carbonate mechanism attribution is not supported"
        )
    elif diagnostics["rank_deficient"]:
        status_reason = (
            "panel is rank deficient, but no unresolved carbonate equivalence "
            "class was identified"
        )
    else:
        status_reason = "no unresolved carbonate equivalence class in this panel"

    return {
        **diagnostics,
        "status": panel_status,
        "status_reason": status_reason,
        "observed_ion_panel": panel,
        "full_ion_order": resolved_ion_order,
        "n_observed_ions": len(panel),
        "unresolved_carbonate_classes": unresolved_carbonate_classes,
        "carbonate_status": (
            "ABSTAIN" if unresolved_carbonate_classes else "NO_UNRESOLVED_CLASS"
        ),
        "panel_claim_guardrail": (
            "ABSTAIN means this measured ion panel does not uniquely identify "
            "an unresolved carbonate mechanism; report an equivalence class or "
            "conditional evidence instead."
        ),
    }


__all__ = [
    "EvidenceLiftedResolution",
    "evidence_lifted_resolution",
    "evidence_score_map",
    "reaction_identifiability_diagnostics",
    "reaction_panel_diagnostics",
    "stoichiometric_equivalence_classes",
]
