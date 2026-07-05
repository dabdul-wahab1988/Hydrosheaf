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
from typing import Iterable, Mapping, Sequence

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

    graph = {str(label): set() for label in labels}
    for i, left in enumerate(labels):
        for j in range(i + 1, len(labels)):
            right = str(labels[j])
            norm = float(np.linalg.norm(arr[i]) * np.linalg.norm(arr[j]))
            cosine = float(np.dot(arr[i], arr[j]) / norm) if norm else 0.0
            similarity = abs(cosine) if signed else cosine
            if similarity >= cosine_threshold:
                graph[str(left)].add(right)
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
    evidence_scores: Mapping[str, float],
    *,
    class_id: str = "",
    score_floor: float = 1e-9,
    weight_transform: str = "odds",
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
    clean_members = tuple(str(member) for member in members)
    if not clean_members:
        raise ValueError("members must contain at least one reaction")

    scores: dict[str, float] = {}
    for member in clean_members:
        value = float(evidence_scores.get(member, 0.0))
        if not math.isfinite(value):
            value = 0.0
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

    if len(clean_members) == 1:
        status = "structurally_unique"
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
    )


def evidence_score_map(
    labels: Iterable[str],
    scores: Iterable[float],
) -> dict[str, float]:
    """Create a reaction-score mapping with float coercion."""
    return {str(label): float(score) for label, score in zip(labels, scores)}
