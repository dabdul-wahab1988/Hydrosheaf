"""Shared definitions and numerical helpers for the M5 benchmark."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import os
import subprocess
import time
from typing import Iterable, Mapping, Sequence

import numpy as np


ION_ORDER = [
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
]

MOLAR_MASS_G_MOL = {
    "Ca": 40.078,
    "Mg": 24.305,
    "Na": 22.989769,
    "K": 39.0983,
    "HCO3": 61.0168,
    "Cl": 35.453,
    "SO4": 96.06,
    "NO3": 62.0049,
    "F": 18.998403,
    "Fe": 55.845,
    "PO4": 94.9714,
}

PHREEQC_TOTALS = {
    "Ca": "Ca",
    "Mg": "Mg",
    "Na": "Na",
    "K": "K",
    "HCO3": "C(4)",
    "Cl": "Cl",
    "SO4": "S(6)",
    "NO3": "N(5)",
    "F": "F",
    "Fe": "Fe",
    "PO4": "P",
}


@dataclass(frozen=True)
class ReactionDefinition:
    label: str
    family: str
    stoich: Mapping[str, float]
    phreeqc_components: tuple[tuple[str, float], ...]
    signed: bool
    si_phase: str | None = None


REACTIONS = [
    ReactionDefinition(
        "calcite",
        "carbonate",
        {"Ca": 1.0, "HCO3": 2.0},
        (("Ca(HCO3)2", 1.0),),
        True,
        "Calcite",
    ),
    ReactionDefinition(
        "dolomite",
        "carbonate",
        {"Ca": 1.0, "Mg": 1.0, "HCO3": 4.0},
        (("CaMg(HCO3)4", 1.0),),
        True,
        "Dolomite",
    ),
    ReactionDefinition(
        "albite",
        "silicate",
        {"Na": 1.0, "HCO3": 1.0},
        (("NaHCO3", 1.0),),
        True,
        "Albite",
    ),
    ReactionDefinition(
        "anorthite",
        "silicate",
        {"Ca": 1.0, "HCO3": 2.0},
        (("Ca(HCO3)2", 1.0),),
        True,
        "Anorthite",
    ),
    ReactionDefinition(
        "k_feldspar",
        "silicate",
        {"K": 1.0, "HCO3": 1.0},
        (("KHCO3", 1.0),),
        True,
        "K-feldspar",
    ),
    ReactionDefinition(
        "gypsum",
        "evaporite",
        {"Ca": 1.0, "SO4": 1.0},
        (("CaSO4:2H2O", 1.0),),
        True,
        "Gypsum",
    ),
    ReactionDefinition(
        "halite",
        "evaporite",
        {"Na": 1.0, "Cl": 1.0},
        (("NaCl", 1.0),),
        True,
        "Halite",
    ),
    ReactionDefinition(
        "fluorite",
        "trace_mineral",
        {"Ca": 1.0, "F": 2.0},
        (("CaF2", 1.0),),
        True,
        "Fluorite",
    ),
    ReactionDefinition(
        "apatite",
        "trace_mineral",
        {"Ca": 5.0, "PO4": 3.0, "F": 1.0},
        (("Ca5(PO4)3F", 1.0),),
        True,
        "Fluorapatite",
    ),
    ReactionDefinition(
        "pyrite_oxidation",
        "redox",
        {"SO4": 2.0, "Fe": 1.0},
        (("Fe(SO4)2", 1.0),),
        False,
        None,
    ),
    ReactionDefinition(
        "NO3src",
        "anthropogenic",
        {"NO3": 1.0},
        (("NaNO3", 1.0), ("NaCl", -1.0), ("HCl", 1.0)),
        False,
        None,
    ),
    ReactionDefinition(
        "denit",
        "redox",
        {"HCO3": 1.0, "NO3": -1.0},
        (("NaHCO3", 1.0), ("NaCl", -1.0), ("HCl", 1.0), ("HNO3", -1.0)),
        False,
        None,
    ),
    ReactionDefinition(
        "CaNa_exch",
        "ion_exchange",
        {"Ca": 1.0, "Na": -2.0},
        (("CaCl2", 1.0), ("NaCl", -2.0)),
        False,
        None,
    ),
    ReactionDefinition(
        "NaCa_exch",
        "ion_exchange",
        {"Ca": -1.0, "Na": 2.0},
        (("CaCl2", -1.0), ("NaCl", 2.0)),
        False,
        None,
    ),
    ReactionDefinition(
        "MgNa_exch",
        "ion_exchange",
        {"Mg": 1.0, "Na": -2.0},
        (("MgCl2", 1.0), ("NaCl", -2.0)),
        False,
        None,
    ),
    ReactionDefinition(
        "NaMg_exch",
        "ion_exchange",
        {"Mg": -1.0, "Na": 2.0},
        (("MgCl2", -1.0), ("NaCl", 2.0)),
        False,
        None,
    ),
]

REACTION_BY_LABEL = {reaction.label: reaction for reaction in REACTIONS}
REACTION_LABELS = [reaction.label for reaction in REACTIONS]


def reaction_matrix(ions: Sequence[str] = ION_ORDER) -> np.ndarray:
    """Return reactions x ions stoichiometric matrix."""
    return np.asarray(
        [[reaction.stoich.get(ion, 0.0) for ion in ions] for reaction in REACTIONS],
        dtype=float,
    )


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


def equivalence_classes(
    ions: Sequence[str] = ION_ORDER,
    cosine_threshold: float = 0.999999,
) -> tuple[dict[str, str], list[dict[str, object]]]:
    """Build exact signed stoichiometric equivalence classes."""
    matrix = reaction_matrix(ions)
    graph = {label: set() for label in REACTION_LABELS}
    for i, left in enumerate(REACTION_LABELS):
        for j in range(i + 1, len(REACTION_LABELS)):
            right = REACTION_LABELS[j]
            norm = np.linalg.norm(matrix[i]) * np.linalg.norm(matrix[j])
            cosine = abs(float(np.dot(matrix[i], matrix[j]) / norm)) if norm else 0.0
            if cosine >= cosine_threshold:
                graph[left].add(right)
                graph[right].add(left)
    class_map: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    for index, members in enumerate(_connected_components(graph), start=1):
        class_id = f"EC{index:02d}"
        for label in members:
            class_map[label] = class_id
        rows.append(
            {
                "class_id": class_id,
                "members": ";".join(members),
                "n_members": len(members),
                "ambiguous": len(members) > 1,
            }
        )
    return class_map, rows


def matrix_diagnostics(ions: Sequence[str]) -> dict[str, float]:
    """Calculate structural identifiability diagnostics for an ion panel."""
    matrix = reaction_matrix(ions)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1e-10))
    positive = singular[singular > 1e-10]
    condition = float(positive.max() / positive.min()) if positive.size else math.inf
    normalized = matrix / np.maximum(np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)
    coherence = 0.0
    for i in range(len(normalized)):
        for j in range(i + 1, len(normalized)):
            coherence = max(coherence, abs(float(np.dot(normalized[i], normalized[j]))))

    class_map, classes = equivalence_classes(ions)
    representatives = []
    for item in classes:
        representatives.append(REACTION_LABELS.index(str(item["members"]).split(";")[0]))
    collapsed = matrix[representatives]
    collapsed_norm = collapsed / np.maximum(
        np.linalg.norm(collapsed, axis=1, keepdims=True), 1e-12
    )
    class_coherence = 0.0
    for i in range(len(collapsed_norm)):
        for j in range(i + 1, len(collapsed_norm)):
            class_coherence = max(
                class_coherence,
                abs(float(np.dot(collapsed_norm[i], collapsed_norm[j]))),
            )
    return {
        "n_ions": len(ions),
        "n_reactions": len(REACTIONS),
        "rank": rank,
        "nullity": len(REACTIONS) - rank,
        "condition_number": condition,
        "mutual_coherence": coherence,
        "class_collapsed_coherence": class_coherence,
        "n_equivalence_classes": len(set(class_map.values())),
        "minimum_singular_value": float(singular.min()) if singular.size else 0.0,
        "maximum_singular_value": float(singular.max()) if singular.size else 0.0,
    }


def class_set(labels: Iterable[str], class_map: Mapping[str, str]) -> set[str]:
    return {class_map[label] for label in labels}


def support_from_extents(
    extents: Sequence[float],
    threshold: float = 0.015,
) -> set[str]:
    return {
        label
        for label, extent in zip(REACTION_LABELS, extents)
        if abs(float(extent)) >= threshold
    }


def precision_recall_f1(
    truth: set[str],
    prediction: set[str],
) -> tuple[float, float, float, float]:
    if not prediction:
        precision = 1.0 if not truth else 0.0
    else:
        precision = len(truth & prediction) / len(prediction)
    recall = len(truth & prediction) / len(truth) if truth else 1.0
    f1 = 0.0 if precision + recall == 0 else 2.0 * precision * recall / (precision + recall)
    false_discovery = 1.0 - precision if prediction else 0.0
    return precision, recall, f1, false_discovery


def jaccard(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 1.0
    return len(left & right) / len(left | right)


def thermodynamic_bounds(
    upstream_si: Mapping[str, float],
    threshold: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert saturation indices to simple dissolution/precipitation bounds."""
    lower = np.zeros(len(REACTIONS), dtype=float)
    upper = np.full(len(REACTIONS), np.inf, dtype=float)
    for index, reaction in enumerate(REACTIONS):
        if reaction.signed:
            lower[index] = -np.inf
        if not reaction.si_phase:
            continue
        value = upstream_si.get(reaction.si_phase, np.nan)
        if not np.isfinite(value):
            continue
        if value > threshold:
            upper[index] = 0.0
        elif value < -threshold:
            lower[index] = 0.0
    return lower, upper


def count_thermodynamic_violations(
    extents: Sequence[float],
    upstream_si: Mapping[str, float],
    threshold: float = 0.1,
    support_threshold: float = 0.015,
) -> int:
    count = 0
    for reaction, extent in zip(REACTIONS, extents):
        if not reaction.si_phase or abs(float(extent)) < support_threshold:
            continue
        value = upstream_si.get(reaction.si_phase, np.nan)
        if not np.isfinite(value):
            continue
        if (extent > 0 and value > threshold) or (extent < 0 and value < -threshold):
            count += 1
    return count


def fit_inverse(
    residual: np.ndarray,
    ions: Sequence[str],
    method: str,
    lambda_l1: float,
    lambda_l2: float,
    upstream_si: Mapping[str, float] | None = None,
    si_threshold: float = 0.1,
    penalty_scales: Sequence[float] | None = None,
    max_iter: int = 250,
) -> dict[str, object]:
    """Fit one inverse model with column normalization and optional SI bounds."""
    matrix = reaction_matrix(ions)
    column_norms = np.maximum(np.linalg.norm(matrix, axis=1), 1e-12)
    scaled = matrix / column_norms[:, None]
    if penalty_scales is None:
        penalty_array = np.ones(len(REACTIONS), dtype=float)
    else:
        penalty_array = np.asarray(penalty_scales, dtype=float)
        if len(penalty_array) != len(REACTIONS):
            raise ValueError(
                "penalty_scales must have one value per M5 reaction "
                f"({len(REACTIONS)} expected, {len(penalty_array)} found)."
            )
        penalty_array = np.clip(penalty_array, 0.05, 20.0)
    signed_mask = [reaction.signed for reaction in REACTIONS]
    lower = np.asarray([-np.inf if signed else 0.0 for signed in signed_mask], dtype=float)
    upper = np.full(len(REACTIONS), np.inf, dtype=float)
    if method.startswith("thermo"):
        lower, upper = thermodynamic_bounds(upstream_si or {}, threshold=si_threshold)

    beta_lower = lower * column_norms
    beta_upper = upper * column_norms
    started = time.perf_counter()
    design = scaled.T
    lipschitz = float(np.linalg.norm(design, ord=2) ** 2 + lambda_l2)
    lipschitz = max(lipschitz, 1e-12)
    beta = np.zeros(len(REACTIONS), dtype=float)
    accelerated = beta.copy()
    momentum = 1.0
    signed_array = np.asarray(signed_mask, dtype=bool)
    converged = False
    for iterations in range(1, max_iter + 1):
        gradient = (
            design.T @ (design @ accelerated - residual)
            + lambda_l2 * accelerated
        )
        candidate = accelerated - gradient / lipschitz
        thresholds = lambda_l1 * penalty_array / lipschitz
        updated = np.sign(candidate) * np.maximum(
            np.abs(candidate) - thresholds, 0.0
        )
        updated[~signed_array] = np.maximum(
            candidate[~signed_array] - thresholds[~signed_array], 0.0
        )
        updated = np.minimum(np.maximum(updated, beta_lower), beta_upper)
        if np.max(np.abs(updated - beta)) <= 1e-6:
            beta = updated
            converged = True
            break
        next_momentum = 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * momentum**2))
        accelerated = updated + (
            (momentum - 1.0) / next_momentum
        ) * (updated - beta)
        beta = updated
        momentum = next_momentum
    extents = beta / column_norms
    prediction = matrix.T @ extents
    post_residual = residual - prediction
    return {
        "extents": extents,
        "prediction": prediction,
        "residual": post_residual,
        "rmse": float(np.sqrt(np.mean(post_residual**2))),
        "l1_norm": float(np.sum(np.abs(extents))),
        "iterations": iterations,
        "converged": converged,
        "runtime_ms": 1000.0 * (time.perf_counter() - started),
        "thermodynamic_violations": count_thermodynamic_violations(
            extents,
            upstream_si or {},
            threshold=si_threshold,
        ),
    }


def find_phreeqc() -> tuple[Path, Path]:
    """Find the local USGS PHREEQC executable and database."""
    executable_candidates = [
        os.environ.get("PHREEQC_EXE"),
        r"C:\Program Files\USGS\phreeqc-3.7.3-15968-x64\bin\Release\phreeqc.exe",
        r"C:\Program Files\USGS\phreeqc-3.7.3-15968-x64\bin\ClrRelease\phreeqc.exe",
        r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\Release\phreeqc.exe",
        r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\ClrRelease\phreeqc.exe",
    ]
    database_candidates = [
        os.environ.get("PHREEQC_DATABASE"),
        r"C:\Program Files\USGS\phreeqc-3.7.3-15968-x64\database\phreeqc.dat",
        r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\database\phreeqc.dat",
    ]
    executable = next(
        (Path(value) for value in executable_candidates if value and Path(value).is_file()),
        None,
    )
    database = next(
        (Path(value) for value in database_candidates if value and Path(value).is_file()),
        None,
    )
    if executable is None or database is None:
        raise FileNotFoundError(
            "PHREEQC was not found. Set PHREEQC_EXE and PHREEQC_DATABASE."
        )
    return executable, database


def run_phreeqc(
    input_path: Path,
    output_path: Path,
    database_path: Path,
    executable_path: Path,
) -> float:
    """Run PHREEQC and return elapsed seconds."""
    started = time.perf_counter()
    completed = subprocess.run(
        [
            str(executable_path),
            input_path.name,
            output_path.name,
            str(database_path),
        ],
        cwd=input_path.parent,
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"PHREEQC failed ({completed.returncode}).\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    output_text = output_path.read_text(encoding="utf-8", errors="replace")
    if "ERROR:" in output_text or "ERRORS:" in output_text:
        raise RuntimeError(f"PHREEQC reported an input error in {output_path}.")
    return elapsed


def finite_or_blank(value: float) -> float | str:
    return float(value) if np.isfinite(value) else ""
