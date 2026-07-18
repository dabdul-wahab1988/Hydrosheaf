"""M6 extended reaction dictionary and self-contained inverse solver.

Extends the M5 reaction set by making SiO2 and Sr *genuine reaction-basis
species*, not hand-imposed evidence gates:
  - silicate weathering (albite, K-feldspar, biotite) releases dissolved SiO2;
  - carbonate / gypsum dissolution releases trace Sr.
Consequently the diagnostic value of Sr/SiO2 becomes a property of the
stoichiometric matrix (rank, null space, coherence) and of ground-truth
reaction recovery — NOT a definitional rule. Removing SiO2/Sr from the measured
panel then genuinely reduces the structural identifiability of the silicate /
carbonate reactions, which the synthetic-field validation confirms against known
truth. This directly answers the circularity critique of the evidence-lift gate.

Nothing here modifies the frozen M5 package; it is a documented M6 extension used
for structural diagnostics and the synthetic-field validation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

# Extended species basis (M5 11 + SiO2 + Sr)
SPECIES = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4",
           "SiO2", "Sr"]


@dataclass(frozen=True)
class Rxn:
    label: str
    family: str
    stoich: Mapping[str, float]
    signed: bool  # True: mineral (dissolve +/precipitate -); False: one-way process


# Stoichiometry as products released per unit dissolution (mmol basis).
# SiO2 on silicates and Sr on carbonate/sulfate are the load-bearing additions.
REACTIONS = [
    Rxn("calcite", "carbonate", {"Ca": 1.0, "HCO3": 2.0, "Sr": 0.015}, True),
    Rxn("dolomite", "carbonate", {"Ca": 1.0, "Mg": 1.0, "HCO3": 4.0, "Sr": 0.020}, True),
    Rxn("albite", "silicate", {"Na": 1.0, "HCO3": 1.0, "SiO2": 2.0}, True),
    Rxn("anorthite", "silicate", {"Ca": 1.0, "HCO3": 2.0, "SiO2": 0.0}, True),
    Rxn("k_feldspar", "silicate", {"K": 1.0, "HCO3": 1.0, "SiO2": 2.0}, True),
    Rxn("biotite", "silicate", {"Mg": 0.5, "K": 1.0, "HCO3": 1.0, "SiO2": 1.5}, True),
    Rxn("gypsum", "evaporite", {"Ca": 1.0, "SO4": 1.0, "Sr": 0.05}, True),
    Rxn("halite", "evaporite", {"Na": 1.0, "Cl": 1.0}, True),
    Rxn("fluorite", "trace_mineral", {"Ca": 1.0, "F": 2.0}, True),
    Rxn("apatite", "trace_mineral", {"Ca": 5.0, "PO4": 3.0, "F": 1.0}, True),
    Rxn("pyrite_oxidation", "redox", {"SO4": 2.0, "Fe": 1.0}, False),
    Rxn("NO3src", "anthropogenic", {"NO3": 1.0}, False),
    Rxn("denit", "redox", {"HCO3": 1.0, "NO3": -1.0}, False),
    Rxn("CaNa_exch", "ion_exchange", {"Ca": 1.0, "Na": -2.0}, False),
    Rxn("NaCa_exch", "ion_exchange", {"Ca": -1.0, "Na": 2.0}, False),
    Rxn("MgNa_exch", "ion_exchange", {"Mg": 1.0, "Na": -2.0}, False),
    Rxn("NaMg_exch", "ion_exchange", {"Mg": -1.0, "Na": 2.0}, False),
]
LABELS = [r.label for r in REACTIONS]
FAMILY_OF = {r.label: r.family for r in REACTIONS}
SIGNED = np.array([r.signed for r in REACTIONS], dtype=bool)


def reaction_matrix(species: Sequence[str] = SPECIES) -> np.ndarray:
    """(n_reactions x n_species) stoichiometric matrix over a measured panel."""
    return np.array([[r.stoich.get(s, 0.0) for s in species] for r in REACTIONS],
                    dtype=float)


# --- structural identifiability diagnostics (gate-independent) ----------------
def _connected(graph):
    seen, comps = set(), []
    for n in graph:
        if n in seen:
            continue
        stack, comp = [n], []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x); comp.append(x)
            stack.extend(graph[x] - seen)
        comps.append(sorted(comp))
    return comps


def equivalence_classes(species: Sequence[str] = SPECIES, thr: float = 0.9995):
    """Exact-ish stoichiometric equivalence classes on a measured panel."""
    M = reaction_matrix(species)
    graph = {l: set() for l in LABELS}
    for i in range(len(LABELS)):
        for j in range(i + 1, len(LABELS)):
            ni, nj = np.linalg.norm(M[i]), np.linalg.norm(M[j])
            if ni == 0 or nj == 0:
                continue
            if abs(float(M[i] @ M[j] / (ni * nj))) >= thr:
                graph[LABELS[i]].add(LABELS[j]); graph[LABELS[j]].add(LABELS[i])
    cmap = {}
    for k, comp in enumerate(_connected(graph), 1):
        for l in comp:
            cmap[l] = f"EC{k:02d}"
    return cmap


def diagnostics(species: Sequence[str]) -> dict:
    M = reaction_matrix(species)
    # drop reactions with no measured species (all-zero rows) from the structural view
    active = np.where(np.any(M != 0, axis=1))[0]
    Ma = M[active]
    sv = np.linalg.svd(Ma, compute_uv=False) if len(Ma) else np.array([0.0])
    rank = int(np.linalg.matrix_rank(Ma, tol=1e-9)) if len(Ma) else 0
    pos = sv[sv > 1e-9]
    cond = float(pos.max() / pos.min()) if pos.size else math.inf
    norm = Ma / np.maximum(np.linalg.norm(Ma, axis=1, keepdims=True), 1e-12)
    coh = 0.0
    for i in range(len(norm)):
        for j in range(i + 1, len(norm)):
            coh = max(coh, abs(float(norm[i] @ norm[j])))
    cmap = equivalence_classes(species)
    return {
        "n_species": len(species), "n_active_reactions": len(active),
        "rank": rank, "nullity": len(active) - rank,
        "condition_number": cond, "mutual_coherence": coh,
        "n_equivalence_classes": len(set(cmap.values())),
        "min_singular_value": float(sv.min()) if sv.size else 0.0,
    }


def silicate_identifiability(species: Sequence[str]) -> float:
    """Max mutual coherence between the SiO2-bearing silicate reactions
    (albite, K-feldspar, biotite — the Na/K silicates whose signature SiO2
    carries) and any other reaction on this panel. 0 = distinguishable,
    1 = collinear/unidentifiable. Gate-independent; SiO2 lowers it. The
    anorthite/calcite (Ca-silicate vs carbonate) collinearity is excluded
    because SiO2 genuinely cannot resolve it."""
    M = reaction_matrix(species)
    norm = M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-12)
    sio2_silicates = [i for i, r in enumerate(REACTIONS)
                      if r.family == "silicate" and r.stoich.get("SiO2", 0) > 0
                      and np.linalg.norm(M[i]) > 0]
    worst = 0.0
    for i in sio2_silicates:
        for j, r in enumerate(REACTIONS):
            if r.family == "silicate" or np.linalg.norm(M[j]) == 0:
                continue  # silicate-vs-NON-silicate discrimination only
            worst = max(worst, abs(float(norm[i] @ norm[j])))
    return worst


# --- self-contained FISTA inverse (matrix-parameterised) ----------------------
def fit(residual: np.ndarray, species: Sequence[str], lambda_l1: float = 0.01,
        lambda_l2: float = 1e-3, penalty_scales: np.ndarray | None = None,
        lower: np.ndarray | None = None, upper: np.ndarray | None = None,
        max_iter: int = 400) -> dict:
    M = reaction_matrix(species)              # (R x S)
    col = np.maximum(np.linalg.norm(M, axis=1), 1e-12)
    scaled = M / col[:, None]
    design = scaled.T                          # (S x R)
    r = np.asarray(residual, float)
    pen = np.ones(len(REACTIONS)) if penalty_scales is None else np.clip(
        np.asarray(penalty_scales, float), 0.05, 50.0)
    lo = (np.where(SIGNED, -np.inf, 0.0) if lower is None else lower) * col
    hi = (np.full(len(REACTIONS), np.inf) if upper is None else upper) * col
    L = max(float(np.linalg.norm(design, 2) ** 2 + lambda_l2), 1e-12)
    beta = np.zeros(len(REACTIONS)); acc = beta.copy(); mom = 1.0
    for _ in range(max_iter):
        grad = design.T @ (design @ acc - r) + lambda_l2 * acc
        cand = acc - grad / L
        thr = lambda_l1 * pen / L
        upd = np.sign(cand) * np.maximum(np.abs(cand) - thr, 0.0)
        upd[~SIGNED] = np.maximum(cand[~SIGNED] - thr[~SIGNED], 0.0)
        upd = np.minimum(np.maximum(upd, lo), hi)
        if np.max(np.abs(upd - beta)) <= 1e-7:
            beta = upd
            break
        m2 = 0.5 * (1 + math.sqrt(1 + 4 * mom ** 2))
        acc = upd + ((mom - 1) / m2) * (upd - beta)
        beta = upd; mom = m2
    extents = beta / col
    pred = M.T @ extents
    return {"extents": extents, "prediction": pred, "rmse": float(np.sqrt(np.mean((r - pred) ** 2))),
            "l1_norm": float(np.abs(extents).sum())}


def support(extents: Sequence[float], thr: float = 0.02) -> set:
    return {LABELS[i] for i, e in enumerate(extents) if abs(float(e)) >= thr}


def prf1(truth: set, pred: set):
    if not pred:
        p = 1.0 if not truth else 0.0
    else:
        p = len(truth & pred) / len(pred)
    r = len(truth & pred) / len(truth) if truth else 1.0
    f1 = 0.0 if p + r == 0 else 2 * p * r / (p + r)
    return p, r, f1


def class_f1(truth: set, pred: set, cmap: Mapping[str, str]):
    tc = {cmap[x] for x in truth if x in cmap}
    pc = {cmap[x] for x in pred if x in cmap}
    return prf1(tc, pc)[2]


def true_resolution_class(truth: set, pred: set, cmap: Mapping[str, str]) -> str:
    """Ground-truth identifiability class (needs known truth) — independent of
    any evidence gate. Mirrors the M5 F1 thresholds."""
    _, _, f1 = prf1(truth, pred)
    cf1 = class_f1(truth, pred, cmap)
    if f1 >= 0.80:
        return "identifiable"
    if cf1 >= 0.80:
        return "equivalence_class"
    if cf1 >= 0.50:
        return "partially_identifiable"
    return "non_identifiable"


# Panels for the synthetic validation (reaction-basis species only)
PANELS = {
    "majors": ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"],
    "plus_F": ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F"],
    "plus_Sr_SiO2": ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "SiO2", "Sr"],
}


if __name__ == "__main__":
    for name, sp in PANELS.items():
        d = diagnostics(sp)
        print(f"{name:14s} rank={d['rank']} nullity={d['nullity']} "
              f"coh={d['mutual_coherence']:.3f} nEC={d['n_equivalence_classes']} "
              f"sil_coh={silicate_identifiability(sp):.3f}")
