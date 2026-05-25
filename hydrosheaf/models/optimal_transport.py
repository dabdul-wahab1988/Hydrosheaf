"""Reaction-aware optimal transport for chemistry-based edge plausibility.

Computes an unbalanced optimal transport cost between upstream and downstream
ion profiles. The cost matrix encodes "how hard" it is to convert one ion species
into another via known geochemical reactions.

Uses ``scipy.optimize.linprog`` (not POT) since the ion count is small.
"""

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linprog

from ..config import Config
from ..log import get_logger

logger = get_logger("models.optimal_transport")

# Species families for cost matrix construction
CONSERVATIVE_SPECIES = {"Cl", "Br", "Na"}
CHARGE_FAMILIES = {
    "Ca": 2, "Mg": 2, "Na": 1, "K": 1,
    "Cl": -1, "SO4": -2, "NO3": -1, "HCO3": -1,
    "F": -1, "Fe": 0, "PO4": -3,
}


def _build_cost_matrix(
    ion_order: List[str],
    config: Config,
) -> np.ndarray:
    """Build the per-species transport cost matrix C[i][j].

    C[i][j] is the cost of converting 1 unit of species i upstream into
    species j downstream.

    Rules:
    - Same species: cost 0 (no conversion needed)
    - Conservative species (Cl, Br, Na): high off-diagonal cost
    - Same charge family: moderate cost
    - Otherwise: high cost
    """
    n = len(ion_order)
    C = np.ones((n, n), dtype=float) * 5.0  # Default: expensive conversion

    for i, ion_i in enumerate(ion_order):
        for j, ion_j in enumerate(ion_order):
            if i == j:
                C[i, j] = 0.0  # Same species: free
            elif ion_i in CONSERVATIVE_SPECIES or ion_j in CONSERVATIVE_SPECIES:
                C[i, j] = 10.0  # Conservative species shouldn't change
            elif (
                ion_i in CHARGE_FAMILIES
                and ion_j in CHARGE_FAMILIES
                and CHARGE_FAMILIES[ion_i] != 0
                and CHARGE_FAMILIES[ion_j] != 0
                and np.sign(CHARGE_FAMILIES[ion_i]) == np.sign(CHARGE_FAMILIES[ion_j])
            ):
                C[i, j] = 1.0  # Same charge sign: plausible ion exchange
            elif (
                ion_i in CHARGE_FAMILIES
                and ion_j in CHARGE_FAMILIES
                and CHARGE_FAMILIES[ion_i] != 0
                and CHARGE_FAMILIES[ion_j] != 0
                and np.sign(CHARGE_FAMILIES[ion_i]) != np.sign(CHARGE_FAMILIES[ion_j])
            ):
                C[i, j] = 2.0  # Opposite charge: needs precipitation/dissolution

    return C


def _reaction_creation_cost(
    ion: str,
    config: Config,
) -> float:
    """Return the creation/destruction penalty for an ion based on known reactions."""
    ion_lower = ion.lower()

    # Check if the ion participates in any active mineral reactions
    from ..data.minerals import get_mineral_stoich

    supported = False
    for mineral_name in config.active_minerals:
        try:
            stoich = get_mineral_stoich(mineral_name)
            if ion in stoich and abs(stoich[ion]) > 1e-9:
                supported = True
                break
        except ValueError:
            continue

    if supported:
        return 0.5  # Reaction-supported: low penalty
    elif ion in CONSERVATIVE_SPECIES:
        return 5.0  # Conservative: strong penalty
    else:
        return 2.0  # Unsupported: moderate penalty


def compute_unbalanced_ot(
    upstream: Sequence[float],
    downstream: Sequence[float],
    ion_order: List[str],
    config: Config,
) -> Dict[str, float]:
    """Compute unbalanced optimal transport between upstream and downstream ions.

    Unbalanced OT allows mass creation/destruction with penalty, making it
    suitable for non-conservative chemical systems where reactions produce
    or consume species.

    Minimises:
        sum(T_ij * C_ij) + create_i * c_create_i + destroy_i * c_destroy_i
    subject to:
        outgoing_i + destroy_i = upstream_i
        incoming_i + create_i = downstream_i
        T_ij, create_i, destroy_i >= 0

    Parameters
    ----------
    upstream : sequence of float
        Ion concentrations at upstream node.
    downstream : sequence of float
        Ion concentrations at downstream node.
    ion_order : list of str
        Ion species names.
    config : Config
        Hydrosheaf configuration.

    Returns
    -------
    dict with keys:
        ot_total_cost : float
            Total unbalanced OT cost.
        ot_balanced_cost : float
            Transport-only cost (T_ij * C_ij term).
        ot_creation_mass : float
            Total mass created.
        ot_destruction_mass : float
            Total mass destroyed.
        ot_conservative_mismatch : float
            Cost attributable to conservative species mismatch.
        ot_reaction_plausibility : float
            Ratio of creation/destruction explained by known reactions (0-1).
    """
    n = len(ion_order)
    upstream_arr = np.array(upstream, dtype=float)
    downstream_arr = np.array(downstream, dtype=float)

    # Validate lengths
    if len(upstream_arr) != n or len(downstream_arr) != n:
        return {
            "ot_total_cost": 0.0,
            "ot_balanced_cost": 0.0,
            "ot_creation_mass": 0.0,
            "ot_destruction_mass": 0.0,
            "ot_conservative_mismatch": 0.0,
            "ot_reaction_plausibility": 0.0,
        }

    total_upstream = float(np.sum(upstream_arr))
    total_downstream = float(np.sum(downstream_arr))

    if total_upstream <= 0 or total_downstream <= 0:
        return {
            "ot_total_cost": 0.0,
            "ot_balanced_cost": 0.0,
            "ot_creation_mass": 0.0,
            "ot_destruction_mass": 0.0,
            "ot_conservative_mismatch": 0.0,
            "ot_reaction_plausibility": 0.0,
        }

    C = _build_cost_matrix(ion_order, config)

    # Creation/destruction penalties
    c_create = np.array([_reaction_creation_cost(ion, config) for ion in ion_order])
    c_destroy = np.array([_reaction_creation_cost(ion, config) for ion in ion_order])

    # Decision variables: T_ij flattened (n*n), then create_i (n), then destroy_i (n)
    n_T = n * n
    n_vars = n_T + 2 * n

    # Objective coefficients
    c_obj = np.zeros(n_vars)
    for i in range(n):
        for j in range(n):
            c_obj[i * n + j] = C[i, j]
    for i in range(n):
        c_obj[n_T + i] = c_create[i]
    for i in range(n):
        c_obj[n_T + n + i] = c_destroy[i]

    # Constraints
    A_rows = []
    b_vals = []

    # Outgoing constraint: sum_j T_ij + destroy_i = upstream_i
    for i in range(n):
        row = np.zeros(n_vars)
        for j in range(n):
            row[i * n + j] = 1.0
        row[n_T + n + i] = 1.0  # destroy_i
        A_rows.append(row)
        b_vals.append(upstream_arr[i])

    # Incoming constraint: sum_i T_ij + create_j = downstream_j
    for j in range(n):
        row = np.zeros(n_vars)
        for i in range(n):
            row[i * n + j] = 1.0
        row[n_T + j] = 1.0  # create_j
        A_rows.append(row)
        b_vals.append(downstream_arr[j])

    A_eq = np.array(A_rows)
    b_eq = np.array(b_vals)

    # Bounds: all variables >= 0
    bounds = [(0, None)] * n_vars

    # Solve
    try:
        result = linprog(
            c_obj,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
            options={"disp": False},
        )
    except Exception:
        # Fallback to interior point if highs not available
        try:
            result = linprog(
                c_obj,
                A_eq=A_eq,
                b_eq=b_eq,
                bounds=bounds,
                method="interior-point",
                options={"disp": False},
            )
        except Exception:
            return {
                "ot_total_cost": 0.0,
                "ot_balanced_cost": 0.0,
                "ot_creation_mass": 0.0,
                "ot_destruction_mass": 0.0,
                "ot_conservative_mismatch": 0.0,
                "ot_reaction_plausibility": 0.0,
            }

    if not result.success:
        logger.debug("OT LP did not converge: %s", result.message)
        return {
            "ot_total_cost": float(np.dot(c_obj, np.zeros(n_vars))),
            "ot_balanced_cost": 0.0,
            "ot_creation_mass": 0.0,
            "ot_destruction_mass": 0.0,
            "ot_conservative_mismatch": 0.0,
            "ot_reaction_plausibility": 0.0,
        }

    x = result.x

    # Extract costs
    transport_cost = 0.0
    for i in range(n):
        for j in range(n):
            transport_cost += x[i * n + j] * C[i, j]

    creation_mass = float(np.sum(x[n_T : n_T + n]))
    destruction_mass = float(np.sum(x[n_T + n : n_T + 2 * n]))

    creation_cost = float(np.dot(x[n_T : n_T + n], c_create))
    destruction_cost = float(np.dot(x[n_T + n : n_T + 2 * n], c_destroy))

    total_cost = transport_cost + creation_cost + destruction_cost

    # Conservative mismatch: cost from conservative species off-diagonal transport
    conservative_mismatch = 0.0
    for i, ion_i in enumerate(ion_order):
        for j, ion_j in enumerate(ion_order):
            if i != j and (ion_i in CONSERVATIVE_SPECIES or ion_j in CONSERVATIVE_SPECIES):
                conservative_mismatch += x[i * n + j] * C[i, j]

    # Reaction plausibility: fraction of off-diagonal transport that moves
    # between ions in the same charge family (capped at 1.0).
    plausible_mass = 0.0
    for i, ion_i in enumerate(ion_order):
        for j, ion_j in enumerate(ion_order):
            if i != j and ion_i in CHARGE_FAMILIES and ion_j in CHARGE_FAMILIES:
                plausible_mass += x[i * n + j]

    total_off_diag = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                total_off_diag += x[i * n + j]

    reaction_plausibility = min(1.0, plausible_mass / max(1e-12, total_off_diag)) if total_off_diag > 1e-12 else 1.0

    # Normalize by total mass for scale-invariant comparison
    scale = max(total_upstream, total_downstream, 1e-12)
    total_cost /= scale
    transport_cost /= scale
    creation_mass /= scale
    destruction_mass /= scale
    conservative_mismatch /= scale

    return {
        "ot_total_cost": total_cost,
        "ot_balanced_cost": transport_cost,
        "ot_creation_mass": creation_mass,
        "ot_destruction_mass": destruction_mass,
        "ot_conservative_mismatch": conservative_mismatch,
        "ot_reaction_plausibility": reaction_plausibility,
    }
