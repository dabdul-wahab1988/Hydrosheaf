"""Reaction dictionary and sparse fitting."""

from dataclasses import dataclass
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple
import math
import numpy as np

from ..config import Config, DEFAULT_ION_ORDER
from ..data.minerals import get_mineral_stoich
from ..log import get_logger

logger = get_logger("models.reactions")

# Mapping of reactions to their mandatory indicator ions (M2-M5 PhD Remediation)
# If honest_modeling is enabled, these reactions are pruned if indicators are missing.
INDICATOR_IONS = {
    "pyrite_oxidation_aerobic": ["Fe", "SO4"],
    "pyrite_oxidation_denit": ["Fe", "SO4"],
    "biotite": ["K", "Mg", "Fe"],
    "chlorite": ["Mg", "Fe"],
    "fluorite": ["F", "Ca"],
    "sylvite": ["K", "Cl"],
}

def _vector_from_coeffs(coeffs: Mapping[str, float], ion_order: Iterable[str]) -> List[float]:
    return [float(coeffs.get(ion, 0.0)) for ion in ion_order]


def build_reaction_dictionary(
    config: Config,
    pre_si_mask: Optional[Mapping[str, float]] = None,
    sample: Optional[Mapping[str, float]] = None,
    dynamic_denit_scale: float = 1.0,
) -> Tuple[List[List[float]], List[str], List[bool], List[float]]:
    ion_order = config.ion_order or DEFAULT_ION_ORDER
    kappa = config.denit_kappa

    # Identify available data (measured ions)
    available = set(config.measured_ions) if config.measured_ions else set(ion_order)

    reactions: List[Tuple[str, Mapping[str, float], bool, float]] = []

    # Geologic Bias Scaling (M2-M5 PhD Remediation & Winner-Takes-All Priority)
    def get_penalty_scale(name: str) -> float:
        n = name.lower().replace(" ", "_")
        if config.geologic_bias == "crystalline":
            if any(s in n for s in ["albite", "anorthite", "feldspar", "biotite", "chlorite", "pyroxene"]):
                return 0.5
            if any(e in n for e in ["gypsum", "anhydrite", "halite", "sylvite"]):
                return 10.0
            # Winner-Takes-All: Massive penalty for carbonates in basement rock unless strongly supported
            if any(c in n for c in ["calcite", "dolomite", "magnesite", "aragonite"]):
                return 50.0 
        elif config.geologic_bias == "sedimentary":
            if any(c in n for c in ["calcite", "dolomite", "magnesite", "aragonite"]):
                return 0.5
            # Winner-Takes-All: Massive penalty for primary silicates in mature sedimentary basins
            if any(s in n for s in ["albite", "anorthite", "feldspar", "biotite", "chlorite", "pyroxene"]):
                return 50.0
        return 1.0

    # 1. Add User-Selected Minerals
    for name in config.active_minerals:
        # Check Technical Remediation: Honest Modeling
        if config.honest_modeling:
            indicators = INDICATOR_IONS.get(name.lower().replace(" ", "_"))
            if indicators:
                missing = [ion for ion in indicators if ion not in available]
                if missing:
                    # Specific Remediation: Pyrite Auto-Substitution (Fixes Flaw 1)
                    if name.lower() == "pyrite_oxidation_aerobic" and "Fe" in missing:
                        logger.warning(f"Honest Modeling: Substituting 'pyrite_net' for '{name}' (Iron not measured).")
                        try:
                            reactions.append(("pyrite_net", get_mineral_stoich("pyrite_net"), True, get_penalty_scale("pyrite_net")))
                            continue
                        except ValueError: pass

                    logger.warning(f"Honest Modeling: Pruning mineral '{name}' because indicators {missing} are not measured.")
                    continue

        # Check Thermodynamic Logic Gate
        if config.use_thermodynamic_logic_gates and pre_si_mask is not None:
            si_val = pre_si_mask.get(name)
            if si_val is not None and si_val > config.si_logic_gate_threshold:
                logger.debug(f"Logic Gate: Pruning mineral '{name}' (SI={si_val:.2f} > {config.si_logic_gate_threshold})")
                continue

        # Check Concentration Logic Gate
        if sample is not None:
            so4_val = float(sample.get("SO4") or 0.0)
            cl_val = float(sample.get("Cl") or 0.0)
            f_val = float(sample.get("F") or 0.0)
            
            if name.lower() in ["gypsum", "anhydrite"]:
                if so4_val < 0.208:
                    logger.debug(f"Logic Gate: Pruning '{name}' (SO4 < 0.208 mmol/L)")
                    continue
                if config.geologic_bias == "crystalline" and so4_val > 0.52:
                    logger.debug(f"Logic Gate: Pruning '{name}' in crystalline setting (SO4 > 0.52 mmol/L)")
                    continue
            elif name.lower() in ["halite", "sylvite"]:
                if cl_val < 0.564:
                    logger.debug(f"Logic Gate: Pruning '{name}' (Cl < 0.564 mmol/L)")
                    continue
            elif name.lower() == "fluorite":
                if f_val < 0.026:
                    logger.debug(f"Logic Gate: Pruning '{name}' (F < 0.026 mmol/L)")
                    continue

        try:
            stoich = get_mineral_stoich(name)
            reactions.append((name, stoich, True, get_penalty_scale(name)))
        except ValueError:
            continue

    # 2. Add Process Reactions
    has_so4_source = any(lbl in ["gypsum", "anhydrite", "pyrite_net", "SO4_input"] for lbl, _, _, _ in reactions)
    if not has_so4_source and "SO4" in available:
        try: reactions.append(("SO4_input", get_mineral_stoich("SO4_input"), False, 1.0))
        except ValueError: pass

    add_no3src = True
    add_denit = True
    no3src_scale = 1.0
    if sample is not None:
        if sample.get("NO3", 0.0) < 0.16:
            add_denit = False
            logger.debug("Logic Gate: Pruning 'denit' (NO3 < 0.16 mmol/L)")
        if sample.get("NO3", 0.0) > 0.8:
            no3src_scale = 0.1
            logger.debug("Logic Gate: Forcing 'NO3src' selection (NO3 > 0.8 mmol/L)")

    if add_no3src:
        reactions.append(("NO3src", {"NO3": 1}, False, no3src_scale))
    if add_denit:
        reactions.append(("denit", {"HCO3": kappa, "NO3": -1}, False, 1.0))

    if config.exchange_enabled:
        # Bidirectional Exchange: Forward = Ca/Mg release (salinization), Reverse = Na release (freshening)
        reactions.append(("CaNa_exch", {"Ca": 1, "Na": -2}, False, 1.0))
        reactions.append(("NaCa_exch", {"Ca": -1, "Na": 2}, False, 1.0))
        reactions.append(("MgNa_exch", {"Mg": 1, "Na": -2}, False, 1.0))
        reactions.append(("NaMg_exch", {"Mg": -1, "Na": 2}, False, 1.0))

    labels = [label for label, _, _, _ in reactions]
    mineral_mask = [is_mineral for _, _, is_mineral, _ in reactions]
    penalty_scales = [scale for _, _, _, scale in reactions]
    matrix = [_vector_from_coeffs(coeffs, ion_order) for _, coeffs, _, _ in reactions]
    return matrix, labels, mineral_mask, penalty_scales


def _dot(a: Iterable[float], b: Iterable[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _combine_reactions(
    matrix: Sequence[Sequence[float]], weights: Sequence[float]
) -> List[float]:
    if not matrix:
        return []
    result = [0.0] * len(matrix[0])
    for reaction, weight in zip(matrix, weights):
        for idx, value in enumerate(reaction):
            result[idx] += value * weight
    return result


def _weighted_vectors(
    reaction_matrix: Sequence[Sequence[float]],
    residual: Sequence[float],
    weights: Sequence[float],
) -> Tuple[List[List[float]], List[float]]:
    if not weights:
        weighted_matrix = [list(row) for row in reaction_matrix]
        weighted_residual = list(residual)
    else:
        weighted_matrix = [
            [value * (weight**0.5) for value, weight in zip(row, weights)]
            for row in reaction_matrix
        ]
        weighted_residual = [
            value * (weight**0.5) for value, weight in zip(residual, weights)
        ]
    return weighted_matrix, weighted_residual


def _weighted_norm_sq(values: Sequence[float], weights: Sequence[float]) -> float:
    if not weights:
        return sum(v * v for v in values)
    return sum(w * v * v for v, w in zip(values, weights))


@dataclass
class ReactionFit:
    extents: List[float]
    residual: List[float]
    residual_norm: float
    l1_norm: float
    iterations: int
    converged: bool
    extents_std: Optional[List[float]] = None
    extents_ci_low: Optional[List[float]] = None
    extents_ci_high: Optional[List[float]] = None
    uncertainty_result: Optional[object] = None

    @property
    def aic(self) -> float:
        """Akaike Information Criterion: AIC = 2k + n*ln(RSS/n)"""
        n = len(self.residual)
        k = len([e for e in self.extents if abs(e) > 1e-6]) + 1 # +1 for variance
        if n == 0 or self.residual_norm <= 0: return float('inf')
        return 2*k + n * math.log(self.residual_norm / n)

    @property
    def aicc(self) -> float:
        """Corrected AIC for small sample sizes: AICc = AIC + 2k(k+1)/(n-k-1)"""
        n = len(self.residual)
        k = len([e for e in self.extents if abs(e) > 1e-6]) + 1
        aic_val = self.aic
        if n - k - 1 <= 0: return float('inf') # Avoid singularity
        return aic_val + (2 * k * (k + 1)) / (n - k - 1)

    @property
    def bic(self) -> float:
        """Bayesian Information Criterion: BIC = ln(n)k + n*ln(RSS/n)"""
        n = len(self.residual)
        k = len([e for e in self.extents if abs(e) > 1e-6]) + 1
        if n == 0 or self.residual_norm <= 0: return float('inf')
        return k * math.log(n) + n * math.log(self.residual_norm / n)


def _soft_threshold(value: float, threshold: float) -> float:
    if value > threshold: return value - threshold
    if value < -threshold: return value + threshold
    return 0.0


def fit_reactions(
    residual: List[float],
    reaction_matrix: List[List[float]],
    weights: List[float],
    lambda_l1: float,
    max_iter: int = 300,
    tol: float = 1e-6,
    signed_mask: Optional[List[bool]] = None,
    lb: Optional[List[float]] = None,
    ub: Optional[List[float]] = None,
    penalty_scales: Optional[List[float]] = None,
    lambda_l2: float = 0.0,
) -> ReactionFit:
    if not reaction_matrix:
        return ReactionFit([], residual, _weighted_norm_sq(residual, weights), 0.0, 0, True)

    weighted_matrix, weighted_residual = _weighted_vectors(reaction_matrix, residual, weights)
    m = len(weighted_matrix)

    gram = [[_dot(weighted_matrix[i], weighted_matrix[j]) for j in range(m)] for i in range(m)]
    s_r = [_dot(weighted_matrix[i], weighted_residual) for i in range(m)]

    if signed_mask is None: signed_mask = [False] * m
    if penalty_scales is None: penalty_scales = [1.0] * m

    max_diag = max(gram[i][i] for i in range(m)) if m > 0 else 0.0
    ridge_epsilon = max_diag * 1e-10 + 1e-20

    z = [0.0] * m
    converged = False

    for iteration in range(1, max_iter + 1):
        max_delta = 0.0
        for j in range(m):
            dot_prod = sum(gram[j][k] * z[k] for k in range(m) if k != j)
            rho = s_r[j] - dot_prod

            # Elastic Net Step: divisor includes L2 penalty (lambda_l2)
            # Objective: 0.5*RSS + lambda_l1*|z| + 0.5*lambda_l2*z^2
            denom = gram[j][j] + lambda_l2 + ridge_epsilon

            # Apply per-reaction penalty scaling
            scaled_lambda = lambda_l1 * penalty_scales[j]
            updated = _soft_threshold(rho, scaled_lambda / 2.0) / denom if denom > 1e-15 else 0.0
            if not signed_mask[j]: updated = max(0.0, updated)
            if lb and lb[j] is not None: updated = max(lb[j], updated)
            if ub and ub[j] is not None: updated = min(ub[j], updated)

            delta = abs(updated - z[j])
            if delta > max_delta: max_delta = delta
            z[j] = updated

        if max_delta <= tol:
            converged = True
            break

    fitted = _combine_reactions(reaction_matrix, z)
    post_residual = [r - f for r, f in zip(residual, fitted)]
    l1_norm = sum(abs(v) * p for v, p in zip(z, penalty_scales))
    return ReactionFit(z, post_residual, _weighted_norm_sq(post_residual, weights), l1_norm, iteration, converged)
