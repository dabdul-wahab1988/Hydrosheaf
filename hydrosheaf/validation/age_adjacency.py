"""Age-based evidence for temporal direction and direct adjacency.

Age ordering and graph adjacency are different estimands.  If the age at a
downstream node is greater than the age at an upstream node, the observation
supports a *forward temporal direction*.  It does not, by itself, establish
that the two nodes are directly adjacent: an indirect path can have the same
ordering and a larger cumulative travel time.

This module keeps those questions separate.  :func:`compute_direction_evidence`
calculates the one-sided probability and cost for downstream water being
older.  :func:`compute_age_adjacency_evidence` additionally scores explicitly
provided direct and indirect travel-time hypotheses with normal likelihoods.
No distance, velocity, path, or intermediate-node values are inferred here;
callers must provide those quantities if they want a direct-versus-indirect
comparison.

The likelihood model is deliberately small and auditable.  It treats the
reported upstream and downstream ages as independent normal observations and
adds an independent normal process/travel-time error.  Thus the uncertainty
of the observed increment is

``sqrt(upstream_sigma_years**2 + downstream_sigma_years**2 +
process_sigma_years**2 - 2*age_covariance_years2)``.

This is an uncertainty-propagation model, not a claim that groundwater age is
itself a direct travel-time measurement.  The returned Bayes factor is valid
only for the supplied travel-time hypotheses and their shared error model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Literal


DirectionStatus = Literal["forward_compatible", "reverse_incompatible", "indeterminate"]
AdjacencyStatus = Literal["scored", "ambiguous", "insufficient_information"]

_MIN_PROBABILITY = 1.0e-12
_SQRT_TWO = math.sqrt(2.0)
_LOG_SQRT_TWO_PI = 0.5 * math.log(2.0 * math.pi)


def _finite(name: str, value: Any) -> float:
    """Convert a scalar to a finite float, rejecting booleans and missing data."""

    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool.")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a finite real number.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _nonnegative(name: str, value: Any) -> float:
    result = _finite(name, value)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative.")
    return result


def _propagated_sigma(
    upstream_sigma_years: float,
    downstream_sigma_years: float,
    process_sigma_years: float,
    age_covariance_years2: float = 0.0,
) -> float:
    """Return endpoint/process uncertainty with an optional endpoint covariance."""

    variance = (
        upstream_sigma_years**2
        + downstream_sigma_years**2
        + process_sigma_years**2
        - 2.0 * age_covariance_years2
    )
    tolerance = 1.0e-12 * max(
        1.0,
        upstream_sigma_years**2,
        downstream_sigma_years**2,
        process_sigma_years**2,
    )
    if variance < -tolerance:
        raise ValueError("age covariance produces a negative variance.")
    result = math.sqrt(max(0.0, variance))
    if not math.isfinite(result):
        # The inputs are finite, but a pathological combination can still
        # overflow on some Python/libm implementations.  Failing explicitly
        # keeps an invalid uncertainty from silently becoming a score.
        raise ValueError("propagated uncertainty is not finite.")
    return result


def _normal_log_likelihood(observed: float, expected: float, sigma: float) -> float:
    """Evaluate a normal log density with overflow-safe tails.

    A very distant finite observation has a log likelihood of negative
    infinity in floating-point arithmetic.  Returning ``-inf`` is preferable
    to overflowing a squared z score or returning NaN; the caller handles the
    resulting comparison explicitly.
    """

    residual = observed - expected
    if not math.isfinite(residual):
        return -math.inf
    z = residual / sigma
    if not math.isfinite(z):
        return -math.inf
    z_squared = z * z
    if not math.isfinite(z_squared):
        return -math.inf
    return -0.5 * z_squared - math.log(sigma) - _LOG_SQRT_TWO_PI


def _direction_probability(increment_years: float, sigma_years: float) -> tuple[float, float | None]:
    """Return P(increment >= 0) and its z score.

    With zero age uncertainty this is the deterministic limit.  The z score
    is ``None`` in that limit so serialised results do not contain artificial
    infinities.
    """

    if sigma_years == 0.0:
        if increment_years > 0.0:
            return 1.0, None
        if increment_years < 0.0:
            return 0.0, None
        return 0.5, None
    z = increment_years / sigma_years
    # erfc is more accurate than 1 + erf in the negative tail.
    probability = 0.5 * math.erfc(-z / _SQRT_TWO)
    return min(1.0, max(0.0, probability)), z


@dataclass(frozen=True)
class DirectionEvidence:
    """One-sided age-ordering evidence for a directed candidate edge.

    ``forward_compatibility_probability`` means the probability that the
    downstream-minus-upstream age increment is non-negative under the stated
    endpoint-error model.  It is a temporal direction statistic, not a
    probability of direct graph adjacency or a probability of any graph path.
    """

    age_increment_years: float
    age_sigma_years: float
    direction_z: float | None
    forward_compatibility_probability: float
    forward_compatibility_cost: float
    status: DirectionStatus
    flags: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly representation for audit tables."""

        result = asdict(self)
        result["flags"] = list(self.flags)
        return result


@dataclass(frozen=True)
class AgeAdjacencyEvidence:
    """Auditable temporal and direct-versus-indirect age evidence.

    ``adjacency_status`` describes information availability, not the truth of
    adjacency:

    * ``scored``: both explicit travel-time hypotheses were scored with
      positive hypothesis-specific uncertainty and a finite log Bayes factor
      is available;
    * ``ambiguous``: only a direct travel-time hypothesis was supplied, so a
      direct score exists but cannot be compared with an indirect hypothesis;
    * ``insufficient_information``: the direct hypothesis or a positive
      hypothesis-specific uncertainty is missing.

    Positive ``log_bayes_factor_direct_vs_indirect`` favours the direct
    hypothesis for the supplied model; negative values favour the indirect
    hypothesis.  This result never labels an edge as direct or indirect.
    """

    upstream_age_years: float
    downstream_age_years: float
    upstream_sigma_years: float
    downstream_sigma_years: float
    age_covariance_years2: float
    process_sigma_years: float
    observed_increment_years: float
    propagated_sigma_years: float
    direction: DirectionEvidence
    direct_travel_years: float | None
    indirect_travel_years: float | None
    direct_sigma_years: float | None
    indirect_sigma_years: float | None
    direct_log_likelihood: float | None
    indirect_log_likelihood: float | None
    direct_standardized_residual: float | None
    indirect_standardized_residual: float | None
    log_bayes_factor_direct_vs_indirect: float | None
    comparison_discriminating: bool | None
    adjacency_status: AdjacencyStatus
    flags: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a nested JSON-friendly representation for audit tables."""

        result = asdict(self)
        result["direction"] = self.direction.to_dict()
        result["flags"] = list(self.flags)
        return result


def compute_direction_evidence(
    upstream_age_years: float,
    downstream_age_years: float,
    *,
    upstream_sigma_years: float = 0.0,
    downstream_sigma_years: float = 0.0,
    age_covariance_years2: float = 0.0,
) -> DirectionEvidence:
    """Compute temporal direction compatibility for a directed edge.

    Ages are interpreted in years, and downstream water is considered
    forward-compatible when its age is greater than or equal to the upstream
    age.  Upstream and downstream uncertainties are propagated in quadrature.
    No travel time, distance, velocity, or graph labels are used.
    """

    upstream = _finite("upstream_age_years", upstream_age_years)
    downstream = _finite("downstream_age_years", downstream_age_years)
    upstream_sigma = _nonnegative("upstream_sigma_years", upstream_sigma_years)
    downstream_sigma = _nonnegative("downstream_sigma_years", downstream_sigma_years)
    covariance = _finite("age_covariance_years2", age_covariance_years2)
    if abs(covariance) > upstream_sigma * downstream_sigma + 1.0e-12:
        raise ValueError(
            "age_covariance_years2 must satisfy |covariance| <= "
            "upstream_sigma_years * downstream_sigma_years."
        )
    age_increment = downstream - upstream
    if not math.isfinite(age_increment):
        raise ValueError("age increment is not finite.")
    age_sigma = _propagated_sigma(
        upstream_sigma,
        downstream_sigma,
        0.0,
        age_covariance_years2=covariance,
    )
    probability, z_score = _direction_probability(age_increment, age_sigma)
    cost = -math.log(max(_MIN_PROBABILITY, probability))
    if age_increment > 0.0:
        status: DirectionStatus = "forward_compatible"
    elif age_increment < 0.0:
        status = "reverse_incompatible"
    else:
        status = "indeterminate"
    flags: list[str] = []
    if age_sigma == 0.0:
        flags.append("zero_age_uncertainty")
    if probability <= _MIN_PROBABILITY:
        flags.append("strong_reverse_age_order")
    elif probability >= 1.0 - _MIN_PROBABILITY:
        flags.append("strong_forward_age_order")
    return DirectionEvidence(
        age_increment_years=float(age_increment),
        age_sigma_years=float(age_sigma),
        direction_z=z_score,
        forward_compatibility_probability=float(probability),
        forward_compatibility_cost=float(cost),
        status=status,
        flags=tuple(flags),
    )


def compute_age_adjacency_evidence(
    upstream_age_years: float,
    downstream_age_years: float,
    *,
    upstream_sigma_years: float = 0.0,
    downstream_sigma_years: float = 0.0,
    age_covariance_years2: float = 0.0,
    process_sigma_years: float = 0.0,
    direct_travel_years: float | None = None,
    indirect_travel_years: float | None = None,
    direct_travel_sigma_years: float | None = None,
    indirect_travel_sigma_years: float | None = None,
) -> AgeAdjacencyEvidence:
    """Score age direction and explicitly supplied travel-time hypotheses.

    ``direct_travel_years`` and ``indirect_travel_years`` are expected travel
    times supplied by the caller, for example from an independently specified
    path/velocity model.  They are not estimated from the age difference.
    Each supplied hypothesis receives a normal log likelihood of the observed
    age increment.  The hypothesis-specific likelihood standard deviation is
    the quadrature propagation of endpoint/process uncertainty and the
    corresponding travel-time dispersion.  ``age_covariance_years2`` is the
    covariance of the endpoint age errors and is subtracted twice because the
    increment is downstream minus upstream.

    A direct-versus-indirect log Bayes factor is returned only when both
    hypotheses are supplied and both hypothesis-specific uncertainties are
    strictly positive.  A direct hypothesis alone gives ``ambiguous`` status;
    without a direct hypothesis (or without a usable uncertainty) the result
    has ``insufficient_information`` status.  The function raises for
    non-finite values and negative uncertainties/travel times rather than
    silently converting them into evidence.
    """

    upstream = _finite("upstream_age_years", upstream_age_years)
    downstream = _finite("downstream_age_years", downstream_age_years)
    upstream_sigma = _nonnegative("upstream_sigma_years", upstream_sigma_years)
    downstream_sigma = _nonnegative("downstream_sigma_years", downstream_sigma_years)
    covariance = _finite("age_covariance_years2", age_covariance_years2)
    if abs(covariance) > upstream_sigma * downstream_sigma + 1.0e-12:
        raise ValueError(
            "age_covariance_years2 must satisfy |covariance| <= "
            "upstream_sigma_years * downstream_sigma_years."
        )
    process_sigma = _nonnegative("process_sigma_years", process_sigma_years)
    direct = (
        None
        if direct_travel_years is None
        else _nonnegative("direct_travel_years", direct_travel_years)
    )
    indirect = (
        None
        if indirect_travel_years is None
        else _nonnegative("indirect_travel_years", indirect_travel_years)
    )
    direct_sigma = (
        0.0
        if direct_travel_sigma_years is None
        else _nonnegative("direct_travel_sigma_years", direct_travel_sigma_years)
    )
    indirect_sigma = (
        0.0
        if indirect_travel_sigma_years is None
        else _nonnegative("indirect_travel_sigma_years", indirect_travel_sigma_years)
    )
    observed_increment = downstream - upstream
    if not math.isfinite(observed_increment):
        raise ValueError("observed age increment is not finite.")

    direction = compute_direction_evidence(
        upstream,
        downstream,
        upstream_sigma_years=upstream_sigma,
        downstream_sigma_years=downstream_sigma,
        age_covariance_years2=covariance,
    )
    propagated_sigma = _propagated_sigma(
        upstream_sigma,
        downstream_sigma,
        process_sigma,
        age_covariance_years2=covariance,
    )

    direct_log_likelihood: float | None = None
    indirect_log_likelihood: float | None = None
    direct_z: float | None = None
    indirect_z: float | None = None
    log_bayes_factor: float | None = None
    comparison_discriminating: bool | None = None
    flags: list[str] = []

    if direct is None:
        flags.append("direct_travel_time_missing")
    if indirect is None:
        flags.append("indirect_travel_time_missing")
    direct_total_sigma = (
        None
        if direct is None
        else math.hypot(propagated_sigma, direct_sigma)
    )
    indirect_total_sigma = (
        None
        if indirect is None
        else math.hypot(propagated_sigma, indirect_sigma)
    )
    if propagated_sigma == 0.0:
        flags.append("zero_propagated_uncertainty")
    if direct is not None and direct_total_sigma and direct_total_sigma > 0.0:
        direct_residual = observed_increment - direct
        direct_z = direct_residual / direct_total_sigma
        direct_log_likelihood = _normal_log_likelihood(
            observed_increment,
            direct,
            direct_total_sigma,
        )
    elif direct is not None:
        flags.append("zero_direct_uncertainty")
    if indirect is not None and indirect_total_sigma and indirect_total_sigma > 0.0:
        indirect_residual = observed_increment - indirect
        indirect_z = indirect_residual / indirect_total_sigma
        indirect_log_likelihood = _normal_log_likelihood(
            observed_increment,
            indirect,
            indirect_total_sigma,
        )
    elif indirect is not None:
        flags.append("zero_indirect_uncertainty")

    if direct is not None and indirect is not None:
        same_hypothesis = math.isclose(
            direct,
            indirect,
            rel_tol=0.0,
            abs_tol=0.0,
        ) and math.isclose(
            direct_sigma,
            indirect_sigma,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        if same_hypothesis:
            comparison_discriminating = False
            flags.append("equal_travel_time_hypotheses")
        elif (
            direct_log_likelihood is not None
            and indirect_log_likelihood is not None
            and not (
                math.isinf(direct_log_likelihood)
                and direct_log_likelihood < 0.0
            )
            and not (
                math.isinf(indirect_log_likelihood)
                and indirect_log_likelihood < 0.0
            )
        ):
            comparison_discriminating = True
            log_bayes_factor = direct_log_likelihood - indirect_log_likelihood
            if not math.isfinite(log_bayes_factor):
                log_bayes_factor = None
                comparison_discriminating = False
                flags.append("likelihood_comparison_undefined")
        else:
            comparison_discriminating = False
            flags.append("likelihood_comparison_undefined")

    if direct is None or direct_log_likelihood is None:
        adjacency_status: AdjacencyStatus = "insufficient_information"
    elif indirect is None:
        adjacency_status = "ambiguous"
    elif indirect_log_likelihood is None or log_bayes_factor is None:
        adjacency_status = "insufficient_information"
    else:
        adjacency_status = "scored"

    return AgeAdjacencyEvidence(
        upstream_age_years=float(upstream),
        downstream_age_years=float(downstream),
        upstream_sigma_years=float(upstream_sigma),
        downstream_sigma_years=float(downstream_sigma),
        age_covariance_years2=float(covariance),
        process_sigma_years=float(process_sigma),
        observed_increment_years=float(observed_increment),
        propagated_sigma_years=float(propagated_sigma),
        direction=direction,
        direct_travel_years=direct,
        indirect_travel_years=indirect,
        direct_sigma_years=(
            None if direct is None else float(direct_sigma)
        ),
        indirect_sigma_years=(
            None if indirect is None else float(indirect_sigma)
        ),
        direct_log_likelihood=direct_log_likelihood,
        indirect_log_likelihood=indirect_log_likelihood,
        direct_standardized_residual=direct_z,
        indirect_standardized_residual=indirect_z,
        log_bayes_factor_direct_vs_indirect=log_bayes_factor,
        comparison_discriminating=comparison_discriminating,
        adjacency_status=adjacency_status,
        flags=tuple(flags),
    )


__all__ = [
    "AdjacencyStatus",
    "AgeAdjacencyEvidence",
    "DirectionEvidence",
    "DirectionStatus",
    "compute_age_adjacency_evidence",
    "compute_direction_evidence",
]
