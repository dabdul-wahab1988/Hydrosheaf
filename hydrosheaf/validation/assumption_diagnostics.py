"""Assumption-audit diagnostics for the five assumptions that had no test.

``docs/hydrosheaf_model_assumptions.md`` and the technical manual's assumption
catalogue enumerate the conditions a HydroSheaf inference relies on. Most of
them are instrumented: a violation surfaces as a named field, a gate code or an
abstention. Five were not, and this module supplies a diagnostic for each:

* ``A4``  missing values are missing at random with respect to the quantity
  inferred;
* ``A7``  elevation and head share a common vertical datum;
* ``A10`` laboratory bias is negligible relative to the effects inferred;
* ``A34`` the tracer is conservative apart from modelled decay;
* ``A37`` observation errors are independent across species.

Design rules
------------
Every diagnostic in this module obeys four rules.

1. **Detection, never confirmation.** A ``HOLDS`` result means *no violation was
   detected*, not that the assumption is true. Each report states the sample
   size and the power limitation. ``UNDETERMINED`` is returned whenever the
   inputs cannot support a test, so that absent evidence cannot be read as
   compliance. A module that reported ``HOLDS`` for an empty input would be
   worse than useless.
2. **No imputation.** Rows lacking a required field are excluded and counted.
   No diagnostic substitutes a default value, in keeping with the framework's
   data contract.
3. **Deterministic.** The only stochastic component is the permutation null in
   ``A37``, driven by ``numpy.random.default_rng`` with a declared default seed
   (``DEFAULT_SEED = 1729`` from :mod:`hydrosheaf.reproducibility.core`).
4. **Unit-free where possible.** ``A34`` and ``A37`` use log-ratios and rank
   statistics, which are invariant to a common scale factor. ``A7`` is
   deliberately *not* unit-free, because a unit mismatch is one of the things it
   is designed to catch.

The reports are written to be appended to a run manifest or serialised into a
manuscript supplement; :meth:`AssumptionReport.to_dict` and
:meth:`AssumptionAudit.to_json` return stable, JSON-compatible payloads.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..data.parsing import sample_list
from ..data.qc import charge_balance_ratio
from ..models.redox import classify_redox

__all__ = [
    "ASSUMPTION_CATALOGUE",
    "DEFAULT_SEED",
    "AssumptionAudit",
    "AssumptionReport",
    "AssumptionStatus",
    "assess_model_assumptions",
    "diagnose_error_independence",
    "diagnose_laboratory_bias",
    "diagnose_missingness_at_random",
    "diagnose_tracer_conservativeness",
    "diagnose_vertical_datum",
]

# --------------------------------------------------------------------------
# Vocabulary
# --------------------------------------------------------------------------

DEFAULT_SEED = 1729
"""Default permutation seed, matching :data:`hydrosheaf.reproducibility.core.DEFAULT_SEED`."""


class AssumptionStatus(str, Enum):
    """Outcome of a single assumption diagnostic.

    ``HOLDS`` is deliberately worded as an absence of detected violation rather
    than as confirmation of the assumption.
    """

    HOLDS = "HOLDS"
    """No violation was detected at the declared significance and effect size."""

    VIOLATED = "VIOLATED"
    """A violation was detected: the data are inconsistent with the assumption."""

    UNDETERMINED = "UNDETERMINED"
    """The inputs cannot support the test (too few rows, no usable fields)."""

    NOT_APPLICABLE = "NOT_APPLICABLE"
    """The assumption is moot for these data (for example, nothing is missing)."""


ASSUMPTION_CATALOGUE: Dict[str, str] = {
    "A4": (
        "Missing values are missing at random with respect to the quantity "
        "inferred."
    ),
    "A7": "Elevation and head share a common vertical datum.",
    "A10": "Laboratory bias is negligible relative to the effects inferred.",
    "A34": "The tracer is conservative apart from modelled decay.",
    "A37": "Observation errors are independent across species.",
}
"""Statements of the five assumptions, keyed as in the manual's Appendix I."""


# --------------------------------------------------------------------------
# Report types
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class AssumptionReport:
    """Auditable outcome of one assumption diagnostic."""

    assumption_id: str
    """Identifier matching :data:`ASSUMPTION_CATALOGUE` (``'A4'`` etc.)."""

    statement: str
    """The assumption that was tested, verbatim from the catalogue."""

    status: AssumptionStatus
    """Outcome; see :class:`AssumptionStatus`."""

    metrics: Mapping[str, float] = field(default_factory=dict)
    """Numeric diagnostics supporting the outcome."""

    findings: Tuple[str, ...] = ()
    """Machine-readable finding codes, e.g. ``'missingness_covariate_association'``."""

    reasons: Tuple[str, ...] = ()
    """Human-readable explanations, one per finding."""

    n_used: int = 0
    """Rows that contributed to the test."""

    n_excluded: int = 0
    """Rows excluded for want of a required field; never imputed."""

    def __post_init__(self) -> None:
        if self.assumption_id not in ASSUMPTION_CATALOGUE:
            raise ValueError(
                f"Unknown assumption id {self.assumption_id!r}; "
                f"expected one of {sorted(ASSUMPTION_CATALOGUE)}."
            )

    @property
    def violated(self) -> bool:
        return self.status is AssumptionStatus.VIOLATED

    def to_dict(self) -> Dict[str, Any]:
        return {
            "assumption_id": self.assumption_id,
            "statement": self.statement,
            "status": self.status.value,
            "metrics": dict(self.metrics),
            "findings": list(self.findings),
            "reasons": list(self.reasons),
            "n_used": self.n_used,
            "n_excluded": self.n_excluded,
        }


@dataclass(frozen=True)
class AssumptionAudit:
    """Aggregate audit over a set of assumption diagnostics."""

    reports: Tuple[AssumptionReport, ...]
    counts: Mapping[str, int] = field(default_factory=dict)
    status: AssumptionStatus = AssumptionStatus.HOLDS
    seed: int = DEFAULT_SEED

    @property
    def by_id(self) -> Dict[str, AssumptionReport]:
        return {r.assumption_id: r for r in self.reports}

    @property
    def violations(self) -> Tuple[AssumptionReport, ...]:
        return tuple(r for r in self.reports if r.violated)

    @property
    def undetermined(self) -> Tuple[AssumptionReport, ...]:
        return tuple(
            r for r in self.reports if r.status is AssumptionStatus.UNDETERMINED
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status.value,
            "counts": dict(self.counts),
            "seed": self.seed,
            "reports": [r.to_dict() for r in self.reports],
        }

    def to_json(self) -> str:
        """Deterministic JSON, suitable for hashing or a run manifest."""
        return json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))


# --------------------------------------------------------------------------
# Small statistics helpers (numpy only; no SciPy dependency for the audit)
# --------------------------------------------------------------------------


def _finite(value: object) -> Optional[float]:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _two_sided_normal_p(z: float) -> float:
    """Two-sided p-value for a standard normal deviate."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def _ranks(values: np.ndarray) -> np.ndarray:
    """Average ranks, handling ties by midpoint ranking."""
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    i = 0
    while i < len(sorted_values):
        j = i
        while j + 1 < len(sorted_values) and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def _rank_sum_test(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    """Wilcoxon rank-sum (Mann-Whitney U) test with tie correction.

    Returns ``(z, p_value, rank_biserial)``. The normal approximation is used,
    which is adequate for ``n1 * n2 >= 8``; callers must check the sample sizes
    and report ``UNDETERMINED`` below that. ``rank_biserial`` is the effect
    size in ``[-1, 1]``: the probability that a random draw from ``a`` exceeds
    a random draw from ``b``, rescaled.
    """
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0, 0.0
    pooled = np.concatenate([np.asarray(a, dtype=float), np.asarray(b, dtype=float)])
    ranks = _ranks(pooled)
    r1 = float(np.sum(ranks[:n1]))

    # Tie correction for the variance of the rank sum.
    _, counts = np.unique(pooled, return_counts=True)
    tie_term = float(np.sum(counts**3 - counts))
    n = n1 + n2
    variance = (n1 * n2 / 12.0) * ((n + 1.0) - tie_term / (n * (n - 1.0))) if n > 1 else 0.0
    if variance <= 0.0:
        return 0.0, 1.0, 0.0

    u1 = r1 - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    # Continuity correction toward the mean.
    correction = 0.5 if u1 > mean_u else (-0.5 if u1 < mean_u else 0.0)
    z = (u1 - mean_u - correction) / math.sqrt(variance)
    # Rank-biserial correlation: 2*U1/(n1*n2) - 1.
    effect = 2.0 * u1 / (n1 * n2) - 1.0
    return float(z), float(_two_sided_normal_p(z)), float(effect)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float(np.dot(xc, xc)) * float(np.dot(yc, yc)))
    if denom <= 0.0:
        return float("nan")
    return float(np.dot(xc, yc) / denom)


def _binomial_tail_at_least(k: int, n: int, p: float) -> float:
    """Exact upper-tail probability P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    total = 0.0
    for i in range(k, n + 1):
        total += math.comb(n, i) * (p**i) * ((1.0 - p) ** (n - i))
    return min(1.0, total)


def _quantile(values: np.ndarray, q: float) -> float:
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, q, method="linear"))


def _iqr(values: np.ndarray) -> float:
    if len(values) == 0:
        return float("nan")
    return _quantile(values, 0.75) - _quantile(values, 0.25)


def _coerce_float(sample: Mapping[str, object], key: str) -> Optional[float]:
    if key not in sample:
        return None
    return _finite(sample.get(key))


# --------------------------------------------------------------------------
# A4 -- missingness at random
# --------------------------------------------------------------------------

_DEFAULT_A4_COVARIATES: Tuple[str, ...] = (
    "head",
    "head_meas",
    "hydraulic_head",
    "elevation",
    "elevation_m",
    "screen_depth",
    "well_depth",
    "latitude",
    "longitude",
    "depth_m",
)


def diagnose_missingness_at_random(
    samples: object,
    *,
    channels: Optional[Iterable[str]] = None,
    covariates: Optional[Iterable[str]] = None,
    alpha: float = 0.01,
    effect_threshold: float = 0.2,
    cluster_threshold: float = 0.6,
    min_pairs: int = 8,
    min_rows: int = 4,
) -> AssumptionReport:
    """Test whether missingness is associated with the sampled quantity.

    Two sub-tests are run.

    **Covariate association.** For every (channel, covariate) pair with enough
    rows in both groups, the covariate distribution is compared between rows
    where the channel is present and rows where it is absent. Two rank-based
    tests are run, because a single location test is blind to the most common
    real pattern:

    * a **location** test (Wilcoxon rank-sum) catches missingness that tracks
      the level of a covariate -- a channel analysed only in deep wells;
    * a **dispersion** test (the same rank-sum applied to absolute deviations
      from the pooled median, the Ansari--Bradley construction) catches
      missingness that tracks the *spread* -- a channel analysed only in a
      middle depth band, which leaves the medians equal and defeats a location
      test entirely.

    A significant association with a non-trivial effect size on either test is a
    violation. A covariate with no variation at all is skipped rather than
    counted, so the reported number of pairs tested is the number that could
    actually discriminate.

    **Missingness clustering.** Pairwise association between the missingness
    indicators is summarised by Cramér's V. Strong clustering means whole panels
    were analysed or omitted together, which is structured (MAR or MNAR)
    missingness rather than MCAR.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    channels : iterable of str, optional
        Channels to test for missingness. Defaults to every key observed in at
        least one row that is absent in at least one row.
    covariates : iterable of str, optional
        Quantities the missingness must be independent of. Defaults to
        :data:`_DEFAULT_A4_COVARIATES`, restricted to keys present in the data.
    alpha : float
        Significance level, Bonferroni-corrected by the number of pairs tested.
    effect_threshold : float
        Minimum absolute rank-biserial effect (``|2U/(n1 n2) - 1|``) before a
        significant association is called a violation.
    cluster_threshold : float
        Cramér's V above which missingness clustering is called a violation.
    min_pairs, min_rows : int
        Sample-size floors; below them the test is ``UNDETERMINED`` rather than
        ``HOLDS``.

    Returns
    -------
    AssumptionReport
        With ``assumption_id='A4'``.

    Notes
    -----
    Failing to detect an association does not establish MCAR: the test has power
    only against associations the available covariates can express. The report
    records the number of pairs actually tested so the reader can judge coverage.
    """
    rows = sample_list(samples)
    statement = ASSUMPTION_CATALOGUE["A4"]

    if len(rows) < min_rows:
        return AssumptionReport(
            assumption_id="A4",
            statement=statement,
            status=AssumptionStatus.UNDETERMINED,
            reasons=(f"only {len(rows)} rows available; need >= {min_rows}",),
            n_used=len(rows),
        )

    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    # A channel is testable when it is present in some rows and absent in others.
    if channels is None:
        candidates = [
            key
            for key in keys
            if any(key in row for row in rows) and any(key not in row for row in rows)
        ]
    else:
        candidates = [c for c in channels if c in keys]

    if not candidates:
        return AssumptionReport(
            assumption_id="A4",
            statement=statement,
            status=AssumptionStatus.NOT_APPLICABLE,
            reasons=("no channel is partially observed; the assumption is moot",),
            n_used=len(rows),
            metrics={"n_rows": float(len(rows)), "n_channels_with_missingness": 0.0},
        )

    present_covariates = [
        c
        for c in (covariates if covariates is not None else _DEFAULT_A4_COVARIATES)
        if c in keys
    ]

    findings: List[str] = []
    reasons: List[str] = []
    metrics: Dict[str, float] = {
        "n_rows": float(len(rows)),
        "n_channels_with_missingness": float(len(candidates)),
        "n_covariates": float(len(present_covariates)),
    }

    # --- covariate association -------------------------------------------
    # Two tests per pair: location and dispersion. A middle-band missingness
    # rule leaves the medians equal, so a location test alone would miss it.
    n_tests = max(1, 2 * len(candidates) * len(present_covariates))
    corrected_alpha = alpha / n_tests
    tested_pairs = 0
    scale_pairs = 0
    degenerate_pairs = 0
    strongest_effect = 0.0
    strongest_scale_effect = 0.0
    smallest_p = 1.0
    smallest_scale_p = 1.0
    for channel in candidates:
        for covariate in present_covariates:
            present: List[float] = []
            absent: List[float] = []
            for row in rows:
                value = _coerce_float(row, covariate)
                if value is None:
                    continue
                (present if channel in row else absent).append(value)
            if len(present) < min_pairs or len(absent) < min_pairs:
                continue
            pooled = np.asarray(present + absent, dtype=float)
            if float(np.ptp(pooled)) == 0.0:
                # A constant covariate cannot support any rank-based test.
                degenerate_pairs += 1
                continue
            tested_pairs += 1

            _, p_value, effect = _rank_sum_test(present, absent)
            smallest_p = min(smallest_p, p_value)
            if abs(effect) > abs(strongest_effect):
                strongest_effect = effect
            if p_value < corrected_alpha and abs(effect) >= effect_threshold:
                findings.append("missingness_covariate_association")
                reasons.append(
                    f"missingness of '{channel}' is associated with the level of "
                    f"'{covariate}' (p={p_value:.3g} < {corrected_alpha:.3g}, "
                    f"rank-biserial={effect:+.3f})"
                )

            # Ansari-Bradley style dispersion test on absolute deviations.
            centre = float(np.median(pooled))
            dev_present = np.abs(np.asarray(present, dtype=float) - centre)
            dev_absent = np.abs(np.asarray(absent, dtype=float) - centre)
            if float(np.ptp(np.concatenate([dev_present, dev_absent]))) > 0.0:
                scale_pairs += 1
                _, p_scale, eff_scale = _rank_sum_test(dev_present, dev_absent)
                smallest_scale_p = min(smallest_scale_p, p_scale)
                if abs(eff_scale) > abs(strongest_scale_effect):
                    strongest_scale_effect = eff_scale
                if p_scale < corrected_alpha and abs(eff_scale) >= effect_threshold:
                    findings.append("missingness_covariate_dispersion")
                    reasons.append(
                        f"missingness of '{channel}' is associated with the spread "
                        f"of '{covariate}' (dispersion test p={p_scale:.3g} < "
                        f"{corrected_alpha:.3g}, rank-biserial={eff_scale:+.3f}); "
                        f"the covariate medians agree, so a location test alone "
                        f"would have missed this"
                    )

    metrics["n_pairs_tested"] = float(tested_pairs)
    metrics["n_dispersion_pairs_tested"] = float(scale_pairs)
    metrics["n_degenerate_pairs_skipped"] = float(degenerate_pairs)
    metrics["smallest_p_value"] = float(smallest_p)
    metrics["smallest_dispersion_p_value"] = float(smallest_scale_p)
    metrics["strongest_rank_biserial"] = float(strongest_effect)
    metrics["strongest_dispersion_rank_biserial"] = float(strongest_scale_effect)
    metrics["bonferroni_alpha"] = float(corrected_alpha)

    # --- missingness clustering -------------------------------------------
    indicator_matrix = np.zeros((len(rows), len(candidates)), dtype=float)
    for j, channel in enumerate(candidates):
        indicator_matrix[:, j] = [1.0 if channel in row else 0.0 for row in rows]

    max_v = 0.0
    n_pairs_compared = 0
    if len(candidates) >= 2:
        for a in range(len(candidates)):
            for b in range(a + 1, len(candidates)):
                va, vb = indicator_matrix[:, a], indicator_matrix[:, b]
                if va.std() == 0.0 or vb.std() == 0.0:
                    continue
                n_pairs_compared += 1
                # Cramer's V for two binary variables equals |phi|.
                phi = float(np.corrcoef(va, vb)[0, 1])
                if not math.isfinite(phi):
                    continue
                max_v = max(max_v, abs(phi))
                if abs(phi) >= cluster_threshold:
                    findings.append("missingness_clustering")
                    reasons.append(
                        f"missingness of '{candidates[a]}' and '{candidates[b]}' "
                        f"co-occur (Cramer's V={abs(phi):.3f} >= {cluster_threshold})"
                    )

    metrics["max_missingness_cramers_v"] = float(max_v)
    metrics["n_missingness_pairs"] = float(n_pairs_compared)

    # --- verdict ----------------------------------------------------------
    if findings:
        status = AssumptionStatus.VIOLATED
    elif tested_pairs == 0 and n_pairs_compared == 0:
        status = AssumptionStatus.UNDETERMINED
        if degenerate_pairs:
            reasons.append(
                f"{degenerate_pairs} channel/covariate pairs were skipped because "
                f"the covariate has no variation, and no other pair was testable"
            )
        else:
            reasons.append(
                "no channel/covariate pair had enough rows in both groups; "
                "the test could not be run"
            )
    else:
        status = AssumptionStatus.HOLDS
        reasons.append(
            f"{tested_pairs} channel/covariate pairs ({scale_pairs} with a usable "
            f"dispersion test) and {n_pairs_compared} missingness pairs tested; no "
            f"association detected at Bonferroni alpha={corrected_alpha:.3g}"
        )
        if degenerate_pairs:
            reasons.append(
                f"{degenerate_pairs} pair(s) skipped: the covariate has no variation"
            )

    return AssumptionReport(
        assumption_id="A4",
        statement=statement,
        status=status,
        metrics=metrics,
        findings=tuple(sorted(set(findings))),
        reasons=tuple(reasons),
        n_used=len(rows),
    )


# --------------------------------------------------------------------------
# A7 -- vertical datum consistency
# --------------------------------------------------------------------------

_UNIT_FACTORS: Dict[str, float] = {
    "metre_to_foot": 1.0 / 0.3048,
    "foot_to_metre": 0.3048,
    "foot_to_metre_us_survey": 1200.0 / 3937.0,
}


def diagnose_vertical_datum(
    samples: object,
    *,
    elevation_key: str = "elevation",
    head_keys: Sequence[str] = ("head", "head_meas", "hydraulic_head"),
    group_key: Optional[str] = "dataset",
    negative_fraction_tolerance: float = 0.02,
    unit_factor_tolerance: float = 0.05,
    min_rows: int = 3,
    min_group_size: int = 3,
    binom_p0: float = 0.02,
) -> AssumptionReport:
    """Test elevation and head for a shared vertical datum.

    Three sub-tests are run.

    **Sign test.** With a common datum and a water table at or below ground, the
    depth ``elevation - head`` should be non-negative. The count of negative
    depths is compared with a binomial null at ``binom_p0`` (default 0.02,
    allowing for artesian wells and measurement error).

    **Unit-factor test.** When the data come from more than one cohort, the ratio
    of median depths between cohorts is compared with unity and with the known
    conversion factors in :data:`_UNIT_FACTORS`. A group ratio within tolerance
    of ``0.3048`` or ``3.280 84`` is a unit or datum mismatch, not
    hydrogeology -- this is the failure the framework's own reproducibility
    notes describe for DEM-derived elevations.

    **Correlation.** Under a single datum, elevation and head are strongly
    positively correlated across a monitoring network. A low pooled correlation
    indicates that at least two of the values are expressed on different
    vertical references.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    elevation_key : str
        Field holding ground elevation.
    head_keys : sequence of str
        Candidate head fields, tried in order; the first present wins.
    group_key : str, optional
        Cohort field used for the unit-factor test (``'dataset'`` from the field
        loader). Pass ``None`` to skip the grouped test.
    negative_fraction_tolerance, unit_factor_tolerance : float
        Tolerances for the sign and unit-factor tests.
    min_rows, min_group_size : int
        Sample-size floors.
    binom_p0 : float
        Null probability of a negative depth under a common datum.

    Returns
    -------
    AssumptionReport
        With ``assumption_id='A7'``.
    """
    rows = sample_list(samples)
    statement = ASSUMPTION_CATALOGUE["A7"]

    def head_of(row: Mapping[str, object]) -> Optional[float]:
        for key in head_keys:
            value = _coerce_float(row, key)
            if value is not None:
                return value
        return None

    usable: List[Tuple[float, float, str]] = []
    for row in rows:
        elevation = _coerce_float(row, elevation_key)
        if elevation is None:
            continue
        head = head_of(row)
        if head is None:
            continue
        group = str(row.get(group_key, "ungrouped")) if group_key else "ungrouped"
        usable.append((elevation, head, group))

    n_excluded = len(rows) - len(usable)
    if not usable:
        elevation_present = any(elevation_key in row for row in rows)
        head_present = any(head_of(row) is not None for row in rows)
        if not elevation_present and not head_present:
            return AssumptionReport(
                assumption_id="A7",
                statement=statement,
                status=AssumptionStatus.NOT_APPLICABLE,
                reasons=(
                    "no elevation and no head field present; the assumption does "
                    "not apply to these data",
                ),
                n_used=0,
                n_excluded=n_excluded,
            )
        return AssumptionReport(
            assumption_id="A7",
            statement=statement,
            status=AssumptionStatus.UNDETERMINED,
            reasons=(
                "elevation and head are not jointly populated on any row; "
                "the datum comparison cannot be made",
            ),
            n_used=0,
            n_excluded=n_excluded,
        )

    depths = np.array([e - h for e, h, _ in usable], dtype=float)
    elevations = np.array([e for e, _, _ in usable], dtype=float)
    heads = np.array([h for _, h, _ in usable], dtype=float)
    groups = [g for _, _, g in usable]

    findings: List[str] = []
    reasons: List[str] = []
    metrics: Dict[str, float] = {
        "n_rows": float(len(rows)),
        "n_rows_used": float(len(usable)),
        "n_rows_excluded": float(n_excluded),
        "depth_median": _quantile(depths, 0.5),
        "depth_p05": _quantile(depths, 0.05),
        "depth_p95": _quantile(depths, 0.95),
        "negative_depth_fraction": float(np.mean(depths < 0.0)),
    }

    too_few = len(usable) < min_rows

    # --- sign test --------------------------------------------------------
    n_negative = int(np.sum(depths < 0.0))
    tail = _binomial_tail_at_least(n_negative, len(usable), binom_p0)
    metrics["negative_depth_count"] = float(n_negative)
    metrics["negative_depth_binomial_p"] = float(tail)
    if not too_few and (
        tail < 0.01 or (metrics["negative_depth_fraction"] > negative_fraction_tolerance)
    ):
        findings.append("head_above_ground")
        reasons.append(
            f"{n_negative}/{len(usable)} rows have head above ground "
            f"(fraction {metrics['negative_depth_fraction']:.3f} > "
            f"tolerance {negative_fraction_tolerance}); check that elevation and "
            f"head share a datum"
        )

    # --- correlation ------------------------------------------------------
    r = _pearson(elevations, heads)
    metrics["elevation_head_pearson_r"] = r
    if not too_few and math.isfinite(r) and r < 0.5:
        findings.append("weak_elevation_head_correlation")
        reasons.append(
            f"elevation-head correlation is only r={r:.3f}; two vertical "
            f"references may be mixed"
        )

    # --- unit-factor test -------------------------------------------------
    if group_key is not None:
        by_group: Dict[str, List[float]] = {}
        for depth, group in zip(depths.tolist(), groups):
            by_group.setdefault(group, []).append(depth)
        eligible = {
            g: float(np.median(np.abs(v)))
            for g, v in by_group.items()
            if len(v) >= min_group_size
        }
        metrics["n_groups_tested"] = float(len(eligible))
        factor_reported = 0.0
        if len(eligible) >= 2:
            medians = sorted(eligible.items(), key=lambda kv: kv[1])
            base = medians[0][1]
            if base > 0.0:
                for name, median in medians[1:]:
                    ratio = median / base
                    for factor_name, factor in _UNIT_FACTORS.items():
                        if abs(ratio - factor) / factor <= unit_factor_tolerance:
                            findings.append("vertical_datum_unit_mismatch")
                            reasons.append(
                                f"median depth of '{name}' is {ratio:.4f} times that "
                                f"of '{medians[0][0]}', within tolerance of the "
                                f"{factor_name} factor {factor:.4f}: a unit or datum "
                                f"mismatch, not hydrogeology"
                            )
                            factor_reported = max(factor_reported, factor)
        metrics["unit_factor_detected"] = float(factor_reported)

    # --- verdict ----------------------------------------------------------
    if findings:
        status = AssumptionStatus.VIOLATED
    elif too_few:
        status = AssumptionStatus.UNDETERMINED
        reasons.append(
            f"only {len(usable)} rows carry both elevation and head; "
            f"need >= {min_rows}"
        )
    else:
        status = AssumptionStatus.HOLDS
        reasons.append(
            f"{len(usable)} rows compared; "
            f"{metrics['negative_depth_fraction']:.3f} negative depths and "
            f"r={r:.3f} are consistent with a shared datum"
        )

    return AssumptionReport(
        assumption_id="A7",
        statement=statement,
        status=status,
        metrics=metrics,
        findings=tuple(sorted(set(findings))),
        reasons=tuple(reasons),
        n_used=len(usable),
        n_excluded=n_excluded,
    )


# --------------------------------------------------------------------------
# A10 -- laboratory bias
# --------------------------------------------------------------------------

_DEFAULT_A10_IONS: Tuple[str, ...] = (
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
)
_REPLICATE_PRECISION_BOUNDS: Tuple[float, float] = (0.001, 0.5)
"""Pooled relative standard deviation bounds, as fractions.

Below the lower bound the replicates are implausibly identical (a copied or
rounded value); above the upper bound they are too dispersed for a single
laboratory method. Both are indicative rather than decisive, and both are
reported as findings rather than as hard failures.
"""


def diagnose_laboratory_bias(
    samples: object,
    *,
    ion_order: Optional[Iterable[str]] = None,
    group_key: Optional[str] = "dataset",
    replicate_keys: Sequence[str] = ("site_id",),
    date_keys: Sequence[str] = ("sample_date", "date", "timestamp"),
    cb_spread_threshold: float = 0.05,
    min_group_size: int = 5,
    min_replicates: int = 3,
    alpha: float = 0.01,
) -> AssumptionReport:
    """Screen for a laboratory or cohort effect in the analytical results.

    Two sub-tests are run.

    **Cohort charge-balance heterogeneity.** Charge balance is an internal
    consistency property of the major-ion analysis. In a well-behaved dataset it
    is a property of the samples, not of which laboratory or campaign produced
    them. A significant difference in charge-balance distribution between
    cohorts -- especially a consistent sign shift -- indicates an analytical
    effect such as one laboratory systematically under-reporting anions, or a
    unit or rounding convention applied in one campaign only.

    **Replicate precision.** Where a site is sampled more than once on the same
    date, the spread between those rows estimates analytical precision. The
    pooled relative standard deviation per ion is compared with
    :data:`_REPLICATE_PRECISION_BOUNDS`. Two failure modes are flagged:
    implausibly tight agreement (suggesting copied or rounded values rather than
    independent analyses) and implausibly loose agreement.

    Both sub-tests assume the ions are expressed in a consistent unit across the
    rows compared; charge balance itself is scale-invariant, so the first test is
    unit-free per cohort.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    ion_order : iterable of str, optional
        Species used for the charge-balance computation. Defaults to
        :data:`_DEFAULT_A10_IONS`.
    group_key : str, optional
        Cohort field. Pass ``None`` to skip the grouped test.
    replicate_keys : sequence of str
        Fields identifying a sampled location.
    date_keys : sequence of str
        Date fields used to separate sampling occasions; if none is present,
        all rows sharing ``replicate_keys`` are treated as one occasion.
    cb_spread_threshold : float
        Absolute between-cohort spread in median charge balance that is treated
        as a violation.
    min_group_size, min_replicates : int
        Sample-size floors.
    alpha : float
        Significance level for the between-cohort rank-sum test.

    Returns
    -------
    AssumptionReport
        With ``assumption_id='A10'``.
    """
    rows = sample_list(samples)
    statement = ASSUMPTION_CATALOGUE["A10"]
    ions = list(ion_order) if ion_order is not None else list(_DEFAULT_A10_IONS)

    findings: List[str] = []
    reasons: List[str] = []
    metrics: Dict[str, float] = {"n_rows": float(len(rows))}

    # --- charge balance per row -------------------------------------------
    cb_rows: List[Tuple[float, str]] = []
    n_cb_excluded = 0
    for row in rows:
        values: List[float] = []
        complete = True
        for ion in ions:
            value = _coerce_float(row, ion)
            if value is None:
                complete = False
                break
            values.append(value)
        if not complete:
            n_cb_excluded += 1
            continue
        cb_rows.append((charge_balance_ratio(values, ions), str(row.get(group_key, "ungrouped")) if group_key else "ungrouped"))

    metrics["n_rows_with_charge_balance"] = float(len(cb_rows))
    metrics["n_rows_excluded_from_cb"] = float(n_cb_excluded)
    tests_run = 0

    if cb_rows:
        metrics["cb_median"] = float(np.median([cb for cb, _ in cb_rows]))
        metrics["cb_iqr"] = _iqr(np.array([cb for cb, _ in cb_rows], dtype=float))

    if group_key is not None and cb_rows:
        by_group: Dict[str, List[float]] = {}
        for cb, group in cb_rows:
            by_group.setdefault(group, []).append(cb)
        eligible = {g: v for g, v in by_group.items() if len(v) >= min_group_size}
        metrics["n_groups_with_charge_balance"] = float(len(eligible))
        if len(eligible) >= 2:
            medians = {g: float(np.median(v)) for g, v in eligible.items()}
            spread = max(medians.values()) - min(medians.values())
            metrics["cb_between_group_median_spread"] = float(spread)
            names = sorted(eligible)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    a, b = eligible[names[i]], eligible[names[j]]
                    _, p_value, effect = _rank_sum_test(a, b)
                    tests_run += 1
                    if spread >= cb_spread_threshold and p_value < alpha:
                        findings.append("cohort_charge_balance_shift")
                        reasons.append(
                            f"charge balance differs between cohorts "
                            f"'{names[i]}' and '{names[j]}' "
                            f"(medians {medians[names[i]]:+.4f} vs "
                            f"{medians[names[j]]:+.4f}, spread {spread:.4f} >= "
                            f"{cb_spread_threshold}, p={p_value:.3g}, "
                            f"rank-biserial={effect:+.3f})"
                        )
        elif len(eligible) == 1:
            reasons.append(
                "only one cohort has enough complete major-ion analyses; the "
                "between-cohort comparison was not possible"
            )

    # --- replicate precision ---------------------------------------------
    replicates: Dict[Tuple[Any, ...], List[Mapping[str, object]]] = {}
    for row in rows:
        if any(key not in row for key in replicate_keys):
            continue
        date_value: Any = None
        for date_key in date_keys:
            if date_key in row:
                date_value = row.get(date_key)
                break
        key = tuple(row.get(k) for k in replicate_keys) + (date_value,)
        replicates.setdefault(key, []).append(row)

    replicate_sets = [v for v in replicates.values() if len(v) >= 2]
    metrics["n_replicate_sets"] = float(len(replicate_sets))
    pooled_rsd: Dict[str, float] = {}
    if replicate_sets:
        for ion in ions:
            spreads: List[float] = []
            for group in replicate_sets:
                values = [
                    v
                    for v in (_coerce_float(row, ion) for row in group)
                    if v is not None
                ]
                if len(values) < min_replicates:
                    continue
                mean = float(np.mean(values))
                if mean <= 0.0:
                    continue
                spreads.append(float(np.std(values, ddof=1)) / mean)
            if spreads:
                pooled_rsd[ion] = float(np.mean(spreads))
        for ion, rsd in pooled_rsd.items():
            low, high = _REPLICATE_PRECISION_BOUNDS
            if rsd < low:
                findings.append("replicate_agreement_implausibly_tight")
                reasons.append(
                    f"replicate relative standard deviation for '{ion}' is "
                    f"{rsd:.5f} < {low}: independent analyses rarely agree this "
                    f"closely, consider copied or rounded values"
                )
            elif rsd > high:
                findings.append("replicate_agreement_implausibly_loose")
                reasons.append(
                    f"replicate relative standard deviation for '{ion}' is "
                    f"{rsd:.3f} > {high}: analytical precision is poor enough to "
                    f"compete with the effects being inferred"
                )
    else:
        reasons.append(
            "no replicate sample sets found (need >= 2 rows sharing "
            f"{list(replicate_keys)} on one date); precision could not be estimated"
        )

    for ion, rsd in pooled_rsd.items():
        metrics[f"replicate_rsd_{ion}"] = rsd
    metrics["n_tests_run"] = float(tests_run)

    # --- verdict ----------------------------------------------------------
    if findings:
        status = AssumptionStatus.VIOLATED
    elif not cb_rows and not replicate_sets:
        status = AssumptionStatus.UNDETERMINED
        reasons.append(
            "neither a complete major-ion panel nor a replicate set was "
            "available; no laboratory-bias screen could be run"
        )
    elif len(cb_rows) < min_group_size and not replicate_sets:
        status = AssumptionStatus.UNDETERMINED
        reasons.append(
            f"only {len(cb_rows)} complete major-ion analyses and no replicates; "
            f"below the floor of {min_group_size}"
        )
    else:
        status = AssumptionStatus.HOLDS
        reasons.append(
            f"{len(cb_rows)} charge-balanced rows and {len(replicate_sets)} "
            f"replicate sets screened; no cohort effect or precision anomaly "
            f"detected"
        )

    return AssumptionReport(
        assumption_id="A10",
        statement=statement,
        status=status,
        metrics=metrics,
        findings=tuple(sorted(set(findings))),
        reasons=tuple(reasons),
        n_used=len(cb_rows) + sum(len(v) for v in replicate_sets),
        n_excluded=n_cb_excluded,
    )


# --------------------------------------------------------------------------
# A34 -- tracer conservativeness
# --------------------------------------------------------------------------


def diagnose_tracer_conservativeness(
    samples: object,
    *,
    tracer: str = "Cl",
    reference: str = "Na",
    ion_order: Optional[Iterable[str]] = None,
    correlation_floor: float = 0.3,
    ratio_dispersion_threshold: float = 0.5,
    alpha: float = 0.01,
    min_pairs: int = 8,
) -> AssumptionReport:
    """Test whether a declared tracer behaves conservatively.

    A conservative tracer is altered only by mixing, evaporation and its own
    modelled decay. Three consequences are testable from a monitoring dataset.

    **Covariation with a reference species.** Conservative species move together
    under evaporation and mixing, so ``log(tracer)`` and ``log(reference)`` are
    strongly positively correlated across a network. A weak or negative
    correlation indicates a tracer-specific source or sink -- halite dissolution
    for chloride, wastewater, or precipitation of a chloride-bearing phase.

    **Redox invariance.** Chloride and the noble-gas tracers are redox-invariant.
    The tracer/reference log-ratio is compared across the redox classes returned
    by :func:`hydrosheaf.models.redox.classify_redox`; a systematic difference
    means the tracer responds to redox chemistry and is not conservative.

    **Ratio stability.** The dispersion (interquartile range) of
    ``log(tracer/reference)`` quantifies how much tracer-specific behaviour
    remains after removing the shared conservative signal. A small, stable ratio
    is the signature of conservativeness.

    Both the correlation and the ratio tests are invariant to a common unit
    factor, so the diagnostic may be run on any consistently expressed dataset.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    tracer : str
        The tracer whose conservativeness is asserted. Defaults to ``'Cl'``,
        matching ``Config.residence_time_tracer``.
    reference : str
        A second species expected to share the conservative signal. Must differ
        from ``tracer``; the diagnostic raises ``ValueError`` otherwise, since a
        tracer cannot validate itself.
    ion_order : iterable of str, optional
        Unused placeholder retained for interface symmetry with the other
        diagnostics; accepted and ignored.
    correlation_floor : float
        Minimum acceptable Pearson correlation of the two log series.
    ratio_dispersion_threshold : float
        Maximum acceptable interquartile range of the log-ratio series.
    alpha : float
        Significance level for the redox-invariance test.
    min_pairs : int
        Sample-size floor; below it the test is ``UNDETERMINED``.

    Returns
    -------
    AssumptionReport
        With ``assumption_id='A34'``.
    """
    rows = sample_list(samples)
    statement = ASSUMPTION_CATALOGUE["A34"]

    if tracer == reference:
        raise ValueError(
            "tracer and reference must differ: a tracer cannot be validated "
            "against itself."
        )

    pairs: List[Tuple[float, float, str]] = []
    n_excluded = 0
    for row in rows:
        t = _coerce_float(row, tracer)
        r = _coerce_float(row, reference)
        if t is None or r is None or t <= 0.0 or r <= 0.0:
            n_excluded += 1
            continue
        pairs.append((t, r, classify_redox(row)))

    if not pairs:
        return AssumptionReport(
            assumption_id="A34",
            statement=statement,
            status=AssumptionStatus.UNDETERMINED,
            reasons=(
                f"no row carries both '{tracer}' and '{reference}' as positive "
                f"values; the conservativeness test needs both",
            ),
            n_used=0,
            n_excluded=n_excluded,
        )

    log_tracer = np.log(np.array([t for t, _, _ in pairs], dtype=float))
    log_reference = np.log(np.array([r for _, r, _ in pairs], dtype=float))
    log_ratio = log_tracer - log_reference
    classes = [c for _, _, c in pairs]

    findings: List[str] = []
    reasons: List[str] = []
    metrics: Dict[str, float] = {
        "n_rows": float(len(rows)),
        "n_rows_used": float(len(pairs)),
        "n_rows_excluded": float(n_excluded),
    }

    r = _pearson(log_tracer, log_reference)
    metrics["log_tracer_log_reference_pearson_r"] = r
    metrics["log_ratio_iqr"] = _iqr(log_ratio)
    metrics["log_ratio_median"] = _quantile(log_ratio, 0.5)

    # A degenerate series carries no information about covariation, so the
    # sub-test is skipped rather than reported as passing.
    tracer_sd = float(np.std(log_tracer, ddof=1)) if len(log_tracer) > 1 else 0.0
    reference_sd = (
        float(np.std(log_reference, ddof=1)) if len(log_reference) > 1 else 0.0
    )
    ratio_sd = float(np.std(log_ratio, ddof=1)) if len(log_ratio) > 1 else 0.0
    metrics["tracer_log_sd"] = tracer_sd
    metrics["reference_log_sd"] = reference_sd
    metrics["log_ratio_sd"] = ratio_sd
    covariation_testable = tracer_sd > 0.0 and reference_sd > 0.0
    dispersion_testable = ratio_sd > 0.0
    sub_tests_run = int(covariation_testable) + int(dispersion_testable)

    if len(pairs) >= min_pairs:
        if covariation_testable and math.isfinite(r) and r < correlation_floor:
            findings.append("tracer_reference_decorrelated")
            reasons.append(
                f"log('{tracer}') and log('{reference}') correlate at only "
                f"r={r:.3f} < {correlation_floor}; '{tracer}' carries a signal "
                f"the conservative reference does not share"
            )
        iqr = metrics["log_ratio_iqr"]
        if dispersion_testable and math.isfinite(iqr) and iqr > ratio_dispersion_threshold:
            findings.append("tracer_reference_ratio_dispersed")
            reasons.append(
                f"IQR of log('{tracer}'/'{reference}') is {iqr:.3f} > "
                f"{ratio_dispersion_threshold}; tracer-specific behaviour is "
                f"large relative to the shared conservative signal"
            )

        # Redox invariance.
        by_class: Dict[str, List[float]] = {}
        for value, klass in zip(log_ratio.tolist(), classes):
            by_class.setdefault(klass, []).append(value)
        eligible = {k: v for k, v in by_class.items() if len(v) >= max(3, min_pairs // 2)}
        metrics["n_redox_classes_tested"] = float(len(eligible))
        names = sorted(eligible)
        if len(eligible) >= 2:
            sub_tests_run += 1
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                _, p_value, effect = _rank_sum_test(eligible[names[i]], eligible[names[j]])
                if p_value < alpha and abs(effect) >= 0.3:
                    findings.append("tracer_redox_dependence")
                    reasons.append(
                        f"log('{tracer}'/'{reference}') differs between redox "
                        f"classes '{names[i]}' and '{names[j]}' "
                        f"(p={p_value:.3g}, rank-biserial={effect:+.3f}); a "
                        f"conservative tracer should be redox-invariant"
                    )
    else:
        reasons.append(
            f"only {len(pairs)} joint observations of '{tracer}' and "
            f"'{reference}'; need >= {min_pairs}"
        )

    metrics["n_sub_tests_run"] = float(sub_tests_run)

    # --- verdict ----------------------------------------------------------
    if findings:
        status = AssumptionStatus.VIOLATED
    elif len(pairs) < min_pairs or sub_tests_run == 0:
        status = AssumptionStatus.UNDETERMINED
        if len(pairs) >= min_pairs and sub_tests_run == 0:
            reasons.append(
                f"'{tracer}' and/or '{reference}' has no variation across the "
                f"{len(pairs)} retained rows (log-sd {tracer_sd:.3g} and "
                f"{reference_sd:.3g}), so neither covariation nor dispersion can "
                f"be assessed; a conservative tracer cannot be confirmed from "
                f"constant data"
            )
    else:
        status = AssumptionStatus.HOLDS
        reasons.append(
            f"{len(pairs)} joint observations; r={r:.3f}, log-ratio IQR="
            f"{metrics['log_ratio_iqr']:.3f}, and no redox dependence detected at "
            f"alpha={alpha}"
        )

    return AssumptionReport(
        assumption_id="A34",
        statement=statement,
        status=status,
        metrics=metrics,
        findings=tuple(sorted(set(findings))),
        reasons=tuple(reasons),
        n_used=len(pairs),
        n_excluded=n_excluded,
    )


# --------------------------------------------------------------------------
# A37 -- cross-species error independence
# --------------------------------------------------------------------------


def diagnose_error_independence(
    samples: object,
    *,
    ion_order: Optional[Iterable[str]] = None,
    n_permutations: int = 200,
    seed: int = DEFAULT_SEED,
    alpha: float = 0.01,
    effect_threshold: float = 0.3,
    min_rows: int = 8,
    min_species: int = 3,
) -> AssumptionReport:
    """Test whether cross-species scatter is consistent with independent errors.

    Raw concentrations are compositional, so their correlations are dominated by
    the constant-sum constraint and would reject independence by construction.
    The diagnostic therefore works in **centred log-ratio** (CLR) coordinates,
    which are free of the closure.

    A common additive-log-ratio choice -- ``log(x_i / x_ref)`` against one
    reference species -- is *not* used, and the reason is worth stating because
    it is a trap. Every such ratio shares the term ``-log(x_ref)``, which induces
    a positive correlation between all ratio coordinates even when the species
    errors are completely independent. With five species of equal relative error
    the induced off-diagonal correlation is about ``+0.5`` (measured on
    independent log-normal draws in this repository's test suite), which exceeds
    any sensible effect threshold and would make the diagnostic fire on clean
    data. CLR avoids this: the induced structure is exactly ``-1/(D-1)`` for
    equal log-variances, and the test below accounts for it rather than assuming
    it away.

    The test statistic is the largest absolute off-diagonal correlation of the
    CLR coordinates. Its null distribution is a **marginal-preserving
    independence surrogate**: per-species log-means and log-standard-deviations
    are estimated from the retained rows, and surrogate datasets are drawn with
    independent Gaussian log-errors having exactly those marginals. Each
    surrogate is pushed through the identical CLR and correlation steps, so the
    closure-induced negative structure is reproduced in the null and cannot be
    mistaken for a violation.

    A rejection means that a cause beyond the modelled closure links the
    species -- a shared analytical bias, a reagent or dilution affecting a
    charge group, or a compositional effect CLR does not remove. The framework's
    weighted least squares treats species errors as independent, so a rejection
    quantifies a genuine mis-specification of the objective.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    ion_order : iterable of str, optional
        Species to test. Defaults to every field with at least ``min_rows``
        numeric values.
    n_permutations : int
        Number of surrogate datasets for the null distribution.
    seed : int
        Seed for the surrogate generator; defaults to :data:`DEFAULT_SEED`.
    alpha : float
        Significance level.
    effect_threshold : float
        Minimum observed maximum absolute correlation before a significant
        result is called a violation.
    min_rows, min_species : int
        Sample-size floors.

    Returns
    -------
    AssumptionReport
        With ``assumption_id='A37'``.
    """
    rows = sample_list(samples)
    statement = ASSUMPTION_CATALOGUE["A37"]

    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)

    if ion_order is not None:
        candidates = [s for s in ion_order if s in keys]
    else:
        candidates = [
            key
            for key in keys
            if sum(1 for row in rows if _coerce_float(row, key) is not None) >= min_rows
        ]

    if len(candidates) < min_species:
        return AssumptionReport(
            assumption_id="A37",
            statement=statement,
            status=AssumptionStatus.UNDETERMINED,
            reasons=(
                f"only {len(candidates)} species have >= {min_rows} usable values; "
                f"need >= {min_species}",
            ),
            n_used=0,
        )

    matrix_rows: List[List[float]] = []
    n_excluded = 0
    for row in rows:
        values: List[float] = []
        complete = True
        for species in candidates:
            value = _coerce_float(row, species)
            if value is None or value <= 0.0:
                complete = False
                break
            values.append(value)
        if not complete:
            n_excluded += 1
            continue
        matrix_rows.append(values)

    metrics: Dict[str, float] = {
        "n_rows": float(len(rows)),
        "n_rows_used": float(len(matrix_rows)),
        "n_rows_excluded": float(n_excluded),
        "n_species_tested": float(len(candidates)),
        "n_permutations": float(n_permutations),
    }
    findings: List[str] = []
    reasons: List[str] = []

    if len(matrix_rows) < min_rows:
        return AssumptionReport(
            assumption_id="A37",
            statement=statement,
            status=AssumptionStatus.UNDETERMINED,
            metrics=metrics,
            reasons=(
                f"only {len(matrix_rows)} rows carry a complete panel across "
                f"{len(candidates)} species; need >= {min_rows}",
            ),
            n_used=len(matrix_rows),
            n_excluded=n_excluded,
        )

    concentrations = np.asarray(matrix_rows, dtype=float)
    log_concentrations = np.log(concentrations)
    # CLR: subtract the row-wise mean of the logs.
    clr_coordinates = log_concentrations - log_concentrations.mean(axis=1, keepdims=True)

    def max_abs_offdiag(x: np.ndarray) -> Tuple[float, float]:
        """Return (max |off-diagonal r|, max positive off-diagonal r)."""
        if x.shape[1] < 2:
            return 0.0, 0.0
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = np.corrcoef(x, rowvar=False)
        corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
        mask = ~np.eye(corr.shape[0], dtype=bool)
        if not mask.any():
            return 0.0, 0.0
        offdiag = corr[mask]
        return float(np.max(np.abs(offdiag))), float(np.max(offdiag))

    observed, observed_positive = max_abs_offdiag(clr_coordinates)
    n_species = len(candidates)
    metrics["max_abs_offdiagonal_correlation"] = observed
    metrics["max_positive_offdiagonal_correlation"] = observed_positive
    metrics["expected_independent_offdiagonal"] = -1.0 / max(1, n_species - 1)
    metrics["clr_trace"] = float(np.trace(np.cov(clr_coordinates, rowvar=False)))

    # Marginal-preserving independence surrogate.
    means = log_concentrations.mean(axis=0)
    sds = log_concentrations.std(axis=0, ddof=1)
    rng = np.random.default_rng(seed)
    n_extreme = 0
    surrogate_maxima: List[float] = []
    for _ in range(max(1, n_permutations)):
        draw = means + rng.normal(0.0, 1.0, size=(len(matrix_rows), n_species)) * sds
        draw = draw - draw.mean(axis=1, keepdims=True)
        value, _ = max_abs_offdiag(draw)
        surrogate_maxima.append(value)
        if value >= observed:
            n_extreme += 1
    p_value = (n_extreme + 1.0) / (len(surrogate_maxima) + 1.0)
    metrics["surrogate_p_value"] = float(p_value)
    metrics["surrogate_null_median"] = (
        float(np.median(surrogate_maxima)) if surrogate_maxima else 0.0
    )
    metrics["surrogate_null_p95"] = (
        float(np.quantile(surrogate_maxima, 0.95)) if surrogate_maxima else 0.0
    )

    if p_value < alpha and observed >= effect_threshold:
        findings.append("cross_species_error_dependence")
        reasons.append(
            f"maximum absolute off-diagonal CLR correlation is {observed:.3f}, "
            f"above the independence surrogate "
            f"(95th percentile {metrics['surrogate_null_p95']:.3f}, "
            f"p={p_value:.4f} < {alpha}); a shared cause links the species, so "
            f"the weighted least-squares independence assumption is mis-specified"
        )
    else:
        reasons.append(
            f"maximum absolute off-diagonal CLR correlation {observed:.3f} against "
            f"an independence surrogate with median "
            f"{metrics['surrogate_null_median']:.3f} and 95th percentile "
            f"{metrics['surrogate_null_p95']:.3f} (p={p_value:.4f}); no dependence "
            f"detected at alpha={alpha}"
        )
        if observed_positive > 0.0:
            reasons.append(
                f"largest positive off-diagonal correlation is "
                f"{observed_positive:.3f}; the closure-induced structure under "
                f"independence is {metrics['expected_independent_offdiagonal']:.3f}"
            )

    status = AssumptionStatus.VIOLATED if findings else AssumptionStatus.HOLDS

    return AssumptionReport(
        assumption_id="A37",
        statement=statement,
        status=status,
        metrics=metrics,
        findings=tuple(sorted(set(findings))),
        reasons=tuple(reasons),
        n_used=len(matrix_rows),
        n_excluded=n_excluded,
    )


# --------------------------------------------------------------------------
# Aggregate audit
# --------------------------------------------------------------------------


def assess_model_assumptions(
    samples: object,
    *,
    ion_order: Optional[Iterable[str]] = None,
    a4_channels: Optional[Iterable[str]] = None,
    a4_covariates: Optional[Iterable[str]] = None,
    a7_elevation_key: str = "elevation",
    a7_head_keys: Sequence[str] = ("head", "head_meas", "hydraulic_head"),
    a10_group_key: Optional[str] = "dataset",
    a10_replicate_keys: Sequence[str] = ("site_id",),
    a34_tracer: str = "Cl",
    a34_reference: str = "Na",
    a37_n_permutations: int = 200,
    seed: int = DEFAULT_SEED,
) -> AssumptionAudit:
    """Run all five assumption diagnostics and aggregate the outcome.

    The aggregate status is the weakest individual status, in the order

    ``VIOLATED`` > ``UNDETERMINED`` > ``NOT_APPLICABLE`` > ``HOLDS``,

    so that a single violation is never averaged away, and an unrun test is
    never reported as compliance.

    Parameters
    ----------
    samples : mapping or sequence of mappings
        Sample records.
    Other parameters are forwarded to the individual diagnostics.

    Returns
    -------
    AssumptionAudit
        With ``reports`` in catalogue order and a ``counts`` mapping.
    """
    reports: List[AssumptionReport] = [
        diagnose_missingness_at_random(
            samples, channels=a4_channels, covariates=a4_covariates
        ),
        diagnose_vertical_datum(
            samples, elevation_key=a7_elevation_key, head_keys=a7_head_keys
        ),
        diagnose_laboratory_bias(
            samples,
            ion_order=ion_order,
            group_key=a10_group_key,
            replicate_keys=a10_replicate_keys,
        ),
        diagnose_tracer_conservativeness(
            samples, tracer=a34_tracer, reference=a34_reference
        ),
        diagnose_error_independence(
            samples, ion_order=ion_order, n_permutations=a37_n_permutations, seed=seed
        ),
    ]

    counts: Dict[str, int] = {status.value: 0 for status in AssumptionStatus}
    for report in reports:
        counts[report.status.value] += 1

    if counts[AssumptionStatus.VIOLATED.value]:
        overall = AssumptionStatus.VIOLATED
    elif counts[AssumptionStatus.UNDETERMINED.value]:
        overall = AssumptionStatus.UNDETERMINED
    elif counts[AssumptionStatus.NOT_APPLICABLE.value] == len(reports):
        overall = AssumptionStatus.NOT_APPLICABLE
    else:
        overall = AssumptionStatus.HOLDS

    return AssumptionAudit(
        reports=tuple(reports), counts=counts, status=overall, seed=seed
    )
