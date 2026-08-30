"""Unit tests for the M3 identified-TTD infeasibility diagnostics.

Every test here is synthetic and hand-checkable.  None of them touches the
national dataset, the frozen protocol, or the solver's behaviour on real
tracer kernels, so the geometry the diagnostic reports can be verified
independently of whether the data are loadable.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import run_m3_infeasibility_diagnostics as diag  # noqa: E402

from hydrosheaf.nuclear.ttd_identified import TracerConstraint  # noqa: E402

# A three-bin response whose achievable range over the simplex is [1, 3].
ENVELOPE_RESPONSE = (1.0, 2.0, 3.0)


def _constraint(
    tracer: str,
    response,
    observed: float,
    sigma: float,
    likelihood: str = "gaussian",
) -> TracerConstraint:
    return TracerConstraint(
        tracer=tracer,
        response=response,
        observed=observed,
        sigma=sigma,
        likelihood=likelihood,
    )


# ---------------------------------------------------------------------------
# Envelope diagnostic
# ---------------------------------------------------------------------------


def test_response_envelope_is_the_simplex_range():
    constraint = _constraint("3H", ENVELOPE_RESPONSE, observed=2.0, sigma=0.1)
    assert diag.response_envelope(constraint) == (1.0, 3.0)


def test_interval_above_the_envelope_is_flagged():
    # [9.8, 10.2] lies entirely above the achievable [1, 3].
    constraint = _constraint("3H", ENVELOPE_RESPONSE, observed=10.0, sigma=0.1)
    assert diag.violates_envelope(constraint, 1.96) is True


def test_interval_below_the_envelope_is_flagged():
    # [-0.2, 0.2] lies entirely below the achievable [1, 3].
    constraint = _constraint("3H", ENVELOPE_RESPONSE, observed=0.0, sigma=0.1)
    assert diag.violates_envelope(constraint, 1.96) is True


def test_interval_inside_the_envelope_is_not_flagged():
    constraint = _constraint("3H", ENVELOPE_RESPONSE, observed=2.0, sigma=0.1)
    assert diag.violates_envelope(constraint, 1.96) is False


def test_widening_the_interval_can_clear_a_violation():
    # observed 4.0 with sigma 0.5: violating at k=1.0 ([3.5, 4.5]), reachable
    # at k=6.0 ([1.0, 7.0]) because the lower end touches the envelope.
    constraint = _constraint("3H", ENVELOPE_RESPONSE, observed=4.0, sigma=0.5)
    assert diag.violates_envelope(constraint, 1.0) is True
    assert diag.violates_envelope(constraint, 6.0) is False


def test_censored_likelihoods_are_one_sided():
    # upper_censored constrains only response @ w <= observed + margin, so an
    # observation far ABOVE the envelope is satisfiable and must not be flagged.
    high = _constraint(
        "CFC11", ENVELOPE_RESPONSE, observed=10.0, sigma=0.1, likelihood="upper_censored"
    )
    assert diag.violates_envelope(high, 1.96) is False
    low = _constraint(
        "CFC11", ENVELOPE_RESPONSE, observed=0.0, sigma=0.1, likelihood="upper_censored"
    )
    assert diag.violates_envelope(low, 1.96) is True

    # lower_censored is the mirror image.
    low_ok = _constraint(
        "CFC12", ENVELOPE_RESPONSE, observed=0.0, sigma=0.1, likelihood="lower_censored"
    )
    assert diag.violates_envelope(low_ok, 1.96) is False
    high_bad = _constraint(
        "CFC12", ENVELOPE_RESPONSE, observed=10.0, sigma=0.1, likelihood="lower_censored"
    )
    assert diag.violates_envelope(high_bad, 1.96) is True


# ---------------------------------------------------------------------------
# Feasibility and minimal infeasible subsets
# ---------------------------------------------------------------------------

# Two bins.  "first" reads w[0], "second" reads w[1], and w[0] + w[1] = 1.
# Demanding w[0] ~ 0.9 and w[1] ~ 0.9 simultaneously is infeasible, while each
# demand alone is satisfiable.  "slack" always evaluates to exactly 1.0.
CONFLICT_K = 1.0


def _first(observed: float, sigma: float = 0.01) -> TracerConstraint:
    return _constraint("first", (1.0, 0.0), observed, sigma)


def _second(observed: float, sigma: float = 0.01) -> TracerConstraint:
    return _constraint("second", (0.0, 1.0), observed, sigma)


def _slack() -> TracerConstraint:
    return _constraint("slack", (1.0, 1.0), 1.0, 0.5)


def test_singletons_of_the_conflicting_pair_are_each_feasible():
    assert diag.is_feasible([_first(0.9)], 2, CONFLICT_K) is True
    assert diag.is_feasible([_second(0.9)], 2, CONFLICT_K) is True


def test_disjoint_half_spaces_are_jointly_infeasible():
    assert diag.is_feasible([_first(0.9), _second(0.9)], 2, CONFLICT_K) is False


def test_minimal_infeasible_subset_of_a_pairwise_conflict_has_size_two():
    panel = [_first(0.9), _second(0.9), _slack()]
    subset = diag.minimal_infeasible_subset(panel, 2, CONFLICT_K)
    assert subset is not None
    assert len(subset) == 2
    assert set(subset) == {"first", "second"}


def test_feasible_panel_yields_no_minimal_infeasible_subset():
    panel = [_first(0.5), _second(0.5), _slack()]
    assert diag.is_feasible(panel, 2, CONFLICT_K) is True
    assert diag.minimal_infeasible_subset(panel, 2, CONFLICT_K) is None


def test_minimal_infeasible_subset_of_an_unreachable_observation_has_size_one():
    # response @ w <= 1.0 for every simplex point, so demanding 10.0 is
    # infeasible on its own: this is the envelope violation seen as an MIS.
    panel = [_constraint("unreachable", (1.0, 1.0), 10.0, 0.01), _first(0.5), _second(0.5)]
    subset = diag.minimal_infeasible_subset(panel, 2, CONFLICT_K)
    assert subset == ("unreachable",)


def test_minimal_infeasible_subset_prefers_the_smallest_size():
    # A size-3-only conflict: each of w[0], w[1], w[2] must be ~0.4, which is
    # pairwise satisfiable but sums above one.  Nothing smaller is infeasible.
    panel = [
        _constraint("a", (1.0, 0.0, 0.0), 0.4, 0.01),
        _constraint("b", (0.0, 1.0, 0.0), 0.4, 0.01),
        _constraint("c", (0.0, 0.0, 1.0), 0.4, 0.01),
    ]
    for pair in ((0, 1), (0, 2), (1, 2)):
        assert diag.is_feasible([panel[i] for i in pair], 3, CONFLICT_K) is True
    subset = diag.minimal_infeasible_subset(panel, 3, CONFLICT_K)
    assert subset is not None
    assert len(subset) == 3


# ---------------------------------------------------------------------------
# Skip accounting
# ---------------------------------------------------------------------------


def test_skip_ledger_counts_reasons_and_exception_types():
    ledger = diag.SkipLedger()
    ledger.record("missing_sample_year")
    ledger.record("missing_sample_year", 2)
    reason = ledger.record_exception("prepare_row_observations", ValueError("bad"))
    assert reason == "prepare_row_observations_error:ValueError"
    assert ledger.total == 4
    assert ledger.as_dict() == {
        "missing_sample_year": 3,
        "prepare_row_observations_error:ValueError": 1,
    }


# ---------------------------------------------------------------------------
# CLI contract
# ---------------------------------------------------------------------------


def test_default_sigma_multipliers_are_the_published_pair():
    args = diag.parse_args([])
    assert args.k_values == [1.96, 6.0]
    assert args.sources == "national"
    assert args.diagnostic == "both"


def test_repeated_k_flags_accumulate():
    args = diag.parse_args(["--k", "3.0", "--k", "4.0"])
    assert args.k_values == [3.0, 4.0]


def test_paths_are_derived_from_the_repository_layout():
    assert diag.REPO_ROOT.is_dir()
    assert (diag.REPO_ROOT / "hydrosheaf").is_dir()
    assert Path(diag.parse_args([]).protocol).name == "identified_ttd_protocol.yaml"


def test_flatten_to_records_renders_counts_and_rates():
    payload = {
        "envelope": {
            "by_k": {
                "k=6": {
                    "sigma_multiplier": 6.0,
                    "folds_examined": 100,
                    "folds_with_envelope_violation": 4,
                    "violating_constraints": 4,
                    "total_constraints": 200,
                    "by_tracer": {"3H": {"violations": 4, "seen": 50}},
                }
            },
            "accounting": {
                "eligible_folds": 101,
                "folds_examined": 100,
                "folds_skipped": 1,
                "fold_skips": {"no_linear_interval_constraints": 1},
                "row_skips": {},
            },
        }
    }
    records = diag.flatten_to_records(payload)
    by_metric = {(r["metric"], r["key"]): r for r in records}
    assert by_metric[("folds_with_envelope_violation", "all")]["rate"] == pytest.approx(0.04)
    assert by_metric[("violating_constraints_by_tracer", "3H")]["count"] == 4
    assert by_metric[("fold_skip", "no_linear_interval_constraints")]["count"] == 1
    assert by_metric[("fold_total", "eligible")]["count"] == 101
