"""Tests for the five assumption-audit diagnostics.

Each diagnostic is exercised on synthetic data constructed so that the
assumption is either satisfied or clearly violated. A diagnostic that cannot
separate those two cases is worthless, so every test here checks the verdict and
not merely that the function runs.

The tests also pin the honesty properties the module claims: an empty or
insufficient input must be ``UNDETERMINED`` rather than ``HOLDS``, and the
permutation test must be reproducible for a fixed seed.
"""

from __future__ import annotations

import json
import math

import pytest

from hydrosheaf.validation.assumption_diagnostics import (
    ASSUMPTION_CATALOGUE,
    AssumptionAudit,
    AssumptionReport,
    AssumptionStatus,
    assess_model_assumptions,
    diagnose_error_independence,
    diagnose_laboratory_bias,
    diagnose_missingness_at_random,
    diagnose_tracer_conservativeness,
    diagnose_vertical_datum,
)

# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------


def _mcar_rows(n: int = 40) -> list[dict]:
    """A panel where Fe is analysed at random, independent of depth."""
    rows = []
    for i in range(n):
        row = {
            "site_id": f"W{i:03d}",
            "screen_depth": 10.0 + (i % 20) * 5.0,
            "Ca": 1.0,
            "Mg": 0.5,
            "Cl": 0.4,
        }
        # Missing by a rule that is independent of screen_depth.
        if i % 2 == 0:
            row["Fe"] = 0.02
        rows.append(row)
    return rows


def _mnar_rows(n: int = 40) -> list[dict]:
    """A panel where Fe is analysed only in the shallowest wells."""
    rows = []
    for i in range(n):
        depth = 10.0 + i * 4.0  # 10 .. 166
        row = {
            "site_id": f"W{i:03d}",
            "screen_depth": depth,
            "Ca": 1.0,
            "Mg": 0.5,
            "Cl": 0.4,
        }
        if depth <= 90.0:
            row["Fe"] = 0.02
        rows.append(row)
    return rows


def _balanced_row(**overrides) -> dict:
    """A row with a modest positive charge balance."""
    row = {
        "Ca": 1.0,
        "Mg": 0.5,
        "Na": 1.0,
        "K": 0.05,
        "HCO3": 2.0,
        "Cl": 0.5,
        "SO4": 0.25,
        "NO3": 0.1,
        "F": 0.01,
    }
    row.update(overrides)
    return row


def _anion_excess_row(scale: float = 1.0) -> dict:
    """A row with too much anion charge, scaled to create within-group spread."""
    return _balanced_row(
        Ca=0.6 * scale,
        Mg=0.3 * scale,
        Na=0.6 * scale,
        HCO3=2.0,
        Cl=0.6 * scale,
        SO4=0.30 * scale,
    )


def _conservative_rows(n: int = 30) -> list[dict]:
    """Cl and Na both driven by a shared evaporation factor."""
    rows = []
    for i in range(n):
        evap = 1.0 + 0.02 * i
        rows.append(
            {
                "site_id": f"W{i:03d}",
                "Cl": 0.4 * evap,
                "Na": 1.0 * evap,
                "NO3": 0.1,
            }
        )
    return rows


def _nonconservative_rows(n: int = 30) -> list[dict]:
    """Cl varies independently of Na: a tracer-specific source or sink."""
    rows = []
    for i in range(n):
        rows.append(
            {
                "site_id": f"W{i:03d}",
                "Cl": 0.1 + 0.05 * ((i * 7) % 13),
                "Na": 1.0 + 0.005 * i,
                "NO3": 0.1,
            }
        )
    return rows


def _shared_factor_rows(n: int = 80, *, dependent: bool) -> list[dict]:
    """Five species with either independent or group-structured log errors.

    ``dependent=False`` draws each species' log-error independently, which is the
    null the diagnostic must not reject. ``dependent=True`` adds a shared
    row-level factor to the anions and a different one to the cations, modelling
    a cause that affects a charge group rather than a single species.
    """
    import numpy as np

    species = ["Ca", "Mg", "Na", "Cl", "SO4"]
    rng = np.random.default_rng(20260101)
    noise = rng.normal(0.0, 0.05, size=(n, len(species)))
    group_factors = rng.normal(0.0, 0.25, size=(n, 2))
    rows = []
    for i in range(n):
        if dependent:
            # Anions share one row factor, cations another. A factor shared by
            # *every* species would cancel in CLR and is correctly invisible.
            anion_factor = float(group_factors[i, 0])
            cation_factor = float(group_factors[i, 1])
        else:
            anion_factor = cation_factor = 0.0
        row = {"site_id": f"W{i:03d}"}
        for k, name in enumerate(species):
            is_anion = name in {"Cl", "SO4"}
            shared = anion_factor if is_anion else cation_factor
            row[name] = math.exp(0.5 * k + float(noise[i, k]) + shared)
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Catalogue and report mechanics
# --------------------------------------------------------------------------


def test_catalogue_covers_the_five_assumptions():
    assert set(ASSUMPTION_CATALOGUE) == {"A4", "A7", "A10", "A34", "A37"}


def test_report_rejects_unknown_assumption_id():
    with pytest.raises(ValueError):
        AssumptionReport(
            assumption_id="A99",
            statement="nonsense",
            status=AssumptionStatus.HOLDS,
        )


def test_report_to_dict_is_json_serialisable():
    report = diagnose_missingness_at_random(_mcar_rows())
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["assumption_id"] == "A4"
    assert payload["status"] in {s.value for s in AssumptionStatus}


def test_audit_json_is_deterministic():
    rows = _balanced_row()
    first = assess_model_assumptions([rows], a37_n_permutations=20, seed=7)
    second = assess_model_assumptions([rows], a37_n_permutations=20, seed=7)
    assert first.to_json() == second.to_json()


def test_audit_counts_match_reports():
    audit = assess_model_assumptions(_mcar_rows(), a37_n_permutations=20)
    assert isinstance(audit, AssumptionAudit)
    assert sum(audit.counts.values()) == len(audit.reports) == 5
    assert {r.assumption_id for r in audit.reports} == set(ASSUMPTION_CATALOGUE)


# --------------------------------------------------------------------------
# A4 -- missingness at random
# --------------------------------------------------------------------------


def test_a4_random_missingness_holds():
    report = diagnose_missingness_at_random(_mcar_rows())
    assert report.assumption_id == "A4"
    assert report.status is AssumptionStatus.HOLDS
    assert report.findings == ()
    assert report.metrics["n_pairs_tested"] >= 1


def test_a4_depth_dependent_missingness_is_violated():
    report = diagnose_missingness_at_random(_mnar_rows())
    assert report.status is AssumptionStatus.VIOLATED
    assert "missingness_covariate_association" in report.findings
    # Depth is lower where Fe was analysed, so the effect is signed negative;
    # what matters for detection is the magnitude.
    assert abs(report.metrics["strongest_rank_biserial"]) > 0.2


def test_a4_no_missingness_is_not_applicable():
    rows = [{"site_id": "W1", "Ca": 1.0, "Cl": 0.5} for _ in range(6)]
    report = diagnose_missingness_at_random(rows)
    assert report.status is AssumptionStatus.NOT_APPLICABLE


def test_a4_insufficient_rows_is_undetermined_not_holds():
    report = diagnose_missingness_at_random([{"site_id": "W1", "Ca": 1.0}])
    assert report.status is AssumptionStatus.UNDETERMINED


def test_a4_middle_band_missingness_is_caught_by_the_dispersion_test():
    """A middle-band rule leaves the covariate medians equal.

    A location test alone has no power against it, so the diagnostic also runs a
    dispersion test (rank-sum on absolute deviations from the pooled median).
    This case was found by adversarial testing and is the reason that test
    exists; if it regresses, the diagnostic silently loses the most common
    real-world pattern of structured missingness.
    """
    rows = []
    for i in range(60):
        depth = 10.0 + 3.0 * i  # 10 .. 187
        row = {"site_id": f"W{i:03d}", "screen_depth": depth, "Ca": 1.0, "Cl": 0.4}
        if 60.0 <= depth <= 120.0:
            row["Fe"] = 0.02
        rows.append(row)

    report = diagnose_missingness_at_random(rows)

    assert report.status is AssumptionStatus.VIOLATED
    assert "missingness_covariate_dispersion" in report.findings
    # The location test on its own genuinely fails to see this pattern.
    assert report.metrics["smallest_p_value"] > 0.01


def test_a4_constant_covariate_is_skipped_and_reports_undetermined():
    """A covariate with no variation cannot discriminate and must not be counted."""
    rows = []
    for i in range(40):
        row = {"site_id": f"W{i:03d}", "screen_depth": 50.0, "Ca": 1.0, "Cl": 0.4}
        if i % 2 == 0:
            row["Fe"] = 0.02
        rows.append(row)

    report = diagnose_missingness_at_random(rows)

    assert report.status is AssumptionStatus.UNDETERMINED
    assert report.metrics["n_pairs_tested"] == 0.0
    assert report.metrics["n_degenerate_pairs_skipped"] >= 1.0


# --------------------------------------------------------------------------
# A7 -- vertical datum
# --------------------------------------------------------------------------


def test_a7_common_datum_holds():
    rows = [
        {
            "site_id": f"W{i:03d}",
            "elevation": 100.0 + 5.0 * i,
            "head": 100.0 + 5.0 * i - 12.0,
        }
        for i in range(10)
    ]
    report = diagnose_vertical_datum(rows)
    assert report.status is AssumptionStatus.HOLDS
    assert report.metrics["negative_depth_fraction"] == 0.0
    assert report.metrics["elevation_head_pearson_r"] > 0.9


def test_a7_head_above_ground_is_violated():
    rows = [
        {
            "site_id": f"W{i:03d}",
            "elevation": 100.0 + 2.0 * i,
            "head": 130.0 + 2.0 * i,  # 30 m above ground everywhere
        }
        for i in range(10)
    ]
    report = diagnose_vertical_datum(rows)
    assert report.status is AssumptionStatus.VIOLATED
    assert "head_above_ground" in report.findings


def test_a7_unit_mismatch_between_cohorts_is_violated():
    rows = []
    for i in range(8):
        depth = 10.0 + 0.5 * i
        rows.append(
            {
                "site_id": f"M{i}",
                "dataset": "metres_cohort",
                "elevation": 200.0 + 3.0 * i,
                "head": 200.0 + 3.0 * i - depth,
            }
        )
    for i in range(8):
        depth_ft = (10.0 + 0.5 * i) / 0.3048  # the same depth expressed in feet
        rows.append(
            {
                "site_id": f"F{i}",
                "dataset": "feet_cohort",
                "elevation": (200.0 + 3.0 * i) / 0.3048,
                "head": (200.0 + 3.0 * i) / 0.3048 - depth_ft,
            }
        )
    report = diagnose_vertical_datum(rows)
    assert report.status is AssumptionStatus.VIOLATED
    assert "vertical_datum_unit_mismatch" in report.findings


def test_a7_missing_fields_is_not_applicable():
    rows = [{"site_id": f"W{i}", "Ca": 1.0, "Cl": 0.5} for i in range(5)]
    report = diagnose_vertical_datum(rows)
    assert report.status is AssumptionStatus.NOT_APPLICABLE


def test_a7_only_one_field_populated_is_undetermined():
    rows = [{"site_id": f"W{i}", "elevation": 100.0 + i} for i in range(5)]
    report = diagnose_vertical_datum(rows)
    assert report.status is AssumptionStatus.UNDETERMINED


# --------------------------------------------------------------------------
# A10 -- laboratory bias
# --------------------------------------------------------------------------


def test_a10_consistent_single_cohort_holds():
    rows = [
        _balanced_row(site_id=f"C{i:02d}", dataset="cohort_one") for i in range(12)
    ]
    report = diagnose_laboratory_bias(rows)
    assert report.status is AssumptionStatus.HOLDS
    assert report.findings == ()


def test_a10_cohort_charge_balance_shift_is_violated():
    rows = []
    for i in range(10):
        rows.append(
            _balanced_row(
                site_id=f"A{i:02d}",
                dataset="cohort_a",
                Cl=0.5 + 0.004 * i,  # small within-cohort spread
            )
        )
    for i in range(10):
        rows.append(
            _anion_excess_row(scale=1.0 + 0.004 * i)
            | {"site_id": f"B{i:02d}", "dataset": "cohort_b"}
        )
    report = diagnose_laboratory_bias(rows)
    assert report.status is AssumptionStatus.VIOLATED
    assert "cohort_charge_balance_shift" in report.findings
    assert report.metrics["cb_between_group_median_spread"] >= 0.05


def test_a10_copied_replicates_are_flagged():
    base = _balanced_row()
    rows = []
    for i in range(4):
        for rep in range(3):
            # Identical analyses across three nominally independent replicates.
            rows.append(dict(base, site_id=f"R{i}", sample_date="2020-01-01"))
    report = diagnose_laboratory_bias(rows, group_key=None)
    assert report.status is AssumptionStatus.VIOLATED
    assert "replicate_agreement_implausibly_tight" in report.findings


def test_a10_no_panel_and_no_replicates_is_undetermined():
    rows = [{"site_id": f"W{i}", "Cl": 0.5} for i in range(6)]
    report = diagnose_laboratory_bias(rows)
    assert report.status is AssumptionStatus.UNDETERMINED


# --------------------------------------------------------------------------
# A34 -- tracer conservativeness
# --------------------------------------------------------------------------


def test_a34_conservative_tracer_holds():
    report = diagnose_tracer_conservativeness(_conservative_rows())
    assert report.status is AssumptionStatus.HOLDS
    assert report.findings == ()
    assert report.metrics["log_tracer_log_reference_pearson_r"] > 0.9


def test_a34_nonconservative_tracer_is_violated():
    report = diagnose_tracer_conservativeness(_nonconservative_rows())
    assert report.status is AssumptionStatus.VIOLATED
    assert report.findings


def test_a34_rejects_self_reference():
    with pytest.raises(ValueError):
        diagnose_tracer_conservativeness(_conservative_rows(), tracer="Cl", reference="Cl")


def test_a34_insufficient_pairs_is_undetermined():
    rows = [{"site_id": "W1", "Cl": 0.4, "Na": 1.0}]
    report = diagnose_tracer_conservativeness(rows)
    assert report.status is AssumptionStatus.UNDETERMINED


def test_a34_constant_series_is_undetermined_not_holds():
    """Constant data cannot demonstrate conservativeness.

    With both series flat the correlation is undefined (nan) and the log-ratio
    IQR is zero, so every sub-test is vacuous. Reporting ``HOLDS`` here would
    launder the absence of evidence into a clean audit, which design rule 2
    forbids.
    """
    rows = [{"site_id": f"W{i}", "Cl": 0.5, "Na": 1.0, "NO3": 0.1} for i in range(20)]

    report = diagnose_tracer_conservativeness(rows)

    assert report.status is AssumptionStatus.UNDETERMINED
    assert report.metrics["n_sub_tests_run"] == 0.0
    assert report.metrics["tracer_log_sd"] == 0.0


# --------------------------------------------------------------------------
# A37 -- cross-species error independence
# --------------------------------------------------------------------------


def test_a37_independent_errors_hold():
    rows = _shared_factor_rows(80, dependent=False)
    report = diagnose_error_independence(
        rows,
        ion_order=["Ca", "Mg", "Na", "Cl", "SO4"],
        n_permutations=200,
    )
    assert report.status is AssumptionStatus.HOLDS
    assert report.findings == ()
    assert report.metrics["surrogate_p_value"] > 0.01


def test_a37_group_structured_errors_are_violated():
    rows = _shared_factor_rows(80, dependent=True)
    report = diagnose_error_independence(
        rows,
        ion_order=["Ca", "Mg", "Na", "Cl", "SO4"],
        n_permutations=200,
    )
    assert report.status is AssumptionStatus.VIOLATED
    assert "cross_species_error_dependence" in report.findings
    assert report.metrics["max_abs_offdiagonal_correlation"] >= 0.3


def test_a37_is_reproducible_for_a_fixed_seed():
    rows = _shared_factor_rows(40, dependent=True)
    kwargs = dict(ion_order=["Ca", "Mg", "Na", "Cl", "SO4"], n_permutations=50, seed=4242)
    first = diagnose_error_independence(rows, **kwargs)
    second = diagnose_error_independence(rows, **kwargs)
    assert first.metrics["surrogate_p_value"] == second.metrics["surrogate_p_value"]


def test_a37_too_few_species_is_undetermined():
    rows = [{"site_id": f"W{i}", "Ca": 1.0 + 0.01 * i} for i in range(10)]
    report = diagnose_error_independence(rows, ion_order=["Ca"])
    assert report.status is AssumptionStatus.UNDETERMINED


# --------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------


def test_aggregate_status_reflects_the_weakest_test():
    """A single violation must not be averaged away."""
    rows = []
    for i in range(20):
        # Depth-dependent missingness (A4 violation) in otherwise clean data.
        row = _balanced_row(site_id=f"W{i:02d}", dataset="one")
        row["elevation"] = 100.0 + i
        row["head"] = 88.0 + i
        row["screen_depth"] = 10.0 + i
        if i < 10:
            row["Fe"] = 0.02
        rows.append(row)
    audit = assess_model_assumptions(rows, a37_n_permutations=20)
    assert audit.status is AssumptionStatus.VIOLATED
    assert audit.by_id["A4"].status is AssumptionStatus.VIOLATED
    assert len(audit.violations) >= 1


def test_empty_input_never_reports_holds():
    """Absent evidence must not be read as compliance."""
    audit = assess_model_assumptions([])
    a4 = audit.by_id["A4"]
    assert a4.status is AssumptionStatus.UNDETERMINED
    assert a4.n_used == 0
    assert audit.status is not AssumptionStatus.HOLDS


# --------------------------------------------------------------------------
# Package surface
# --------------------------------------------------------------------------


def test_functions_are_exported_from_the_package_root():
    import hydrosheaf

    for name in (
        "assess_model_assumptions",
        "diagnose_missingness_at_random",
        "diagnose_vertical_datum",
        "diagnose_laboratory_bias",
        "diagnose_tracer_conservativeness",
        "diagnose_error_independence",
        "AssumptionStatus",
        "AssumptionReport",
        "AssumptionAudit",
    ):
        assert name in hydrosheaf.__all__, name
        assert getattr(hydrosheaf, name) is not None
