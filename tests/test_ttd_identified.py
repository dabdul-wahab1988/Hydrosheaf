import numpy as np
import pytest

from hydrosheaf.nuclear.joint_lpm import (
    TracerFitObservation,
    predict_lpm_tracers,
    tracer_response_kernel,
    transit_time_pdf,
    _integration_grid,
)
from hydrosheaf.nuclear.ttd_identified import (
    AgeFunctional,
    TracerConstraint,
    assess_held_out_tracer,
    infer_ttd_evidence,
    solve_ttd_identified_set,
    standard_age_functionals,
)


def _constraint(name, response, observed, sigma=1.0):
    return TracerConstraint(
        tracer=name,
        response=response,
        observed=observed,
        sigma=sigma,
    )


def test_full_rank_constraints_identify_age_bin_masses():
    ages = np.asarray([0.0, 10.0, 100.0])
    constraints = [
        _constraint("first", [1.0, 0.0, 0.0], 0.2),
        _constraint("second", [0.0, 1.0, 0.0], 0.3),
    ]
    functionals = [
        AgeFunctional("first_fraction", [1.0, 0.0, 0.0], 1.0e-10),
        AgeFunctional("old_fraction", [0.0, 0.0, 1.0], 1.0e-10),
    ]

    report = solve_ttd_identified_set(
        ages, constraints, functionals, sigma_multiplier=0.0
    )

    assert report.status == "IDENTIFIED"
    assert report.nullity == 0
    assert report.bounds["first_fraction"].lower == pytest.approx(0.2)
    assert report.bounds["old_fraction"].upper == pytest.approx(0.5)
    assert report.bounds["old_fraction"].status == "IDENTIFIED"


def test_age_fraction_cutoffs_are_explicit_inputs():
    functionals = standard_age_functionals(
        [0.0, 10.0, 20.0, 100.0], cutoffs_years=[10.0, 20.0]
    )
    coefficients = {item.name: item.coefficients.tolist() for item in functionals}

    assert coefficients["anthropocene"] == [1.0, 1.0, 0.0, 0.0]
    assert coefficients["holocene"] == [0.0, 0.0, 1.0, 0.0]
    assert coefficients["pleistocene"] == [0.0, 0.0, 0.0, 1.0]


def test_rank_deficiency_reports_identified_and_unresolved_functionals_separately():
    ages = np.asarray([0.0, 10.0, 100.0])
    report = solve_ttd_identified_set(
        ages,
        [_constraint("old-sensitive", [0.0, 0.0, 1.0], 0.5)],
        [
            AgeFunctional("old", [0.0, 0.0, 1.0], 1.0e-10),
            AgeFunctional("young", [1.0, 0.0, 0.0], 0.1),
        ],
        sigma_multiplier=0.0,
    )

    assert report.status == "PARTIALLY_IDENTIFIED"
    assert report.nullity == 1
    assert report.bounds["old"].status == "IDENTIFIED"
    assert report.bounds["old"].lower == pytest.approx(0.5)
    assert report.bounds["young"].status == "UNRESOLVED"
    assert report.bounds["young"].lower == pytest.approx(0.0)
    assert report.bounds["young"].upper == pytest.approx(0.5)


def test_inconsistent_tracers_abstain_instead_of_returning_a_point_estimate():
    report = solve_ttd_identified_set(
        [0.0, 10.0],
        [
            _constraint("a", [1.0, 0.0], 0.2),
            _constraint("b", [1.0, 0.0], 0.8),
        ],
        [AgeFunctional("young", [1.0, 0.0], 0.1)],
        sigma_multiplier=0.0,
    )

    assert report.status == "ABSTAIN"
    assert report.bounds == {}
    assert report.feasible_witness is None
    assert "No non-negative unit-mass TTD" in report.message


def test_constraint_order_and_affine_unit_conversion_do_not_change_bounds():
    ages = [0.0, 10.0, 100.0]
    original = [
        _constraint("a", [1.0, 0.0, 0.0], 0.2, 0.02),
        _constraint("b", [0.0, 1.0, 0.0], 0.3, 0.03),
    ]
    converted = [
        _constraint("b", [0.0, 1000.0, 0.0], 300.0, 30.0),
        _constraint("a", [1000.0, 0.0, 0.0], 200.0, 20.0),
    ]
    functional = AgeFunctional("old", [0.0, 0.0, 1.0])
    first = solve_ttd_identified_set(ages, original, [functional])
    second = solve_ttd_identified_set(ages, converted, [functional])

    assert first.bounds["old"].lower == pytest.approx(second.bounds["old"].lower)
    assert first.bounds["old"].upper == pytest.approx(second.bounds["old"].upper)


def test_adding_a_constraint_cannot_widen_a_feasible_bound():
    ages = [0.0, 10.0, 100.0]
    functional = AgeFunctional("old", [0.0, 0.0, 1.0])
    unconstrained = solve_ttd_identified_set(ages, [], [functional])
    constrained = solve_ttd_identified_set(
        ages,
        [_constraint("old-sensitive", [0.0, 0.0, 1.0], 0.4, 0.05)],
        [functional],
        sigma_multiplier=1.0,
    )

    assert constrained.bounds["old"].width <= unconstrained.bounds["old"].width
    assert constrained.bounds["old"].lower >= unconstrained.bounds["old"].lower
    assert constrained.bounds["old"].upper <= unconstrained.bounds["old"].upper


def test_no_functional_meets_a_gate_abstains_but_retains_feasible_bounds():
    report = solve_ttd_identified_set(
        [0.0, 10.0, 100.0],
        [_constraint("old-sensitive", [0.0, 0.0, 1.0], 0.4, 0.05)],
        [AgeFunctional("young", [1.0, 0.0, 0.0], 0.1)],
        sigma_multiplier=1.0,
    )

    assert report.status == "ABSTAIN"
    assert report.feasible_witness is not None
    assert report.bounds["young"].status == "UNRESOLVED"
    assert report.bounds["young"].width > 0.1


def test_exact_duplicate_constraint_does_not_create_information():
    ages = [0.0, 10.0, 100.0]
    constraint = _constraint("old-sensitive", [0.0, 0.0, 1.0], 0.4, 0.05)
    functional = AgeFunctional("young", [1.0, 0.0, 0.0])
    once = solve_ttd_identified_set(ages, [constraint], [functional])
    duplicated = solve_ttd_identified_set(ages, [constraint, constraint], [functional])

    assert duplicated.response_rank == once.response_rank
    assert duplicated.bounds["young"].lower == pytest.approx(once.bounds["young"].lower)
    assert duplicated.bounds["young"].upper == pytest.approx(once.bounds["young"].upper)


def test_upper_censored_observation_creates_only_an_upper_constraint():
    report = solve_ttd_identified_set(
        [0.0, 100.0],
        [
            TracerConstraint(
                "censored",
                [0.0, 1.0],
                observed=0.4,
                sigma=0.1,
                likelihood="upper_censored",
            )
        ],
        [AgeFunctional("old", [0.0, 1.0])],
        sigma_multiplier=1.0,
    )

    assert report.bounds["old"].lower == pytest.approx(0.0)
    assert report.bounds["old"].upper == pytest.approx(0.5)


def test_returned_extremal_witnesses_satisfy_mass_and_tracer_constraints():
    constraint = _constraint("old-sensitive", [0.0, 0.0, 1.0], 0.4, 0.05)
    report = solve_ttd_identified_set(
        [0.0, 10.0, 100.0],
        [constraint],
        [AgeFunctional("young", [1.0, 0.0, 0.0])],
        sigma_multiplier=1.0,
    )

    for witness in (
        report.bounds["young"].lower_witness,
        report.bounds["young"].upper_witness,
    ):
        assert witness.sum() == pytest.approx(1.0)
        assert np.all(witness >= -report.feasibility_tolerance)
        prediction = float(np.asarray(constraint.response) @ witness)
        assert 0.35 - 1.0e-8 <= prediction <= 0.45 + 1.0e-8


def test_grid_refinement_preserves_a_locked_linear_synthetic_bound():
    def solve(ages):
        response = np.asarray(ages) / 100.0
        return solve_ttd_identified_set(
            ages,
            [_constraint("linear-clock", response, 0.3)],
            [AgeFunctional("mean_age_years", ages, 1.0e-8, "years")],
            sigma_multiplier=0.0,
        ).bounds["mean_age_years"]

    coarse = solve([0.0, 50.0, 100.0])
    refined = solve([0.0, 25.0, 50.0, 75.0, 100.0])

    assert coarse.lower == pytest.approx(30.0)
    assert coarse.upper == pytest.approx(30.0)
    assert refined.lower == pytest.approx(coarse.lower)
    assert refined.upper == pytest.approx(coarse.upper)


@pytest.mark.parametrize(
    "tracer",
    ["3H", "3H/3He", "14C", "39Ar", "SF6", "CFC12", "85Kr", "4He"],
)
def test_response_kernel_matches_existing_lpm_forward_model(tracer):
    parameters = {"mean_age_years": 35.0, "dispersion": 0.1}
    ages, quadrature = _integration_grid(35.0, 500.0)
    density = transit_time_pdf("DM", parameters, ages, quadrature)
    mass = density * quadrature
    kernel_prediction = float(tracer_response_kernel(tracer, ages, 2020.0) @ mass)
    existing_prediction = predict_lpm_tracers(
        "DM", parameters, 2020.0, [tracer], max_age_years=500.0
    )[tracer]

    assert kernel_prediction == pytest.approx(existing_prediction, rel=1.0e-12)


def test_contaminated_observation_is_explicitly_excluded():
    report = infer_ttd_evidence(
        [
            TracerFitObservation(
                tracer="SF6",
                value=10.0,
                sigma=1.0,
                units="pptv",
                likelihood="contaminated_mixture",
            )
        ],
        sample_year=2020.0,
        age_grid_years=[0.0, 10.0, 100.0],
    )

    assert report.status == "ABSTAIN"
    assert report.excluded_observations[0]["tracer"] == "SF6"
    assert report.excluded_observations[0]["reason"] == (
        "no_declared_linear_interval_semantics"
    )


def test_held_out_check_does_not_condition_on_held_out_concentration():
    ages = [0.0, 10000.0]
    c14_response = tracer_response_kernel("14C", ages, 2020.0)
    report = solve_ttd_identified_set(
        ages,
        [_constraint("14C", c14_response, c14_response[0], 0.01)],
        [AgeFunctional("young", [1.0, 0.0], 0.01)],
        sigma_multiplier=0.0,
    )
    before = len(report.constraints)
    check = assess_held_out_tracer(
        report,
        TracerFitObservation("39Ar", 100.0, 1.0, "pmc"),
        sample_year=2020.0,
    )

    assert check["status"] == "COMPATIBLE"
    assert check["conditioned_on_held_out_observation"] is False
    assert len(report.constraints) == before


def test_public_root_exports_are_available():
    import hydrosheaf

    assert hydrosheaf.AgeFunctional is AgeFunctional
    assert hydrosheaf.solve_ttd_identified_set is solve_ttd_identified_set
