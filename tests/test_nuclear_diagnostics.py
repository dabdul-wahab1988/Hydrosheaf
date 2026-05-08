from __future__ import annotations

import math

from hydrosheaf.nuclear import (
    audit_graph_age_coherence,
    diagnose_lpm_identifiability,
    diagnose_tracer_disagreement,
    fit_lpm_models,
    infer_multi_tracer_age,
    tracer_removal_sensitivity,
)


def test_tracer_disagreement_flags_young_gas_old_carbon14():
    result = infer_multi_tracer_age(
        {
            "c14_pmc": 55.0,
            "sf6_pptv": 10.0,
        },
        sample_year=2025.0,
    )

    disagreement = result["disagreement"]
    assert disagreement["flags"]["uncorroborated_young_gas_with_old_14c"]
    assert disagreement["flags"]["young_old_contradiction"]
    assert disagreement["flags"]["possible_gas_contamination"]
    assert result["skipped_estimates"]


def test_diagnose_tracer_disagreement_reports_influence():
    report = diagnose_tracer_disagreement(
        {
            "age_years": 100.0,
            "n_tracers": 3,
            "estimates": [
                {
                    "tracer": "3H/3He",
                    "age_years": 12.0,
                    "ci_low_years": 8.0,
                    "ci_high_years": 18.0,
                    "method": "closed_system_ingrowth",
                    "observed_value": 1.0,
                },
                {
                    "tracer": "14C",
                    "age_years": 5000.0,
                    "ci_low_years": 3000.0,
                    "ci_high_years": 8000.0,
                    "method": "PFM",
                    "observed_value": 55.0,
                },
                {
                    "tracer": "4He",
                    "age_years": 7000.0,
                    "ci_low_years": 1000.0,
                    "ci_high_years": 20000.0,
                    "method": "screening",
                    "observed_value": 1.0e-7,
                },
            ],
            "skipped_estimates": [],
        }
    )

    assert report["conflict_score_log10_age"] > 1.0
    assert report["flags"]["tracer_conflict_high"]
    assert report["per_tracer_influence"]


def test_lpm_identifiability_flags_multiple_acceptable_families():
    observations = {
        "tritium_TU": 2.0,
        "he3_trit_TU": 1.0,
        "sf6_pptv": 8.0,
    }
    fits = fit_lpm_models(
        observations,
        sample_year=2025.0,
        models=["PFM", "EM", "DM"],
        max_age_years=120.0,
        age_steps=25,
    )

    diagnostic = diagnose_lpm_identifiability(fits, aic_delta_threshold=10.0)
    assert diagnostic["best_model"] in {"PFM", "EM", "DM"}
    assert diagnostic["acceptable_models"]
    assert "not_identifiable_from_available_tracers" in diagnostic["flags"]


def test_tracer_removal_sensitivity_reports_each_available_tracer():
    observations = {
        "tritium_TU": 2.5,
        "he3_trit_TU": 1.2,
        "sf6_pptv": 7.0,
    }

    sensitivity = tracer_removal_sensitivity(
        observations,
        sample_year=2025.0,
        models=["PFM", "EM"],
        max_age_years=120.0,
        age_steps=20,
    )

    removed = {row["removed_tracer"] for row in sensitivity["removal_sensitivity"]}
    assert {"3H", "3H/3HE", "SF6"}.issubset(removed)
    assert math.isfinite(sensitivity["full_fit"]["best_mean_age_years"])


def test_audit_graph_age_coherence_detects_downstream_younger_edge():
    audit = audit_graph_age_coherence(
        [("A", "B"), ("B", "C")],
        {
            "A": {"age_years": 50.0, "ci_low_years": 40.0, "ci_high_years": 60.0},
            "B": {
                "age_years": 20.0,
                "ci_low_years": 18.0,
                "ci_high_years": 22.0,
                "flags": {"possible_gas_contamination": True},
            },
            "C": {"age_years": 90.0, "ci_low_years": 70.0, "ci_high_years": 120.0},
        },
    )

    assert audit["n_violations"] == 1
    edge = audit["edges"][0]
    assert edge["violation"]
    assert edge["severe"]
    assert "contamination" in edge["explanation"]
