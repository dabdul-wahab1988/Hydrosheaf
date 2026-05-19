from __future__ import annotations

import math

import pandas as pd

from hydrosheaf.nuclear.old_groundwater import (
    OldGroundwaterPrior,
    aggregate_c14_correction_candidates,
    apply_he4_uncertainty_mode,
    build_old_groundwater_priors,
    diagnose_old_groundwater_constraints,
    prepare_c14_observation,
    _lookup_oldwater_prior,
)


def test_aggregate_c14_candidates_collects_plausible_models():
    df = pd.DataFrame(
        [
            {"SampleID": "A", "Model": "PHREEQC", "Corrected_Ao_pmC": 100.0, "Corrected_14C_sample_pmC": 60.0},
            {"SampleID": "A", "Model": "REV", "Corrected_Ao_pmC": 95.0, "Corrected_14C_sample_pmC": 58.0},
            {"SampleID": "A", "Model": "No 14C", "Corrected_Ao_pmC": math.nan, "Corrected_14C_sample_pmC": math.nan},
        ]
    )
    out = aggregate_c14_correction_candidates(df)
    assert out.loc[0, "c14_candidate_count"] == 2
    assert "PHREEQC" in out.loc[0, "c14_candidate_models_json"]


def test_prepare_c14_observation_ensemble_uses_candidate_median():
    obs = {"c14_pmc": 75.0, "corrected_c14_pmc": 62.0, "corrected_a0_pmc": 100.0}
    out, initial, diag = prepare_c14_observation(
        obs,
        mode="ensemble",
        candidate_corrected_pmcs="[60.0, 54.0, 66.0]",
        candidate_initial_pmcs="[100.0, 95.0, 105.0]",
        candidate_models='["PHREEQC", "REV", "F&G"]',
    )
    assert out["c14_pmc"] == 60.0
    assert initial == 100.0
    assert diag["c14_candidate_count"] == 3
    assert diag["c14_effective_source"] == "ensemble_weighted"


def test_apply_he4_uncertainty_mode_inflates_sigma_with_calibration_fields():
    out, diag = apply_he4_uncertainty_mode(
        {
            "he4_ccpg": 1.0e-6,
            "he4_sigma_ccpg": 1.0e-9,
            "he4_background_ccpg": 0.0,
            "he4_u_ppm": 2.5,
            "he4_th_ppm": 9.0,
            "he4_porosity": 0.25,
            "he4_bulk_density": 2.1,
        },
        mode="calibrated_uncertainty",
    )
    assert out["he4_sigma_ccpg"] > 1.0e-9
    assert abs(diag["he4_rate_uncertainty_fraction"] - 0.25) < 1.0e-12


def test_old_groundwater_diagnostic_distinguishes_conflict_and_single_tracer():
    conflict = diagnose_old_groundwater_constraints(
        {
            "c14_pmc": 80.0,
            "c14_initial_pmc": 100.0,
            "he4_ccpg": 2.0e-6,
            "he4_background_ccpg": 0.0,
            "he4_accumulation_rate_ccpg_per_year": 1.0e-11,
        }
    )
    assert conflict["old_groundwater_case"] == "14C+4He"
    assert conflict["old_groundwater_status"] == "conflict"

    single = diagnose_old_groundwater_constraints(
        {
            "he4_ccpg": 5.0e-7,
            "he4_background_ccpg": 0.0,
            "he4_accumulation_rate_ccpg_per_year": 1.0e-11,
        }
    )
    assert single["old_groundwater_case"] == "4He_only"
    assert single["old_groundwater_status"] == "single_tracer"


# --- Phase 5: Hierarchical Old-Water Priors ---

def test_build_old_groundwater_priors_groups_by_study_unit():
    df = pd.DataFrame(
        [
            {"site_id": "A", "StudyUnit": "SU1", "AqGroup": "AQ1", "Corrected_Ao_pmC": 100.0, "he4_background_ccpg": 1e-7, "he4_accumulation_rate_ccpg_per_year": 1e-11},
            {"site_id": "B", "StudyUnit": "SU1", "AqGroup": "AQ1", "Corrected_Ao_pmC": 102.0, "he4_background_ccpg": 1.2e-7, "he4_accumulation_rate_ccpg_per_year": 1.1e-11},
            {"site_id": "C", "StudyUnit": "SU2", "AqGroup": "AQ2", "Corrected_Ao_pmC": 95.0, "he4_background_ccpg": 2e-7, "he4_accumulation_rate_ccpg_per_year": 2e-11},
        ]
    )
    priors = build_old_groundwater_priors(df)
    assert len(priors) >= 2
    # Study-unit/aquifer prior should exist
    su1_aq1 = priors.get("SU1|AQ1")
    assert su1_aq1 is not None
    assert su1_aq1.n_support >= 2


def test_lookup_oldwater_prior_prefers_study_unit_over_aquifer():
    priors = {
        "SU1|AQ1": OldGroundwaterPrior(
            study_unit="SU1", aquifer_group="AQ1",
            a0_pmc_mean=100.0, a0_pmc_sigma=5.0,
            he4_background_mean=1e-7, he4_background_sigma=1e-8,
            he4_rate_mean=1e-11, he4_rate_sigma=1e-12,
            n_support=2,
        ),
        "|AQ1": OldGroundwaterPrior(
            study_unit="", aquifer_group="AQ1",
            a0_pmc_mean=90.0, a0_pmc_sigma=10.0,
            he4_background_mean=2e-7, he4_background_sigma=2e-8,
            he4_rate_mean=2e-11, he4_rate_sigma=2e-12,
            n_support=1,
        ),
    }
    prior, scope = _lookup_oldwater_prior(priors, "SU1", "AQ1")
    assert prior is not None
    assert scope == "study_unit_aquifer"
    assert prior.a0_pmc_mean == 100.0


def test_lookup_oldwater_prior_falls_back_to_global():
    global_prior = OldGroundwaterPrior(
        study_unit="", aquifer_group="",
        a0_pmc_mean=85.0, a0_pmc_sigma=10.0,
        he4_background_mean=1e-7, he4_background_sigma=1e-8,
        he4_rate_mean=1e-11, he4_rate_sigma=1e-12,
        n_support=10,
    )
    priors = {"global|fallback": global_prior}
    prior, scope = _lookup_oldwater_prior(priors, "UNKNOWN", "UNKNOWN")
    assert prior is not None
    assert scope == "global_fallback"
    assert prior.a0_pmc_mean == 85.0


def test_prepare_c14_observation_uses_prior_when_hierarchical():
    obs = {"c14_pmc": 75.0, "corrected_c14_pmc": 62.0, "corrected_a0_pmc": 100.0}
    prior = OldGroundwaterPrior(
        study_unit="SU1", aquifer_group="AQ1",
        a0_pmc_mean=95.0, a0_pmc_sigma=3.0,
        he4_background_mean=1e-7, he4_background_sigma=1e-8,
        he4_rate_mean=1e-11, he4_rate_sigma=1e-12,
        n_support=5,
    )
    out, initial, diag = prepare_c14_observation(obs, mode="hierarchical", prior=prior)
    assert diag["c14_effective_source"] == "hierarchical_prior"
    assert initial == 98.75
    assert out["c14_initial_pmc"] == 98.75
    assert out["c14_pmc"] == 62.0


def test_apply_he4_uses_prior_when_hierarchical():
    obs = {
        "he4_ccpg": 1.0e-6,
        "he4_sigma_ccpg": 1.0e-9,
        "he4_background_ccpg": 0.0,
        "he4_u_ppm": 2.5,
        "he4_th_ppm": 9.0,
        "he4_porosity": 0.25,
        "he4_bulk_density": 2.1,
    }
    prior = OldGroundwaterPrior(
        study_unit="SU1", aquifer_group="AQ1",
        a0_pmc_mean=95.0, a0_pmc_sigma=3.0,
        he4_background_mean=1e-7, he4_background_sigma=1e-8,
        he4_rate_mean=1e-11, he4_rate_sigma=1e-12,
        n_support=5,
    )
    out, diag = apply_he4_uncertainty_mode(obs, mode="hierarchical", prior=prior)
    assert diag["he4_uncertainty_mode"] == "hierarchical"
    assert out["he4_background_ccpg"] == 1e-7
