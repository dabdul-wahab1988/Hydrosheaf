from __future__ import annotations

import math

import pandas as pd

from hydrosheaf.nuclear.old_groundwater import (
    aggregate_c14_correction_candidates,
    apply_he4_uncertainty_mode,
    diagnose_old_groundwater_constraints,
    prepare_c14_observation,
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
    assert diag["c14_effective_source"] == "ensemble_median"


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
