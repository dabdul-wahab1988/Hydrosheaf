from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "M3" / "m3_age_benchmark" / "scripts" / "run_m3_usgs_benchmark.py"
DESIGN_SCRIPT_PATH = REPO_ROOT / "M3" / "m3_age_benchmark" / "scripts" / "run_m3_design_matrix.py"
GAS_AUDIT_SCRIPT_PATH = REPO_ROOT / "M3" / "m3_age_benchmark" / "scripts" / "audit_m3_gas_corrections.py"
DESIGN_CONFIG_PATH = REPO_ROOT / "M3" / "m3_age_benchmark" / "configs" / "design_matrix.yaml"
USGS_INPUT = (
    REPO_ROOT
    / "M2"
    / "m2_benchmark"
    / "external"
    / "usgs_age"
    / "input"
    / "DataForNationalGroundwaterAge_1_1"
)


def _load_m3_module():
    spec = importlib.util.spec_from_file_location("run_m3_usgs_benchmark", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_design_module():
    spec = importlib.util.spec_from_file_location("run_m3_design_matrix", DESIGN_SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_gas_audit_module():
    spec = importlib.util.spec_from_file_location("audit_m3_gas_corrections", GAS_AUDIT_SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_m3_usgs_loader_preserves_one_row_per_site():
    if not USGS_INPUT.exists():
        pytest.skip("USGS source tables are not available in this checkout.")
    module = _load_m3_module()
    df = module.load_usgs_national_dataset()
    assert len(df) == df["site_id"].nunique()
    assert df["site_id"].duplicated().sum() == 0
    assert len(df) >= 1000


def test_m3_fit_age_uses_joint_lpm_canonical_age():
    module = _load_m3_module()

    class Fit:
        age_years = 12.5
        parameters = {"mean_age_years": 99.0}

    assert module._fit_age(Fit()) == 12.5


def test_m3_loader_adds_corrected_gases_and_helium_calibration():
    if not USGS_INPUT.exists():
        pytest.skip("USGS source tables are not available in this checkout.")
    module = _load_m3_module()
    df = module.load_usgs_national_dataset()
    assert "dissolved_gas_correction" in df.columns
    assert "he4_accumulation_rate_ccpg_per_year" in df.columns
    assert df["he4_accumulation_rate_ccpg_per_year"].notna().sum() > 0
    assert df["he4_source"].notna().sum() > 0
    assert "raw_sf6_pptv" in df.columns
    assert "raw_tritium_TU" in df.columns


def test_m3_c14_selection_prefers_geochemical_correction():
    module = _load_m3_module()
    df = module.pd.DataFrame(
        [
            {
                "SampleID": "A",
                "Model": "Use scaling factor",
                "Corrected_Ao_pmC": module.np.nan,
                "Corrected_14C_sample_pmC": module.np.nan,
            },
            {
                "SampleID": "A",
                "Model": "PHREEQC",
                "Corrected_Ao_pmC": 55.0,
                "Corrected_14C_sample_pmC": 80.0,
            },
        ]
    )
    selected = module._deduplicate_c14_rows(df)
    assert selected.loc[0, "Model"] == "PHREEQC"


def test_m3_nonidentifiability_reason_for_3he_without_tritium():
    module = _load_m3_module()
    reason = module._nonidentifiability_reason(
        {
            "tritium_TU": module.np.nan,
            "he3_trit_TU": 4.0,
        },
        "3He(trit)",
    )
    assert "3H/3He requires finite positive tritium" in reason


def test_m3_age_class_and_log_error_helpers():
    module = _load_m3_module()
    assert module._age_class(25.0) == "modern_le_50"
    assert module._age_class(500.0) == "intermediate_50_1000"
    assert module._age_class(5000.0) == "old_1000_30000"
    assert module._age_class(50000.0) == "very_old_gt_30000"
    assert abs(module._log10_abs_error(200.0, 100.0) - module.math.log10(2.0)) < 1e-12


def test_m3_design_factors_restore_raw_gases_and_disable_he4():
    module = _load_m3_module()
    row = module.pd.Series(
        {
            "sf6_pptv": 8.0,
            "raw_sf6_pptv": 3.0,
            "sf6_sigma_pptv": 0.8,
            "raw_sf6_sigma_pptv": 0.3,
            "he4_ccpg": 1.0e-5,
            "he4_accumulation_rate_ccpg_per_year": 2.0e-11,
        }
    )
    obs = {"sf6_pptv": 8.0, "sf6_sigma_pptv": 0.8, "he4_ccpg": 1.0e-5}
    adjusted, factors = module._apply_design_factors(
        obs,
        row,
        {"gas_correction_mode": "raw", "he4_mode": "disabled"},
    )
    assert adjusted["sf6_pptv"] == 3.0
    assert adjusted["sf6_sigma_pptv"] == 0.3
    assert adjusted["he4_ccpg"] is None
    assert factors["factor_gas_correction_mode"] == "raw"


def test_m3_screened_gas_selector_prefers_raw_for_pathological_corrected_fit():
    module = _load_m3_module()
    selected, reason = module._choose_screened_gas_result(
        {
            "est_age_multi": 42.6533,
            "young_gas_proxy_coherence": 0.45,
            "fit_aic": 20.0,
            "fit_objective": 300.0,
        },
        {
            "est_age_multi": 8.0,
            "young_gas_proxy_coherence": 0.28,
            "fit_aic": 8.0,
            "fit_objective": 20.0,
        },
    )
    assert selected == "raw"
    assert "raw" in reason


def test_m3_screened_gas_selector_preserves_corrected_when_proxy_coherence_is_better():
    module = _load_m3_module()
    selected, reason = module._choose_screened_gas_result(
        {
            "est_age_multi": 22.4,
            "young_gas_proxy_coherence": 0.09,
            "fit_aic": 17.0,
            "fit_objective": 260.0,
        },
        {
            "est_age_multi": 6.2,
            "young_gas_proxy_coherence": 0.44,
            "fit_aic": 8.0,
            "fit_objective": 12.0,
        },
    )
    assert selected == "usgs_dgm"
    assert "corrected" in reason


def test_m3_young_gas_masking_flags_supersaturated_sf6():
    module = _load_m3_module()
    masked, diag = module.calculate_tracer_reliability_weights(
        {
            "tritium_TU": 1.0,
            "he3_trit_TU": 2.0,
            "sf6_pptv": 100.0,
            "sf6_sigma_pptv": 0.2,
        },
        2010.0,
        12.0,
    )
    assert abs(masked["sf6_weight"] - 0.005) < 1e-6
    assert "SF6" in diag["young_gas_masked_tracers"]
    assert "above_historical_max_downweighted" in diag["young_gas_masked_reason"]


def test_m3_young_gas_masking_keeps_plausible_sf6():
    module = _load_m3_module()
    masked, diag = module.calculate_tracer_reliability_weights(
        {
            "tritium_TU": 1.5,
            "he3_trit_TU": 1.2,
            "sf6_pptv": 5.0,
            "sf6_sigma_pptv": 0.2,
        },
        2010.0,
        10.0,
    )
    assert masked["sf6_pptv"] == 5.0
    assert diag["young_gas_masked_count"] == 0


def test_m3_detects_screenable_gas_difference():
    module = _load_m3_module()
    row = module.pd.Series(
        {
            "sf6_pptv": 8.0,
            "raw_sf6_pptv": 4.0,
            "tritium_TU": 1.0,
            "raw_tritium_TU": 1.0,
        }
    )
    assert module._has_screenable_gas_difference(row)


def test_m3_skips_dual_screen_when_no_gas_difference():
    module = _load_m3_module()
    row = module.pd.Series(
        {
            "site_id": "A",
            "sample_year": 2010.0,
            "reference_age_years": 20.0,
            "lat": 0.0,
            "lon": 0.0,
            "tritium_TU": 1.0,
            "tritium_sigma_TU": 0.1,
            "raw_tritium_TU": 1.0,
            "raw_tritium_sigma_TU": 0.1,
            "LPM_TracersMod": "3H",
            "LPM_Name": "EM",
        }
    )

    original = module._fit_prepared_benchmark_row
    calls = []

    def fake_fit(*args, **kwargs):
        calls.append(1)
        return {
            "factor_gas_correction_mode": "usgs_dgm",
            "est_age_multi": 10.0,
            "est_age_3h": 10.0,
            "fit_aic": 1.0,
            "fit_objective": 1.0,
            "young_gas_proxy_coherence": 0.1,
            "site_id": "A",
            "ref_age": 20.0,
            "age_class": "modern_le_50",
            "modern_fraction": 1.0,
            "modern_age": module.np.nan,
            "old_age": module.np.nan,
            "reported_model": "EM",
            "multi_model": "EM",
            "model_strategy": "reported",
            "tracer_mode": "3H",
            "n_tracers": 1,
            "c14_initial_pmc": 100.0,
            "c14_correction_mode": "selected",
            "c14_correction_model": None,
            "c14_correction_strategy": "selected",
            "c14_candidate_count": 0,
            "c14_candidate_models": "",
            "c14_candidate_spread_pmc": 0.0,
            "c14_effective_source": "raw_fallback",
            "c14_effective_pmc": module.np.nan,
            "he4_calibrated": False,
            "he4_uncertainty_mode": "calibrated",
            "he4_accumulation_rate": module.np.nan,
            "he4_rate_uncertainty_fraction": module.np.nan,
            "he4_sigma_effective_ccpg": module.np.nan,
            "he4_source": None,
            "dissolved_gas_correction": "",
            "dgm_name": "",
            "young_gas_masking_mode": "screen_contaminated",
            "young_gas_masked_tracers": "",
            "young_gas_masked_reason": "",
            "young_gas_masked_count": 0,
            "young_gas_proxy_names": "3H_only",
            "young_gas_proxy_count": 1,
            "young_gas_proxy_values": "3H_only:10.0",
            "young_gas_selected_mode": "usgs_dgm",
            "young_gas_screen_reason": "",
            "young_gas_corrected_proxy_coherence": module.np.nan,
            "young_gas_raw_proxy_coherence": module.np.nan,
            "young_gas_corrected_fit_aic": module.np.nan,
            "young_gas_raw_fit_aic": module.np.nan,
            "young_gas_corrected_fit_objective": module.np.nan,
            "young_gas_raw_fit_objective": module.np.nan,
            "young_gas_delta_aic_raw_minus_corrected": module.np.nan,
            "young_gas_delta_objective_raw_minus_corrected": module.np.nan,
            "old_groundwater_case": "none",
            "old_groundwater_active_constraints": "",
            "old_groundwater_constraint_count": 0,
            "old_groundwater_status": "none",
            "old_groundwater_apparent_c14_age": module.np.nan,
            "old_groundwater_apparent_he4_age": module.np.nan,
            "old_groundwater_constraint_gap_log10": module.np.nan,
            "fit_converged": True,
            "error_multi": 10.0,
            "log10_error": 0.3,
            "within_factor_2": True,
            "within_factor_10": True,
            "failure_reason": "",
        }

    module._fit_prepared_benchmark_row = fake_fit
    try:
        result = module.fit_benchmark_row(row, age_steps=8, model_strategy="reported", factors={"gas_correction_mode": "screened_dgm"})
    finally:
        module._fit_prepared_benchmark_row = original

    assert len(calls) == 1
    assert result["young_gas_screen_reason"] == "no_screenable_difference"


def test_m3_design_matrix_config_loads():
    module = _load_design_module()
    config = module.load_design_matrix(DESIGN_CONFIG_PATH)
    scenario_ids = {scenario["scenario_id"] for scenario in config["scenarios"]}
    assert config["defaults"]["age_steps"] == module.usgs.M3_DEFAULT_AGE_STEPS
    assert "parity_reported_corrected" in scenario_ids
    assert "ablation_raw_gases" in scenario_ids
    assert "screened_dgm_gases" in scenario_ids
    assert "oldwater_c14_ensemble" in scenario_ids
    assert "oldwater_he4_uncertainty" in scenario_ids
    assert "hydrosheaf_selection_corrected" in scenario_ids


def test_m3_design_matrix_rejects_coarse_full_age_grid():
    module = _load_design_module()
    config = {"defaults": {"age_steps": 8}, "scenarios": [{"scenario_id": "parity_reported_corrected"}]}
    with pytest.raises(ValueError, match="age_steps"):
        module.resolve_age_steps(config, None, max_rows=None)
    assert module.resolve_age_steps(
        config,
        None,
        max_rows=None,
        allow_coarse_full_grid=True,
    ) == 8


def test_m3_design_matrix_pairwise_deltas():
    module = _load_design_module()
    df = module.pd.DataFrame(
        [
            {
                "scenario_id": "parity_reported_corrected",
                "site_id": "A",
                "log10_error": 0.3,
                "est_age_multi": 10.0,
                "within_factor_2": False,
                "within_factor_10": True,
            },
            {
                "scenario_id": "ablation_raw_gases",
                "site_id": "A",
                "log10_error": 0.1,
                "est_age_multi": 12.0,
                "within_factor_2": True,
                "within_factor_10": True,
            },
        ]
    )
    summary = module.summarize_pairwise_deltas(df)
    assert summary.loc[0, "scenario_id"] == "ablation_raw_gases"
    assert summary.loc[0, "median_delta_log10_error"] < 0
    assert summary.loc[0, "gained_factor_2_rows"] == 1


def test_m3_design_matrix_companion_paths_follow_custom_output():
    module = _load_design_module()
    output = REPO_ROOT / "M3" / "m3_age_benchmark" / "results" / "phase3_oldwater.csv"
    paths = module.companion_output_paths(output)
    assert paths["summary"].name == "phase3_oldwater_summary.csv"
    assert paths["pairwise"].name == "phase3_oldwater_pairwise_deltas.csv"
    assert paths["qa"].name == "phase3_oldwater_qa.md"


def test_m3_gas_audit_tracks_changed_values_and_effects():
    module = _load_gas_audit_module()
    design = module.pd.DataFrame(
        [
            {
                "scenario_id": "parity_reported_corrected",
                "site_id": "A",
                "log10_error": 0.4,
                "est_age_multi": 20.0,
                "within_factor_2": False,
            },
            {
                "scenario_id": "ablation_raw_gases",
                "site_id": "A",
                "log10_error": 0.1,
                "est_age_multi": 12.0,
                "within_factor_2": True,
            },
        ]
    )
    source = module.pd.DataFrame(
        [
            {
                "site_id": "A",
                "reference_age_years": 10.0,
                "sf6_pptv": 8.0,
                "raw_sf6_pptv": 4.0,
                "tritium_TU": 3.0,
                "raw_tritium_TU": 3.0,
                "he3_trit_TU": 3.0,
                "raw_he3_trit_TU": 1.0,
                "dissolved_gas_correction": "dgm_sf6",
            }
        ]
    )
    audit = module.build_gas_audit(design, source)
    assert audit.loc[0, "delta_log10_error_raw_minus_corrected"] < 0
    assert audit.loc[0, "sf6_changed_by_dgm"]
    assert "sf6" in audit.loc[0, "changed_gas_fields"]
    assert audit.loc[0, "corrected_3h3he_apparent_age"] > audit.loc[0, "raw_3h3he_apparent_age"]
    summary = module.summarize_audit(audit)
    assert summary.loc[0, "raw_improved_fraction"] == 1.0
