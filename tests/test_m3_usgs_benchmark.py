from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
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


def test_m3_fit_identifiability_rejects_underdetermined_misfit_and_broad_profiles():
    module = _load_m3_module()

    class Fit:
        converged = True
        n_tracers = 3
        effective_n_params = 2
        rmse_standardized = 1.0
        age_profile_log10_span = 0.2

    assert module._fit_identifiability(Fit()) == (True, "")
    Fit.n_tracers = 2
    assert "residual_degrees" in module._fit_identifiability(Fit())[1]
    Fit.n_tracers = 3
    Fit.rmse_standardized = 2.1
    assert "gross_standardized_misfit" in module._fit_identifiability(Fit())[1]
    Fit.rmse_standardized = 1.0
    Fit.age_profile_log10_span = 0.6
    assert "factor_3" in module._fit_identifiability(Fit())[1]
    Fit.age_profile_log10_span = 0.2
    supported, reason = module._fit_identifiability(Fit(), reported_model_supported=False)
    assert supported is False
    assert "reported_model_not_supported" in reason


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
    assert "raw_sf6_dissolved" in df.columns
    assert "dgm_sf6_pptv" in df.columns
    assert not df["raw_gas_atmospheric_equivalent_available"].any()
    assert df["he4_source"].eq("usgs_lpm_site_solution_rate").any()


def test_m3_loader_preserves_reported_lpm_input_scales_and_configuration():
    if not USGS_INPUT.exists():
        pytest.skip("USGS source tables are not available in this checkout.")
    module = _load_m3_module()
    df = module.load_usgs_national_dataset().set_index("site_id")
    scales = module._reported_prediction_scale_factors(df.loc["VRPDPAS1-06"])

    assert scales["3H"] == pytest.approx(1.0)
    assert scales["3H/3He"] == pytest.approx(1.0)
    assert scales["14C"] == pytest.approx(0.3)
    assert df.loc["VRPDPAS1-06", "reported_lpm_optimization_parameters"] == "Mean Age, Fraction"
    assert df.loc["VRPDPAS1-06", "reported_lpm_initial_age_c2_years"] == pytest.approx(75.0)


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
    assert factors.get("scenario_supported", True)


def test_m3_raw_gas_mode_is_unsupported_without_comparable_atmospheric_values():
    module = _load_m3_module()
    row = module.pd.Series({"sf6_pptv": 8.0, "raw_sf6_dissolved": 2.0})
    adjusted, factors = module._apply_design_factors(
        {"sf6_pptv": 8.0, "sf6_sigma_pptv": 0.8},
        row,
        {"gas_correction_mode": "raw"},
    )
    assert adjusted["sf6_pptv"] is None
    assert factors["gas_correction_comparison_supported"] is False
    assert factors["scenario_supported"] is False


def test_m3_corrected_mode_uses_usgs_table6_dgm_values():
    module = _load_m3_module()
    row = module.pd.Series(
        {
            "sf6_pptv": 8.0,
            "sf6_sigma_pptv": 0.8,
            "dgm_sf6_pptv": 6.0,
            "dgm_sf6_sigma_pptv": 0.6,
        }
    )
    adjusted, factors = module._apply_design_factors(
        {"sf6_pptv": 8.0, "sf6_sigma_pptv": 0.8},
        row,
        {"gas_correction_mode": "usgs_dgm"},
    )
    assert adjusted["sf6_pptv"] == 6.0
    assert adjusted["sf6_sigma_pptv"] == 0.6
    assert factors["gas_correction_input_source"] == "usgs_table6_dgm"


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
    )
    assert masked["sf6_pptv"] == 5.0
    assert diag["young_gas_masked_count"] == 0


def test_m3_young_gas_conflict_marks_contaminated_mixture():
    module = _load_m3_module()
    masked, diag = module.calculate_tracer_reliability_weights(
        {
            "tritium_TU": 1.0,
            "he3_trit_TU": 300.0,
            "sf6_pptv": 7.0,
            "sf6_sigma_pptv": 0.2,
        },
        2010.0,
    )
    assert masked["sf6_likelihood"] == "contaminated_mixture"
    assert "SF6:contaminated_mixture" in diag["young_gas_likelihood_assignments"]


def test_m3_reliability_weighting_api_cannot_accept_reference_age():
    import inspect

    module = _load_m3_module()
    signature = inspect.signature(module.calculate_tracer_reliability_weights)
    assert "reference_age" not in signature.parameters


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
    scenarios = {scenario["scenario_id"]: scenario for scenario in config["scenarios"]}
    assert scenarios["tracerlpm_strict_parity"]["reported_input_scaling"] is True
    assert scenarios["tracerlpm_strict_parity"]["reported_parameter_configuration"] is True
    assert scenarios["tracerlpm_parity_hier_oldwater"]["enabled"] is False
    assert scenarios["tracerlpm_parity_hier_oldwater"]["withdrawal_reason"]


def test_m3_design_matrix_excludes_withdrawn_scenarios_by_default():
    module = _load_design_module()
    config = {
        "defaults": {"age_steps": 8},
        "scenarios": [{"scenario_id": "withdrawn", "enabled": False}],
    }
    with pytest.raises(ValueError, match="No design-matrix scenarios selected"):
        module.run_design_matrix(pd.DataFrame(), config)


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


def test_m3_design_summaries_exclude_unsupported_gas_rows():
    module = _load_design_module()
    results = module.pd.DataFrame(
        [
            {
                "scenario_id": "parity_reported_corrected",
                "site_id": "A",
                "scenario_supported": True,
                "ref_age": 10.0,
                "est_age_multi": 12.0,
                "log10_error": 0.1,
                "within_factor_2": True,
                "within_factor_10": True,
                "he4_calibrated": False,
            },
            {
                "scenario_id": "ablation_raw_gases",
                "site_id": "A",
                "scenario_supported": False,
                "ref_age": 10.0,
                "est_age_multi": 9.0,
                "log10_error": 0.05,
                "within_factor_2": True,
                "within_factor_10": True,
                "he4_calibrated": False,
            },
        ]
    )
    summary = module.summarize_results(results)
    raw = summary[summary["scenario_id"] == "ablation_raw_gases"].iloc[0]
    assert raw["supported_rows"] == 0
    assert module.summarize_pairwise_deltas(results).empty


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


def test_reported_tracer_mask_keeps_only_usgs_reported_tracers():
    module = _load_m3_module()
    obs = {
        "tritium_TU": 1.0,
        "tritium_sigma_TU": 0.1,
        "he3_trit_TU": 2.0,
        "he3_trit_sigma_TU": 0.2,
        "sf6_pptv": 5.0,
        "sf6_sigma_pptv": 0.5,
        "cfc12_pptv": 400.0,
        "cfc12_sigma_pptv": 20.0,
        "c14_pmc": 80.0,
        "c14_sigma_pmc": 5.0,
    }
    masked = module._apply_reported_tracer_mask(obs, "3H, 3He(trit)")
    assert masked["tritium_TU"] == 1.0
    assert masked["he3_trit_TU"] == 2.0
    assert masked["sf6_pptv"] is None
    assert masked["sf6_sigma_pptv"] is None
    assert masked["cfc12_pptv"] is None
    assert masked["cfc12_sigma_pptv"] is None
    assert masked["c14_pmc"] is None
    assert masked["c14_sigma_pmc"] is None


def test_strict_parity_uses_reported_lpm_name():
    module = _load_m3_module()
    assert module._supported_reported_model("DM") == "DM"
    assert module._supported_reported_model("dm") == "DM"
    assert module._supported_reported_model("GA") == "GA"
    assert module._supported_reported_model("BMM-DM-DM") == "BMM-DM-DM"


def test_reported_bmm_configuration_uses_initial_not_final_parameters():
    module = _load_m3_module()
    row = {
        "reported_lpm_optimization_parameters": "Mean Age, Fraction",
        "reported_lpm_initial_age_c1_years": 25.0,
        "reported_lpm_initial_model_param1_c1": 0.01,
        "reported_lpm_initial_fraction_c1": 0.5,
        "reported_lpm_initial_age_c2_years": 75.0,
        "reported_lpm_initial_model_param1_c2": 0.01,
        # This final value must never enter the configuration helper.
        "reference_age_years": 40.6,
    }
    template, free = module._reported_fit_configuration(row, "BMM-DM-DM")

    assert template["mean_age_1_years"] == 25.0
    assert template["mean_age_2_years"] == 75.0
    assert template["dispersion"] == 0.01
    assert template["dispersion_2"] == 0.01
    assert free == ("mean_age_1_years", "binary_fraction")


def test_unsupported_reported_model_records_fallback():
    module = _load_m3_module()

    class FakeFit:
        model = "GA"
        parameters = {"mean_age_years": 10.0}
        objective = 1.0
        rmse_standardized = 0.1
        aic = 10.0
        bic = 12.0
        n_tracers = 1
        residuals = []
        converged = True
        note = ""
        age_years = 10.0

    original_fit_lpm_model = module.fit_lpm_model
    module.fit_lpm_model = lambda *args, **kwargs: FakeFit()
    try:
        row = module.pd.Series(
            {
                "site_id": "A",
                "sample_year": 2010.0,
                "reference_age_years": 20.0,
                "lat": 0.0,
                "lon": 0.0,
                "tritium_TU": 1.0,
                "tritium_sigma_TU": 0.1,
                "LPM_Name": "UNKNOWN_MODEL",
                "LPM_TracersMod": "3H",
            }
        )
        result = module._fit_prepared_benchmark_row(
            row,
            age_steps=8,
            factors={"lpm_strategy": "reported", "tracer_set": "reported"},
        )
        assert result["reported_model_supported"] is False
        assert "Unsupported model" in result["reported_model_fallback_reason"]
        assert result["multi_model"] == "GA"
    finally:
        module.fit_lpm_model = original_fit_lpm_model


def test_reported_uztt_missing_returns_zero():
    module = _load_m3_module()
    assert module._reported_uztt_years({}) == 0.0
    assert module._reported_uztt_years({"Rpt_UZtt_yrs": None}) == 0.0
    assert module._reported_uztt_years({"Rpt_UZtt_yrs": float("nan")}) == 0.0
    assert module._reported_uztt_years({"Rpt_UZtt_yrs": -5.0}) == 0.0


def test_apply_age_target_add_reported():
    module = _load_m3_module()
    row = {"Rpt_UZtt_yrs": 5.0, "Rpt_TotAge_yrs": 30.0}
    factors = {"age_target_mode": "reported_total", "uztt_mode": "add_reported"}
    est_age, diag = module._apply_age_target_mode(20.0, row, factors)
    assert est_age == 25.0
    assert diag["est_age_total_years"] == 25.0
    assert diag["reference_age_years"] == 30.0


def test_apply_age_target_saturated_only_never_below_0_1():
    module = _load_m3_module()
    row = {"Rpt_UZtt_yrs": 5.0, "Rpt_TotAge_yrs": 5.0}
    factors = {"age_target_mode": "saturated_only", "uztt_mode": "ignore"}
    est_age, diag = module._apply_age_target_mode(2.0, row, factors)
    assert est_age == 2.0
    assert diag["reference_age_years"] == 0.1


# --- Phase 3: Site-Specific Input Histories ---

def test_make_site_context_extracts_latitude_and_study_unit():
    module = _load_m3_module()
    ctx = module._make_site_context(
        {"site_id": "A", "sample_year": 2010.0, "lat": 45.0, "lon": -120.0, "StudyUnit": "SU1", "AqGroup": "AQ1", "recharge_temperature_c": 12.0}
    )
    assert ctx.latitude == 45.0
    assert ctx.study_unit == "SU1"
    assert ctx.recharge_temperature_c == 12.0


def test_get_site_histories_returns_mapping_and_meta():
    module = _load_m3_module()
    row = {"site_id": "A", "sample_year": 2010.0, "lat": 45.0, "lon": -120.0, "StudyUnit": "SU1"}
    histories, meta = module._get_site_histories(row)
    assert isinstance(histories, dict)
    assert "input_history_mode" in meta
    assert "input_history_region" in meta


# --- Phase 4: Censored Gas Likelihoods ---

def test_gas_likelihood_counts_zero_when_no_gases():
    module = _load_m3_module()
    counts = module._gas_likelihood_counts({"tritium_TU": 1.0, "tritium_sigma_TU": 0.1})
    assert counts["gaussian"] == 0
    assert counts["upper_censored"] == 0
    assert counts["contaminated_mixture"] == 0


def test_gas_likelihood_counts_tracks_contaminated_gas():
    module = _load_m3_module()
    counts = module._gas_likelihood_counts(
        {
            "sample_year": 2010.0,
            "sf6_pptv": 7.0,
            "sf6_sigma_pptv": 0.2,
            "sf6_likelihood": "contaminated_mixture",
        }
    )
    assert counts["contaminated_mixture"] == 1


# --- Phase 5: Hierarchical Old-Water Priors ---

def test_lookup_oldwater_prior_prefers_study_unit():
    module = _load_m3_module()
    from hydrosheaf.nuclear.old_groundwater import OldGroundwaterPrior

    priors = {
        "SU1|AQ1": OldGroundwaterPrior(
            study_unit="SU1", aquifer_group="AQ1",
            a0_pmc_mean=100.0, a0_pmc_sigma=5.0,
            he4_background_mean=1e-7, he4_background_sigma=1e-8,
            he4_rate_mean=1e-11, he4_rate_sigma=1e-12,
            n_support=2,
        ),
    }
    prior, scope = module._lookup_oldwater_prior(priors, "SU1", "AQ1")
    assert scope == "study_unit_aquifer"
    assert prior.a0_pmc_mean == 100.0


# --- Phase 6: Age Fraction Constraints ---

def test_fit_benchmark_row_includes_age_fraction_columns_when_enabled():
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
            "LPM_Name": "EM",
            "LPM_TracersMod": "3H",
            "FracAnthropocene": 0.8,
            "FracHolocene": 0.15,
            "FracPleistocene": 0.05,
        }
    )
    result = module.fit_benchmark_row(
        row,
        age_steps=8,
        model_strategy="reported",
        factors={"tracer_set": "reported", "use_age_fractions": True, "age_target_mode": "reported_total", "uztt_mode": "add_reported"},
    )
    assert "pred_frac_anthropocene" in result
    assert "age_fraction_loss" in result


# --- Phase 7: Gated BMA ---

def test_selection_strategy_records_bma_columns():
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
            "sf6_pptv": 5.0,
            "sf6_sigma_pptv": 0.5,
            "LPM_Name": "EM",
            "LPM_TracersMod": "3H, SF6",
        }
    )
    result = module.fit_benchmark_row(
        row,
        age_steps=8,
        model_strategy="selection",
        factors={"tracer_set": "reported", "age_target_mode": "reported_total", "uztt_mode": "add_reported"},
    )
    assert "bma_used" in result
    assert "bma_skip_reason" in result
    assert "bma_n_models" in result


# --- Phase 8: BMM Refinement Diagnostics ---

def test_bmm_refinement_diagnostics_exist():
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
            "LPM_Name": "BMM-DM-DM",
            "LPM_TracersMod": "3H",
        }
    )
    result = module.fit_benchmark_row(
        row,
        age_steps=8,
        model_strategy="reported",
        factors={"tracer_set": "reported"},
    )
    assert "refinement_attempted" in result
    assert "refinement_success" in result
    assert "refinement_message" in result


# --- Phase 10: Figure/Table Loaders ---

def test_figure_primary_results_prefers_strict_parity_scenario():
    import importlib.util
    fig_script = (
        REPO_ROOT
        / "M3"
        / "m3_age_benchmark"
        / "scripts"
        / "make_publication_figures.py"
    )
    if not fig_script.exists():
        pytest.skip("Figure script not found")
    spec = importlib.util.spec_from_file_location("make_publication_figures", fig_script)
    fig_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fig_module)
    # _mode_label_from_df should recognize strict parity
    df = pd.DataFrame({"scenario_id": ["tracerlpm_strict_parity"], "model_strategy": ["reported"]})
    assert fig_module._mode_label_from_df(df) == "Reported-Configuration Emulation"
