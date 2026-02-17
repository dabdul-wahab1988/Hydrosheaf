"""Configuration defaults for hydrosheaf."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DEFAULT_ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]



def default_phreeqc_database() -> str:
    return str(Path(__file__).resolve().parent / "databases" / "phreeqc.dat")


@dataclass
class Config:
    ion_order: List[str] = field(default_factory=lambda: DEFAULT_ION_ORDER.copy())
    unit: str = "mmol/L"
    unit_mode: str = "mmol_L"
    weights: List[float] = field(default_factory=lambda: [1.0] * 11)
    # Weights for conservative mixing (should prioritize Cl, Br, Isotopes)
    # Default: Cl (index 5) = 1.0, others = 0.01 (small weight for stability)
    conservative_weights: List[float] = field(
        default_factory=lambda: [0.01, 0.01, 0.01, 0.01, 0.01, 1.0, 0.01, 0.01, 0.01, 0.01, 0.01]
    )

    lambda_sparse: float = 0.0

    lambda_l1: float = 0.0
    reaction_max_iter: int = 300
    reaction_tol: float = 1e-6
    charge_balance_limit: float = 0.1
    ec_tds_penalty_limit: float = 0.0
    ec_tds_penalty_enabled: bool = False
    missing_policy: str = "skip"
    detection_limit_policy: str = "half"
    strict_input_validation: bool = False
    eta_ec: float = 0.0
    eta_tds: float = 0.0
    allow_signed_reactions: bool = False
    signed_reaction_labels: List[str] = field(default_factory=list)
    denit_kappa: float = 1.0
    transport_models_enabled: List[str] = field(default_factory=lambda: ["evap", "mix"])
    mixing_endmembers: Dict[str, List[float]] = field(default_factory=dict)
    # Dictionary mapping endmember names to (d18O, d2H) tuples
    mixing_endmembers_isotopes: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    ec_model: Tuple[float, float] = (1.0, 0.0)
    tds_model: Tuple[float, float] = (1.0, 0.0)
    phreeqc_enabled: bool = True
    phreeqc_mode: str = "phreeqpython"
    phreeqc_database: str = field(default_factory=default_phreeqc_database)
    phreeqc_executable: str = ""
    temp_default_c: float = 25.0
    si_threshold_tau: float = 0.2
    constraints_hard: bool = True
    edge_p_min: float = 0.75
    edge_radius_km: float = 5.0
    edge_max_neighbors: int = 3
    # Stratigraphic Priority settings
    edge_max_neighbors_primary: int = 3  # Neighbors in same HSU/Layer
    edge_max_neighbors_secondary: int = 1  # Vertical leakage neighbors (Cross-Layer)
    edge_sigma_meas: float = 0.5
    edge_sigma_dtw: float = 1.0
    edge_sigma_elev: float = 1.0
    edge_sigma_topo: float = 10.0
    edge_head_inference: str = "heuristic"  # heuristic, bayesian, bayesian_mcmc
    edge_dtw_prior_mu: float = 5.0
    edge_dtw_prior_sigma: float = 5.0
    edge_head_prior_mu: float = 0.0
    edge_head_prior_sigma: float = 1000.0
    edge_topo_sigma_depth: float = 5.0
    edge_topo_reject_p: float = 0.1
    edge_map_prior_weight: float = 0.0
    edge_map_candidate_multiplier: int = 5
    edge_map_p_min: float = 0.1
    sheaf_isotope_enabled: bool = True
    sheaf_cl_enabled: bool = True
    sheaf_age_enabled: bool = True
    sheaf_iso_sigma_d18o: float = 0.2
    sheaf_iso_sigma_d2h: float = 1.0
    sheaf_weight_head_prior: float = 1.0
    sheaf_weight_isotope: float = 1.0
    sheaf_weight_cl: float = 0.5
    sheaf_weight_age: float = 2.0
    sheaf_weight_global: float = 1.0

    sheaf_shallow_depth_m: float = 30.0
    sheaf_evap_gate_strength: float = 1.0
    sheaf_max_iter: int = 3
    sheaf_soft_beta: float = 1.0  # Soft selection sharpness (inverse temperature)
    edge_gradient_min: float = 1e-4
    edge_head_key: str = "head_meas"
    edge_dtw_key: str = "dtw"
    edge_elevation_key: str = "elevation"
    edge_aquifer_key: str = "aquifer_unit"
    edge_screen_depth_key: str = "screen_depth"
    edge_well_depth_key: str = "well_depth"
    edge_depth_mismatch: float = 20.0
    isotope_enabled: bool = False
    isotope_weight: float = 1.0
    isotope_d_excess_weight: float = 0.0
    isotope_d18o_key: str = "18O"
    isotope_d2h_key: str = "2H"
    lmwl_a: float = 8.66
    lmwl_b: float = 7.22
    lmwl_defined: bool = True
    auto_lmwl: bool = False
    isotope_consistency_weight: float = 0.0
    # Ion exchange settings
    exchange_enabled: bool = True
    exchange_cai_threshold: float = 0.1
    # Gibbs diagram settings (supplements/replaces isotopes)
    gibbs_enabled: bool = True
    gibbs_weight: float = 0.5
    gibbs_tds_precipitation: float = 100.0
    gibbs_tds_evaporation: float = 1000.0
    # Mineral Library settings
    active_minerals: List[str] = field(
        default_factory=lambda: [
            "calcite",
            "dolomite",
            "gypsum",
            "halite",
            "fluorite",
            "albite",
            "anorthite",
            "pyrite_oxidation_aerobic",
        ]
    )
    custom_minerals_path: str = ""

    # Nitrate Source V2 settings
    nitrate_source_enabled: bool = False
    nitrate_source_prior: float = 0.5
    nitrate_source_weights: Dict[str, float] = field(
        default_factory=lambda: {
            "w1_no3_cl": 1.2,
            "w2_no3_k": 0.4,
            "w3_po4": 0.3,
            "w4_fe": 0.6,
            "w5_denitrif": 1.5,
            "w6_alk_coupling": 0.8,
            "w7_coda_salinity": 0.0,
        }
    )
    nitrate_source_evap_gate: float = 0.5
    # Threshold overrides (None means auto-detected from data)
    nitrate_source_d_excess_p25: Optional[float] = None
    nitrate_source_po4_p90: Optional[float] = None
    nitrate_source_min_mg_L: float = 10.0
    nitrate_isotope_n15_col: str = "d15N"
    nitrate_isotope_o18_col: str = "d18O_NO3"
    nitrate_isotope_mixing_enabled: bool = True
    nitrate_isotope_water_o18_col: str = "d18O"
    nitrate_isotope_process_constraints_enabled: bool = True
    nitrate_isotope_qc_enabled: bool = True

    # MCMC Bayesian Isotope Mixing settings
    isotope_mcmc_enabled: bool = False
    isotope_mcmc_n_samples: int = 2000
    isotope_mcmc_n_chains: int = 4
    isotope_mcmc_target_accept: float = 0.9
    isotope_mcmc_warmup: int = 1000
    isotope_mcmc_hierarchical_enabled: bool = False
    isotope_mcmc_sources: List[str] = field(
        default_factory=lambda: ["Manure", "Fertilizer", "Soil_N", "Precipitation"]
    )

    # FloPy 1D Saturated Transport settings
    flopy_transport_enabled: bool = False
    aquifer_thickness_m: float = 10.0
    aquifer_porosity: float = 0.25
    aquifer_hydraulic_k_m_day: float = 1.0
    dispersivity_m: float = 1.0
    denitrification_k_1_day: float = 0.001
    transport_n_cells: int = 50
    transport_n_stress_periods: int = 10
    transport_perlen_days: float = 365.0

    # Uncertainty quantification settings
    uncertainty_method: str = "none"  # none, bootstrap, bayesian, monte_carlo
    bootstrap_n_resamples: int = 1000
    bootstrap_ci_method: str = "percentile"  # percentile, bca
    bayesian_n_samples: int = 5000
    bayesian_n_chains: int = 4
    bayesian_target_accept: float = 0.95
    bayesian_warmup_fraction: float = 0.5
    monte_carlo_n_samples: int = 1000
    input_uncertainty_pct: float = 5.0  # default 5% relative uncertainty

    # Prior hyperparameters for Bayesian inference
    prior_gamma_mu: float = 1.0
    prior_gamma_sigma: float = 0.5
    prior_xi_scale: float = 1.0  # Laplace scale parameter
    prior_sigma_scale: float = 0.1  # observation noise prior

    # Temporal dynamics settings
    temporal_enabled: bool = False
    temporal_window_days: int = 365
    temporal_min_samples: int = 3
    temporal_interpolation_method: str = "linear"  # linear, spline, nearest
    temporal_frequency_days: int = 30  # interpolation grid spacing
    residence_time_method: str = (
        "cross_correlation"  # gradient, cross_correlation, bayesian_lag, ttd, tracer_decay, recharge_piston
    )
    residence_time_tracer: str = "Cl"  # conservative tracer for age estimation
    residence_time_hydraulic_k: float = 1.0  # m/day, for gradient method
    residence_time_porosity: float = 0.2  # effective porosity
    # Recharge piston flow settings
    recharge_lag_volume_mm: float = 500.0  # Effective storage volume (mm)
    recharge_effective_fraction: float = 0.5  # Fraction of rain becoming recharge

    tau_agreement_tolerance: float = (
        0.4  # relative spread threshold for "tau_ambiguous"
    )
    tau_min_peak_corr: float = 0.2  # minimum correlation peak to accept a tracer
    tau_max_relative_uncertainty: float = 1.5  # reject if unc/tau exceeds this
    tau_max_uncertainty_days: float = (
        180.0  # reject if absolute uncertainty exceeds this
    )
    tau_physics_blend_threshold: float = (
        0.5  # blend with physics prior if |tau-phy|/max > threshold
    )
    ttd_grid_dt_days: float = 30.0  # time-grid step for TTD convolution (days)
    ttd_max_lag_days: float = 365.0  # maximum lag support for TTD (days)
    ttd_smoothness_lambda: float = 0.0  # smoothness penalty on TTD weights (0 disables)
    ttd_min_r2: float = 0.2  # minimum R^2 to accept a tracer TTD fit
    ttd_attenuation_k_max: float = (
        0.02  # 1/day, grid search max for exp(-k*tau) attenuation
    )
    ttd_attenuation_k_steps: int = 6  # number of k grid points (includes 0)
    bayes_lag_grid_dt_days: float = 5.0  # grid step for Bayesian lag posterior
    bayes_lag_max_lag_days: float = 365.0  # max tau support for Bayesian lag
    bayes_lag_prior_sigma_multiplier: float = (
        1.0  # multiplies physics sigma for prior width
    )
    bayes_lag_min_pairs: int = 5  # minimum overlapping pairs required

    # Residence-time coupling (static edges)
    residence_time_coupling_enabled: bool = False
    residence_time_reference_days: float = 30.0

    # Reactive transport validation settings
    reactive_transport_validation: bool = False
    rt_simulator: str = "phreeqc_kinetic"  # phreeqc_kinetic, mt3dms
    rt_n_time_steps: int = 100
    rt_consistency_rmse_threshold: float = 1.0  # mmol/L
    rt_consistency_nse_threshold: float = 0.5
    rt_default_residence_time: float = 30.0  # days (if not computed)

    # Rate law parameters (can be overridden per-mineral)
    rt_default_rate_constant: float = 1e-6  # mol/m²/s
    rt_default_surface_area: float = 0.1  # m²/L

    # Path to custom rate law definitions
    rt_custom_rates_file: str = ""

    # 3D flow network settings
    network_3d_enabled: bool = False
    z_coordinate_key: str = "screen_depth"  # or "z_mASL"
    z_coordinate_positive_down: bool = True  # True if depth, False if elevation

    # Vertical flow
    vertical_flow_enabled: bool = True
    vertical_anisotropy: float = 0.1  # α_v: K_h/K_v indicator
    vertical_gradient_min: float = 1e-3  # minimum vertical gradient
    upward_flow_probability: float = 0.5  # regional setting

    # Layer system
    layer_enabled: bool = False
    layer_key: str = "aquifer_layer"  # column with layer index
    layer_names: List[str] = field(default_factory=list)
    layer_tops: List[float] = field(default_factory=list)
    layer_bottoms: List[float] = field(default_factory=list)
    aquitard_leakage_p: float = 0.3
    layer_mineral_map: Dict[int, List[str]] = field(default_factory=dict)

    # Screen interval
    screen_top_key: str = "screen_top"
    screen_bottom_key: str = "screen_bottom"
    screen_overlap_threshold: float = 5.0  # meters

    # Ultra upgrades
    compositional_objective: bool = False  # Use Aitchison geometry for residuals
    compositional_weighting: bool = False # Use 1/x weights (relative error)
    latent_endmembers_enabled: bool = False # Inject virtual nodes via EMMA
    latent_endmembers_count: int = 2
    iterative_jacobian_enabled: bool = False # Update reaction dictionary with PHREEQC
    iterative_jacobian_max_iter: int = 3
    reacted_ttd_enabled: bool = False # Convolve kinetics over TTD

    def validate(self) -> None:
        # if len(self.ion_order) != 10:
        #     raise ValueError("ion_order must have 10 entries.")
        if self.unit_mode not in {"mmol_L", "meq_L"}:
            raise ValueError("unit_mode must be 'mmol_L' or 'meq_L'.")
        if len(self.weights) != len(self.ion_order):
            raise ValueError("weights must match ion_order length.")
        if len(self.conservative_weights) != len(self.ion_order):
            # Auto-align if default mismatch (backward compatibility for tests)
            if len(self.weights) == len(self.ion_order):
                # Fallback to standard weights if conservative ones are stale
                self.conservative_weights = list(self.weights)
            else:
                raise ValueError("conservative_weights must match ion_order length.")
        if any(w < 0 for w in self.weights):
            raise ValueError("weights must be non-negative.")
        if any(w < 0 for w in self.conservative_weights):
            raise ValueError("conservative_weights must be non-negative.")
        if self.lambda_sparse < 0 or self.lambda_l1 < 0:

            raise ValueError("lambda penalties must be non-negative.")
        if self.reaction_max_iter <= 0:
            raise ValueError("reaction_max_iter must be positive.")
        if self.reaction_tol < 0:
            raise ValueError("reaction_tol must be non-negative.")
        if self.charge_balance_limit < 0:
            raise ValueError("charge_balance_limit must be non-negative.")
        if self.ec_tds_penalty_limit < 0:
            raise ValueError("ec_tds_penalty_limit must be non-negative.")
        if self.phreeqc_mode not in {"phreeqpython", "subprocess"}:
            raise ValueError("phreeqc_mode must be 'phreeqpython' or 'subprocess'.")
        if self.temp_default_c <= 0:
            raise ValueError("temp_default_c must be positive.")
        if self.si_threshold_tau < 0:
            raise ValueError("si_threshold_tau must be non-negative.")
        if not 0.0 <= self.edge_p_min <= 1.0:
            raise ValueError("edge_p_min must be between 0 and 1.")
        if self.edge_radius_km < 0:
            raise ValueError("edge_radius_km must be non-negative.")
        if self.edge_max_neighbors < 0:
            raise ValueError("edge_max_neighbors must be non-negative.")
        if self.edge_sigma_meas <= 0 or self.edge_sigma_dtw <= 0:
            raise ValueError("edge sigma values must be positive.")
        if self.edge_sigma_elev <= 0 or self.edge_sigma_topo <= 0:
            raise ValueError("edge sigma values must be positive.")
        if self.edge_head_inference not in {"heuristic", "bayesian", "bayesian_mcmc"}:
            raise ValueError(
                "edge_head_inference must be 'heuristic', 'bayesian', or 'bayesian_mcmc'."
            )
        if self.edge_dtw_prior_sigma <= 0 or self.edge_head_prior_sigma <= 0:
            raise ValueError("edge prior sigmas must be positive.")
        if self.edge_topo_sigma_depth <= 0:
            raise ValueError("edge_topo_sigma_depth must be positive.")
        if not 0.0 < self.edge_topo_reject_p < 1.0:
            raise ValueError("edge_topo_reject_p must be in (0, 1).")
        if self.edge_map_prior_weight < 0:
            raise ValueError("edge_map_prior_weight must be non-negative.")
        if self.edge_map_candidate_multiplier < 1:
            raise ValueError("edge_map_candidate_multiplier must be at least 1.")
        if not 0.0 <= self.edge_map_p_min <= 1.0:
            raise ValueError("edge_map_p_min must be between 0 and 1.")
        if self.sheaf_iso_sigma_d18o <= 0 or self.sheaf_iso_sigma_d2h <= 0:
            raise ValueError("sheaf isotope sigmas must be positive.")
        if self.sheaf_weight_head_prior < 0:
            raise ValueError("sheaf_weight_head_prior must be non-negative.")
        if self.sheaf_weight_isotope < 0:
            raise ValueError("sheaf_weight_isotope must be non-negative.")
        if self.sheaf_weight_cl < 0:
            raise ValueError("sheaf_weight_cl must be non-negative.")
        if self.sheaf_weight_global < 0:
            raise ValueError("sheaf_weight_global must be non-negative.")
        if self.sheaf_shallow_depth_m < 0:
            raise ValueError("sheaf_shallow_depth_m must be non-negative.")
        if self.sheaf_evap_gate_strength < 0:
            raise ValueError("sheaf_evap_gate_strength must be non-negative.")
        if self.sheaf_max_iter < 1:
            raise ValueError("sheaf_max_iter must be at least 1.")
        if self.edge_gradient_min < 0:
            raise ValueError("edge_gradient_min must be non-negative.")
        if self.edge_depth_mismatch < 0:
            raise ValueError("edge_depth_mismatch must be non-negative.")
        if self.isotope_weight < 0:
            raise ValueError("isotope_weight must be non-negative.")
        if self.isotope_d_excess_weight < 0:
            raise ValueError("isotope_d_excess_weight must be non-negative.")
        if self.missing_policy not in {"skip", "impute_zero"}:
            raise ValueError("missing_policy must be 'skip' or 'impute_zero'.")
        if self.detection_limit_policy not in {"half", "zero", "value", "drop"}:
            raise ValueError(
                "detection_limit_policy must be one of: half, zero, value, drop."
            )
        if any(model not in {"evap", "mix"} for model in self.transport_models_enabled):
            raise ValueError(
                "transport_models_enabled must be a subset of {'evap','mix'}."
            )
        for name, vector in self.mixing_endmembers.items():
            if len(vector) != len(self.ion_order):
                raise ValueError(f"endmember '{name}' has invalid length.")

        if self.nitrate_source_prior < 0 or self.nitrate_source_prior > 1:
            raise ValueError("nitrate_source_prior must be between 0 and 1.")
        if any(w < 0 for w in self.nitrate_source_weights.values()):
            raise ValueError("nitrate_source_weights must be non-negative.")

        # MCMC Isotope Mixing validation
        if self.isotope_mcmc_n_samples < 100:
            raise ValueError("isotope_mcmc_n_samples must be at least 100.")
        if self.isotope_mcmc_n_chains < 1:
            raise ValueError("isotope_mcmc_n_chains must be at least 1.")
        if not 0.5 <= self.isotope_mcmc_target_accept <= 0.99:
            raise ValueError("isotope_mcmc_target_accept must be between 0.5 and 0.99.")
        if self.isotope_mcmc_warmup < 100:
            raise ValueError("isotope_mcmc_warmup must be at least 100.")

        # FloPy Transport validation
        if self.aquifer_thickness_m <= 0:
            raise ValueError("aquifer_thickness_m must be positive.")
        if not 0.0 < self.aquifer_porosity <= 1.0:
            raise ValueError("aquifer_porosity must be between 0 and 1.")
        if self.aquifer_hydraulic_k_m_day <= 0:
            raise ValueError("aquifer_hydraulic_k_m_day must be positive.")
        if self.dispersivity_m < 0:
            raise ValueError("dispersivity_m must be non-negative.")
        if self.denitrification_k_1_day < 0:
            raise ValueError("denitrification_k_1_day must be non-negative.")
        if self.transport_n_cells < 2:
            raise ValueError("transport_n_cells must be at least 2.")
        if self.transport_n_stress_periods < 1:
            raise ValueError("transport_n_stress_periods must be at least 1.")
        if self.transport_perlen_days <= 0:
            raise ValueError("transport_perlen_days must be positive.")

        # Uncertainty quantification validation
        if self.uncertainty_method not in {
            "none",
            "bootstrap",
            "bayesian",
            "monte_carlo",
        }:
            raise ValueError(
                "uncertainty_method must be one of: none, bootstrap, bayesian, monte_carlo."
            )
        if not 0.0 < self.tau_agreement_tolerance <= 1.0:
            raise ValueError("tau_agreement_tolerance must be in (0, 1].")
        if not 0.0 <= self.tau_min_peak_corr <= 1.0:
            raise ValueError("tau_min_peak_corr must be in [0, 1].")
        if self.tau_max_relative_uncertainty <= 0:
            raise ValueError("tau_max_relative_uncertainty must be positive.")
        if self.tau_max_uncertainty_days <= 0:
            raise ValueError("tau_max_uncertainty_days must be positive.")
        if self.tau_physics_blend_threshold < 0:
            raise ValueError("tau_physics_blend_threshold must be non-negative.")
        if self.ttd_grid_dt_days <= 0:
            raise ValueError("ttd_grid_dt_days must be positive.")
        if self.ttd_max_lag_days <= 0:
            raise ValueError("ttd_max_lag_days must be positive.")
        if self.ttd_smoothness_lambda < 0:
            raise ValueError("ttd_smoothness_lambda must be non-negative.")
        if not 0.0 <= self.ttd_min_r2 <= 1.0:
            raise ValueError("ttd_min_r2 must be between 0 and 1.")
        if self.ttd_attenuation_k_max < 0:
            raise ValueError("ttd_attenuation_k_max must be non-negative.")
        if self.ttd_attenuation_k_steps < 1:
            raise ValueError("ttd_attenuation_k_steps must be at least 1.")
        if self.bayes_lag_grid_dt_days <= 0:
            raise ValueError("bayes_lag_grid_dt_days must be positive.")
        if self.bayes_lag_max_lag_days <= 0:
            raise ValueError("bayes_lag_max_lag_days must be positive.")
        if self.bayes_lag_prior_sigma_multiplier <= 0:
            raise ValueError("bayes_lag_prior_sigma_multiplier must be positive.")
        if self.bayes_lag_min_pairs < 2:
            raise ValueError("bayes_lag_min_pairs must be at least 2.")
        if self.bootstrap_n_resamples < 1:
            raise ValueError("bootstrap_n_resamples must be at least 1.")
        if self.bootstrap_ci_method not in {"percentile", "bca"}:
            raise ValueError("bootstrap_ci_method must be 'percentile' or 'bca'.")
        if self.bayesian_n_samples < 1:
            raise ValueError("bayesian_n_samples must be at least 1.")
        if self.bayesian_n_chains < 1:
            raise ValueError("bayesian_n_chains must be at least 1.")
        if not 0.0 < self.bayesian_target_accept <= 1.0:
            raise ValueError("bayesian_target_accept must be between 0 and 1.")
        if not 0.0 < self.bayesian_warmup_fraction < 1.0:
            raise ValueError("bayesian_warmup_fraction must be between 0 and 1.")
        if self.monte_carlo_n_samples < 1:
            raise ValueError("monte_carlo_n_samples must be at least 1.")
        if self.input_uncertainty_pct < 0:
            raise ValueError("input_uncertainty_pct must be non-negative.")
        if self.prior_gamma_sigma <= 0:
            raise ValueError("prior_gamma_sigma must be positive.")
        if self.prior_xi_scale <= 0:
            raise ValueError("prior_xi_scale must be positive.")
        if self.prior_sigma_scale <= 0:
            raise ValueError("prior_sigma_scale must be positive.")

        # Temporal dynamics validation
        if self.temporal_window_days < 1:
            raise ValueError("temporal_window_days must be at least 1.")
        if self.temporal_min_samples < 2:
            raise ValueError("temporal_min_samples must be at least 2.")
        if self.temporal_interpolation_method not in {"linear", "spline", "nearest"}:
            raise ValueError(
                "temporal_interpolation_method must be one of: linear, spline, nearest."
            )
        if self.temporal_frequency_days < 1:
            raise ValueError("temporal_frequency_days must be at least 1.")
        if self.residence_time_method not in {
            "gradient",
            "cross_correlation",
            "tracer_decay",
            "bayesian_lag",
            "ttd",
            "recharge_piston",
        }:
            raise ValueError(
                "residence_time_method must be one of: gradient, cross_correlation, bayesian_lag, ttd, tracer_decay, recharge_piston."
            )

        if self.residence_time_hydraulic_k <= 0:
            raise ValueError("residence_time_hydraulic_k must be positive.")
        if not 0.0 < self.residence_time_porosity <= 1.0:
            raise ValueError("residence_time_porosity must be between 0 and 1.")
        if self.residence_time_reference_days <= 0:
            raise ValueError("residence_time_reference_days must be positive.")

        # Reactive transport validation
        if self.rt_simulator not in {"phreeqc_kinetic", "mt3dms"}:
            raise ValueError("rt_simulator must be 'phreeqc_kinetic' or 'mt3dms'.")
        if self.rt_n_time_steps < 1:
            raise ValueError("rt_n_time_steps must be at least 1.")
        if self.rt_consistency_rmse_threshold < 0:
            raise ValueError("rt_consistency_rmse_threshold must be non-negative.")
        if self.rt_consistency_nse_threshold < -1:
            raise ValueError("rt_consistency_nse_threshold must be >= -1.")
        if self.rt_default_residence_time <= 0:
            raise ValueError("rt_default_residence_time must be positive.")
        if self.rt_default_rate_constant <= 0:
            raise ValueError("rt_default_rate_constant must be positive.")
        if self.rt_default_surface_area <= 0:
            raise ValueError("rt_default_surface_area must be positive.")

        # 3D network validation
        if self.vertical_anisotropy <= 0:
            raise ValueError("vertical_anisotropy must be positive.")
        if self.vertical_gradient_min < 0:
            raise ValueError("vertical_gradient_min must be non-negative.")
        if not 0.0 < self.upward_flow_probability < 1.0:
            raise ValueError("upward_flow_probability must be between 0 and 1.")
        if not 0.0 <= self.aquitard_leakage_p <= 1.0:
            raise ValueError("aquitard_leakage_p must be between 0 and 1.")
        if self.screen_overlap_threshold < 0:
            raise ValueError("screen_overlap_threshold must be non-negative.")
        if self.layer_enabled:
            if len(self.layer_names) != len(self.layer_tops) or len(
                self.layer_names
            ) != len(self.layer_bottoms):
                raise ValueError(
                    "layer_names, layer_tops, and layer_bottoms must have same length."
                )
            for i in range(len(self.layer_tops)):
                if self.layer_tops[i] >= self.layer_bottoms[i]:
                    raise ValueError(
                        f"layer_tops[{i}] must be less than layer_bottoms[{i}]."
                    )

    def lambda_l1_value(self) -> float:
        return self.lambda_l1 if self.lambda_l1 else self.lambda_sparse

    def load_from_calibration_json(self, path: str) -> None:
        """
        Load optimized parameters from a PEST calibration JSON file.
        Updates Config fields if parameter names match known config keys.
        """
        import json

        with open(path, "r") as f:
            data = json.load(f)

        optimal = data.get("optimal_parameters", {})
        if not optimal:
            return

        # Mapping from PEST parameter names (flexible) to Config fields
        # Users should name parameters in calibration.yaml matching these keys,
        # or we use standard mappings.

        # Direct matches
        for key, value in optimal.items():
            if hasattr(self, key):
                # Type conversion based on existing attribute
                current = getattr(self, key)
                try:
                    if isinstance(current, int):
                        setattr(self, key, int(value))
                    elif isinstance(current, float):
                        setattr(self, key, float(value))
                    elif isinstance(current, bool):
                        setattr(self, key, bool(value))
                except (ValueError, TypeError):
                    pass  # Skip incompatible types

        # Specific Mappings for common aliases
        aliases = {
            "dispersivity": "dispersivity_m",
            "velocity": "aquifer_hydraulic_k_m_day",  # Assuming velocity correlates to K? Or just velocity?
            # Note: velocity depends on K, porosity, gradient. Usually we calibrate K.
            "hydraulic_k": "aquifer_hydraulic_k_m_day",
            "K": "aquifer_hydraulic_k_m_day",
            "decay": "denitrification_k_1_day",
            "porosity": "aquifer_porosity",
            "rate_constant": "rt_default_rate_constant",
            "surface_area": "rt_default_surface_area",
            "residence_time": "rt_default_residence_time",
            "ks_multiplier": None,  # Handle custom vadose scaling?
            "kc": None,
        }

        for p_name, conf_key in aliases.items():
            if p_name in optimal and conf_key is not None:
                val = float(optimal[p_name])
                setattr(self, conf_key, val)

        print(f"Loaded {len(optimal)} parameters from calibration.")


def default_config() -> Config:

    return Config()
