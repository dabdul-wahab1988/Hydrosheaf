"""Configuration defaults for hydrosheaf."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DEFAULT_ION_ORDER = ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]


def default_phreeqc_database() -> str:
    return str(Path(__file__).resolve().parent / "databases" / "phreeqc.dat")


@dataclass
class Config:
    ion_order: List[str] = field(default_factory=lambda: DEFAULT_ION_ORDER.copy())
    unit: str = "mmol/L"
    unit_mode: str = "mmol_L"
    weights: List[float] = field(default_factory=lambda: [1.0] * 10)
    lambda_sparse: float = 0.0
    lambda_l1: float = 0.0
    reaction_max_iter: int = 300
    reaction_tol: float = 1e-6
    charge_balance_limit: float = 0.1
    ec_tds_penalty_limit: float = 0.0
    ec_tds_penalty_enabled: bool = False
    missing_policy: str = "skip"
    detection_limit_policy: str = "half"
    eta_ec: float = 0.0
    eta_tds: float = 0.0
    allow_signed_reactions: bool = False
    signed_reaction_labels: List[str] = field(default_factory=list)
    denit_kappa: float = 1.0
    transport_models_enabled: List[str] = field(default_factory=lambda: ["evap", "mix"])
    mixing_endmembers: Dict[str, List[float]] = field(default_factory=dict)
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
    edge_sigma_meas: float = 0.5
    edge_sigma_dtw: float = 1.0
    edge_sigma_elev: float = 1.0
    edge_sigma_topo: float = 10.0
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
    residence_time_method: str = "cross_correlation"  # gradient, cross_correlation, tracer_decay
    residence_time_tracer: str = "Cl"  # conservative tracer for age estimation
    residence_time_hydraulic_k: float = 1.0  # m/day, for gradient method
    residence_time_porosity: float = 0.2  # effective porosity

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

    # Screen interval
    screen_top_key: str = "screen_top"
    screen_bottom_key: str = "screen_bottom"
    screen_overlap_threshold: float = 5.0  # meters

    def validate(self) -> None:
        if len(self.ion_order) != 10:
            raise ValueError("ion_order must have 10 entries.")
        if self.unit_mode not in {"mmol_L", "meq_L"}:
            raise ValueError("unit_mode must be 'mmol_L' or 'meq_L'.")
        if len(self.weights) != len(self.ion_order):
            raise ValueError("weights must match ion_order length.")
        if any(w < 0 for w in self.weights):
            raise ValueError("weights must be non-negative.")
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
            raise ValueError("transport_models_enabled must be a subset of {'evap','mix'}.")
        for name, vector in self.mixing_endmembers.items():
            if len(vector) != len(self.ion_order):
                raise ValueError(f"endmember '{name}' has invalid length.")

        if self.nitrate_source_prior < 0 or self.nitrate_source_prior > 1:
            raise ValueError("nitrate_source_prior must be between 0 and 1.")
        if any(w < 0 for w in self.nitrate_source_weights.values()):
            raise ValueError("nitrate_source_weights must be non-negative.")

        # Uncertainty quantification validation
        if self.uncertainty_method not in {"none", "bootstrap", "bayesian", "monte_carlo"}:
            raise ValueError("uncertainty_method must be one of: none, bootstrap, bayesian, monte_carlo.")
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
            raise ValueError("temporal_interpolation_method must be one of: linear, spline, nearest.")
        if self.temporal_frequency_days < 1:
            raise ValueError("temporal_frequency_days must be at least 1.")
        if self.residence_time_method not in {"gradient", "cross_correlation", "tracer_decay"}:
            raise ValueError("residence_time_method must be one of: gradient, cross_correlation, tracer_decay.")
        if self.residence_time_hydraulic_k <= 0:
            raise ValueError("residence_time_hydraulic_k must be positive.")
        if not 0.0 < self.residence_time_porosity <= 1.0:
            raise ValueError("residence_time_porosity must be between 0 and 1.")

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
            if len(self.layer_names) != len(self.layer_tops) or len(self.layer_names) != len(self.layer_bottoms):
                raise ValueError("layer_names, layer_tops, and layer_bottoms must have same length.")
            for i in range(len(self.layer_tops)):
                if self.layer_tops[i] >= self.layer_bottoms[i]:
                    raise ValueError(f"layer_tops[{i}] must be less than layer_bottoms[{i}].")

    def lambda_l1_value(self) -> float:
        return self.lambda_l1 if self.lambda_l1 else self.lambda_sparse


def default_config() -> Config:
    return Config()
