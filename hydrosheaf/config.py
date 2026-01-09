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

    def lambda_l1_value(self) -> float:
        return self.lambda_l1 if self.lambda_l1 else self.lambda_sparse


def default_config() -> Config:
    return Config()
