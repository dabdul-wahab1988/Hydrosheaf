"""
Calibration Configuration Parser.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
import yaml
import pandas as pd

from .definitions import AdjustableParameter, Observation


@dataclass
class CalibrationConfig:
    problem_type: str  # "kinetic", "transport", "vadose", "nitrate", "age", "topology", "composite"
    n_workers: int = 1
    max_nfev: int = 50
    output_dir: str = "calibration_results"
    
    # Engine configuration
    engine: str = "internal"  # "internal" (PESTGLM), "pestpp-glm", "pestpp-ies", "pestpp-sen", "pestpp-swp", "pestpp-mou", "pestpp-opt", "pestpp-da"
    work_dir: Optional[str] = None # Workspace for PEST++ files
    pestpp_version: Optional[str] = None
    pestpp_options: Dict[str, Any] = field(default_factory=dict)
    ies_settings: Dict[str, Any] = field(default_factory=dict)
    sen_settings: Dict[str, Any] = field(default_factory=dict)
    swp_settings: Dict[str, Any] = field(default_factory=dict)
    opt_settings: Dict[str, Any] = field(default_factory=dict)
    da_settings: Dict[str, Any] = field(default_factory=dict)
    loss: str = "linear"  # Robust loss options: linear, huber, soft_l1, cauchy

    # Generic Parameter Definitions (for simple/generic problems)
    parameters: List[AdjustableParameter] = field(default_factory=list)

    # File paths for specific adapters
    observations_file: Optional[str] = None
    model_config_file: Optional[str] = None  # e.g. vadose_profile.json
    
    # Model specific settings (single model)
    adapter_settings: Dict[str, Any] = field(default_factory=dict)
    
    # Composite configuration (list of model definitions)
    sub_models: List[Dict[str, Any]] = field(default_factory=list)


def load_calibration_config(path: str) -> CalibrationConfig:
    with open(path, "r") as f:
        data = yaml.safe_load(f)

    cal_section = data.get("calibration", data)
    settings = cal_section.get("settings", {})

    config = CalibrationConfig(
        problem_type=cal_section.get("type", "generic"),
        n_workers=settings.get("n_workers", 1),
        max_nfev=settings.get("max_iterations", 50),
        output_dir=settings.get("output_dir", "calibration_results"),
        engine=settings.get("engine", "internal"),
        work_dir=settings.get("work_dir"),
        pestpp_version=settings.get("pestpp_version"),
        pestpp_options=settings.get("pestpp_options", {}),
        ies_settings=settings.get("ies", {}),
        sen_settings=settings.get("sen", {}),
        swp_settings=settings.get("swp", {}),
        opt_settings=settings.get("opt", {}),
        da_settings=settings.get("da", {}),
        loss=settings.get("loss", "linear"),
        observations_file=cal_section.get("observations", {}).get("file") if isinstance(cal_section.get("observations"), dict) else None,
        model_config_file=cal_section.get("model", {}).get("config_file") if isinstance(cal_section.get("model"), dict) else None,
        adapter_settings=cal_section.get("model", {}) if isinstance(cal_section.get("model"), dict) else {},
        sub_models=cal_section.get("sub_models", [])
    )

    # Parse Parameters (Global/Generic)
    for p_def in cal_section.get("parameters", []):
        config.parameters.append(
            AdjustableParameter(
                name=p_def["name"],
                value=float(p_def["initial"]),
                lower_bound=float(p_def["bounds"][0]),
                upper_bound=float(p_def["bounds"][1]),
                log_transform=p_def.get("log", False),
                prior_mean=p_def.get("prior_mean"),
                prior_sigma=p_def.get("prior_sigma"),
                fixed=p_def.get("fixed", False),
                tied_to=p_def.get("tied_to", None),
            )
        )

    return config



def load_observations_from_csv(path: str) -> List[Observation]:
    """
    Load observations from a simple CSV: id, value, weight
    """
    df = pd.read_csv(path)
    obs = []
    for _, row in df.iterrows():
        name = str(row.get("id", f"obs_{_}"))
        val = float(row.get("value", 0.0))
        weight = float(row.get("weight", 1.0))
        obs.append(Observation(name, val, weight))
    return obs
