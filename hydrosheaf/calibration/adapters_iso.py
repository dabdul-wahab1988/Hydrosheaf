"""
Calibration Adapter for Water Isotopes (d18O, d2H).
Implements mixing models with softmax parameterization.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional
import numpy as np

from .definitions import AdjustableParameter, Observation
from .problem import CalibrationProblem
from ..log import get_logger

logger = get_logger("calibration.adapters_iso")


@dataclass
class WaterEndmember:
    """Definition of a water source endmember."""
    name: str
    d18O: float
    d2H: float
    d18O_sigma: float = 0.5
    d2H_sigma: float = 2.0


class WaterIsotopeMixingAdapter(CalibrationProblem):
    """
    Calibrates mixing fractions of water sources to match observed d18O and d2H.
    Uses softmax parameterization to ensure fractions are positive and sum to 1.

    Parameters
    ----------
    endmembers : List[WaterEndmember]
        List of potential water sources.
    observations : Dict[str, Dict[str, float]]
        Map of sample_id -> {"18O": val, "2H": val}
    group_map : Dict[str, str]
        Map of sample_id -> group_id. Fractions are estimated per group.
    weights : Dict[str, float]
        Weights for observations (default 1.0/sigma). e.g. {"18O": 5.0, "2H": 1.0}
    """

    def __init__(
        self,
        endmembers: List[WaterEndmember],
        observations: Dict[str, Dict[str, float]],
        group_map: Dict[str, str],
        weights: Optional[Dict[str, float]] = None,
    ):
        self.endmembers = endmembers
        self.observations = observations
        self.group_map = group_map
        self.weights = weights or {"18O": 5.0, "2H": 1.0}

        # Verify groups
        self.groups = sorted(list(set(group_map.values())))

        # We need N endmembers
        self.n_sources = len(endmembers)
        if self.n_sources < 2:
            raise ValueError("At least 2 endmembers required for mixing.")

    def get_parameters(self) -> List[AdjustableParameter]:
        """
        Generate parameters: theta_{group}_{source}
        One theta per source per group.

        To fix identifiability of softmax (since softmax is invariant to adding const),
        we can fix the last theta to 0.0 or just let regularization handle it.
        Fixing one is cleaner. Let's fix the last source's theta to 0.
        """
        params = []

        # For each group
        for group in self.groups:
            # We create params for first N-1 sources
            # The last source (N-th) has implicit theta=0
            for i in range(self.n_sources - 1):
                src = self.endmembers[i]
                p_name = f"theta_{group}_{src.name}"

                # Theta is unbounded in principle, but logits usually -10 to 10
                params.append(
                    AdjustableParameter(
                        name=p_name,
                        value=0.0,  # Start at equal mixing (if all are 0)
                        lower_bound=-10.0,
                        upper_bound=10.0,
                        prior_mean=0.0,
                        prior_sigma=2.0  # Weak prior to keep logits reasonable
                    )
                )
        return params

    def get_observations(self) -> List[Observation]:
        """Generate d18O and d2H observations for all samples."""
        obs_list = []
        for sample_id, vals in self.observations.items():
            if "18O" in vals:
                obs_list.append(
                    Observation(
                        name=f"{sample_id}_18O",
                        value=vals["18O"],
                        weight=self.weights.get("18O", 1.0)
                    )
                )
            if "2H" in vals:
                obs_list.append(
                    Observation(
                        name=f"{sample_id}_2H",
                        value=vals["2H"],
                        weight=self.weights.get("2H", 1.0)
                    )
                )
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        """
        Predict isotopes based on mixing fractions derived from thetas.
        """
        results = {}

        # 1. Resolve fractions for each group
        group_fractions = {}  # group -> [f1, f2, ..., fN]

        for group in self.groups:
            # Gather logits
            logits = []
            for i in range(self.n_sources - 1):
                src_name = self.endmembers[i].name
                p_name = f"theta_{group}_{src_name}"
                val = param_values.get(p_name, 0.0)
                logits.append(val)

            # Last source has logit 0
            logits.append(0.0)

            # Softmax
            exp_logits = np.exp(np.array(logits))
            sum_exp = np.sum(exp_logits)
            fractions = exp_logits / sum_exp

            group_fractions[group] = fractions

        # 2. Predict for each sample
        # Pre-compute endmember vectors
        em_d18O = np.array([e.d18O for e in self.endmembers])
        em_d2H = np.array([e.d2H for e in self.endmembers])

        for sample_id, _ in self.observations.items():
            group = self.group_map.get(sample_id)
            if not group:
                logger.warning(f"Sample {sample_id} not found in group_map. Skipping prediction.")
                continue

            # Ensure group is a string key for lookup (in case group_map has mixed types or unexpected)
            # though type hint says Dict[str, str], robust runtime check is safer if needed.
            # But we trust the type hint here for now.

            if group not in group_fractions:
                # Should not happen if self.groups is built from group_map values
                continue

            fracs = group_fractions[group]

            pred_18O = np.dot(fracs, em_d18O)
            pred_2H = np.dot(fracs, em_d2H)

            results[f"{sample_id}_18O"] = float(pred_18O)
            results[f"{sample_id}_2H"] = float(pred_2H)

        return results
