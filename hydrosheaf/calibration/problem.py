"""
Abstract Base Class for Calibration Problems.
"""

from abc import ABC, abstractmethod
from typing import Dict, List
from .definitions import AdjustableParameter, Observation


class CalibrationProblem(ABC):
    """
    Abstract base class for defining a calibration problem.
    Adapts a Hydrosheaf model (or sub-model) to the PESTGLM engine.
    """

    @abstractmethod
    def get_parameters(self) -> List[AdjustableParameter]:
        """Return the list of adjustable parameters."""
        pass

    @abstractmethod
    def get_observations(self) -> List[Observation]:
        """Return the list of observations to match."""
        pass

    @abstractmethod
    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        """
        Run the model with the given parameter values.

        Parameters
        ----------
        param_values : Dict[str, float]
            Dictionary of {parameter_name: value} from the optimizer.

        Returns
        -------
        Dict[str, float]
            Dictionary of {observation_name: simulated_value}.
        """
        pass
