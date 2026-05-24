"""
Calibration Adapters for different Hydrosheaf modules.
"""

from dataclasses import dataclass
import math
import re
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union, cast
import numpy as np

from .definitions import AdjustableParameter, Observation
from .problem import CalibrationProblem
from ..reactive_transport import KineticParameters
from ..reactive_transport.kinetic_phreeqc import (
    run_phreeqc_kinetic,
    build_kinetic_block,
)
from ..data.units import mgL_to_mmolL
from ..log import get_logger

logger = get_logger("calibration.adapters")


@dataclass
class KineticExperiment:
    """
    Data for a single kinetic experiment (e.g. one batch reactor or one flow path).
    """

    id: str
    initial_solution: Dict[str, float]
    residence_time_days: float
    reaction_extents: List[float]  # Estimated extents from inverse model (or 0 if calibration drives reaction)
    reaction_labels: List[str]  # Names of reactions corresponding to extents
    observations: Dict[str, float]  # {ion_name: concentration}
    temperature_c: float = 25.0
    weight: float = 1.0
    # Hydrosheaf internal convention is mmol/L. Use mg/L explicitly when providing
    # field measurements directly.
    units: str = "mmol/L"  # "mg/L" or "mmol/L"
    edge_id: Optional[str] = None  # Optional edge identifier for parameter scoping
    default_extent_mmol_L: Optional[float] = None  # If set, activates reactions with missing/zero extents
    geological_layer: Optional[str] = None
    aquifer_unit: Optional[str] = None
    site_group: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class KineticCalibrationAdapter(CalibrationProblem):
    """
    Adapts the Kinetic PHREEQC module for PEST calibration.
    Calibrates rate constants and surface areas.
    """

    def __init__(
        self,
        base_params: Dict[str, KineticParameters],
        experiments: List[KineticExperiment],
        config: Any,  # Hydrosheaf Config object
        params_to_fit: Optional[List[str]] = None,  # List of "reaction:param" or "reaction:param:scope"
        strict_units: bool = True,
    ):
        self.base_params = base_params
        self.experiments = experiments
        self.config = config
        self.strict_units = strict_units

        # Validate and convert units immediately
        self._normalize_experiments()

        # Default to fitting k for all reactions in base_params if not specified
        if params_to_fit is None:
            self.params_to_fit = [f"{r}:k" for r in base_params.keys()]
        else:
            self.params_to_fit = params_to_fit

    def _normalize_experiments(self):
        """Convert all mg/L inputs to mmol/L for internal consistency."""
        for exp in self.experiments:
            if exp.units == "mg/L":
                # Convert initial solution
                new_init = {}
                for ion, val in exp.initial_solution.items():
                    if ion in ["pH", "temp_c", "pe", "temperature"]:
                        new_init[ion] = val
                    else:
                        try:
                            # Use mgL_to_mmolL
                            new_init[ion] = mgL_to_mmolL(val, ion)
                        except KeyError:
                            if self.strict_units:
                                raise ValueError(
                                    f"Unknown species '{ion}' in initial_solution with units mg/L. "
                                    "Either add a molar mass mapping (hydrosheaf/data/units.py) "
                                    "or provide this experiment in mmol/L."
                                )
                            new_init[ion] = val
                exp.initial_solution = new_init

                # Convert observations
                new_obs = {}
                for ion, val in exp.observations.items():
                    try:
                        new_obs[ion] = mgL_to_mmolL(val, ion)
                    except KeyError:
                        if self.strict_units:
                            raise ValueError(
                                f"Unknown species '{ion}' in observations with units mg/L. "
                                "Either add a molar mass mapping (hydrosheaf/data/units.py) "
                                "or provide this experiment in mmol/L."
                            )
                        new_obs[ion] = val
                exp.observations = new_obs
                
                exp.units = "mmol/L"

    @staticmethod
    def _safe_scope_id(value: Any) -> str:
        return re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_") or "unknown"

    @staticmethod
    def _normalize_param_type(p_type: str) -> str:
        aliases = {
            "rate": "k",
            "rate_constant": "k",
            "surface_area": "A",
            "area": "A",
            "extent": "extent",
            "xi": "extent",
        }
        return aliases.get(p_type.strip(), p_type.strip())

    def _scope_values(self, scope: str) -> List[str]:
        scope = scope.strip()
        if scope == "global":
            return [""]
        if scope == "per_edge":
            return sorted(list(set(
                self._safe_scope_id(exp.edge_id)
                for exp in self.experiments
                if exp.edge_id
            )))
        if scope in {"per_layer", "per_geological_layer"}:
            return sorted(list(set(
                self._safe_scope_id(exp.geological_layer)
                for exp in self.experiments
                if exp.geological_layer
            )))
        if scope in {"per_aquifer", "per_aquifer_unit"}:
            return sorted(list(set(
                self._safe_scope_id(exp.aquifer_unit)
                for exp in self.experiments
                if exp.aquifer_unit
            )))
        if scope in {"per_group", "per_site_group"}:
            return sorted(list(set(
                self._safe_scope_id(exp.site_group)
                for exp in self.experiments
                if exp.site_group
            )))
        return [self._safe_scope_id(scope)]

    def _parameter_name(self, reaction: str, p_type: str, scope_value: str = "") -> str:
        prefix = f"extent_{reaction}" if p_type == "extent" else f"log_{p_type}_{reaction}"
        return f"{prefix}__{scope_value}" if scope_value else prefix

    def get_parameters(self) -> List[AdjustableParameter]:
        """
        Generate PEST parameters based on requested fit configuration.
        Supports global, edge-specific, and simple hierarchical scoping.
        """
        pest_params = []
        seen_params = set()

        for p_spec in self.params_to_fit:
            parts = p_spec.split(":")
            if len(parts) < 2:
                raise ValueError(
                    f"Invalid kinetic calibration parameter spec '{p_spec}'. "
                    "Expected 'reaction:parameter[:scope]'."
                )
            reaction = parts[0]
            p_type = self._normalize_param_type(parts[1])
            scope = parts[2] if len(parts) > 2 else "global"

            if reaction not in self.base_params and p_type != "extent":
                continue

            scope_values = self._scope_values(scope)
            if scope != "global" and not scope_values:
                raise ValueError(
                    f"No experiments contain metadata for kinetic calibration scope "
                    f"'{scope}' in '{p_spec}'."
                )

            for scope_value in scope_values:
                param_name = self._parameter_name(reaction, p_type, scope_value)
                if param_name in seen_params:
                    continue
                seen_params.add(param_name)

                if p_type == "extent":
                    val = self._default_extent_for_reaction(reaction)
                else:
                    base_p = self.base_params[reaction]
                    val = base_p.rate_constant if p_type == "k" else base_p.surface_area
                pest_params.append(self._make_param(param_name, val, p_type))

        return pest_params

    def _make_param(self, name: str, val: float, p_type: str) -> AdjustableParameter:
        p_type = self._normalize_param_type(p_type)
        if p_type == "k":
            return AdjustableParameter(
                name=name,
                value=val,
                lower_bound=val * 1e-4,
                upper_bound=val * 1e4,
                log_transform=True,
                prior_mean=val,
                prior_sigma=1.0,
            )
        if p_type == "extent":
            start = max(float(val), 1e-9)
            upper = max(start * 100.0, self._max_configured_extent() * 10.0, 1.0)
            return AdjustableParameter(
                name=name,
                value=start,
                lower_bound=0.0,
                upper_bound=upper,
                log_transform=False,
                prior_mean=start,
                prior_sigma=max(start, 1.0),
            )
        else: # A
            return AdjustableParameter(
                name=name,
                value=val,
                lower_bound=val * 0.1,
                upper_bound=val * 10.0,
                log_transform=True,
                prior_mean=val,
                prior_sigma=0.5,
            )

    def _default_extent_for_reaction(self, reaction: str) -> float:
        values: List[float] = []
        for exp in self.experiments:
            for label, extent in zip(exp.reaction_labels, exp.reaction_extents):
                if label == reaction and float(extent) > 0.0:
                    values.append(float(extent))
            if exp.default_extent_mmol_L is not None and reaction in exp.reaction_labels:
                values.append(float(exp.default_extent_mmol_L))
        return float(np.median(values)) if values else 1.0

    def _max_configured_extent(self) -> float:
        values: List[float] = []
        for exp in self.experiments:
            values.extend(abs(float(v)) for v in exp.reaction_extents)
            if exp.default_extent_mmol_L is not None:
                values.append(abs(float(exp.default_extent_mmol_L)))
        return max(values) if values else 1.0

    def get_observations(self) -> List[Observation]:
        """Generate PEST observations from experiments."""
        obs_list = []
        for exp in self.experiments:
            for ion, val in exp.observations.items():
                obs_name = f"{exp.id}_{ion}"
                obs_list.append(
                    Observation(name=obs_name, value=val, weight=exp.weight)
                )
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        """
        Run PHREEQC for all experiments with updated kinetic parameters.
        """
        # 1. Run Experiments
        results = {}

        for exp in self.experiments:
            # Resolve parameters for this experiment
            current_params = self._resolve_parameters(exp, param_values)

            if len(exp.reaction_labels) != len(exp.reaction_extents):
                raise ValueError(
                    f"KineticExperiment '{exp.id}' has {len(exp.reaction_labels)} reaction_labels but "
                    f"{len(exp.reaction_extents)} reaction_extents. These must match."
                )
            
            # Ensure dissolution happens even if extents are 0
            # If extents are 0, build_kinetic_block might skip them.
            # We force them to have a small positive extent if they are being calibrated
            # OR we pass a flag to build_kinetic_block to force inclusion.
            # Let's modify reaction_extents temporarily if they are 0 but reaction is active.
            
            active_extents = self._resolve_extents(exp, param_values)
            if exp.default_extent_mmol_L is not None and exp.default_extent_mmol_L > 0:
                for i, label in enumerate(exp.reaction_labels):
                    if i >= len(active_extents):
                        continue
                    if active_extents[i] != 0.0:
                        continue
                    if not self._reaction_is_calibrated(label, exp, param_values):
                        continue
                    active_extents[i] = float(exp.default_extent_mmol_L)
            
            # Build input
            kinetics_block = build_kinetic_block(
                reaction_labels=exp.reaction_labels,
                extents=active_extents,
                residence_time_days=exp.residence_time_days,
                kinetic_params=current_params,
                temperature_c=exp.temperature_c,
            )

            # Run
            out = run_phreeqc_kinetic(
                initial_solution=exp.initial_solution,
                kinetics_block=kinetics_block,
                residence_time_days=exp.residence_time_days,
                config=self.config,
            )

            if out["success"]:
                # Map outputs to observations
                # final_composition is a list in ion_order
                ion_order = self.config.ion_order
                # Explicit cast to List[float] to satisfy mypy
                raw_vals = out.get("final_composition", [])
                if isinstance(raw_vals, list):
                    final_vals = cast(List[float], raw_vals)
                else:
                    final_vals = []

                # Create a dict for easy lookup
                sim_map = {}
                for i, ion in enumerate(ion_order):
                    if i < len(final_vals):
                        sim_map[ion] = final_vals[i]

                # Add to results
                for ion in exp.observations:
                    obs_name = f"{exp.id}_{ion}"
                    if ion in sim_map:
                        results[obs_name] = sim_map[ion]
                    else:
                        results[obs_name] = -9999.0
                        logger.warning(
                            f"Missing simulation result for {ion} in experiment {exp.id}"
                        )
            else:
                # Failure penalty
                logger.error(
                    f"PHREEQC run failed for {exp.id}: {out.get('error_message')}"
                )
                for ion in exp.observations:
                    results[f"{exp.id}_{ion}"] = -1e6

        return results

    def _reaction_is_calibrated(
        self, reaction: str, exp: KineticExperiment, param_values: Dict[str, float]
    ) -> bool:
        if f"log_k_{reaction}" in param_values or f"log_A_{reaction}" in param_values:
            return True
        if f"extent_{reaction}" in param_values:
            return True
        if exp.edge_id:
            if self._has_scoped_reaction_param(
                reaction, self._safe_scope_id(exp.edge_id), param_values
            ):
                return True
        for value in (exp.geological_layer, exp.aquifer_unit, exp.site_group):
            if value and self._has_scoped_reaction_param(
                reaction, self._safe_scope_id(value), param_values
            ):
                return True
        if exp.metadata:
            for value in exp.metadata.values():
                if value and self._has_scoped_reaction_param(
                    reaction, self._safe_scope_id(value), param_values
                ):
                    return True
        return False

    def _has_scoped_reaction_param(
        self, reaction: str, scope_value: str, param_values: Dict[str, float]
    ) -> bool:
        return any(
            name in param_values
            for name in (
                f"log_k_{reaction}__{scope_value}",
                f"log_A_{reaction}__{scope_value}",
                f"extent_{reaction}__{scope_value}",
            )
        )

    def _resolve_scoped_value(
        self,
        exp: KineticExperiment,
        reaction: str,
        p_type: str,
        param_values: Dict[str, float],
        default: float,
    ) -> float:
        prefix = f"extent_{reaction}" if p_type == "extent" else f"log_{p_type}_{reaction}"
        value = float(param_values.get(prefix, default))
        suffixes: List[Any] = [
            exp.edge_id,
            exp.geological_layer,
            exp.aquifer_unit,
            exp.site_group,
        ]
        if exp.metadata:
            suffixes.extend(v for v in exp.metadata.values() if v not in (None, ""))

        for suffix in suffixes:
            if not suffix:
                continue
            param_name = f"{prefix}__{self._safe_scope_id(suffix)}"
            if param_name in param_values:
                value = float(param_values[param_name])
        return value

    def _resolve_extents(
        self, exp: KineticExperiment, param_values: Dict[str, float]
    ) -> List[float]:
        active_extents = list(exp.reaction_extents)
        for i, label in enumerate(exp.reaction_labels):
            if i >= len(active_extents):
                continue
            active_extents[i] = self._resolve_scoped_value(
                exp, label, "extent", param_values, float(active_extents[i])
            )
        return active_extents

    def _resolve_parameters(self, exp: KineticExperiment, param_values: Dict[str, float]) -> Dict[str, KineticParameters]:
        """
        Resolve final KineticParameters for an experiment based on global/edge scopes.
        """
        resolved = {}
        for name, base_p in self.base_params.items():
            # Start with base
            k = base_p.rate_constant
            A = base_p.surface_area
            
            # NOTE: parameter names contain "log_" for readability, but param_values are
            # real-space values after transform handling in the calibration engine.
            k = self._resolve_scoped_value(exp, name, "k", param_values, k)
            A = self._resolve_scoped_value(exp, name, "A", param_values, A)
            
            # Create new params object
            from ..reactive_transport import KineticParameters as KP
            resolved[name] = KP(
                reaction_name=base_p.reaction_name,
                rate_constant=k,
                surface_area=A,
                reaction_order=base_p.reaction_order,
                activation_energy=base_p.activation_energy,
                reference_temp_k=base_p.reference_temp_k
            )
        return resolved


class CompositeCalibrationAdapter(CalibrationProblem):
    """
    Combines multiple sub-problems into a single calibration problem.
    Allows simultaneous calibration of different model components (e.g. Vadose + Transport)
    that may share parameters or be independent.
    """

    def __init__(self, sub_problems: List[CalibrationProblem]):
        self.sub_problems = sub_problems
        self._parameter_map = {} # param_name -> [problem_indices]
        self._parameters = []
        self._observations = []
        
        self._initialize()

    def _initialize(self):
        """Aggregate parameters and observations."""
        seen_params = {}
        
        for i, prob in enumerate(self.sub_problems):
            # Parameters
            params = prob.get_parameters()
            for p in params:
                if p.name not in seen_params:
                    seen_params[p.name] = p
                    self._parameters.append(p)
                    self._parameter_map[p.name] = [i]
                else:
                    # Parameter shared across models
                    # Verify consistency? (bounds etc). For now assume first def wins.
                    self._parameter_map[p.name].append(i)
            
            # Observations
            # We assume observation names are unique globally.
            # If not, we might need to prefix them.
            obs = prob.get_observations()
            self._observations.extend(obs)

    def get_parameters(self) -> List[AdjustableParameter]:
        return self._parameters

    def get_observations(self) -> List[Observation]:
        return self._observations

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        """
        Run all sub-models and aggregate results.
        """
        aggregated_results = {}
        
        # In serial execution (Python side), we just loop.
        # This might be slow if models are heavy. 
        # But if running inside PEST++ parallel agents, each agent runs this sequentially.
        
        for prob in self.sub_problems:
            # Filter params for this problem?
            # Or just pass all? Sub-problems usually ignore unknown kwargs/dict keys
            # if using .get(). If strict, we need to filter.
            # Our adapters (Generic, Kinetic, Transport) use .get() or loop over expected.
            # So passing the full dict is safe and easier.
            
            # Run sub-model
            res = prob.run_model(param_values)
            aggregated_results.update(res)
            
        return aggregated_results


class GenericFunctionAdapter(CalibrationProblem):
    """
    Adapts a generic Python function for testing PEST logic.
    """

    def __init__(
        self,
        func: Callable[[Dict[str, float]], Dict[str, float]],
        parameters: List[AdjustableParameter],
        observations: List[Observation],
    ):
        self.func = func
        self._parameters = parameters
        self._observations = observations

    def get_parameters(self) -> List[AdjustableParameter]:
        return self._parameters

    def get_observations(self) -> List[Observation]:
        return self._observations

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        return self.func(param_values)


# --- New Adapters for Vadose and Transport ---


@dataclass
class TransportExperiment:
    """
    Data for a 1D transport experiment (column test or field tracer test).
    """

    id: str
    times: List[float]  # days
    observed_concentrations: List[float]  # mg/L or mmol/L
    distance_m: float  # distance from source
    weights: Optional[List[float]] = None  # optional weights per point

    # Fixed parameters (not calibrated)
    velocity_m_day: float = 0.1
    source_concentration: float = 1.0


class TransportCalibrationAdapter(CalibrationProblem):
    """
    Adapts 1D Transport models (Analytical or FloPy) for PEST calibration.
    Calibrates: dispersivity, decay_rate, effective_porosity, retardation.
    """

    def __init__(
        self,
        experiments: List[TransportExperiment],
        use_analytical: bool = True,
        params_to_fit: List[str] = [
            "dispersivity",
            "decay",
        ],  # ["dispersivity", "decay", "velocity"]
        base_dispersivity: float = 1.0,
        base_decay: float = 0.0,
        base_velocity: float = 0.1,
    ):
        self.experiments = experiments
        self.use_analytical = use_analytical
        self.params_to_fit = params_to_fit

        self.base_vals = {
            "dispersivity": base_dispersivity,
            "decay": base_decay,
            "velocity": base_velocity,
        }

    def get_parameters(self) -> List[AdjustableParameter]:
        pest_params = []

        if "dispersivity" in self.params_to_fit:
            val = self.base_vals["dispersivity"]
            pest_params.append(
                AdjustableParameter(
                    "dispersivity",
                    val,
                    0.01,
                    100.0,
                    log_transform=True,
                    prior_mean=val,
                    prior_sigma=2.0,
                )
            )

        if "decay" in self.params_to_fit:
            val = self.base_vals["decay"]
            # Decay can be 0, so maybe not log-transform if it can be 0.
            # But usually we fit k > 0. If base is 0, start at 1e-5.
            start = val if val > 0 else 1e-4
            pest_params.append(
                AdjustableParameter(
                    "decay",
                    start,
                    1e-6,
                    1.0,
                    log_transform=True,
                    prior_mean=start,
                    prior_sigma=2.0,
                )
            )

        if "velocity" in self.params_to_fit:
            val = self.base_vals["velocity"]
            pest_params.append(
                AdjustableParameter(
                    "velocity", val, val * 0.1, val * 10.0, log_transform=True
                )
            )

        return pest_params

    def get_observations(self) -> List[Observation]:
        obs_list = []
        for exp in self.experiments:
            for i, (t, val) in enumerate(zip(exp.times, exp.observed_concentrations)):
                w = exp.weights[i] if exp.weights else 1.0
                obs_list.append(Observation(f"{exp.id}_{i}", val, weight=w))
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        # Update parameters
        disp = param_values.get("dispersivity", self.base_vals["dispersivity"])
        decay = param_values.get("decay", self.base_vals["decay"])
        # Global velocity override? Or per-experiment?
        # If "velocity" is calibrated, it overrides experiment defaults.
        global_vel = param_values.get("velocity", None)

        from ..transport.flopy_1d import run_analytical_1d_transport

        results = {}
        for exp in self.experiments:
            vel = global_vel if global_vel is not None else exp.velocity_m_day

            # Run simulation
            if self.use_analytical:
                # We use the analytical runner which returns TransportResult
                res = run_analytical_1d_transport(
                    times=np.array(exp.times),
                    distance_m=exp.distance_m,
                    velocity_m_day=vel,
                    dispersivity_m=disp,
                    decay_rate_1_day=decay,
                    source_concentration=exp.source_concentration,
                )

                # Map to observations
                # res.concentrations corresponds to exp.times
                for i, val in enumerate(res.concentrations):
                    results[f"{exp.id}_{i}"] = val
            else:
                # FloPy runner logic
                from ..transport.flopy_1d import build_1d_transport_model, run_1d_transport
                
                # Build model
                mf, mt, model_params = build_1d_transport_model(
                    model_name=f"trans_{exp.id}",
                    # Use temp workspace managed by runner or tempfile inside build
                    workspace=None, 
                    aquifer_length_m=exp.distance_m,
                    # Assume unit thickness/width if not specified in experiment
                    # Ideally experiment would hold geometry info. 
                    # For 1D calibration, we assume standard column or flowpath.
                    dispersivity_m=disp,
                    decay_rate_1_day=decay,
                    velocity_m_day=vel, # Used to calc gradient/K inside builder
                    source_concentration=exp.source_concentration,
                    # Time discretization: match experiment times?
                    # FloPy model needs max time.
                    perlen_days=max(exp.times) * 1.1 if exp.times else 10.0,
                    n_time_steps=50 # Fixed resolution
                )
                
                # Run
                res = run_1d_transport(mf, mt, model_params, cleanup=True)
                
                if res.success:
                    # Interpolate to observation times
                    # res.times vs exp.times
                    if len(res.times) > 1:
                        interp_vals = np.interp(exp.times, res.times, res.concentrations)
                        for i, val in enumerate(interp_vals):
                            results[f"{exp.id}_{i}"] = val
                    else:
                        # Fallback
                        for i in range(len(exp.times)):
                            results[f"{exp.id}_{i}"] = res.concentrations[0]
                else:
                    for i in range(len(exp.times)):
                        results[f"{exp.id}_{i}"] = -999.0 # Penalty
class GenericFunctionAdapter(CalibrationProblem):
    """
    Adapts a generic Python function for testing PEST logic.
    """

    def __init__(
        self,
        func: Callable[[Dict[str, float]], Dict[str, float]],
        parameters: List[AdjustableParameter],
        observations: List[Observation],
    ):
        self.func = func
        self._parameters = parameters
        self._observations = observations

    def get_parameters(self) -> List[AdjustableParameter]:
        return self._parameters

    def get_observations(self) -> List[Observation]:
        return self._observations

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        return self.func(param_values)


# --- New Adapters for Vadose and Transport ---


@dataclass
class TransportExperiment:
    """
    Data for a 1D transport experiment (column test or field tracer test).
    """

    id: str
    times: List[float]  # days
    observed_concentrations: List[float]  # mg/L or mmol/L
    distance_m: float  # distance from source
    weights: Optional[List[float]] = None  # optional weights per point

    # Fixed parameters (not calibrated)
    velocity_m_day: float = 0.1
    source_concentration: float = 1.0


class TransportCalibrationAdapter(CalibrationProblem):
    """
    Adapts 1D Transport models (Analytical or FloPy) for PEST calibration.
    Calibrates: dispersivity, decay_rate, effective_porosity, retardation.
    """

    def __init__(
        self,
        experiments: List[TransportExperiment],
        use_analytical: bool = True,
        params_to_fit: List[str] = [
            "dispersivity",
            "decay",
        ],  # ["dispersivity", "decay", "velocity"]
        base_dispersivity: float = 1.0,
        base_decay: float = 0.0,
        base_velocity: float = 0.1,
    ):
        self.experiments = experiments
        self.use_analytical = use_analytical
        self.params_to_fit = params_to_fit

        self.base_vals = {
            "dispersivity": base_dispersivity,
            "decay": base_decay,
            "velocity": base_velocity,
        }

    def get_parameters(self) -> List[AdjustableParameter]:
        pest_params = []

        if "dispersivity" in self.params_to_fit:
            val = self.base_vals["dispersivity"]
            pest_params.append(
                AdjustableParameter(
                    "dispersivity",
                    val,
                    0.01,
                    100.0,
                    log_transform=True,
                    prior_mean=val,
                    prior_sigma=2.0,
                )
            )

        if "decay" in self.params_to_fit:
            val = self.base_vals["decay"]
            # Decay can be 0, so maybe not log-transform if it can be 0.
            # But usually we fit k > 0. If base is 0, start at 1e-5.
            start = val if val > 0 else 1e-4
            pest_params.append(
                AdjustableParameter(
                    "decay",
                    start,
                    1e-6,
                    1.0,
                    log_transform=True,
                    prior_mean=start,
                    prior_sigma=2.0,
                )
            )

        if "velocity" in self.params_to_fit:
            val = self.base_vals["velocity"]
            pest_params.append(
                AdjustableParameter(
                    "velocity", val, val * 0.1, val * 10.0, log_transform=True
                )
            )

        return pest_params

    def get_observations(self) -> List[Observation]:
        obs_list = []
        for exp in self.experiments:
            for i, (t, val) in enumerate(zip(exp.times, exp.observed_concentrations)):
                w = exp.weights[i] if exp.weights else 1.0
                obs_list.append(Observation(f"{exp.id}_{i}", val, weight=w))
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        # Update parameters
        disp = param_values.get("dispersivity", self.base_vals["dispersivity"])
        decay = param_values.get("decay", self.base_vals["decay"])
        # Global velocity override? Or per-experiment?
        # If "velocity" is calibrated, it overrides experiment defaults.
        global_vel = param_values.get("velocity", None)

        from ..transport.flopy_1d import run_analytical_1d_transport

        results = {}
        for exp in self.experiments:
            vel = global_vel if global_vel is not None else exp.velocity_m_day

            # Run simulation
            if self.use_analytical:
                # We use the analytical runner which returns TransportResult
                res = run_analytical_1d_transport(
                    times=np.array(exp.times),
                    distance_m=exp.distance_m,
                    velocity_m_day=vel,
                    dispersivity_m=disp,
                    decay_rate_1_day=decay,
                    source_concentration=exp.source_concentration,
                )

                # Map to observations
                # res.concentrations corresponds to exp.times
                for i, val in enumerate(res.concentrations):
                    results[f"{exp.id}_{i}"] = val
            else:
                # FloPy runner logic
                from ..transport.flopy_1d import build_1d_transport_model, run_1d_transport
                
                # Build model
                mf, mt, model_params = build_1d_transport_model(
                    model_name=f"trans_{exp.id}",
                    # Use temp workspace managed by runner or tempfile inside build
                    workspace=None, 
                    aquifer_length_m=exp.distance_m,
                    # Assume unit thickness/width if not specified in experiment
                    # Ideally experiment would hold geometry info. 
                    # For 1D calibration, we assume standard column or flowpath.
                    dispersivity_m=disp,
                    decay_rate_1_day=decay,
                    velocity_m_day=vel, # Used to calc gradient/K inside builder
                    source_concentration=exp.source_concentration,
                    # Time discretization: match experiment times?
                    # FloPy model needs max time.
                    perlen_days=max(exp.times) * 1.1 if exp.times else 10.0,
                    n_time_steps=50 # Fixed resolution
                )
                
                # Run
                res = run_1d_transport(mf, mt, model_params, cleanup=True)
                
                if res.success:
                    # Interpolate to observation times
                    # res.times vs exp.times
                    if len(res.times) > 1:
                        interp_vals = np.interp(exp.times, res.times, res.concentrations)
                        for i, val in enumerate(interp_vals):
                            results[f"{exp.id}_{i}"] = val
                    else:
                        # Fallback
                        for i in range(len(exp.times)):
                            results[f"{exp.id}_{i}"] = res.concentrations[0]
                else:
                    for i in range(len(exp.times)):
                        results[f"{exp.id}_{i}"] = -999.0 # Penalty

        return results


class VadoseCalibrationAdapter(CalibrationProblem):
    """
    Adapts Vadose Zone Richards Model for PEST calibration.
    Calibrates: saturated conductivity (Ks), alpha, n, kc, preferential flow, and TTD CV.
    """

    def __init__(
        self,
        profile: Any,  # VadoseProfile
        forcing: Any,  # List[VadoseForcingSample]
        observations: Union[Dict[str, Any], List[Dict[str, Any]]],  # {time_idx: {"theta": val, "depth": d}} OR List of dicts
        config: Any = None,
        layers_to_fit: List[int] = [0],  # Indices of layers to calibrate
        params_to_fit: Optional[List[str]] = None,
    ):
        self.profile = profile
        self.forcing = forcing
        self.observations = observations
        self.config = config
        self.layers_to_fit = layers_to_fit
        self.params_to_fit = params_to_fit or []

    def get_parameters(self) -> List[AdjustableParameter]:
        pest_params = []
        fit_set = set(self.params_to_fit) if self.params_to_fit else set()

        for layer_idx in self.layers_to_fit:
            layer = self.profile.layers[layer_idx]
            if layer.ks_m_day is not None:
                ks = float(layer.ks_m_day)
                if not fit_set or f"ks_L{layer_idx}" in fit_set:
                    pest_params.append(
                        AdjustableParameter(
                            f"ks_L{layer_idx}", ks, ks * 1e-2, ks * 1e2, log_transform=True,
                            description=f"Saturated hydraulic conductivity for layer {layer_idx}"
                        )
                    )
            if layer.alpha_1_m is not None:
                alpha = float(layer.alpha_1_m)
                if not fit_set or f"alpha_L{layer_idx}" in fit_set:
                    pest_params.append(
                        AdjustableParameter(
                            f"alpha_L{layer_idx}", alpha, alpha * 0.1, alpha * 10.0, log_transform=True,
                            description=f"van Genuchten alpha for layer {layer_idx}"
                        )
                    )
            if layer.n is not None:
                n_val = float(layer.n)
                if f"n_L{layer_idx}" in fit_set:
                    pest_params.append(
                        AdjustableParameter(
                            f"n_L{layer_idx}", n_val, 1.01, 10.0, log_transform=False,
                            description=f"van Genuchten n for layer {layer_idx}"
                        )
                    )

        if "kc" in fit_set or "Kc" in fit_set:
            kc_val = float(self.config.kc if self.config and hasattr(self.config, "kc") else 1.0)
            pest_params.append(
                AdjustableParameter(
                    "kc", kc_val, 0.1, 5.0, log_transform=True,
                    description="Evapotranspiration crop coefficient Kc"
                )
            )

        if "preferential_flow_fraction" in fit_set or "pref_flow_frac" in fit_set:
            pff = float(self.config.preferential_flow_fraction if self.config and hasattr(self.config, "preferential_flow_fraction") else 0.0)
            pest_params.append(
                AdjustableParameter(
                    "preferential_flow_fraction", max(pff, 0.0), 0.0, 1.0, log_transform=False,
                    description="Preferential flow fraction"
                )
            )

        if "ttd_cv" in fit_set or "ttd_default_cv" in fit_set:
            cv = float(self.config.ttd_default_cv if self.config and hasattr(self.config, "ttd_default_cv") else 0.7)
            pest_params.append(
                AdjustableParameter(
                    "ttd_cv", cv, 0.01, 5.0, log_transform=True,
                    description="TTD coefficient of variation"
                )
            )

        return pest_params

    def get_observations(self) -> List[Observation]:
        obs_list = []
        if isinstance(self.observations, list):
            for i, o in enumerate(self.observations):
                o_dict = cast(Dict[str, Any], o)
                name = f"obs_{i}_{o_dict['type']}"
                obs_list.append(Observation(name, float(o_dict["value"]), weight=1.0))
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        from ..vadose.richards1d import run_richards_column

        new_layers = list(self.profile.layers)

        for layer_idx in self.layers_to_fit:
            original = new_layers[layer_idx]
            ks = param_values.get(f"ks_L{layer_idx}", original.ks_m_day)
            alpha = param_values.get(f"alpha_L{layer_idx}", original.alpha_1_m)
            n_val = param_values.get(f"n_L{layer_idx}", original.n)

            try:
                from dataclasses import replace
                new_layer = replace(original, ks_m_day=ks, alpha_1_m=alpha, n=n_val)
                new_layers[layer_idx] = new_layer
            except Exception:
                pass

        from dataclasses import replace
        new_profile = replace(self.profile, layers=new_layers)

        kc = param_values.get("kc", param_values.get("Kc", self.config.kc if self.config else 1.0))
        pff = param_values.get("preferential_flow_fraction", 
                               param_values.get("pref_flow_frac", 
                                                self.config.preferential_flow_fraction if self.config else 0.0))
        ttd_cv = param_values.get("ttd_cv", 
                                  param_values.get("ttd_default_cv", 
                                                   self.config.ttd_default_cv if self.config else 0.7))

        new_config = self.config
        if self.config:
            try:
                new_config = replace(self.config, kc=kc, preferential_flow_fraction=pff, ttd_default_cv=ttd_cv)
            except Exception:
                pass

        # 2. Run Model
        sim = run_richards_column(new_profile, self.forcing, config=new_config)

        # 3. Extract Observations
        results: Dict[str, float] = {}
        if isinstance(self.observations, list):
            for i, o in enumerate(self.observations):
                o_dict = cast(Dict[str, Any], o)
                t_idx = int(o_dict["time_idx"])
                target_z = float(o_dict["depth_m"])

                if hasattr(sim, "theta") and t_idx < len(sim.theta):
                    if o_dict["type"] == "theta":
                        profile_vals = sim.theta[t_idx]
                    else:
                        profile_vals = sim.psi_m[t_idx]
                    
                    z_arr = np.array(sim.z_m)
                    val_arr = np.array(profile_vals)
                    val = float(np.interp(target_z, z_arr, val_arr))
                    results[f"obs_{i}_{o_dict['type']}"] = val
                else:
                    results[f"obs_{i}_{o_dict['type']}"] = -9999.0

        return results


@dataclass
class NitrateSourceCalibrationObservation:
    node_id: str
    value: float
    target: str = "p_manure"  # "p_manure" or "p_fertilizer"
    weight: float = 1.0


class NitrateSourceCalibrationAdapter(CalibrationProblem):
    """
    Calibrates nitrate source-discrimination weights, priors, and QC thresholds.
    """

    DEFAULT_WEIGHTS = {
        "w1_no3_cl": 1.2,
        "w2_no3_k": 0.4,
        "w3_po4": 0.3,
        "w4_fe": 0.6,
        "w5_denitrif": 1.0,
        "w6_alk_coupling": 0.8,
    }

    def __init__(
        self,
        nodes_df: Any,
        observations: Sequence[NitrateSourceCalibrationObservation],
        edge_results: Optional[List[Any]] = None,
        config: Any = None,
        params_to_fit: Optional[List[str]] = None,
        base_overrides: Optional[Dict[str, Any]] = None,
    ):
        self.nodes_df = nodes_df
        self.observations = list(observations)
        self.edge_results = edge_results or []
        self.config = config
        self.params_to_fit = params_to_fit or [
            "prior_logit",
            "w1_no3_cl",
            "w2_no3_k",
            "w3_po4",
            "w4_fe",
            "w5_denitrif",
            "w6_alk_coupling",
            "evap_gate_factor",
            "denitrification_min_extent",
            "min_top_probability",
            "min_top_margin",
        ]
        self.base_overrides = base_overrides or {}

    @staticmethod
    def _sigmoid(value: float) -> float:
        value = max(-60.0, min(60.0, float(value)))
        return 1.0 / (1.0 + math.exp(-value))

    @staticmethod
    def _logit(probability: float) -> float:
        p = min(1.0 - 1e-9, max(1e-9, float(probability)))
        return math.log(p / (1.0 - p))

    def get_parameters(self) -> List[AdjustableParameter]:
        params: List[AdjustableParameter] = []
        for name in self.params_to_fit:
            if name == "prior_logit":
                prior_p = float(self.base_overrides.get("prior_prob", 0.5))
                params.append(
                    AdjustableParameter(
                        name=name,
                        value=self._logit(prior_p),
                        lower_bound=-8.0,
                        upper_bound=8.0,
                        prior_mean=0.0,
                        prior_sigma=2.0,
                    )
                )
            elif name in self.DEFAULT_WEIGHTS:
                value = float(
                    self.base_overrides.get("weights", {}).get(
                        name, self.DEFAULT_WEIGHTS[name]
                    )
                )
                params.append(
                    AdjustableParameter(
                        name=name,
                        value=value,
                        lower_bound=1e-6,
                        upper_bound=10.0,
                        log_transform=True,
                        prior_mean=value,
                        prior_sigma=1.0,
                    )
                )
            elif name == "evap_gate_factor":
                value = float(self.base_overrides.get(name, 0.5))
                params.append(
                    AdjustableParameter(name, value, 0.0, 1.0, prior_mean=value, prior_sigma=0.25)
                )
            elif name == "nitrate_source_min_mg_L":
                value = float(self.base_overrides.get(name, 10.0))
                params.append(
                    AdjustableParameter(
                        name,
                        value,
                        0.0,
                        250.0,
                        log_transform=False,
                        prior_mean=value,
                        prior_sigma=max(value, 1.0),
                    )
                )
            elif name in {"min_top_probability", "min_top_margin", "max_normalized_entropy", 
                          "min_tail_probability", "max_sensitivity_tv", "sensitivity_delta", 
                          "min_top_ci_gap"}:
                defaults = {
                    "min_top_probability": 0.55,
                    "min_top_margin": 0.12,
                    "max_normalized_entropy": 0.90,
                    "min_tail_probability": 0.02,
                    "max_sensitivity_tv": 0.35,
                    "sensitivity_delta": 0.5,
                    "min_top_ci_gap": 0.0,
                }
                value = float(self.base_overrides.get("isotope_qc", {}).get(name, defaults[name]))
                bounds = {
                    "min_top_probability": (0.0, 1.0),
                    "min_top_margin": (0.0, 1.0),
                    "max_normalized_entropy": (0.01, 1.0),
                    "min_tail_probability": (0.001, 0.5),
                    "max_sensitivity_tv": (0.01, 1.0),
                    "sensitivity_delta": (1e-4, 1.0),
                    "min_top_ci_gap": (-1.0, 1.0),
                }
                params.append(
                    AdjustableParameter(name, value, bounds[name][0], bounds[name][1], prior_mean=value, prior_sigma=0.2)
                )
            elif name in {"denitrification_min_extent", "denitrification_strength", "denitrification_slope",
                          "denitrification_slope_sigma", "denitrification_perp_sigma", "nitrification_strength",
                          "nitrification_sigma", "nitrification_o2_o18"}:
                defaults = {
                    "denitrification_min_extent": 0.05,
                    "denitrification_strength": 1.5,
                    "denitrification_slope": 0.65,
                    "denitrification_slope_sigma": 0.25,
                    "denitrification_perp_sigma": 4.0,
                    "nitrification_strength": 1.0,
                    "nitrification_sigma": 4.0,
                    "nitrification_o2_o18": 23.5,
                }
                value = float(self.base_overrides.get("isotope_process_constraints", {}).get(name, defaults[name]))
                bounds = {
                    "denitrification_min_extent": (0.0, 10.0),
                    "denitrification_strength": (0.0, 20.0),
                    "denitrification_slope": (0.1, 2.0),
                    "denitrification_slope_sigma": (0.01, 1.0),
                    "denitrification_perp_sigma": (0.1, 20.0),
                    "nitrification_strength": (0.0, 20.0),
                    "nitrification_sigma": (0.1, 20.0),
                    "nitrification_o2_o18": (0.0, 50.0),
                }
                params.append(
                    AdjustableParameter(name, value, bounds[name][0], bounds[name][1], prior_mean=value, prior_sigma=0.25)
                )
        return params

    def get_observations(self) -> List[Observation]:
        return [
            Observation(
                name=self._obs_name(obs),
                value=float(obs.value),
                weight=float(obs.weight),
            )
            for obs in self.observations
        ]

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        from ..nitrate_source_v2 import infer_node_posteriors

        overrides = dict(self.base_overrides)
        weights = dict(overrides.get("weights", {}))
        for key in self.DEFAULT_WEIGHTS:
            if key in param_values:
                weights[key] = float(param_values[key])
        if weights:
            overrides["weights"] = weights

        if "prior_logit" in param_values:
            overrides["prior_prob"] = self._sigmoid(param_values["prior_logit"])
        if "evap_gate_factor" in param_values:
            overrides["evap_gate_factor"] = min(1.0, max(0.0, param_values["evap_gate_factor"]))
        if "nitrate_source_min_mg_L" in param_values:
            overrides["nitrate_source_min_mg_L"] = max(
                0.0, float(param_values["nitrate_source_min_mg_L"])
            )
            
        process = dict(overrides.get("isotope_process_constraints", {}))
        for key in ("denitrification_min_extent", "denitrification_strength", "denitrification_slope",
                    "denitrification_slope_sigma", "denitrification_perp_sigma", "nitrification_strength",
                    "nitrification_sigma", "nitrification_o2_o18"):
            if key in param_values:
                process[key] = float(param_values[key])
        if process:
            overrides["isotope_process_constraints"] = process

        qc = dict(overrides.get("isotope_qc", {}))
        for key in ("min_top_probability", "min_top_margin", "max_normalized_entropy", 
                    "min_tail_probability", "max_sensitivity_tv", "sensitivity_delta", 
                    "min_top_ci_gap"):
            if key in param_values:
                qc[key] = float(param_values[key])
        if qc:
            overrides["isotope_qc"] = qc

        inferred = infer_node_posteriors(
            self.nodes_df,
            self.edge_results,
            config_overrides=overrides,
            config=self.config,
        )
        results: Dict[str, float] = {}
        for obs in self.observations:
            result = inferred.get(str(obs.node_id))
            value = -9999.0
            if result is not None:
                if obs.target == "p_fertilizer":
                    value = (
                        float(result.p_fertilizer)
                        if result.p_fertilizer is not None
                        else 1.0
                    )
                else:
                    value = float(result.p_manure) if result.p_manure is not None else 0.0
            results[self._obs_name(obs)] = value
        return results

    @staticmethod
    def _obs_name(obs: NitrateSourceCalibrationObservation) -> str:
        safe_node = re.sub(r"[^A-Za-z0-9_]+", "_", str(obs.node_id)).strip("_")
        return f"nitrate_{safe_node}_{obs.target}"


@dataclass
class AgeTemporalExperiment:
    id: str
    observed_age_days: Optional[float] = None
    weight: float = 1.0
    default_tau_days: float = 30.0
    distance_m: Optional[float] = None
    tracer_initial: Optional[float] = None
    tracer_observed: Optional[float] = None
    
    # For time-series node pairs:
    node_u_times: Optional[List[float]] = None
    node_u_values: Optional[List[float]] = None
    node_v_times: Optional[List[float]] = None
    node_v_values: Optional[List[float]] = None
    tracer_ion: Optional[str] = "Cl"


class AgeTemporalCalibrationAdapter(CalibrationProblem):
    """
    Calibrates groundwater age and temporal residence-time parameters.
    """

    def __init__(
        self,
        experiments: Sequence[AgeTemporalExperiment],
        params_to_fit: Optional[List[str]] = None,
        base_decay_rate_1_day: float = 1e-3,
        base_velocity_m_day: float = 0.1,
        base_porosity: float = 0.2,
        base_ttd_cv: float = 0.7,
    ):
        self.experiments = list(experiments)
        self.params_to_fit = params_to_fit or ["tau"]
        self.base_decay_rate_1_day = float(base_decay_rate_1_day)
        self.base_velocity_m_day = float(base_velocity_m_day)
        self.base_porosity = float(base_porosity)
        self.base_ttd_cv = float(base_ttd_cv)

    def get_parameters(self) -> List[AdjustableParameter]:
        params: List[AdjustableParameter] = []
        seen = set()
        for spec in self.params_to_fit:
            parts = spec.split(":")
            name = parts[0]
            scope = parts[1] if len(parts) > 1 else "global"
            if name in {"tau", "tau_days"}:
                if scope == "per_experiment":
                    for exp in self.experiments:
                        p_name = f"tau_days__{exp.id}"
                        if p_name in seen:
                            continue
                        seen.add(p_name)
                        params.append(
                            AdjustableParameter(
                                p_name,
                                max(float(exp.default_tau_days), 1e-6),
                                1e-6,
                                1e6,
                                log_transform=True,
                                prior_mean=max(float(exp.default_tau_days), 1e-6),
                                prior_sigma=1.0,
                            )
                        )
                elif "tau_days" not in seen:
                    seen.add("tau_days")
                    # Calculate default_tau from all default_tau_days
                    default_tau = float(np.median([e.default_tau_days for e in self.experiments]))
                    params.append(
                        AdjustableParameter(
                            "tau_days",
                            max(default_tau, 1e-6),
                            1e-6,
                            1e6,
                            log_transform=True,
                            prior_mean=max(default_tau, 1e-6),
                            prior_sigma=1.0,
                        )
                    )
            elif name in {"decay", "decay_rate", "decay_rate_1_day"} and "decay_rate_1_day" not in seen:
                seen.add("decay_rate_1_day")
                params.append(
                    AdjustableParameter(
                        "decay_rate_1_day",
                        max(self.base_decay_rate_1_day, 1e-12),
                        1e-12,
                        10.0,
                        log_transform=True,
                        prior_mean=max(self.base_decay_rate_1_day, 1e-12),
                        prior_sigma=2.0,
                    )
                )
            elif name in {"velocity", "velocity_m_day"} and "velocity_m_day" not in seen:
                seen.add("velocity_m_day")
                params.append(
                    AdjustableParameter(
                        "velocity_m_day",
                        max(self.base_velocity_m_day, 1e-9),
                        1e-9,
                        1e4,
                        log_transform=True,
                        prior_mean=max(self.base_velocity_m_day, 1e-9),
                        prior_sigma=2.0,
                    )
                )
            elif name == "porosity" and "porosity" not in seen:
                seen.add("porosity")
                params.append(
                    AdjustableParameter("porosity", self.base_porosity, 0.01, 0.9, prior_mean=self.base_porosity, prior_sigma=0.2)
                )
            elif name in {"ttd_cv", "cv"} and "ttd_cv" not in seen:
                seen.add("ttd_cv")
                params.append(
                    AdjustableParameter("ttd_cv", self.base_ttd_cv, 0.01, 5.0, prior_mean=self.base_ttd_cv, prior_sigma=0.5)
                )
        return params

    def get_observations(self) -> List[Observation]:
        obs = []
        for exp in self.experiments:
            if exp.node_v_values is not None:
                for idx, val in enumerate(exp.node_v_values):
                    t = exp.node_v_times[idx] if exp.node_v_times else idx
                    obs.append(
                        Observation(
                            name=f"{exp.id}_v_{idx}",
                            value=float(val),
                            weight=float(exp.weight),
                            time=float(t)
                        )
                    )
            else:
                obs.append(
                    Observation(
                        name=f"{exp.id}_age_days",
                        value=float(exp.observed_age_days) if exp.observed_age_days is not None else 0.0,
                        weight=float(exp.weight),
                    )
                )
        return obs

    @staticmethod
    def _simulate_ttd_point(
        t: float,
        u_times: List[float],
        u_vals: List[float],
        tau: float,
        cv: float,
        decay: float
    ) -> float:
        if cv < 0.05 or tau <= 0.0:
            t_lag = t - tau
            c_in = float(np.interp(t_lag, u_times, u_vals))
            return c_in * math.exp(-decay * tau)
        
        # Gamma TTD PDF integration
        shape = 1.0 / (cv * cv)
        scale = tau * (cv * cv)
        
        from scipy.stats import gamma
        
        max_lag = max(10.0, tau * 4.0)
        n_steps = 40
        lags = np.linspace(0.0, max_lag, n_steps)
        
        pdf_vals = gamma.pdf(lags, a=shape, scale=scale)
        w_sum = np.sum(pdf_vals)
        if w_sum > 0:
            weights = pdf_vals / w_sum
        else:
            weights = np.zeros_like(pdf_vals)
            weights[0] = 1.0
            
        c_out = 0.0
        for w, lag in zip(weights, lags):
            t_lag = t - lag
            c_in = float(np.interp(t_lag, u_times, u_vals))
            c_out += w * c_in * math.exp(-decay * lag)
        return c_out

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        results: Dict[str, float] = {}
        decay = max(float(param_values.get("decay_rate_1_day", self.base_decay_rate_1_day)), 1e-12)
        velocity = max(float(param_values.get("velocity_m_day", self.base_velocity_m_day)), 1e-12)
        porosity = min(0.99, max(1e-6, float(param_values.get("porosity", self.base_porosity))))
        cv = max(float(param_values.get("ttd_cv", self.base_ttd_cv)), 0.01)

        for exp in self.experiments:
            tau = float(param_values.get(f"tau_days__{exp.id}", param_values.get("tau_days", exp.default_tau_days)))
            if exp.distance_m is not None and "velocity_m_day" in param_values:
                tau = max(0.0, float(exp.distance_m) * porosity / velocity)
            
            if exp.node_u_times is not None and exp.node_v_times is not None:
                for idx, vt in enumerate(exp.node_v_times):
                    pred_val = self._simulate_ttd_point(
                        vt, exp.node_u_times, exp.node_u_values, tau, cv, decay
                    )
                    results[f"{exp.id}_v_{idx}"] = pred_val
            else:
                if exp.tracer_initial is not None and exp.tracer_observed is not None:
                    c0 = max(float(exp.tracer_initial), 1e-30)
                    c = max(float(exp.tracer_observed), 1e-30)
                    tau = max(0.0, -math.log(min(c / c0, 1.0)) / decay)
                results[f"{exp.id}_age_days"] = float(tau)
        return results


@dataclass
class TopologyCalibrationObservation:
    edge_id: str
    observed_present: float
    weight: float = 1.0


class TopologyCalibrationAdapter(CalibrationProblem):
    """
    Calibrates soft edge-presence probabilities for candidate flow-network topology.
    """

    def __init__(
        self,
        candidate_edges: Sequence[Any],
        observations: Sequence[TopologyCalibrationObservation],
        prior_sigma: float = 2.0,
        normalize_by_upstream: bool = False,
    ):
        self.candidate_edges = list(candidate_edges)
        self.observations = list(observations)
        self.prior_sigma = float(prior_sigma)
        self.normalize_by_upstream = bool(normalize_by_upstream)
        self._edge_by_id = {str(edge.edge_id): edge for edge in self.candidate_edges}

    @staticmethod
    def _safe_edge_id(edge_id: Any) -> str:
        return re.sub(r"[^A-Za-z0-9_]+", "_", str(edge_id)).strip("_") or "edge"

    @classmethod
    def _param_name(cls, edge_id: Any) -> str:
        return f"edge_logit__{cls._safe_edge_id(edge_id)}"

    @staticmethod
    def _logit(probability: float) -> float:
        p = min(1.0 - 1e-9, max(1e-9, float(probability)))
        return math.log(p / (1.0 - p))

    @staticmethod
    def _sigmoid(value: float) -> float:
        value = max(-60.0, min(60.0, float(value)))
        return 1.0 / (1.0 + math.exp(-value))

    def get_parameters(self) -> List[AdjustableParameter]:
        params: List[AdjustableParameter] = []
        for edge in self.candidate_edges:
            attrs = getattr(edge, "attrs", {}) or {}
            prior = float(attrs.get("edge_confidence", attrs.get("p_uv", 0.5)) or 0.5)
            prior_logit = self._logit(prior)
            params.append(
                AdjustableParameter(
                    name=self._param_name(edge.edge_id),
                    value=prior_logit,
                    lower_bound=-12.0,
                    upper_bound=12.0,
                    prior_mean=prior_logit,
                    prior_sigma=self.prior_sigma,
                )
            )
        return params

    def get_observations(self) -> List[Observation]:
        return [
            Observation(
                name=f"topology_{self._safe_edge_id(obs.edge_id)}",
                value=float(obs.observed_present),
                weight=float(obs.weight),
            )
            for obs in self.observations
        ]

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        probabilities = self.edge_probabilities(param_values)
        return {
            f"topology_{self._safe_edge_id(obs.edge_id)}": float(
                probabilities.get(str(obs.edge_id), 0.0)
            )
            for obs in self.observations
        }

    def edge_probabilities(self, param_values: Dict[str, float]) -> Dict[str, float]:
        logits = {
            str(edge.edge_id): float(param_values.get(self._param_name(edge.edge_id), 0.0))
            for edge in self.candidate_edges
        }
        if not self.normalize_by_upstream:
            return {edge_id: self._sigmoid(logit) for edge_id, logit in logits.items()}

        grouped: Dict[str, List[str]] = {}
        for edge in self.candidate_edges:
            grouped.setdefault(str(edge.u), []).append(str(edge.edge_id))

        probs: Dict[str, float] = {}
        for edge_ids in grouped.values():
            vals = np.array([logits[edge_id] for edge_id in edge_ids], dtype=float)
            vals = vals - float(np.max(vals))
            exp_vals = np.exp(vals)
            denom = float(np.sum(exp_vals)) or 1.0
            for edge_id, val in zip(edge_ids, exp_vals):
                probs[edge_id] = float(val / denom)
        return probs

    def selected_edges(
        self, param_values: Dict[str, float], threshold: float = 0.5
    ) -> List[Any]:
        probabilities = self.edge_probabilities(param_values)
        return [
            edge
            for edge in self.candidate_edges
            if probabilities.get(str(edge.edge_id), 0.0) >= threshold
        ]


class TopologyOuterLoopCalibrator:
    """
    Implements discrete outer-loop search to optimize network flow topology.
    Alternates between probability calibration, edge selection, network chemical fit,
    and add/remove/swap operations under constraints.
    """

    def __init__(
        self,
        samples_df: Any,
        candidate_edges: Sequence[Any],
        observations: Sequence[TopologyCalibrationObservation],
        config: Any = None,
        max_iterations: int = 5,
        max_neighbors: int = 2,
        head_key: str = "hydraulic_head",
        elevation_key: str = "elevation",
        layer_key: str = "aquifer_layer",
    ):
        import pandas as pd
        self.samples_df = samples_df
        self.sample_map = {}
        if isinstance(samples_df, pd.DataFrame):
            for _, row in samples_df.iterrows():
                site_id = str(row.get("site_id", row.get("SampleID", row.get("Code", ""))))
                if site_id:
                    self.sample_map[site_id] = row.to_dict()
        elif isinstance(samples_df, (list, tuple)):
            for row in samples_df:
                site_id = str(row.get("site_id", row.get("SampleID", row.get("Code", ""))))
                if site_id:
                    self.sample_map[site_id] = dict(row)
        elif isinstance(samples_df, dict):
            self.sample_map = {str(k): dict(v) for k, v in samples_df.items()}

        self.candidate_edges = list(candidate_edges)
        self.observations = list(observations)
        from ..config import Config
        self.config = config or Config()
        self.max_iterations = max_iterations
        self.max_neighbors = max_neighbors
        self.head_key = head_key
        self.elevation_key = elevation_key
        self.layer_key = layer_key

    def check_constraints(self, edges: List[Any]) -> bool:
        """
        Validates topology constraints: acyclicity, max neighbors, head direction, layer compatibility.
        """
        # 1. Acyclicity (directed graph check)
        adj = {}
        for edge in edges:
            adj.setdefault(edge.u, []).append(edge.v)
        visited = {}

        def dfs(node):
            if visited.get(node) == 1:
                return False  # cycle
            if visited.get(node) == 2:
                return True
            visited[node] = 1
            for neighbor in adj.get(node, []):
                if not dfs(neighbor):
                    return False
            visited[node] = 2
            return True

        for edge in edges:
            if visited.get(edge.u) is None:
                if not dfs(edge.u):
                    return False

        # 2. Max downstream/upstream neighbors
        downstream = {}
        upstream = {}
        for edge in edges:
            downstream[edge.u] = downstream.get(edge.u, 0) + 1
            upstream[edge.v] = upstream.get(edge.v, 0) + 1
            if downstream[edge.u] > self.max_neighbors or upstream[edge.v] > self.max_neighbors:
                return False

        # 3. Head direction & layer compatibility
        for edge in edges:
            u_sample = self.sample_map.get(edge.u)
            v_sample = self.sample_map.get(edge.v)
            if u_sample and v_sample:
                # Flow from higher head to lower head
                head_u = u_sample.get(self.head_key)
                head_v = v_sample.get(self.head_key)
                if head_u is not None and head_v is not None:
                    try:
                        if float(head_u) < float(head_v) - 0.01:
                            return False  # uphill flow
                    except (ValueError, TypeError):
                        pass

                # Layer compatibility (screens/aquifers)
                layer_u = u_sample.get(self.layer_key)
                layer_v = v_sample.get(self.layer_key)
                if layer_u is not None and layer_v is not None:
                    try:
                        if float(layer_u) > float(layer_v) + 0.1:
                            if not getattr(self.config, "upward_flow_allowed", True):
                                return False
                    except (ValueError, TypeError):
                        pass
        return True

    def calculate_score(self, edges: List[Any]) -> Tuple[float, Dict[str, Any]]:
        """
        Runs network chemistry fit on the active edges and calculates AIC/BIC.
        """
        if not edges:
            return 1e9, {"phi": 1e9, "aic": 1e9, "bic": 1e9}

        from ..inference.network_fit import fit_network
        try:
            results = fit_network(
                samples=list(self.sample_map.values()),
                edges=edges,
                config=self.config,
            )
            phi = sum(r.objective_score ** 2 for r in results)
        except Exception as e:
            logger.warning(f"fit_network failed during topology search: {e}")
            phi = 1e9

        n_obs = len(self.observations)
        n_par = len(edges)

        # Calculate AIC/BIC
        if n_obs > 0:
            mse = phi / n_obs
            log_lik = -0.5 * n_obs * (math.log(2.0 * math.pi * max(mse, 1e-12)) + 1.0)
            aic = 2.0 * n_par - 2.0 * log_lik
            bic = math.log(n_obs) * n_par - 2.0 * log_lik
        else:
            aic = phi + 2 * n_par
            bic = phi + math.log(max(1, n_par)) * n_par

        return aic, {"phi": phi, "aic": aic, "bic": bic, "results": results}

    def search(self) -> Dict[str, Any]:
        """
        Executes discrete outer-loop search over the topology space.
        """
        sorted_candidates = sorted(
            self.candidate_edges,
            key=lambda e: float(e.attrs.get("edge_confidence", e.attrs.get("p_uv", 0.5))),
            reverse=True,
        )

        current_edges = []
        for edge in sorted_candidates:
            test_edges = current_edges + [edge]
            if self.check_constraints(test_edges):
                current_edges.append(edge)

        best_edges = list(current_edges)
        best_score, best_meta = self.calculate_score(best_edges)

        logger.info(f"Initial topology selected {len(best_edges)} edges. Score (AIC) = {best_score:.4f}")

        # Iterate: add/remove/swap edge operations
        improved = True
        iteration = 0
        while improved and iteration < self.max_iterations:
            improved = False
            iteration += 1
            logger.info(f"Topology Search Iteration {iteration}...")

            # Try Add operation
            for edge in self.candidate_edges:
                if edge not in current_edges:
                    test_edges = current_edges + [edge]
                    if self.check_constraints(test_edges):
                        score, meta = self.calculate_score(test_edges)
                        if score < best_score:
                            best_score = score
                            best_edges = list(test_edges)
                            best_meta = meta
                            current_edges = list(test_edges)
                            improved = True
                            logger.info(f"  [ADD] Edge {edge.edge_id} improved score to {best_score:.4f}")
                            break

            if improved:
                continue

            # Try Remove operation
            for edge in current_edges:
                test_edges = [e for e in current_edges if e != edge]
                if self.check_constraints(test_edges):
                    score, meta = self.calculate_score(test_edges)
                    if score < best_score:
                        best_score = score
                        best_edges = list(test_edges)
                        best_meta = meta
                        current_edges = list(test_edges)
                        improved = True
                        logger.info(f"  [REMOVE] Edge {edge.edge_id} improved score to {best_score:.4f}")
                        break

            if improved:
                continue

            # Try Swap operation (swap an active edge with an inactive candidate)
            for active_edge in current_edges:
                for inactive_edge in self.candidate_edges:
                    if inactive_edge not in current_edges:
                        test_edges = [e for e in current_edges if e != active_edge] + [inactive_edge]
                        if self.check_constraints(test_edges):
                            score, meta = self.calculate_score(test_edges)
                            if score < best_score:
                                best_score = score
                                best_edges = list(test_edges)
                                best_meta = meta
                                current_edges = list(test_edges)
                                improved = True
                                logger.info(f"  [SWAP] Edge {active_edge.edge_id} -> {inactive_edge.edge_id} improved score to {best_score:.4f}")
                                break
                if improved:
                    break

        # Compile final ensemble / uncertainty if requested
        ensemble = [best_edges]
        for edge in self.candidate_edges:
            if edge not in best_edges:
                test_edges = best_edges + [edge]
                if self.check_constraints(test_edges):
                    score, _ = self.calculate_score(test_edges)
                    if abs(score - best_score) < 5.0:
                        ensemble.append(test_edges)
            else:
                test_edges = [e for e in best_edges if e != edge]
                if self.check_constraints(test_edges):
                    score, _ = self.calculate_score(test_edges)
                    if abs(score - best_score) < 5.0:
                        ensemble.append(test_edges)

        edge_probs = {}
        for edge in self.candidate_edges:
            count = sum(1 for ens in ensemble if edge in ens)
            edge_probs[edge.edge_id] = count / len(ensemble)

        return {
            "success": True,
            "optimal_parameters": {f"edge_presence__{e.edge_id}": (1.0 if e in best_edges else 0.0) for e in self.candidate_edges},
            "selected_edges": best_edges,
            "edge_probabilities": edge_probs,
            "phi": best_meta.get("phi", 0.0),
            "aic": best_meta.get("aic", 0.0),
            "bic": best_meta.get("bic", 0.0),
            "n_iterations": iteration,
            "ensemble_size": len(ensemble),
        }
