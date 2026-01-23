"""
Calibration Adapters for different Hydrosheaf modules.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Union, cast
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

    def get_parameters(self) -> List[AdjustableParameter]:
        """
        Generate PEST parameters based on requested fit configuration.
        Supports global and edge-specific scoping.
        """
        pest_params = []
        seen_params = set()

        # Gather all required parameters from all experiments + global definitions
        # If params_to_fit contains "calcite:k", it implies global.
        # If it contains "calcite:k:edge_A_B", it implies specific.
        
        # Actually, let's auto-generate required parameters based on experiments and fit spec.
        # Strategy:
        # 1. Parse params_to_fit rules
        # 2. Iterate experiments to see which rules apply (e.g. edge-specific vs global)
        
        # Simplified: params_to_fit lists the *templates* of what we want to fit.
        # e.g. "calcite:k" -> fit k for calcite everywhere.
        # If we want edge-specific, we'd need to know if we are doing that.
        
        # Let's assume params_to_fit defines the SCOPE too.
        # "calcite:k:global" (default) or "calcite:k:per_edge"
        
        # But wait, self.params_to_fit is usually a list of parameter NAMES or keys.
        # Let's stick to a robust convention:
        # If the string has 2 parts "calcite:k" -> Global parameter "log_k_calcite"
        # If the string has 3 parts "calcite:k:edge1" -> Edge parameter "log_k_calcite__edge1"
        
        # However, usually we configure "fit k for calcite per edge".
        # That generates many parameters.
        
        # Let's interpret self.params_to_fit as a list of "active calibration targets".
        # We'll scan experiments to instantiate them.
        
        # For now, implemented GLOBAL only as per previous logic, but allowing extension.
        # If a param name contains "__", it's scoped.
        
        for p_spec in self.params_to_fit:
            parts = p_spec.split(":")
            reaction = parts[0]
            p_type = parts[1]
            scope = parts[2] if len(parts) > 2 else "global"
            
            if reaction not in self.base_params:
                continue
                
            base_p = self.base_params[reaction]
            
            if scope == "global":
                param_name = f"log_{p_type}_{reaction}"
                if param_name in seen_params: continue
                seen_params.add(param_name)
                
                val = base_p.rate_constant if p_type == "k" else base_p.surface_area
                # Add parameter
                pest_params.append(self._make_param(param_name, val, p_type))
            
            elif scope == "per_edge":
                # Find all edges used in experiments
                edges = {exp.edge_id for exp in self.experiments if exp.edge_id}
                for edge in edges:
                    param_name = f"log_{p_type}_{reaction}__{edge}"
                    if param_name in seen_params: continue
                    seen_params.add(param_name)
                    
                    val = base_p.rate_constant if p_type == "k" else base_p.surface_area
                    pest_params.append(self._make_param(param_name, val, p_type))

            else:
                raise ValueError(
                    f"Unsupported kinetic calibration scope '{scope}' in '{p_spec}'. "
                    "Use 'global' or 'per_edge'."
                )

        return pest_params

    def _make_param(self, name: str, val: float, p_type: str) -> AdjustableParameter:
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
            
            active_extents = list(exp.reaction_extents)
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
        if exp.edge_id:
            if (
                f"log_k_{reaction}__{exp.edge_id}" in param_values
                or f"log_A_{reaction}__{exp.edge_id}" in param_values
            ):
                return True
        return False

    def _resolve_parameters(self, exp: KineticExperiment, param_values: Dict[str, float]) -> Dict[str, KineticParameters]:
        """
        Resolve final KineticParameters for an experiment based on global/edge scopes.
        """
        resolved = {}
        for name, base_p in self.base_params.items():
            # Start with base
            k = base_p.rate_constant
            A = base_p.surface_area
            
            # NOTE: parameter names contain "log_" for readability, but param_values are real-space
            # values when using PEST++ with PARTRANS=log (PEST optimizes log, but writes real values).
            if f"log_k_{name}" in param_values:
                k = float(param_values[f"log_k_{name}"])
                
            if f"log_A_{name}" in param_values:
                A = float(param_values[f"log_A_{name}"])

            # Check Edge overrides (higher priority)
            if exp.edge_id:
                if f"log_k_{name}__{exp.edge_id}" in param_values:
                    k = float(param_values[f"log_k_{name}__{exp.edge_id}"])
                if f"log_A_{name}__{exp.edge_id}" in param_values:
                    A = float(param_values[f"log_A_{name}__{exp.edge_id}"])
            
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
                # FloPy runner logic would go here (slower)
                # Not implemented for this adapter yet to keep it simple,
                # but easily added using run_1d_transport
                pass

        return results


class VadoseCalibrationAdapter(CalibrationProblem):
    """
    Adapts Vadose Zone Richards Model for PEST calibration.
    Calibrates: saturated conductivity (Ks), alpha, n.
    """

    def __init__(
        self,
        profile: Any,  # VadoseProfile
        forcing: Any,  # List[VadoseForcingSample]
        observations: Union[Dict[str, Any], List[Dict[str, Any]]],  # {time_idx: {"theta": val, "depth": d}} OR List of dicts
        config: Any = None,
        layers_to_fit: List[int] = [0],  # Indices of layers to calibrate
    ):
        self.profile = profile
        self.forcing = forcing
        self.observations = observations
        self.config = config
        self.layers_to_fit = layers_to_fit

    def get_parameters(self) -> List[AdjustableParameter]:
        pest_params = []

        for layer_idx in self.layers_to_fit:
            layer = self.profile.layers[layer_idx]
            # Assume we fit Ks and Alpha
            if layer.ks_m_day is not None:
                ks = float(layer.ks_m_day)
                pest_params.append(
                    AdjustableParameter(
                        f"ks_L{layer_idx}", ks, ks * 1e-2, ks * 1e2, log_transform=True
                    )
                )
            if layer.alpha_1_m is not None:
                alpha = float(layer.alpha_1_m)
                pest_params.append(
                    AdjustableParameter(
                        f"alpha_L{layer_idx}",
                        alpha,
                        alpha * 0.5,
                        alpha * 2.0,
                        log_transform=True,
                    )
                )

        return pest_params

    def get_observations(self) -> List[Observation]:
        obs_list = []
        # Observations format: List of (time_index, depth, value, type="theta"|"psi")
        # Simplified for demo: We assume observations is a list of dicts
        # [{"time_idx": 10, "depth_m": 0.5, "value": 0.3, "type": "theta"}]

        if isinstance(self.observations, list):
            for i, o in enumerate(self.observations):
                # Cast o to dict to satisfy mypy
                o_dict = cast(Dict[str, Any], o)
                name = f"obs_{i}_{o_dict['type']}"
                obs_list.append(Observation(name, float(o_dict["value"]), weight=1.0))
        return obs_list

    def run_model(self, param_values: Dict[str, float]) -> Dict[str, float]:
        # 1. Update Profile with new params
        # Deep copy profile structure?
        # For efficiency, we might just modify a copy of the layer objects

        from ..vadose.richards1d import run_richards_column

        # Create a new profile with updated layers
        new_layers = list(self.profile.layers)  # Shallow copy of list

        for layer_idx in self.layers_to_fit:
            # We need to construct a new layer object because they might be frozen/dataclasses
            # Assuming VadoseLayer is a dataclass
            original = new_layers[layer_idx]

            # Get new values
            ks = param_values.get(f"ks_L{layer_idx}", original.ks_m_day)
            alpha = param_values.get(f"alpha_L{layer_idx}", original.alpha_1_m)

            # Update (assuming dataclass with replace or constructor)
            # We don't have easy access to the class definition here to check if it's frozen.
            # Assuming we can instantiate a new one with same fields.
            # Let's rely on type(original)(...) if generic, or specific import.

            # Hack: modify if not frozen, or use replace
            try:
                from dataclasses import replace

                new_layer = replace(original, ks_m_day=ks, alpha_1_m=alpha)
                new_layers[layer_idx] = new_layer
            except Exception:
                pass  # Fallback if not dataclass

        # Reconstruct profile
        from dataclasses import replace

        new_profile = replace(self.profile, layers=new_layers)

        # 2. Run Model
        sim = run_richards_column(new_profile, self.forcing, config=self.config)

        # 3. Extract Observations
        results: Dict[str, float] = {}
        if isinstance(self.observations, list):
            for i, o in enumerate(self.observations):
                o_dict = cast(Dict[str, Any], o)
                # Find time index
                t_idx = int(o_dict["time_idx"])
                target_z = float(o_dict["depth_m"])

                # Interpolate depth
                # sim.z_m is grid. sim.theta[t_idx] is profile.
                # Assuming sim has these attributes (it is Any)
                if hasattr(sim, "theta") and t_idx < len(sim.theta):
                    if o_dict["type"] == "theta":
                        profile_vals = sim.theta[t_idx]
                    else:
                        profile_vals = sim.psi_m[t_idx]
                    
                    # np.interp expects arrays
                    z_arr = np.array(sim.z_m)
                    val_arr = np.array(profile_vals)
                    
                    val = float(np.interp(target_z, z_arr, val_arr))
                    results[f"obs_{i}_{o_dict['type']}"] = val
                else:
                    results[f"obs_{i}_{o_dict['type']}"] = -9999.0

        return results
