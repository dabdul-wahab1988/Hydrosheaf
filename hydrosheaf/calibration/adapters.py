"""
Calibration Adapters for different Hydrosheaf modules.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Callable, Any
import numpy as np

from .definitions import AdjustableParameter, Observation
from .problem import CalibrationProblem
from ..reactive_transport import KineticParameters
from ..reactive_transport.kinetic_phreeqc import (
    run_phreeqc_kinetic,
    build_kinetic_block,
)
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
    reaction_extents: List[float]  # Estimated extents from inverse model
    reaction_labels: List[str]  # Names of reactions corresponding to extents
    observations: Dict[str, float]  # {ion_name: concentration_mg_L}
    temperature_c: float = 25.0
    weight: float = 1.0


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
        params_to_fit: List[
            str
        ] = None,  # List of "reaction:param" strings, e.g. ["calcite:k", "calcite:A"]
    ):
        self.base_params = base_params
        self.experiments = experiments
        self.config = config

        # Default to fitting k for all reactions in base_params if not specified
        if params_to_fit is None:
            self.params_to_fit = [f"{r}:k" for r in base_params.keys()]
        else:
            self.params_to_fit = params_to_fit

    def get_parameters(self) -> List[AdjustableParameter]:
        """Generate PEST parameters based on requested fit configuration."""
        pest_params = []

        for p_spec in self.params_to_fit:
            reaction, p_type = p_spec.split(":")
            if reaction not in self.base_params:
                continue

            base_p = self.base_params[reaction]

            if p_type == "k":  # Rate constant
                # Usually log-transformed
                val = base_p.rate_constant
                pest_params.append(
                    AdjustableParameter(
                        name=f"log_k_{reaction}",
                        value=val,
                        lower_bound=val * 1e-4,
                        upper_bound=val * 1e4,
                        log_transform=True,
                        prior_mean=val,
                        prior_sigma=1.0,  # 1 order of magnitude uncertainty
                    )
                )
            elif p_type == "A":  # Surface Area
                val = base_p.surface_area
                pest_params.append(
                    AdjustableParameter(
                        name=f"log_A_{reaction}",
                        value=val,
                        lower_bound=val * 0.1,
                        upper_bound=val * 10.0,
                        log_transform=True,
                        prior_mean=val,
                        prior_sigma=0.5,
                    )
                )

        return pest_params

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
        # 1. Update Parameters
        current_params = self.base_params.copy()

        for p_name, p_val in param_values.items():
            # Parse name: log_k_calcite -> reaction="calcite", type="k"
            if p_name.startswith("log_k_"):
                reaction = p_name[6:]
                if reaction in current_params:
                    # Update a COPY of the dataclass
                    # (But here current_params is a shallow copy of the dict, objects are shared.
                    # We must copy the object to avoid side effects across runs if we were parallel,
                    # though here it's serial. Safer to replace.)
                    orig = current_params[reaction]
                    # Create new instance with updated value
                    from ..reactive_transport import KineticParameters as KP

                    new_p = KP(
                        reaction_name=orig.reaction_name,
                        rate_constant=p_val,  # PEST handles transform, so this is real value
                        surface_area=orig.surface_area,
                        reaction_order=orig.reaction_order,
                        activation_energy=orig.activation_energy,
                        reference_temp_k=orig.reference_temp_k,
                    )
                    current_params[reaction] = new_p

            elif p_name.startswith("log_A_"):
                reaction = p_name[6:]
                if reaction in current_params:
                    orig = current_params[reaction]
                    from ..reactive_transport import KineticParameters as KP

                    new_p = KP(
                        reaction_name=orig.reaction_name,
                        rate_constant=orig.rate_constant,
                        surface_area=p_val,
                        reaction_order=orig.reaction_order,
                        activation_energy=orig.activation_energy,
                        reference_temp_k=orig.reference_temp_k,
                    )
                    current_params[reaction] = new_p

        # 2. Run Experiments
        results = {}

        for exp in self.experiments:
            # Build input
            kinetics_block = build_kinetic_block(
                reaction_labels=exp.reaction_labels,
                extents=exp.reaction_extents,
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
                final_vals = out["final_composition"]

                # Create a dict for easy lookup
                sim_map = {}
                for i, ion in enumerate(ion_order):
                    if i < len(final_vals):
                        sim_map[ion] = final_vals[i]  # mmol/L?
                        # run_phreeqc_kinetic returns final_composition in mmol/L?
                        # Let's check the source code of run_phreeqc_kinetic...
                        # Yes: "final_comp.append(float(final_row...)*1000.0)" -> mol to mmol.

                # Add to results
                for ion in exp.observations:
                    obs_name = f"{exp.id}_{ion}"
                    # Check units! Adapter assumes exp.observations are in mmol/L
                    # because internal units are mmol/L.
                    # If user provides mg/L, they should convert before creating KineticExperiment.
                    if ion in sim_map:
                        results[obs_name] = sim_map[ion]
                    else:
                        results[obs_name] = -9999.0  # Missing
                        logger.science(
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
        observations: Dict[str, Any],  # {time_idx: {"theta": val, "depth": d}}
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
                name = f"obs_{i}_{o['type']}"
                obs_list.append(Observation(name, o["value"], weight=1.0))
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
        results = {}
        if isinstance(self.observations, list):
            for i, o in enumerate(self.observations):
                # Find time index
                t_idx = o["time_idx"]
                target_z = o["depth_m"]

                # Interpolate depth
                # sim.z_m is grid. sim.theta[t_idx] is profile.
                if t_idx < len(sim.theta):
                    profile_vals = (
                        sim.theta[t_idx] if o["type"] == "theta" else sim.psi_m[t_idx]
                    )
                    val = np.interp(target_z, sim.z_m, profile_vals)
                    results[f"obs_{i}_{o['type']}"] = float(val)
                else:
                    results[f"obs_{i}_{o['type']}"] = -9999.0

        return results
