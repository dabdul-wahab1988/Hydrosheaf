"""Robust certified measurement design for bounded reaction and mixing models.

This adapter compiles an explicit *linear, conditional* hydrochemical model into
ACMD's paired-state minimax loop.  It does not run PHREEQC or infer isotope
response coefficients from generic mineral names.  Saturation index can inform
reaction direction, but finite extent bounds and tracer responses must be
supplied by the caller in consistent units.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import itertools
import json
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np

from ..acmd import ACMD, ACMDAction, ACMDLoop, ACMDMode
from ..acmd_polytope import PolyhedralStateSpace
from ..nuclear.ttd_identified import AgeFunctional


def _nonnegative_concentrations(values: Mapping[str, float], species: Sequence[str], name: str) -> dict[str, float]:
    if set(values) != set(species):
        raise ValueError(f"{name} must contain exactly the declared observed species.")
    result = {key: float(values[key]) for key in species}
    if any(not math.isfinite(value) or value < 0.0 for value in result.values()):
        raise ValueError(f"{name} concentrations must be finite and non-negative.")
    return result


def reaction_stoichiometry_from_config(config: Any, **dictionary_kwargs: Any) -> dict[str, dict[str, float]]:
    """Orient the existing reaction dictionary as ``reaction -> ion -> coefficient``.

    The existing dictionary's rows are reactions. Its penalty scales and
    context gates are not converted into hard constraints or extent bounds.
    """
    from ..config import DEFAULT_ION_ORDER
    from ..models.reactions import build_reaction_dictionary

    matrix, labels, _, _ = build_reaction_dictionary(config, **dictionary_kwargs)
    ions = tuple(config.ion_order or DEFAULT_ION_ORDER)
    if len(matrix) != len(labels):
        raise ValueError("Reaction dictionary matrix and labels have different lengths.")
    if len(set(labels)) != len(labels):
        raise ValueError("Reaction dictionary contains duplicate labels.")
    result: dict[str, dict[str, float]] = {}
    for label, row in zip(labels, matrix):
        if len(row) != len(ions):
            raise ValueError(f"Reaction {label!r} row does not match ion_order.")
        result[label] = {ion: float(value) for ion, value in zip(ions, row) if value != 0}
    return result


@dataclass(frozen=True)
class ChemicalPolytope:
    """Compiled compatible set and its concentration-model provenance."""

    state_space: PolyhedralStateSpace
    species: tuple[str, ...]
    reaction_names: tuple[str, ...]
    endmember_names: tuple[str, ...]
    concentration_unit: str
    reaction_stoichiometry: Mapping[str, Mapping[str, float]]
    endmember_concentrations: Mapping[str, Mapping[str, float]]

    @property
    def state_names(self) -> tuple[str, ...]:
        return self.state_space.state_names


def compile_geochemical_polytope(
    major_ions_u: Mapping[str, float],
    major_ions_v: Mapping[str, float],
    reaction_stoichiometry: Mapping[str, Mapping[str, float]],
    extent_bounds: Mapping[str, tuple[float, float]],
    ion_error_bounds: Mapping[str, float],
    *,
    endmembers: Optional[Mapping[str, Mapping[str, float]]] = None,
    absolute_extent_cap: Optional[tuple[Sequence[str], float]] = None,
    concentration_unit: str = "mmol/L",
) -> ChemicalPolytope:
    """Compile ``c_v = S ξ + C_end λ`` with bounded extents and mixing.

    ``major_ions_u`` is always the upstream endmember.  Optional additional
    endmembers are mixed with it; their fractions are non-negative and sum to
    one.  Every listed concentration and extent uses ``concentration_unit``.
    ``ion_error_bounds`` are absolute intervals, not Gaussian standard errors.
    ``absolute_extent_cap`` is a cap on the sum of absolute extents for a small
    named group, in the same units; it is not a raw CEC value.
    """
    species = tuple(major_ions_v)
    reactions = tuple(reaction_stoichiometry)
    if not species or not reactions:
        raise ValueError("At least one observed species and reaction are required.")
    if any(not str(name).strip() for name in (*species, *reactions)):
        raise ValueError("Species and reaction names must be non-empty.")
    if not concentration_unit.strip():
        raise ValueError("concentration_unit must be declared.")
    upstream = _nonnegative_concentrations(major_ions_u, species, "major_ions_u")
    downstream = _nonnegative_concentrations(major_ions_v, species, "major_ions_v")
    if set(ion_error_bounds) != set(species):
        raise ValueError("ion_error_bounds must contain exactly the observed species.")
    errors = {ion: float(ion_error_bounds[ion]) for ion in species}
    if any(not math.isfinite(value) or value < 0 for value in errors.values()):
        raise ValueError("ion_error_bounds must be finite and non-negative.")
    if set(extent_bounds) != set(reactions):
        raise ValueError("extent_bounds must contain exactly the declared reactions.")

    endmember_data: dict[str, dict[str, float]] = {"upstream": upstream}
    for name, composition in (endmembers or {}).items():
        if not str(name).strip() or name == "upstream":
            raise ValueError("Additional endmember names must be non-empty and differ from 'upstream'.")
        endmember_data[name] = _nonnegative_concentrations(composition, species, f"endmember {name!r}")
    endmember_names = tuple(endmember_data)

    stoich: dict[str, dict[str, float]] = {}
    for reaction in reactions:
        coefficients = {ion: float(value) for ion, value in reaction_stoichiometry[reaction].items()}
        if any(not math.isfinite(value) for value in coefficients.values()):
            raise ValueError(f"Reaction {reaction!r} has non-finite stoichiometry.")
        stoich[reaction] = coefficients

    names = tuple(f"extent:{reaction}" for reaction in reactions) + tuple(
        f"mix:{name}" for name in endmember_names
    )
    n_reactions = len(reactions)
    n_states = len(names)
    bounds: list[tuple[float, float]] = []
    for reaction in reactions:
        pair = extent_bounds[reaction]
        if len(pair) != 2:
            raise ValueError(f"Extent bound for {reaction!r} must be a pair.")
        bounds.append((float(pair[0]), float(pair[1])))
    bounds.extend((0.0, 1.0) for _ in endmember_names)
    if len(endmember_names) == 1:
        bounds[-1] = (1.0, 1.0)

    rows: list[np.ndarray] = []
    limits: list[float] = []
    for ion in species:
        row = np.asarray(
            [stoich[reaction].get(ion, 0.0) for reaction in reactions]
            + [endmember_data[name][ion] for name in endmember_names],
            dtype=float,
        )
        rows.extend((row, -row))
        limits.extend((downstream[ion] + errors[ion], -downstream[ion] + errors[ion]))

    if absolute_extent_cap is not None:
        group, cap = absolute_extent_cap
        group = tuple(group)
        if not group or len(set(group)) != len(group) or len(group) > 12:
            raise ValueError("absolute_extent_cap needs 1 to 12 distinct reactions.")
        if any(name not in reactions for name in group):
            raise ValueError("absolute_extent_cap names an unknown reaction.")
        cap = float(cap)
        if not math.isfinite(cap) or cap < 0:
            raise ValueError("absolute_extent_cap must be finite and non-negative.")
        for signs in itertools.product((-1.0, 1.0), repeat=len(group)):
            row = np.zeros(n_states)
            for reaction, sign in zip(group, signs):
                row[reactions.index(reaction)] = sign
            rows.append(row)
            limits.append(cap)

    mixing_sum = np.zeros((1, n_states))
    mixing_sum[0, n_reactions:] = 1.0
    state_space = PolyhedralStateSpace(
        state_names=names,
        bounds=bounds,
        a_ub=rows,
        b_ub=limits,
        a_eq=mixing_sum,
        b_eq=[1.0],
    )
    return ChemicalPolytope(
        state_space=state_space,
        species=species,
        reaction_names=reactions,
        endmember_names=endmember_names,
        concentration_unit=concentration_unit,
        reaction_stoichiometry=stoich,
        endmember_concentrations=endmember_data,
    )


@dataclass(frozen=True)
class ChemicalAction:
    """An explicitly calibrated linear chemical or isotope-mass measurement."""

    action_id: str
    measurement_type: str
    target_well: str
    response: Mapping[str, float]
    error_bound: float
    lab_cost: float
    field_travel_cost: float = 0.0
    isotope_system: Optional[str] = None
    element: Optional[str] = None
    units: str = ""
    accessibility: float = 1.0
    feasible: bool = True
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        response = {str(name): float(value) for name, value in self.response.items()}
        if not response or any(not math.isfinite(value) for value in response.values()):
            raise ValueError("ChemicalAction requires an explicit response with finite linear coefficients.")
        if not math.isfinite(float(self.error_bound)) or self.error_bound <= 0:
            raise ValueError("ChemicalAction error_bound must be finite and positive.")
        object.__setattr__(self, "response", response)
        object.__setattr__(self, "metadata", dict(self.metadata))

    def to_acmd_action(self, state_names: Sequence[str]) -> ACMDAction:
        unknown = set(self.response) - set(state_names)
        if unknown:
            raise ValueError(f"Action {self.action_id!r} refers to unknown states: {sorted(unknown)}")
        return ACMDAction(
            action_id=self.action_id,
            measurement_type=self.measurement_type,
            target_id=self.target_well,
            response=[self.response.get(name, 0.0) for name in state_names],
            error_bound=self.error_bound,
            standard_deviation=self.error_bound / 1.96,
            cost=self.lab_cost,
            travel_cost=self.field_travel_cost,
            accessibility=self.accessibility,
            feasible=self.feasible,
            metadata={
                **self.metadata,
                "isotope_system": self.isotope_system,
                "element": self.element,
                "units": self.units,
                "response_semantics": "declared_linear_observation",
            },
        )


def build_isotope_mass_action(
    polytope: ChemicalPolytope,
    *,
    action_id: str,
    target_well: str,
    isotope_system: str,
    element: str,
    isotope_atoms_per_species: float,
    source_atom_fractions: Mapping[str, float],
    error_bound: float,
    lab_cost: float,
    field_travel_cost: float = 0.0,
    accessibility: float = 1.0,
) -> ChemicalAction:
    """Build an isotope-*mass* row from declared source atom fractions.

    The observation is heavy-isotope concentration, not a raw delta or ratio.
    Fractions must be calibrated for every source that contributes the element.
    The caller declares the number of isotope-bearing atoms per aqueous species
    (for example, three oxygen atoms per nitrate molecule).
    Negative element stoichiometry is refused because its isotope composition
    generally requires a separate fractionation or sink model.
    """
    atom_count = float(isotope_atoms_per_species)
    if not math.isfinite(atom_count) or atom_count <= 0.0 or not atom_count.is_integer():
        raise ValueError("isotope_atoms_per_species must be a positive integer.")
    response: dict[str, float] = {}
    needed: set[str] = set()
    for reaction in polytope.reaction_names:
        amount = float(polytope.reaction_stoichiometry[reaction].get(element, 0.0))
        if amount < 0:
            raise ValueError("Isotope-mass helper does not model isotope removal or fractionation.")
        if amount > 0:
            state = f"extent:{reaction}"
            needed.add(state)
            response[state] = amount
    for endmember in polytope.endmember_names:
        amount = float(polytope.endmember_concentrations[endmember].get(element, 0.0))
        if amount > 0:
            state = f"mix:{endmember}"
            needed.add(state)
            response[state] = amount
    if not needed or set(source_atom_fractions) != needed:
        raise ValueError(f"source_atom_fractions must name exactly these contributing states: {sorted(needed)}")
    fractions = {name: float(source_atom_fractions[name]) for name in needed}
    if any(not math.isfinite(value) or not 0.0 <= value <= 1.0 for value in fractions.values()):
        raise ValueError("Source atom fractions must be finite values between zero and one.")
    response = {state: atom_count * amount * fractions[state] for state, amount in response.items()}
    return ChemicalAction(
        action_id=action_id,
        measurement_type=f"isotope_mass_{isotope_system}",
        target_well=target_well,
        response=response,
        error_bound=error_bound,
        lab_cost=lab_cost,
        field_travel_cost=field_travel_cost,
        accessibility=accessibility,
        isotope_system=isotope_system,
        element=element,
        units=f"{polytope.concentration_unit} heavy-isotope equivalent",
        metadata={
            "source_atom_fractions": dict(sorted(fractions.items())),
            "isotope_atoms_per_species": atom_count,
            "assumption": "conservative_source_isotope_mass_no_fractionation",
        },
    )


class ChemACMD:
    """Chemical facade for ACMD's robust, set-based active design loop."""

    @staticmethod
    def create_loop(
        polytope: ChemicalPolytope,
        candidates: Sequence[ChemicalAction],
        target_coefficients: Mapping[str, float],
        target_tolerance: float,
        *,
        target_name: str,
        target_units: Optional[str] = None,
        budget: Optional[float] = None,
        cost_exponent: float = 1.0,
        feasibility_tolerance: float = 1.0e-7,
    ) -> ACMDLoop:
        names = polytope.state_names
        unknown = set(target_coefficients) - set(names)
        if unknown or not target_coefficients:
            raise ValueError(f"Target must name known states; unknown: {sorted(unknown)}")
        q = [float(target_coefficients.get(name, 0.0)) for name in names]
        if not all(math.isfinite(value) for value in q):
            raise ValueError("Target coefficients must be finite.")
        actions = [candidate.to_acmd_action(names) for candidate in candidates]
        if len({action.action_id for action in actions}) != len(actions):
            raise ValueError("Chemical action IDs must be unique.")
        payload = {
            "model": "linear_reaction_mixing_v1",
            "state_space": polytope.state_space.to_dict(),
            "concentration_unit": polytope.concentration_unit,
            "actions": [
                {
                    "action_id": action.action_id,
                    "response": action.response.tolist(),
                    "error_bound": action.error_bound,
                    "cost": action.cost,
                    "travel_cost": action.travel_cost,
                    "accessibility": action.accessibility,
                    "target_id": action.target_id,
                    "feasible": action.feasible,
                    "measurement_type": action.measurement_type,
                    "units": action.metadata.get("units", ""),
                    "isotope_system": action.metadata.get("isotope_system"),
                    "element": action.metadata.get("element"),
                    "source_atom_fractions": action.metadata.get("source_atom_fractions"),
                    "isotope_atoms_per_species": action.metadata.get("isotope_atoms_per_species"),
                    "assumption": action.metadata.get("assumption"),
                }
                for action in actions
            ],
            "target": q,
            "target_name": target_name,
            "target_units": target_units or polytope.concentration_unit,
            "target_tolerance": float(target_tolerance),
            "budget": budget,
            "cost_exponent": cost_exponent,
            "feasibility_tolerance": feasibility_tolerance,
        }
        problem_hash = hashlib.sha256(
            json.dumps(payload, sort_keys=True, allow_nan=False).encode("utf-8")
        ).hexdigest()
        return ACMD.create_loop(
            age_grid_years=np.arange(len(names), dtype=float),
            initial_constraints=(),
            candidates=actions,
            target_functional=AgeFunctional(
                name=target_name,
                coefficients=q,
                units=target_units or polytope.concentration_unit,
            ),
            target_tolerance=target_tolerance,
            budget=budget,
            cost_exponent=cost_exponent,
            mode=ACMDMode.ROBUST_MINIMAX,
            feasibility_tolerance=feasibility_tolerance,
            state_space=polytope.state_space,
            metadata={
                "problem_hash": problem_hash,
                "model": "linear_reaction_mixing_v1",
                "state_names": names,
                "concentration_unit": polytope.concentration_unit,
                "claim_scope": "conditional_on_declared_stoichiometry_bounds_and_linear_measurement_rows",
            },
        )
