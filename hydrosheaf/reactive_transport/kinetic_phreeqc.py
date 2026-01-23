"""
PHREEQC kinetics integration for forward reactive transport validation.
"""

from typing import Dict, List, Optional

from ..config import Config
from . import KineticParameters
from .rate_laws import (
    DEFAULT_KINETIC_PARAMS,
    RATE_LAW_TEMPLATES,
    apply_temperature_correction,
)


def build_kinetic_block(
    reaction_labels: List[str],
    extents: List[float],
    residence_time_days: float,
    kinetic_params: Optional[Dict[str, KineticParameters]] = None,
    temperature_c: float = 25.0,
) -> str:
    """
    Build PHREEQC KINETICS and RATES blocks from inverse results.

    Parameters
    ----------
    reaction_labels : List[str]
        Names of reactions from Hydrosheaf inverse model
    extents : List[float]
        Reaction extents in mmol/L
    residence_time_days : float
        Travel time for the edge
    kinetic_params : Optional[Dict]
        Override rate constants and surface areas
    temperature_c : float
        Temperature in Celsius for Arrhenius correction

    Returns
    -------
    str
        PHREEQC input string with KINETICS and RATES blocks

    Mathematical Implementation
    ---------------------------
    For each reaction k with extent ξ_k:

    1. Convert residence time to seconds:
       τ_s = residence_time_days * 86400

    2. Estimate average rate:
       r̄_k = ξ_k / τ_s  (mmol/L/s)

    3. Build RATES block with TST rate law:
       rate = k_k * A_k * (1 - 10^SI_k)

    4. Set initial moles m0 such that:
       - For dissolution (ξ_k > 0): m0 > ξ_k (sufficient mineral)
       - For precipitation (ξ_k < 0): m0 = 0 (mineral forms)

    5. Set -steps τ_s in 100 steps
    """
    if kinetic_params is None:
        kinetic_params = {}

    tau_seconds = residence_time_days * 86400.0

    # Build RATES block
    rates_block = "RATES\n"
    kinetics_entries = []

    for label, extent in zip(reaction_labels, extents):
        # Skip if extent is negligible
        if abs(extent) < 1e-6:
            continue

        # Get kinetic parameters
        if label in kinetic_params:
            params = kinetic_params[label]
        elif label in DEFAULT_KINETIC_PARAMS:
            params = DEFAULT_KINETIC_PARAMS[label]
        else:
            # Skip reactions without kinetic parameters
            continue

        # Apply temperature correction
        k_T = apply_temperature_correction(params, temperature_c)

        # Get rate law template
        if label in RATE_LAW_TEMPLATES:
            rates_block += RATE_LAW_TEMPLATES[label]
            rates_block += "\n"

            # Build KINETICS entry
            # Initial moles: for dissolution, need sufficient mineral
            # For precipitation, start with zero
            if extent > 0:
                # Dissolution: set m0 = extent * 10 (excess)
                m0 = extent * 10.0 / 1000.0  # Convert mmol/L to mol/L
            else:
                # Precipitation: start with zero
                m0 = 0.0

            kinetics_entry = f"""    {label}
        -m0 {m0:.6e}
        -m {m0:.6e}
        -parms {k_T:.6e} {params.surface_area:.6e}
        -formula {_get_mineral_formula(label)}
"""
            kinetics_entries.append(kinetics_entry)

    # Build KINETICS block
    if kinetics_entries:
        kinetics_block = "KINETICS 1\n"
        kinetics_block += "".join(kinetics_entries)

        # Set time steps
        n_steps = 100
        kinetics_block += f"    -steps {tau_seconds:.6e} in {n_steps} steps\n"
        kinetics_block += "    -step_divide 10\n"
    else:
        kinetics_block = ""

    return rates_block + "\n" + kinetics_block


def _get_mineral_formula(reaction_name: str) -> str:
    """
    Get chemical formula for mineral reactions.

    Parameters
    ----------
    reaction_name : str
        Reaction name

    Returns
    -------
    str
        Chemical formula
    """
    formulas = {
        "calcite": "CaCO3",
        "dolomite": "CaMg(CO3)2",
        "gypsum": "CaSO4:2H2O",
        "halite": "NaCl",
        "fluorite": "CaF2",
        "albite": "NaAlSi3O8",
        "anorthite": "CaAl2Si2O8",
        "pyrite_oxidation_aerobic": "FeS2",
    }
    return formulas.get(reaction_name, "")


def run_phreeqc_kinetic(
    initial_solution: Dict[str, float],
    kinetics_block: str,
    residence_time_days: float,
    config: Config,
    n_output_steps: int = 100,
) -> Dict[str, object]:
    """
    Run PHREEQC kinetic simulation.

    Parameters
    ----------
    initial_solution : Dict[str, float]
        Upstream water composition (mmol/L for each ion)
    kinetics_block : str
        KINETICS and RATES blocks from build_kinetic_block
    residence_time_days : float
        Simulation time
    config : Config
        Hydrosheaf configuration
    n_output_steps : int
        Number of output time steps

    Returns
    -------
    Dict[str, object]
        {
            "final_composition": List[float],  # x_v predicted
            "time_series": Dict[str, List[float]],  # concentration vs time
            "si_series": Dict[str, List[float]],  # SI vs time
            "success": bool,
            "error_message": Optional[str]
        }

    Implementation
    --------------
    1. Build SOLUTION block from initial_solution (use existing build_solution_block)

    2. Append kinetics_block

    3. Add SELECTED_OUTPUT with:
       -time true
       -totals Ca Mg Na K Cl S(6) N(5) Alkalinity Fe P
       -saturation_indices Calcite Dolomite Gypsum ...

    4. Run via phreeqpython or subprocess (based on config.phreeqc_mode)

    5. Parse output time series

    6. Return final state and trajectories
    """
    from ..phreeqc.templates import build_solution_block

    result = {
        "final_composition": [],
        "time_series": {},
        "si_series": {},
        "success": False,
        "error_message": None,
    }

    try:
        # Build PHREEQC input
        input_str = ""

        # SOLUTION block
        solution_dict = {
            "Ca": initial_solution.get("Ca", 0.0),
            "Mg": initial_solution.get("Mg", 0.0),
            "Na": initial_solution.get("Na", 0.0),
            "K": initial_solution.get("K", 0.0),
            "Cl": initial_solution.get("Cl", 0.0),
            "SO4": initial_solution.get("SO4", 0.0),
            "NO3": initial_solution.get("NO3", 0.0),
            "HCO3": initial_solution.get("HCO3", 0.0),
            "F": initial_solution.get("F", 0.0),
            "Fe": initial_solution.get("Fe", 0.0),
            "PO4": initial_solution.get("PO4", 0.0),
        }

        # Add pH and temp
        solution_dict["pH"] = initial_solution.get("pH", 7.0)
        solution_dict["temp_c"] = initial_solution.get(
            "temp_c",
            initial_solution.get("temperature", config.temp_default_c),
        )

        input_str += build_solution_block(solution_dict, 1, config.temp_default_c)
        input_str += "\n"

        # Add kinetics
        input_str += kinetics_block
        input_str += "\n"

        # Add SELECTED_OUTPUT
        # In subprocess mode we must explicitly write to a known file.
        selected_output_file = None
        if config.phreeqc_mode == "subprocess":
            selected_output_file = "selected_output.csv"

        if selected_output_file:
            input_str += f"""
SELECTED_OUTPUT 1
    -reset false
    -time true
    -csv true
    -file {selected_output_file}
    -totals Ca Mg Na K Cl S(6) N(5) Alkalinity F Fe P
    -si Calcite Dolomite Gypsum Halite Fluorite Albite Anorthite
    -step true
"""
        else:
            input_str += """
SELECTED_OUTPUT 1
    -reset false
    -time true
    -totals Ca Mg Na K Cl S(6) N(5) Alkalinity F Fe P
    -si Calcite Dolomite Gypsum Halite Fluorite Albite Anorthite
    -step true
"""

        # Run PHREEQC
        if config.phreeqc_mode == "phreeqpython":
            try:
                import phreeqpython as pp  # type: ignore[import-untyped]

                phreeqc = pp.PhreeqPython(database=config.phreeqc_database)
                output = phreeqc.run_string(input_str)

                # Extract final composition
                final_sol = output.selected_output()
                if final_sol and len(final_sol) > 0:
                    final_row = final_sol.iloc[-1]

                    # Extract final concentrations
                    ion_order = config.ion_order
                    final_comp = []

                    for ion in ion_order:
                        if ion == "HCO3":
                            final_comp.append(
                                float(final_row.get("Alk(mol/kgw)", 0.0)) * 1000.0
                            )
                        elif ion == "SO4":
                            final_comp.append(
                                float(final_row.get("S(6)(mol/kgw)", 0.0)) * 1000.0
                            )
                        elif ion == "NO3":
                            final_comp.append(
                                float(final_row.get("N(5)(mol/kgw)", 0.0)) * 1000.0
                            )
                        else:
                            final_comp.append(
                                float(final_row.get(f"{ion}(mol/kgw)", 0.0)) * 1000.0
                            )

                    result["final_composition"] = final_comp

                    # Extract SI series
                    si_minerals = [
                        "Calcite",
                        "Dolomite",
                        "Gypsum",
                        "Halite",
                        "Fluorite",
                        "Albite",
                        "Anorthite",
                    ]
                    for mineral in si_minerals:
                        si_col = f"si_{mineral}"
                        if si_col in final_sol.columns:
                            result["si_series"][mineral] = final_sol[si_col].tolist()

                    result["success"] = True

            except ImportError:
                result["error_message"] = "phreeqpython not available"
                return result
            except Exception as e:
                result["error_message"] = f"PHREEQC error: {str(e)}"
                return result

        else:
            # Subprocess mode for kinetic simulation
            import csv
            import subprocess
            import tempfile
            from pathlib import Path

            if not config.phreeqc_executable:
                result["error_message"] = "phreeqc_executable not configured"
                return result

            if not selected_output_file:
                result["error_message"] = "Internal error: selected_output_file not set for subprocess mode"
                return result

            try:
                with tempfile.TemporaryDirectory() as temp_dir:
                    temp_dir_path = Path(temp_dir)
                    input_path = temp_dir_path / "input.pqi"
                    output_path = temp_dir_path / "output.txt"
                    selected_path = temp_dir_path / selected_output_file

                    input_path.write_text(input_str, encoding="utf-8")

                    subprocess.run(
                        [config.phreeqc_executable, str(input_path), str(output_path)],
                        cwd=str(temp_dir_path),
                        check=True,
                        capture_output=True,
                    )

                    if not selected_path.exists():
                        result["error_message"] = f"Selected output not generated: {selected_output_file}"
                        return result

                    with open(selected_path, "r", encoding="utf-8", newline="") as handle:
                        reader = csv.DictReader(handle)
                        rows = list(reader)

                    if not rows:
                        result["error_message"] = "Selected output file is empty"
                        return result

                    final_row = rows[-1]
                    ion_order = config.ion_order

                    mapping = {
                        "Ca": "Ca(mol/kgw)",
                        "Mg": "Mg(mol/kgw)",
                        "Na": "Na(mol/kgw)",
                        "K": "K(mol/kgw)",
                        "Cl": "Cl(mol/kgw)",
                        "SO4": "S(6)(mol/kgw)",
                        "NO3": "N(5)(mol/kgw)",
                        "HCO3": "Alk(mol/kgw)",
                        "F": "F(mol/kgw)",
                        "Fe": "Fe(mol/kgw)",
                        "PO4": "P(mol/kgw)",
                    }

                    final_comp = []
                    for ion in ion_order:
                        col = mapping.get(ion, f"{ion}(mol/kgw)")
                        raw = final_row.get(col)
                        try:
                            val_mol = float(raw) if raw not in (None, "") else 0.0
                        except ValueError:
                            val_mol = 0.0
                        final_comp.append(val_mol * 1000.0)

                    result["final_composition"] = final_comp

                    # SI columns are typically named like si_Calcite
                    for key, raw in final_row.items():
                        if not key.startswith("si_"):
                            continue
                        mineral = key[3:]
                        try:
                            result["si_series"].setdefault(mineral, []).append(float(raw))
                        except (TypeError, ValueError):
                            continue

                    result["success"] = True
                    return result

            except FileNotFoundError:
                result["error_message"] = f"PHREEQC executable not found: {config.phreeqc_executable}"
                return result
            except subprocess.CalledProcessError as exc:
                stderr = exc.stderr.decode("utf-8", errors="ignore") if exc.stderr else ""
                result["error_message"] = f"PHREEQC subprocess failed: {stderr.strip()}"
                return result
            except Exception as e:
                result["error_message"] = f"Subprocess error: {str(e)}"
                return result

    except Exception as e:

        result["error_message"] = f"Error building kinetic input: {str(e)}"

    return result
