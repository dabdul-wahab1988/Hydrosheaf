"""
PEST++ Runner Module.
Handles execution of PEST++ binaries for calibration.
"""

import os
import sys
import shutil
import pickle
import subprocess
import platform
import zipfile
import urllib.request
import socket
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
import re
import pandas as pd
import numpy as np

from ..definitions import AdjustableParameter, Observation
from ...log import get_logger

logger = get_logger("calibration.pestpp")

PESTPP_VERSION = "5.2.25"
PESTPP_BASE_URL = "https://github.com/usgs/pestpp/releases/download"

def get_bin_dir(version: Optional[str] = None) -> Path:
    """Return the local bin directory, optionally version-specific."""
    if version:
        return Path("bin") / version
    return Path("bin")

def get_executable_suffix() -> str:
    return ".exe" if os.name == "nt" else ""

def get_executable_path(name: str, version: Optional[str] = None) -> str:
    """
    Find the path to a PEST++ executable.
    Downloads it if not found in local bin/.
    """
    v = version or PESTPP_VERSION
    bin_dir = get_bin_dir(v)
    suffix = get_executable_suffix()
    
    # Check if name already has suffix
    if not name.lower().endswith(suffix.lower()) and suffix:
        name += suffix
        
    # Check local bin
    local_path = bin_dir / name
    if local_path.exists():
        return str(local_path.absolute())
        
    # Check PATH
    path_exe = shutil.which(name)
    if path_exe:
        return path_exe
        
    # Not found -> Download
    logger.info(f"Executable {name} not found. Attempting to download PEST++ binaries for version {v}...")
    try:
        download_pestpp_binaries(bin_dir, v)
    except Exception as e:
        logger.error(f"Failed to download PEST++ binaries: {e}")
        # Continue to check if it exists anyway (e.g. partial success)

    if local_path.exists():
        return str(local_path.absolute())
    
    return name

def download_pestpp_binaries(target_dir: Path, version: str):
    """Download and extract PEST++ binaries."""
    target_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine URL based on platform
    system = platform.system()
    if system == "Windows":
        filename = f"pestpp-{version}-win.zip"
        zip_path_internal = "bin"
    elif system == "Linux":
        filename = f"pestpp-{version}-linux.zip"
        zip_path_internal = "bin"
    elif system == "Darwin":
        filename = f"pestpp-{version}-mac.zip"
        zip_path_internal = "bin"
    else:
        raise OSError(f"Unsupported platform for auto-download: {system}")

    url = f"{PESTPP_BASE_URL}/{version}/{filename}"
    temp_zip = target_dir / "temp_pestpp.zip"
    
    logger.info(f"Downloading from {url}...")
    try:
        urllib.request.urlretrieve(url, temp_zip)
    except Exception as e:
        raise RuntimeError(f"Download failed: {e}")
    
    logger.info("Extracting binaries...")
    try:
        with zipfile.ZipFile(temp_zip, 'r') as zip_ref:
            # Filter files to extract (only executables in bin/)
            for member in zip_ref.namelist():
                # Skip directories
                if member.endswith("/"):
                    continue
                    
                if member.startswith(zip_path_internal) and (member.endswith(".exe") or "." not in Path(member).name):
                     # Extract to target_dir, flattening structure
                     try:
                         source = zip_ref.open(member)
                         target_name = Path(member).name
                         if not target_name: continue
                         
                         with open(target_dir / target_name, "wb") as target:
                             shutil.copyfileobj(source, target)
                         
                         # Make executable on Unix
                         if system != "Windows":
                             os.chmod(target_dir / target_name, 0o755)
                     except Exception as e:
                         logger.warning(f"Failed to extract {member}: {e}")
    except Exception as e:
        raise RuntimeError(f"Extraction failed: {e}")
    finally:
        # cleanup
        if temp_zip.exists():
            try:
                os.remove(temp_zip)
            except Exception as e:
                logger.warning(f"Could not remove temp zip: {e}")

    logger.info("PEST++ binaries installed successfully.")



class PestControlData:
    """Represents the '* control data' section of a PEST control file."""

    def __init__(
        self,
        npar: int = 0,
        nobs: int = 0,
        npargp: int = 1,
        nprior: int = 0,
        nobsgp: int = 1,
    ):
        self.restart_estimation: str = "restart estimation"
        self.npar = npar
        self.nobs = nobs
        self.npargp = npargp
        self.nprior = nprior
        self.nobsgp = nobsgp
        self.ntplfle = 1
        self.ninsfle = 1
        self.precis = "single"
        self.dpoint = "point"
        self.numcom = 1
        self.jacfile = 0
        self.messfile = 0
        self.rlambda1 = 10.0
        self.rlamfac = -3.0
        self.phiratsuf = 0.3
        self.phiredlam = 0.01
        self.numlam = 10
        self.relparmax = 10.0
        self.facparmax = 10.0
        self.facorig = 0.001
        self.phiredstp = 0.1
        self.nrelpar = 5
        self.nphistp = 5
        self.nphinorel = 5
        self.noptmax = 100
        self.phiredstp_opt = 0.0
        self.nphistp_opt = 0
        self.nphinorel_opt = 0
        self.phiredphistp = 0.0
        self.nphistpphistp = 0

    def format(self) -> str:
        lines = [
            "pcf",
            "* control data",
            f"{self.restart_estimation}",
            f"{self.npar} {self.nobs} {self.npargp} {self.nprior} {self.nobsgp}",
            f"{self.ntplfle} {self.ninsfle} {self.precis} {self.dpoint} {self.numcom} {self.jacfile} {self.messfile}",
            f"{self.rlambda1} {self.rlamfac} {self.phiratsuf} {self.phiredlam} {self.numlam}",
            f"{self.relparmax} {self.facparmax} {self.facorig}",
            f"{self.phiredstp} {self.nrelpar} {self.nphistp} {self.nphinorel}",
            f"{self.noptmax} {self.phiredstp_opt} {self.nphistp_opt} {self.nphinorel_opt} {self.phiredphistp} {self.nphistpphistp}",
        ]
        return "\n".join(lines) + "\n"


class PestParameter:
    """Represents a single parameter entry in PEST parameter data."""

    def __init__(
        self,
        name: str,
        trans: str,
        change_limit: str,
        value: float,
        lower_bound: float,
        upper_bound: float,
        group: str = "pargp",
        scale: float = 1.0,
        offset: float = 0.0,
        prior_mean: Optional[float] = None,
        prior_sigma: Optional[float] = None,
        tied_to: Optional[str] = None,
    ):
        self.name = name
        self.trans = trans
        self.change_limit = change_limit
        self.value = value
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.group = group
        self.scale = scale
        self.offset = offset
        self.prior_mean = prior_mean
        self.prior_sigma = prior_sigma
        self.tied_to = tied_to

    @classmethod
    def from_adjustable_parameter(cls, p: AdjustableParameter) -> "PestParameter":
        if getattr(p, "fixed", False):
            trans = "fixed"
        elif getattr(p, "tied_to", None):
            trans = "tied"
        else:
            trans = "log" if p.log_transform else "none"
        group = p.group if (p.group and p.group != "default") else "pargp"
        return cls(
            name=p.name,
            trans=trans,
            change_limit="factor",
            value=p.value,
            lower_bound=p.lower_bound,
            upper_bound=p.upper_bound,
            group=group,
            prior_mean=p.prior_mean,
            prior_sigma=p.prior_sigma,
            tied_to=p.tied_to,
        )

    def format(self) -> str:
        base = f"{self.name} {self.trans} {self.change_limit} {self.value} {self.lower_bound} {self.upper_bound} {self.group} {self.scale} {self.offset}"
        if self.trans == "tied" and self.tied_to:
            base += f" {self.tied_to}"
        return base


class PestObservation:
    """Represents a single observation entry in PEST observation data."""

    def __init__(
        self, name: str, value: float, weight: float = 1.0, group: str = "obsgp"
    ):
        self.name = name
        self.value = value
        self.weight = weight
        self.group = group

    @classmethod
    def from_observation(cls, o: Observation) -> "PestObservation":
        group = o.group if (o.group and o.group != "default") else "obsgp"
        return cls(
            name=o.name,
            value=o.value,
            weight=o.weight,
            group=group,
        )

    def format(self) -> str:
        return f"{self.name} {self.value} {self.weight} {self.group}"


class PestControlFile:
    """Manages programmatic building, modifying, writing, and loading of PEST control files."""

    def __init__(self):
        self.control_data = PestControlData()
        self.parameter_groups: List[str] = [
            "pargp relative 0.01 0.0 switch 2.0 parabolic"
        ]
        self.parameters: List[PestParameter] = []
        self.observation_groups: List[str] = ["obsgp"]
        self.observations: List[PestObservation] = []
        self.model_command_line: str = "python worker.py"
        
        # Support multiple template and instruction files
        self.tpl_files: List[str] = ["case.tpl"]
        self.ins_files: List[str] = ["case.ins"]
        self.model_input_files: List[str] = ["params.csv"]
        self.model_output_files: List[str] = ["results.dat"]
        self.pestpp_options: Dict[str, Any] = {}

    @property
    def tpl_file(self) -> str:
        return self.tpl_files[0] if self.tpl_files else ""
    @tpl_file.setter
    def tpl_file(self, val: str):
        if self.tpl_files:
            self.tpl_files[0] = val
        else:
            self.tpl_files.append(val)

    @property
    def ins_file(self) -> str:
        return self.ins_files[0] if self.ins_files else ""
    @ins_file.setter
    def ins_file(self, val: str):
        if self.ins_files:
            self.ins_files[0] = val
        else:
            self.ins_files.append(val)

    @property
    def model_input_file(self) -> str:
        return self.model_input_files[0] if self.model_input_files else ""
    @model_input_file.setter
    def model_input_file(self, val: str):
        if self.model_input_files:
            self.model_input_files[0] = val
        else:
            self.model_input_files.append(val)

    @property
    def model_output_file(self) -> str:
        return self.model_output_files[0] if self.model_output_files else ""
    @model_output_file.setter
    def model_output_file(self, val: str):
        if self.model_output_files:
            self.model_output_files[0] = val
        else:
            self.model_output_files.append(val)

    def add_parameter(self, param: PestParameter):
        self.parameters.append(param)
        self.control_data.npar = len(self.parameters)

    def add_observation(self, obs: PestObservation):
        self.observations.append(obs)
        self.control_data.nobs = len(self.observations)

    def write(self, filepath: Path):
        """Write the complete control file to disk."""
        # Auto-collect parameter groups
        used_pargps = set(pg.split()[0].lower() for pg in self.parameter_groups)
        for p in self.parameters:
            if p.group.lower() not in used_pargps:
                self.parameter_groups.append(f"{p.group} relative 0.01 0.0 switch 2.0 parabolic")
                used_pargps.add(p.group.lower())

        # Auto-collect observation groups
        used_obsgps = set(og.lower() for og in self.observation_groups)
        for o in self.observations:
            if o.group.lower() not in used_obsgps:
                self.observation_groups.append(o.group)
                used_obsgps.add(o.group.lower())

        # Generate Prior Information equations from parameters
        prior_lines = []
        for p in self.parameters:
            if p.prior_mean is not None and p.prior_sigma is not None:
                pgp = "prior"
                if pgp.lower() not in used_obsgps:
                    self.observation_groups.append(pgp)
                    used_obsgps.add(pgp.lower())
                
                if p.trans == "log":
                    val = np.log10(max(1e-12, p.prior_mean))
                    eq = f"1.0 * log({p.name})"
                else:
                    val = p.prior_mean
                    eq = f"1.0 * {p.name}"
                
                weight = 1.0 / max(1e-12, p.prior_sigma)
                pi_name = f"pi_{p.name}"[:20]
                prior_lines.append(f"{pi_name} {eq} = {val} {weight} {pgp}")

        # Sync control data dimensions
        self.control_data.npar = len(self.parameters)
        self.control_data.nobs = len(self.observations)
        self.control_data.ntplfle = len(self.tpl_files)
        self.control_data.ninsfle = len(self.ins_files)
        self.control_data.npargp = len(self.parameter_groups)
        self.control_data.nobsgp = len(self.observation_groups)
        self.control_data.nprior = len(prior_lines)

        with open(filepath, "w") as f:
            f.write(self.control_data.format())

            f.write("* parameter groups\n")
            for pg in self.parameter_groups:
                f.write(f"{pg}\n")

            f.write("* parameter data\n")
            for p in self.parameters:
                f.write(f"{p.format()}\n")

            f.write("* observation groups\n")
            for og in self.observation_groups:
                f.write(f"{og}\n")

            f.write("* observation data\n")
            for o in self.observations:
                f.write(f"{o.format()}\n")

            if prior_lines:
                f.write("* prior information\n")
                for pl in prior_lines:
                    f.write(f"{pl}\n")

            f.write("* model command line\n")
            f.write(f"{self.model_command_line}\n")

            f.write("* model input/output\n")
            for tpl, inp in zip(self.tpl_files, self.model_input_files):
                f.write(f"{tpl} {inp}\n")
            for ins, out in zip(self.ins_files, self.model_output_files):
                f.write(f"{ins} {out}\n")

            f.write("* pestpp options\n")
            for k, v in self.pestpp_options.items():
                key = k if k.startswith("++") else f"++{k}"
                f.write(f"{key}({v})\n")

    def write_tpl(self, filepath: Path):
        """Write the PEST template file (.tpl)."""
        field_width = 20
        with open(filepath, "w") as f:
            f.write("ptf ~\n")
            for p in self.parameters:
                placeholder = f"{p.name:<{field_width}}"
                f.write(f"{p.name},~{placeholder}~\n")

    def write_ins(self, filepath: Path):
        """Write the PEST instruction file (.ins)."""
        with open(filepath, "w") as f:
            f.write("pif @\n")
            for o in self.observations:
                f.write(f"l1 w !{o.name}!\n")


def write_pst_file(
    filepath: Path,
    params: List[AdjustableParameter],
    obs: List[Observation],
    tpl_file: str,
    ins_file: str,
    model_input_file: str,
    model_output_file: str,
    pestpp_options: Dict[str, Any],
    noptmax: int,
):
    """Write PEST control file (.pst) using PestControlFile class."""
    pst = PestControlFile()
    pst.tpl_file = tpl_file
    pst.ins_file = ins_file
    pst.model_input_file = model_input_file
    pst.model_output_file = model_output_file
    pst.pestpp_options = pestpp_options
    pst.control_data.noptmax = noptmax

    for p in params:
        pst.add_parameter(PestParameter.from_adjustable_parameter(p))
    for o in obs:
        pst.add_observation(PestObservation.from_observation(o))

    pst.write(filepath)


def write_tpl_file(
    filepath: Path, params: List[AdjustableParameter], model_input_file: str
):
    """Write template file mapping params to model input CSV using PestControlFile class."""
    pst = PestControlFile()
    pst.model_input_file = model_input_file
    for p in params:
        pst.add_parameter(PestParameter.from_adjustable_parameter(p))
    pst.write_tpl(filepath)


def write_ins_file(filepath: Path, obs: List[Observation]):
    """Write instruction file reading model output CSV using PestControlFile class."""
    pst = PestControlFile()
    for o in obs:
        pst.add_observation(PestObservation.from_observation(o))
    pst.write_ins(filepath)


def generate_worker_script(filepath: Path):
    """Generate the worker.py script."""
    script_content = r'''
import sys
import os
import pickle
from pathlib import Path
import pandas as pd
import numpy as np

def find_repo_root(start_dir: str) -> Path:
    current = Path(start_dir).resolve()
    for _ in range(6):
        if (current / "hydrosheaf").is_dir():
            return current
        if current.parent == current:
            break
        current = current.parent
    return Path(start_dir).resolve()

repo_root = find_repo_root(os.getcwd())
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
try:
    import hydrosheaf
except ImportError:
    pass

def run():
    def load_obs_names() -> list:
        names = []
        names_path = Path("obs_names.txt")
        if names_path.exists():
            try:
                for line in names_path.read_text().splitlines():
                    name = line.strip()
                    if name:
                        names.append(name)
            except Exception:
                return []
        return names

    expected_obs = load_obs_names()
    problem = None
    # 1. Load Problem
    try:
        problem_path = Path("problem.pkl")
        if not problem_path.exists():
            for parent in Path.cwd().resolve().parents:
                candidate = parent / "problem.pkl"
                if candidate.exists():
                    problem_path = candidate
                    break
        with open(problem_path, "rb") as f:
            problem = pickle.load(f)
        expected_obs = [obs.name for obs in problem.get_observations()]
    except Exception as e:
        print(f"Error loading problem: {e}")

    # 2. Read Parameters
    # Format: name,value
    params = {}
    parse_errors = []
    try:
        with open("params.csv", "r") as f:
            for line in f:
                if "," not in line:
                    continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    name = parts[0].strip()
                    value = parts[1].strip()
                    try:
                        params[name] = float(value)
                    except ValueError:
                        parse_errors.append(name)
    except Exception as e:
        print(f"Error reading params: {e}")
        params = {}

    # 3. Run Model
    results = {}
    if problem is not None:
        try:
            param_defs = {p.name: p for p in problem.get_parameters()}
            for name, pdef in param_defs.items():
                value = params.get(name, pdef.value)
                if value < pdef.lower_bound:
                    value = pdef.lower_bound
                if value > pdef.upper_bound:
                    value = pdef.upper_bound
                params[name] = value
            if parse_errors:
                print(f"Warning: using defaults for {', '.join(parse_errors)}")
            results = problem.run_model(params)
        except Exception as e:
            print(f"Error running model: {e}")
            results = {}

    # 4. Write Results
    # Format: name value
    # Must match the order in .ins file if using "l1" instructions sequentially
    # The .ins file iterates over observations.
    # To be safe, we should write them in the order of problem.get_observations()
    # OR the runner.py ensures .ins order matches this output.
    
    # Let's get the expected order from the problem definition
    # (Since problem is loaded, we can ask it)
    try:
        with open("results.dat", "w") as f:
            for obs_name in expected_obs:
                val = results.get(obs_name, -9999.0)
                f.write(f"{obs_name} {val}\n")
    except Exception as e:
        print(f"Error writing results: {e}")

if __name__ == "__main__":
    try:
        run()
    except BaseException as e:
        print(f"Worker fatal error: {e}")
        try:
            names = []
            names_path = Path("obs_names.txt")
            if names_path.exists():
                for line in names_path.read_text().splitlines():
                    name = line.strip()
                    if name:
                        names.append(name)
            with open("results.dat", "w") as f:
                for obs_name in names:
                    f.write(f"{obs_name} -9999.0\n")
        except Exception as write_exc:
            print(f"Worker failed to write results: {write_exc}")
'''
    with open(filepath, "w") as f:
        f.write(script_content)


def write_parameter_ensemble(
    filepath: Path, params: List[AdjustableParameter], n_reals: int
) -> None:
    rng = np.random.default_rng(123)
    header = "real_name," + ",".join(p.name for p in params)
    with open(filepath, "w") as f:
        f.write(f"{header}\n")
        for idx in range(n_reals):
            values = []
            for p in params:
                low = p.lower_bound
                high = p.upper_bound
                if p.log_transform:
                    low = max(low, 1e-12)
                    high = max(high, low * 10.0)
                    log_low = np.log10(low)
                    log_high = np.log10(high)
                    value = 10**rng.uniform(log_low, log_high)
                else:
                    value = rng.uniform(low, high)
                values.append(value)
            values_str = ",".join(f"{value}" for value in values)
            f.write(f"{idx},{values_str}\n")


def write_observation_ensemble(
    filepath: Path, obs: List[Observation], n_reals: int
) -> None:
    rng = np.random.default_rng(123)
    header = "real_name," + ",".join(o.name for o in obs)
    with open(filepath, "w") as f:
        f.write(f"{header}\n")
        for idx in range(n_reals):
            values = []
            for o in obs:
                sigma = 1.0 / max(o.weight, 1e-6)
                if o.weight == 0:
                    sigma = 0.0
                val = rng.normal(o.value, sigma)
                values.append(val)
            values_str = ",".join(f"{value}" for value in values)
            f.write(f"{idx},{values_str}\n")


def generate_sweep_csv(
    filepath: Path,
    params: List[AdjustableParameter],
    n_runs: int = 50,
    method: str = "latin_hypercube"
) -> None:
    rng = np.random.default_rng(123)
    header = "run_id," + ",".join(p.name for p in params)
    
    if method == "latin_hypercube":
        samples = np.zeros((n_runs, len(params)))
        for j in range(len(params)):
            grid = np.linspace(0.0, 1.0, n_runs + 1)
            pts = grid[:-1] + rng.uniform(0.0, 1.0 / n_runs, n_runs)
            rng.shuffle(pts)
            samples[:, j] = pts
    elif method == "grid":
        samples = np.zeros((n_runs, len(params)))
        for j in range(len(params)):
            pts = np.linspace(0.0, 1.0, n_runs)
            samples[:, j] = pts
    else:  # random prior
        samples = rng.uniform(0.0, 1.0, size=(n_runs, len(params)))
        
    with open(filepath, "w") as f:
        f.write(f"{header}\n")
        for idx in range(n_runs):
            values = []
            for j, p in enumerate(params):
                val_01 = samples[idx, j]
                low = p.lower_bound
                high = p.upper_bound
                if p.log_transform:
                    low = max(low, 1e-12)
                    high = max(high, low * 10.0)
                    log_low = np.log10(low)
                    log_high = np.log10(high)
                    val = 10 ** (log_low + val_01 * (log_high - log_low))
                else:
                    val = low + val_01 * (high - low)
                values.append(val)
            values_str = ",".join(f"{v}" for v in values)
            f.write(f"{idx},{values_str}\n")


def _write_da_state_ensemble(
    filepath: Path,
    params: List[AdjustableParameter],
    state_names: List[str],
    state_initial_values: List[float],
    n_reals: int,
) -> None:
    """
    Generate a combined parameter + state variable ensemble CSV for DA.
    Parameters are sampled from bounds; states are initialized from state_initial_values
    with small perturbation.
    """
    rng = np.random.default_rng(456)
    all_names = [p.name for p in params] + list(state_names)
    header = "real_name," + ",".join(all_names)

    # Default state values if not enough provided
    defaults = list(state_initial_values) if state_initial_values else [0.0] * len(state_names)
    while len(defaults) < len(state_names):
        defaults.append(0.0)

    with open(filepath, "w") as f:
        f.write(f"{header}\n")
        for idx in range(n_reals):
            values = []
            # Sample parameters
            for p in params:
                if p.log_transform:
                    low = max(p.lower_bound, 1e-12)
                    log_low = np.log10(low)
                    log_high = np.log10(max(p.upper_bound, low * 10.0))
                    val = 10 ** rng.uniform(log_low, log_high)
                else:
                    val = rng.uniform(p.lower_bound, p.upper_bound)
                values.append(val)
            # Sample state variables (Gaussian noise around initial)
            for s_val in defaults:
                val = max(s_val + rng.normal(0.0, abs(s_val) * 0.1 + 0.01), 0.0)
                values.append(val)
            values_str = ",".join(f"{v:.6f}" for v in values)
            f.write(f"{idx},{values_str}\n")


def parse_pest_residuals(res_path: Path) -> tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
    residuals = {}
    simulated = {}
    phi_by_group = {}
    if not res_path.exists():
        return residuals, simulated, phi_by_group

    with open(res_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                name = parts[0]
                group = parts[1]
                try:
                    meas = float(parts[2])
                    sim = float(parts[3])
                    res = float(parts[4])
                    weight = float(parts[5])
                    residuals[name] = res
                    simulated[name] = sim
                    contrib = (res * weight) ** 2
                    phi_by_group[group] = phi_by_group.get(group, 0.0) + contrib
                except ValueError:
                    continue
    return residuals, simulated, phi_by_group


def parse_pest_sensitivities(sen_path: Path) -> Dict[str, float]:
    sens = {}
    if not sen_path.exists():
        return sens
    with open(sen_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                name = parts[0]
                try:
                    val = float(parts[1])
                    sens[name] = val
                except ValueError:
                    continue
    return sens


def parse_sen_msn(msn_path: Path) -> Dict[str, Dict[str, float]]:
    """Parse PESTPP-SEN Method of Morris summary files (*.msn)."""
    if not msn_path.exists():
        return {}
    try:
        df = pd.read_csv(msn_path)
    except Exception as e:
        logger.warning(f"Failed to parse SEN msn file {msn_path}: {e}")
        return {}

    table: Dict[str, Dict[str, float]] = {}
    for _, row in df.iterrows():
        name = str(row.get("parameter_name", row.get("parameter", row.get("name", "")))).strip()
        if not name:
            continue
        metrics: Dict[str, float] = {}
        for key in ("n_samples", "sen_mean", "sen_mean_abs", "sen_std_dev"):
            if key in row and not pd.isna(row[key]):
                try:
                    metrics[key] = float(row[key])
                except (TypeError, ValueError):
                    continue
        table[name] = metrics
    return table


def parse_pest_identifiability(id_path: Path) -> Dict[str, float]:
    ident = {}
    if not id_path.exists():
        return ident
    with open(id_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                name = parts[0]
                try:
                    val = float(parts[1])
                    ident[name] = val
                except ValueError:
                    continue
    return ident


def parse_ies_csv(csv_path: Path) -> Optional[Dict[str, List[float]]]:
    if not csv_path.exists():
        return None
    try:
        df = pd.read_csv(csv_path)
        if "real_name" in df.columns:
            df = df.drop(columns=["real_name"])
        elif "realname" in df.columns:
            df = df.drop(columns=["realname"])
        return df.to_dict(orient="list")
    except Exception as e:
        logger.warning(f"Failed to parse IES CSV {csv_path}: {e}")
        return None


def _parse_csv_records(csv_path: Path) -> Optional[List[Dict[str, Any]]]:
    if not csv_path.exists():
        return None
    try:
        return pd.read_csv(csv_path).to_dict(orient="records")
    except Exception as e:
        logger.warning(f"Failed to parse CSV records {csv_path}: {e}")
        return None


def _parse_phi_history_csv(csv_path: Path) -> List[float]:
    if not csv_path.exists():
        return []
    try:
        df = pd.read_csv(csv_path)
        real_cols = [
            c for c in df.columns
            if c.lower() not in ("iteration", "iter", "cycle", "real_name", "realname")
        ]
        history = []
        for _, row in df.iterrows():
            vals = pd.to_numeric(row[real_cols], errors="coerce").dropna()
            if not vals.empty:
                history.append(float(vals.mean()))
        return history
    except Exception as e:
        logger.warning(f"Failed to parse phi history {csv_path}: {e}")
        return []


def write_parcov_matrix(filepath: Path, params: List[AdjustableParameter]) -> None:
    """Write a diagonal parameter covariance matrix from parameter bounds."""
    header = "name," + ",".join(p.name for p in params)
    lines = [header]
    for i, p in enumerate(params):
        row = [p.name]
        for j, _ in enumerate(params):
            if i == j:
                # variance = (range / par_sigma_range) ** 2   ; par_sigma_range default = 4
                rng_val = p.upper_bound - p.lower_bound
                if p.log_transform:
                    rng_val = np.log10(p.upper_bound) - np.log10(max(p.lower_bound, 1e-12))
                var = (rng_val / 4.0) ** 2
                row.append(f"{var:.6e}")
            else:
                row.append("0.0")
        lines.append(",".join(row))
    filepath.write_text("\n".join(lines) + "\n")


def write_obscov_matrix(filepath: Path, obs: List[Observation]) -> None:
    """Write a diagonal observation covariance matrix from observation weights."""
    header = "name," + ",".join(o.name for o in obs)
    lines = [header]
    for i, o in enumerate(obs):
        row = [o.name]
        variance = (1.0 / max(o.weight, 1e-6)) ** 2
        for j, _ in enumerate(obs):
            row.append(f"{variance:.6e}" if i == j else "0.0")
        lines.append(",".join(row))
    filepath.write_text("\n".join(lines) + "\n")


def write_localizer_matrix(filepath: Path, params: List[AdjustableParameter], obs: List[Observation]) -> None:
    """Write a default identity localizer matrix (all parameters influence all observations)."""
    header = "name," + ",".join(p.name for p in params)
    lines = [header]
    for o in obs:
        row = [o.name] + ["1.0"] * len(params)
        lines.append(",".join(row))
    filepath.write_text("\n".join(lines) + "\n")


def _pop_option(options: Dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in options:
            return options.pop(name)
        prefixed = name if name.startswith("++") else f"++{name}"
        if prefixed in options:
            return options.pop(prefixed)
    return None


def _set_if_value(options: Dict[str, Any], name: str, value: Any) -> None:
    if value is not None and name not in options and f"++{name}" not in options:
        options[name] = value


def _split_option_tokens(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if str(v).strip()]
    return [tok for tok in re.split(r"[\s,]+", str(value).strip()) if tok]


def _normalize_opt_dec_var_groups(
    raw_value: Any,
    params: List[AdjustableParameter],
) -> Optional[str]:
    tokens = _split_option_tokens(raw_value)
    if tokens and tokens[0].lower() in ("dec_var", "decision_variables", "decision_vars"):
        tokens = tokens[1:]
    if not tokens:
        return None

    group_by_param = {}
    for p in params:
        group = p.group if (p.group and p.group != "default") else "pargp"
        group_by_param[p.name.lower()] = group

    groups: List[str] = []
    seen = set()
    for token in tokens:
        mapped = group_by_param.get(token.lower(), token)
        key = mapped.lower()
        if key in seen:
            continue
        seen.add(key)
        groups.append(mapped)
    return ",".join(groups) if groups else None


def _normalize_constraint_groups(
    raw_value: Any,
    obs: List[Observation],
) -> Optional[str]:
    tokens = _split_option_tokens(raw_value)
    if not tokens:
        return None

    group_by_obs = {}
    for o in obs:
        group_by_obs[o.name.lower()] = o.group if (o.group and o.group != "default") else "obsgp"

    groups: List[str] = []
    seen = set()
    for token in tokens:
        group = group_by_obs.get(token.lower(), token)
        group_l = group.lower()
        if not group_l.startswith(("l_", "less_", "less_than", "g_", "greater_", "greater_than")):
            logger.warning(
                f"Ignoring constraint group '{group}' because PEST++ expects "
                "constraint groups to start with l_, less_, g_, or greater_."
            )
            continue
        if group_l in seen:
            continue
        seen.add(group_l)
        groups.append(group)
    return ",".join(groups) if groups else None


def _normalize_legacy_pestpp_options(
    options: Dict[str, Any],
    engine: str,
    params: List[AdjustableParameter],
    obs: List[Observation],
) -> Optional[int]:
    """
    Translate Hydrosheaf's older convenience option names into real PEST++ names.

    Returns an optional NOPTMAX override when a removed iteration-count alias is used.
    """
    noptmax_override = None

    if engine == "pestpp-sen":
        sen_method = _pop_option(options, "sen_method")
        gsa_method = _pop_option(options, "gsa_method")
        method = gsa_method or sen_method
        if method:
            _set_if_value(options, "gsa_method", method)

        sen_num_samples = _pop_option(options, "sen_num_samples")
        if sen_num_samples is not None:
            target_method = str(method or "morris").lower()
            if target_method == "sobol":
                _set_if_value(options, "gsa_sobol_samples", sen_num_samples)
            else:
                _set_if_value(options, "gsa_morris_r", sen_num_samples)

        _set_if_value(
            options,
            "gsa_morris_delta",
            _pop_option(options, "sen_morris_delta", "morris_delta"),
        )
        _set_if_value(
            options,
            "gsa_sobol_samples",
            _pop_option(options, "sen_sobol_samples"),
        )
        sobol_dist = _pop_option(options, "sen_sobol_par_dist")
        if sobol_dist is not None:
            dist = str(sobol_dist).strip().lower()
            if dist == "uniform":
                dist = "unif"
            elif dist == "normal":
                dist = "norm"
            _set_if_value(options, "gsa_sobol_par_dist", dist)

    if engine in ("pestpp-mou", "pestpp-opt"):
        legacy_risk = _pop_option(options, "risk")
        _set_if_value(options, "opt_risk", legacy_risk)

        dec_groups = _pop_option(options, "dec_var_groups")
        normalized_dec_groups = _normalize_opt_dec_var_groups(dec_groups, params)
        _set_if_value(options, "opt_dec_var_groups", normalized_dec_groups)

        constraint_groups = _pop_option(options, "mou_constraints", "constraint_groups")
        normalized_constraints = _normalize_constraint_groups(constraint_groups, obs)
        _set_if_value(options, "opt_constraint_groups", normalized_constraints)

    if engine == "pestpp-mou":
        mou_max_generations = _pop_option(options, "mou_max_generations")
        if mou_max_generations is not None:
            try:
                noptmax_override = int(mou_max_generations)
            except (TypeError, ValueError):
                logger.warning(f"Ignoring invalid mou_max_generations value: {mou_max_generations}")

    if engine == "pestpp-da":
        da_num_cycles = _pop_option(options, "da_num_cycles")
        if da_num_cycles is not None:
            logger.warning(
                "Ignoring da_num_cycles; PESTPP-DA cycle control is defined by "
                "cycle-aware control-file sections and da_stop_cycle/da_noptmax_schedule."
            )
        da_restart = _pop_option(options, "da_restart_cycle")
        _set_if_value(options, "da_hotstart_cycle", da_restart)

    return noptmax_override


class PestRunner:
    """Wraps PEST++ execution pipeline programmatically."""

    def __init__(
        self,
        problem,
        engine: str = "pestpp-ies",
        work_dir: str = "pest_workspace",
        case_name: str = "case1",
        max_nfev: int = 100,
        n_workers: int = 1,
        pestpp_options: Optional[Dict[str, Any]] = None,
        pestpp_version: Optional[str] = None,
    ):
        self.problem = problem
        self.engine = engine
        self.work_dir = Path(work_dir)
        self.case_name = case_name
        self.max_nfev = max_nfev
        self.n_workers = n_workers
        self.pestpp_options = (pestpp_options or {}).copy()
        self._runtime_options: Dict[str, Any] = {}
        for runtime_key in ("panther_timeout_secs",):
            value = _pop_option(self.pestpp_options, runtime_key)
            if value is not None:
                self._runtime_options[runtime_key] = value
        self.pestpp_version = pestpp_version
        self._agents: List[subprocess.Popen] = []

    def prepare(self) -> PestControlFile:
        """Create workspace, serialize problem, generate scripts, and construct PestControlFile."""
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # 1. Serialize Problem
        with open(self.work_dir / "problem.pkl", "wb") as f:
            pickle.dump(self.problem, f)

        # 2. Generate Worker
        generate_worker_script(self.work_dir / "worker.py")

        # 3. Write observation names
        params = self.problem.get_parameters()
        obs = self.problem.get_observations()
        with open(self.work_dir / "obs_names.txt", "w") as f:
            for obs_def in obs:
                f.write(f"{obs_def.name}\n")

        # 4. Initialize PestControlFile
        pst = PestControlFile()
        pst.tpl_file = f"{self.case_name}.tpl"
        pst.ins_file = f"{self.case_name}.ins"
        pst.model_input_file = "params.csv"
        pst.model_output_file = "results.dat"
        pst.pestpp_options = self.pestpp_options.copy()
        pst.pestpp_options.setdefault("max_run_fail", 1)
        noptmax_override = _normalize_legacy_pestpp_options(
            pst.pestpp_options,
            self.engine,
            params,
            obs,
        )

        if self.engine == "pestpp-glm":
            pst.pestpp_options.setdefault("lambdas", "1.0")

        # Config IES files
        if self.engine == "pestpp-ies":
            n_reals = int(pst.pestpp_options.get("ies_num_reals", 100))
            pst.pestpp_options.setdefault("ies_num_reals", n_reals)
            
            # Prior ensemble
            if "ies_parameter_ensemble" not in pst.pestpp_options:
                par_ens_file = self.work_dir / f"{self.case_name}.0.par.csv"
                write_parameter_ensemble(par_ens_file, params, n_reals)
                pst.pestpp_options["ies_parameter_ensemble"] = par_ens_file.name
                
            # Observation ensemble
            if "ies_observation_ensemble" not in pst.pestpp_options:
                obs_ens_file = self.work_dir / f"{self.case_name}.obs.csv"
                write_observation_ensemble(obs_ens_file, obs, n_reals)
                pst.pestpp_options["ies_observation_ensemble"] = obs_ens_file.name

            # Parameter covariance (parcov)
            parcov_opt = pst.pestpp_options.get("parcov")
            if parcov_opt is True:  # auto-generate from bounds
                parcov_file = self.work_dir / f"{self.case_name}.parcov.csv"
                write_parcov_matrix(parcov_file, params)
                pst.pestpp_options["parcov"] = parcov_file.name
            elif isinstance(parcov_opt, (str, Path)):
                src = Path(parcov_opt)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["parcov"] = src.name

            # Observation covariance (obscov)
            obscov_opt = pst.pestpp_options.get("obscov")
            if obscov_opt is True:
                obscov_file = self.work_dir / f"{self.case_name}.obscov.csv"
                write_obscov_matrix(obscov_file, obs)
                pst.pestpp_options["obscov"] = obscov_file.name
            elif isinstance(obscov_opt, (str, Path)):
                src = Path(obscov_opt)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["obscov"] = src.name

            # Localizer
            loc_opt = pst.pestpp_options.get("ies_localizer")
            if loc_opt is True:
                loc_file = self.work_dir / f"{self.case_name}.localizer.csv"
                write_localizer_matrix(loc_file, params, obs)
                pst.pestpp_options["ies_localizer"] = loc_file.name
            elif isinstance(loc_opt, (str, Path)):
                src = Path(loc_opt)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["ies_localizer"] = src.name

            # Restart ensembles
            restart_par = pst.pestpp_options.get("ies_restart_parameter_ensemble")
            restart_obs = pst.pestpp_options.get("ies_restart_observation_ensemble")
            if restart_par and isinstance(restart_par, (str, Path)):
                src = Path(restart_par)
                if not src.is_absolute():
                    src = Path.cwd() / src
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["ies_restart_parameter_ensemble"] = src.name
                else:
                    logger.warning(f"Restart parameter ensemble not found: {restart_par}")
            if restart_obs and isinstance(restart_obs, (str, Path)):
                src = Path(restart_obs)
                if not src.is_absolute():
                    src = Path.cwd() / src
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["ies_restart_observation_ensemble"] = src.name
                else:
                    logger.warning(f"Restart observation ensemble not found: {restart_obs}")

        # Config Sweep files
        if self.engine == "pestpp-swp":
            sweep_file = self.work_dir / f"{self.case_name}.swp.in.csv"
            n_runs = int(_pop_option(pst.pestpp_options, "sweep_n_runs") or 50)
            method = _pop_option(pst.pestpp_options, "sweep_sampler") or "latin_hypercube"
            generate_sweep_csv(sweep_file, params, n_runs, method)
            pst.pestpp_options.setdefault("sweep_parameter_csv_file", sweep_file.name)
            # Additional SWP options
            if "sweep_output_csv_file" not in pst.pestpp_options:
                pst.pestpp_options["sweep_output_csv_file"] = f"{self.case_name}.swp.out.csv"
            if "sweep_forgive" in self.pestpp_options:
                # passthrough - already copied above, no action needed
                pass

        # Config MOU/OPT files
        if self.engine in ("pestpp-mou", "pestpp-opt"):
            if self.engine == "pestpp-mou":
                pst.pestpp_options.setdefault("mou_population_size", 50)
                pst.pestpp_options.setdefault("mou_generator", "de")

                # Auto-detect objectives from observation names/groups
                obj_names = pst.pestpp_options.get("mou_objectives")
                if not obj_names:
                    # Objectives must use PEST++ direction-aware observation groups.
                    obj_obs = [
                        o.name for o in obs
                        if o.name.startswith(("obj_", "mou_obj_"))
                        or o.group.lower().startswith(("l_", "less_", "less_than", "g_", "greater_", "greater_than"))
                    ]
                    if obj_obs:
                        pst.pestpp_options.setdefault(
                            "mou_objectives",
                            ",".join(obj_obs)
                        )

            # Constraint groups for OPT
            if self.engine == "pestpp-opt":
                con_groups = set(
                    o.group for o in obs
                    if o.group.lower().startswith(("l_", "less_", "less_than", "g_", "greater_", "greater_than"))
                )
                if con_groups:
                    pst.pestpp_options.setdefault(
                        "opt_constraint_groups",
                        ",".join(sorted(con_groups))
                    )

            # Warm-start Jacobian
            base_jac = pst.pestpp_options.get("base_jacobian")
            if not base_jac:
                hot_jac = self.work_dir / f"{self.case_name}.jcb"
                if hot_jac.exists():
                    pst.pestpp_options.setdefault("base_jacobian", hot_jac.name)

            # Risk level for chance constraints
            pst.pestpp_options.setdefault("opt_risk", 0.95)

        # Config DA files
        if self.engine == "pestpp-da":
            da_restart = pst.pestpp_options.get("da_hotstart_cycle", 0)
            if isinstance(da_restart, int) and da_restart > 0:
                pass  # passthrough

            # Cycle tables for time-variant observations/weights
            da_cycle_table = pst.pestpp_options.get("da_cycle_table")
            if da_cycle_table and isinstance(da_cycle_table, (str, Path)):
                src = Path(da_cycle_table)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["da_cycle_table"] = src.name

            da_obs_table = pst.pestpp_options.get("da_observation_cycle_table")
            if da_obs_table and isinstance(da_obs_table, (str, Path)):
                src = Path(da_obs_table)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["da_observation_cycle_table"] = src.name

            # State-parameter ensemble handling
            state_names_raw = _pop_option(pst.pestpp_options, "da_state_names")
            state_names = _split_option_tokens(state_names_raw)
            if state_names:
                # If states are provided, generate initial state-parameter ensemble
                da_state_ens = self.work_dir / f"{self.case_name}.state.0.par.csv"
                if not da_state_ens.exists():
                    # Combine parameter values with state variable initial values
                    n_reals = int(pst.pestpp_options.get("ies_num_reals", 50))
                    state_vals_raw = _pop_option(pst.pestpp_options, "da_state_initial_values")
                    state_vals = []
                    for val in _split_option_tokens(state_vals_raw):
                        try:
                            state_vals.append(float(val))
                        except ValueError:
                            logger.warning(f"Ignoring non-numeric DA state initial value: {val}")
                    _write_da_state_ensemble(da_state_ens, params, state_names, state_vals, n_reals)
                if da_state_ens.exists():
                    pst.pestpp_options.setdefault(
                        "da_state_parameter_ensemble", da_state_ens.name
                    )

            # Observation ensemble per cycle
            da_ens = pst.pestpp_options.get("da_parameter_ensemble")
            if da_ens and isinstance(da_ens, (str, Path)):
                src = Path(da_ens)
                if src.exists():
                    dst = self.work_dir / src.name
                    if not dst.exists():
                        shutil.copy2(str(src), str(dst))
                    pst.pestpp_options["da_parameter_ensemble"] = src.name

        pst.control_data.noptmax = noptmax_override if noptmax_override is not None else self.max_nfev

        for p in params:
            pst.add_parameter(PestParameter.from_adjustable_parameter(p))
        for o in obs:
            pst.add_observation(PestObservation.from_observation(o))

        pst.write_tpl(self.work_dir / pst.tpl_file)
        pst.write_ins(self.work_dir / pst.ins_file)
        pst.write(self.work_dir / f"{self.case_name}.pst")

        return pst

    @staticmethod
    def _find_free_port() -> int:
        """Find an ephemeral TCP port available on localhost."""
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("", 0))
            return int(s.getsockname()[1])

    def _terminate_agents(self) -> None:
        """Clean up any spawned agent (worker) processes."""
        for agent in list(self._agents):
            self._terminate_process(agent)
        self._agents = []

    @staticmethod
    def _terminate_process(proc: Optional[subprocess.Popen]) -> None:
        """Terminate a spawned PEST++ process and close captured pipes."""
        if proc is None:
            return
        try:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait(timeout=5)
        except Exception as e:
            logger.warning(f"Failed to terminate process {getattr(proc, 'pid', '?')}: {e}")
        for stream in (getattr(proc, "stdout", None), getattr(proc, "stderr", None)):
            try:
                if stream:
                    stream.close()
            except Exception:
                pass

    def _has_expected_outputs(self) -> bool:
        """Return True when the selected engine wrote parseable output files."""
        case = self.case_name
        if self.engine == "pestpp-glm":
            return any((self.work_dir / f"{case}{suffix}").exists() for suffix in (".par", ".rei", ".res"))
        if self.engine == "pestpp-ies":
            return (
                (self.work_dir / f"{case}.phi.actual.csv").exists()
                or bool(list(self.work_dir.glob(f"{case}.*.par.csv")))
            )
        if self.engine == "pestpp-sen":
            return any(
                (self.work_dir / f"{case}{suffix}").exists()
                for suffix in (".msn", ".sbl", ".sobol.si.csv", ".raw.csv")
            )
        if self.engine == "pestpp-swp":
            return (self.work_dir / f"{case}.swp.out.csv").exists()
        if self.engine == "pestpp-mou":
            return bool(list(self.work_dir.glob(f"{case}*.dv_pop.csv"))) or any(
                (self.work_dir / f"{case}{suffix}").exists()
                for suffix in (".pareto.summary.csv", ".pareto.archive.summary.csv", ".pareto.csv")
            )
        if self.engine == "pestpp-opt":
            return any((self.work_dir / f"{case}{suffix}").exists() for suffix in (".par", ".rei", ".res"))
        if self.engine == "pestpp-da":
            return (
                (self.work_dir / f"{case}.global.phi.actual.csv").exists()
                or bool(list(self.work_dir.glob(f"{case}.*.*.par.csv")))
            )
        return False

    def _accept_return_code(self, returncode: int, cmd: List[str]) -> None:
        if returncode == 0:
            return
        if self._has_expected_outputs():
            logger.warning(
                f"{self.engine} exited with code {returncode} after writing expected outputs; "
                "continuing with output parsing."
            )
            return
        raise subprocess.CalledProcessError(returncode, cmd)

    def run(self) -> Dict[str, Any]:
        """Execute PEST++ binary and parse the outputs."""
        pst = self.prepare()
        pst_file = f"{self.case_name}.pst"

        # Find Executable
        exe_path = get_executable_path(self.engine, self.pestpp_version)

        try:
            if self.n_workers <= 1:
                # Single-process (serial) execution
                logger.info(f"Starting {self.engine} (serial)...")
                cmd = [exe_path, pst_file]
                proc = subprocess.run(cmd, cwd=str(self.work_dir), check=False)
                self._accept_return_code(getattr(proc, "returncode", 0), cmd)
            else:
                # Parallel manager / agent (PANTHER) execution
                port = self._find_free_port()
                logger.info(
                    f"Starting {self.engine} manager on port {port} "
                    f"with {self.n_workers} agents ..."
                )

                # Default parallel options if not already set
                self.pestpp_options.setdefault("max_run_fail", 3)
                self.pestpp_options.setdefault("panther_agent_restart_on_error", True)
                self.pestpp_options.setdefault("panther_agent_no_ping_timeout_secs", 300)
                self.pestpp_options.setdefault("overdue_resched_fac", 1.15)
                self.pestpp_options.setdefault("overdue_giveup_fac", 2.0)
                # Re-write .pst with updated options
                pst.pestpp_options = self.pestpp_options.copy()
                pst.write(self.work_dir / f"{self.case_name}.pst")

                # Launch manager (non-blocking)
                manager_cmd = [exe_path, pst_file, "/h", f":{port}"]
                manager = subprocess.Popen(
                    manager_cmd, cwd=str(self.work_dir)
                )

                # Allow manager time to bind and start listening
                time.sleep(2)

                # Launch agents
                for i in range(self.n_workers):
                    agent_cmd = [exe_path, pst_file, "/h", f"localhost:{port}"]
                    agent = subprocess.Popen(
                        agent_cmd, cwd=str(self.work_dir)
                    )
                    self._agents.append(agent)
                    logger.debug(
                        f"Spawned agent {i + 1}/{self.n_workers} (pid={agent.pid})"
                    )

                # Wait for manager to finish
                timeout_opt = self._runtime_options.get("panther_timeout_secs")
                manager_timeout = None
                if timeout_opt is not None:
                    try:
                        manager_timeout = float(timeout_opt)
                    except (TypeError, ValueError):
                        logger.warning(f"Ignoring invalid panther_timeout_secs value: {timeout_opt}")
                try:
                    if manager_timeout is None:
                        manager.wait()
                    else:
                        manager.wait(timeout=manager_timeout)
                except subprocess.TimeoutExpired:
                    if self._has_expected_outputs():
                        logger.warning(
                            f"{self.engine} PANTHER manager timed out after "
                            f"{manager_timeout} seconds after writing expected outputs; "
                            "terminating manager and parsing outputs."
                        )
                        self._terminate_process(manager)
                    else:
                        self._terminate_process(manager)
                        raise subprocess.CalledProcessError(-1, manager_cmd)
                else:
                    self._accept_return_code(manager.returncode, manager_cmd)
        except subprocess.CalledProcessError as e:
            logger.error(f"PEST++ failed: {e}")
            return {"success": False, "phi": -1}
        except FileNotFoundError:
            logger.error(f"Executable {exe_path} not found.")
            return {"success": False, "phi": -1}
        finally:
            self._terminate_agents()

        # Parse Results
        result = {"success": True, "phi": 0.0}

        # GLM parser
        if self.engine == "pestpp-glm":
            optimal_params = {}
            par_file = self.work_dir / f"{self.case_name}.par"
            if par_file.exists():
                try:
                    with open(par_file, "r") as f:
                        next(f)
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                optimal_params[parts[0]] = float(parts[1])
                except Exception as e:
                    logger.warning(f"Failed to read par file: {e}")
            
            res_file = self.work_dir / f"{self.case_name}.rei"
            if not res_file.exists():
                res_file = self.work_dir / f"{self.case_name}.res"
            residuals, simulated, phi_by_group = parse_pest_residuals(res_file)
            
            sen_file = self.work_dir / f"{self.case_name}.sen"
            sensitivities = parse_pest_sensitivities(sen_file)
            
            id_file = self.work_dir / f"{self.case_name}.id"
            ident = parse_pest_identifiability(id_file)
            
            cov_path = self.work_dir / f"{self.case_name}.cov"
            cov_str = str(cov_path.absolute()) if cov_path.exists() else None
            
            cor_path = self.work_dir / f"{self.case_name}.cor"
            cor_str = str(cor_path.absolute()) if cor_path.exists() else None

            phi_total = sum(phi_by_group.values()) if phi_by_group else sum(r**2 for r in residuals.values())

            result.update({
                "phi": phi_total,
                "optimal_parameters": optimal_params,
                "residuals": residuals,
                "simulated_observations": simulated,
                "phi_by_group": phi_by_group,
                "covariance_matrix_path": cov_str,
                "correlation_matrix_path": cor_str,
                "sensitivity_table": sensitivities,
                "identifiability": ident
            })

        # IES parser
        elif self.engine == "pestpp-ies":
            prior_par = parse_ies_csv(self.work_dir / f"{self.case_name}.0.par.csv")
            
            post_par = None
            for i in sorted(range(self.max_nfev + 1), reverse=True):
                path = self.work_dir / f"{self.case_name}.{i}.par.csv"
                if path.exists():
                    post_par = parse_ies_csv(path)
                    break
            
            prior_obs = parse_ies_csv(self.work_dir / f"{self.case_name}.0.obs.csv")
            
            post_obs = None
            for i in sorted(range(self.max_nfev + 1), reverse=True):
                path = self.work_dir / f"{self.case_name}.{i}.obs.csv"
                if path.exists():
                    post_obs = parse_ies_csv(path)
                    break

            phi_history = []
            phi_file = self.work_dir / f"{self.case_name}.phi.actual.csv"
            if phi_file.exists():
                phi_history = _parse_phi_history_csv(phi_file)

            optimal_params = {}
            if post_par:
                for k, vals in post_par.items():
                    optimal_params[k] = float(np.mean(vals))

            # Posterior forecast summaries (if ++forecasts option set)
            forecasts = []
            for key in ("forecasts", "++forecasts"):
                fc = self.pestpp_options.get(key, "")
                if fc:
                    forecasts = [f.strip() for f in str(fc).split(",") if f.strip()]
                    break
            posterior_forecast_summaries = {}
            if post_obs and forecasts:
                for fc_name in forecasts:
                    if fc_name in post_obs:
                        vals = post_obs[fc_name]
                        posterior_forecast_summaries[fc_name] = {
                            "mean": float(np.mean(vals)),
                            "std": float(np.std(vals)),
                            "min": float(np.min(vals)),
                            "max": float(np.max(vals)),
                            "median": float(np.median(vals)),
                        }

            # Track generated covariance / localizer file paths
            parcov_path = self.pestpp_options.get("parcov")
            obscov_path = self.pestpp_options.get("obscov")
            localizer_path = self.pestpp_options.get("ies_localizer")
            restart_from = self.pestpp_options.get("ies_restart_parameter_ensemble")

            result.update({
                "phi": phi_history[-1] if phi_history else 0.0,
                "optimal_parameters": optimal_params,
                "prior_parameters": prior_par,
                "posterior_parameters": post_par,
                "prior_observations": prior_obs,
                "simulated_observations": post_obs,
                "phi_history": phi_history,
                "posterior_forecast_summaries": posterior_forecast_summaries if posterior_forecast_summaries else None,
                "parcov_path": parcov_path,
                "obscov_path": obscov_path,
                "localizer_path": localizer_path,
                "restart_from": restart_from,
            })

        # SEN parser
        elif self.engine == "pestpp-sen":
            sobol_file = self.work_dir / f"{self.case_name}.sobol.csv"
            first_order_sobol = {}
            total_sobol = {}
            if sobol_file.exists():
                try:
                    df = pd.read_csv(sobol_file)
                    for _, row in df.iterrows():
                        p_name = str(row.get("parameter", row.get("name", "")))
                        if p_name:
                            first_order_sobol[p_name] = float(row.get("first_order", row.get("si", 0.0)))
                            total_sobol[p_name] = float(row.get("total_order", row.get("sti", 0.0)))
                except Exception as e:
                    logger.warning(f"Failed to parse Sobol csv: {e}")

            sobol_si_file = self.work_dir / f"{self.case_name}.sobol.si.csv"
            sobol_sti_file = self.work_dir / f"{self.case_name}.sobol.sti.csv"
            if sobol_si_file.exists():
                try:
                    df = pd.read_csv(sobol_si_file)
                    par_col = next((c for c in df.columns if c.lower() in ("parameter", "parnme", "name", "parameter_name")), df.columns[0])
                    value_cols = [c for c in df.columns if c != par_col]
                    for _, row in df.iterrows():
                        p_name = str(row[par_col]).strip()
                        vals = pd.to_numeric(row[value_cols], errors="coerce").dropna()
                        if p_name and not vals.empty:
                            first_order_sobol[p_name] = float(vals.mean())
                except Exception as e:
                    logger.warning(f"Failed to parse Sobol first-order csv: {e}")
            if sobol_sti_file.exists():
                try:
                    df = pd.read_csv(sobol_sti_file)
                    par_col = next((c for c in df.columns if c.lower() in ("parameter", "parnme", "name", "parameter_name")), df.columns[0])
                    value_cols = [c for c in df.columns if c != par_col]
                    for _, row in df.iterrows():
                        p_name = str(row[par_col]).strip()
                        vals = pd.to_numeric(row[value_cols], errors="coerce").dropna()
                        if p_name and not vals.empty:
                            total_sobol[p_name] = float(vals.mean())
                except Exception as e:
                    logger.warning(f"Failed to parse Sobol total-order csv: {e}")
                    
            morris_file = self.work_dir / f"{self.case_name}.morris.csv"
            morris_effects = {}
            if morris_file.exists():
                try:
                    df = pd.read_csv(morris_file)
                    for _, row in df.iterrows():
                        p_name = str(row.get("parameter", row.get("name", "")))
                        if p_name:
                            morris_effects[p_name] = [
                                float(row.get("mean", row.get("mu_star", 0.0))),
                                float(row.get("std", row.get("sigma", 0.0)))
                            ]
                except Exception as e:
                    logger.warning(f"Failed to parse Morris csv: {e}")

            morris_summary = parse_sen_msn(self.work_dir / f"{self.case_name}.msn")
            if morris_summary:
                for p_name, metrics in morris_summary.items():
                    morris_effects[p_name] = [
                        float(metrics.get("sen_mean_abs", metrics.get("sen_mean", 0.0))),
                        float(metrics.get("sen_std_dev", 0.0)),
                    ]

            ranked_importance = []
            if total_sobol:
                ranked_importance = sorted(total_sobol.keys(), key=lambda k: total_sobol[k], reverse=True)
            elif morris_effects:
                ranked_importance = sorted(morris_effects.keys(), key=lambda k: abs(morris_effects[k][0]), reverse=True)
            
            recommended_subset = []
            if total_sobol:
                recommended_subset = [k for k, v in total_sobol.items() if v > 0.1]
            elif morris_effects:
                recommended_subset = [k for k, v in morris_effects.items() if abs(v[0]) > 0.5]

            result.update({
                "first_order_sobol": first_order_sobol,
                "total_sobol": total_sobol,
                "morris_elementary_effects": morris_effects,
                "sensitivity_table": morris_summary,
                "ranked_importance": ranked_importance,
                "recommended_calibratable_subset": recommended_subset
            })

        # SWP parser
        elif self.engine == "pestpp-swp":
            swp_file = self.work_dir / f"{self.case_name}.swp.out.csv"
            sweep_table = None
            if swp_file.exists():
                try:
                    df = pd.read_csv(swp_file)
                    sweep_table = df.to_dict(orient="list")
                except Exception as e:
                    logger.warning(f"Failed to parse sweep csv: {e}")
            result.update({
                "sweep_table": sweep_table
            })

        # OPT / MOU parser
        elif self.engine in ("pestpp-mou", "pestpp-opt"):
            pareto_candidates = [
                self.work_dir / f"{self.case_name}.pareto.summary.csv",
                self.work_dir / f"{self.case_name}.pareto.archive.summary.csv",
                self.work_dir / f"{self.case_name}.pareto.csv",
                self.work_dir / f"{self.case_name}.pareto.out.csv",
            ]
            pareto_file = next((p for p in pareto_candidates if p.exists()), pareto_candidates[0])
            pareto_front = None
            if pareto_file.exists():
                try:
                    df = pd.read_csv(pareto_file)
                    pareto_front = df.to_dict(orient="records")
                except Exception as e:
                    logger.warning(f"Failed to parse Pareto CSV: {e}")

            dv_population = None
            obs_population = None
            if self.engine == "pestpp-mou":
                dv_candidates = (
                    sorted(self.work_dir.glob(f"{self.case_name}.*.dv_pop.csv"))
                    + [self.work_dir / f"{self.case_name}.dv_pop.csv"]
                    + sorted(self.work_dir.glob(f"{self.case_name}*.archive*.dv_pop.csv"))
                )
                obs_candidates = (
                    sorted(self.work_dir.glob(f"{self.case_name}.*.obs_pop.csv"))
                    + [self.work_dir / f"{self.case_name}.obs_pop.csv"]
                    + sorted(self.work_dir.glob(f"{self.case_name}*.archive*.obs_pop.csv"))
                )
                for candidate in reversed([p for p in dv_candidates if p.exists()]):
                    dv_population = _parse_csv_records(candidate)
                    if dv_population is not None:
                        break
                for candidate in reversed([p for p in obs_candidates if p.exists()]):
                    obs_population = _parse_csv_records(candidate)
                    if obs_population is not None:
                        break
                    
            opt_par_file = self.work_dir / f"{self.case_name}.par"
            optimal_dec_vars = {}
            if opt_par_file.exists():
                try:
                    with open(opt_par_file, "r") as f:
                        next(f)
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                optimal_dec_vars[parts[0]] = float(parts[1])
                except Exception as e:
                    logger.warning(f"Failed to read decision variables: {e}")

            result.update({
                "pareto_front": pareto_front,
                "decision_variable_population": dv_population,
                "objective_population": obs_population,
                "optimal_decision_variables": optimal_dec_vars,
                "optimal_parameters": optimal_dec_vars
            })

        # DA parser
        elif self.engine == "pestpp-da":
            cycles = []
            prior_par = None
            posterior_par = None
            state_summary = {}

            cycle_files = []
            for p_file in sorted(self.work_dir.glob(f"{self.case_name}.da.cycle_*.par.csv")):
                try:
                    cycle_idx = int(p_file.name.split("cycle_")[1].split(".")[0])
                    cycle_files.append((cycle_idx, 0, p_file, self.work_dir / f"{self.case_name}.da.cycle_{cycle_idx}.obs.csv"))
                except Exception:
                    continue
            pattern = re.compile(rf"^{re.escape(self.case_name)}\.(\d+)\.(\d+)\.par\.csv$", re.IGNORECASE)
            for p_file in sorted(self.work_dir.glob(f"{self.case_name}.*.*.par.csv")):
                match = pattern.match(p_file.name)
                if not match:
                    continue
                cycle_idx = int(match.group(1))
                iter_idx = int(match.group(2))
                cycle_files.append((cycle_idx, iter_idx, p_file, self.work_dir / f"{self.case_name}.{cycle_idx}.{iter_idx}.obs.csv"))

            for cycle_idx, iter_idx, p_file, o_file in cycle_files:
                try:
                    p_data = parse_ies_csv(p_file)
                    o_data = parse_ies_csv(o_file) if o_file.exists() else None

                    # Track state-parameter separation
                    state_names = self.pestpp_options.get("da_state_names", [])
                    state_ens = {}
                    param_ens = {}
                    if p_data:
                        for key, vals in p_data.items():
                            if any(key.startswith(s) for s in state_names):
                                state_ens[key] = vals
                            else:
                                param_ens[key] = vals

                    if cycle_idx == 0:
                        prior_par = param_ens if param_ens else p_data
                    # Last cycle posterior
                    posterior_par = param_ens if param_ens else p_data

                    # Summarize state uncertainties per cycle
                    if state_ens:
                        state_summary[cycle_idx] = {
                            k: {"mean": float(np.mean(v)), "std": float(np.std(v))}
                            for k, v in state_ens.items()
                        }

                    cycles.append({
                        "cycle": cycle_idx,
                        "iteration": iter_idx,
                        "parameters": param_ens,
                        "observations": o_data,
                        "state_ensemble": state_ens,
                    })
                except Exception:
                    continue

            global_prior = parse_ies_csv(self.work_dir / f"{self.case_name}.global.prior.pe.csv")
            if global_prior:
                prior_par = global_prior
            global_pe_files = sorted(self.work_dir.glob(f"{self.case_name}.global.*.pe.csv"))
            if global_pe_files:
                latest_global = parse_ies_csv(global_pe_files[-1])
                if latest_global:
                    posterior_par = latest_global

            phi_history = _parse_phi_history_csv(self.work_dir / f"{self.case_name}.global.phi.actual.csv")
            if not phi_history:
                phi_history = _parse_phi_history_csv(self.work_dir / f"{self.case_name}.phi.actual.csv")

            result.update({
                "phi": phi_history[-1] if phi_history else 0.0,
                "phi_history": phi_history,
                "cycles": cycles,
                "prior_parameters": prior_par,
                "posterior_parameters": posterior_par,
                "state_summaries": state_summary if state_summary else None,
            })

        return result


def run_pestpp(
    problem,
    engine: str = "pestpp-ies",
    work_dir: str = "pest_workspace",
    case_name: str = "case1",
    max_nfev: int = 100,
    n_workers: int = 1,
    pestpp_options: Optional[Dict[str, Any]] = None,
    pestpp_version: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run PEST++ calibration.
    """
    runner = PestRunner(
        problem=problem,
        engine=engine,
        work_dir=work_dir,
        case_name=case_name,
        max_nfev=max_nfev,
        n_workers=n_workers,
        pestpp_options=pestpp_options,
        pestpp_version=pestpp_version,
    )
    return runner.run()
