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
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np

from ..definitions import AdjustableParameter, Observation

logger = logging.getLogger("calibration.pestpp")

PESTPP_VERSION = "5.2.25"
PESTPP_BASE_URL = "https://github.com/usgs/pestpp/releases/download"

def get_bin_dir() -> Path:
    """Return the local bin directory."""
    # Assuming running from project root
    return Path("bin")

def get_executable_suffix() -> str:
    return ".exe" if platform.system() == "Windows" else ""

def get_executable_path(name: str) -> str:
    """
    Find the path to a PEST++ executable.
    Downloads it if not found in local bin/.
    """
    bin_dir = get_bin_dir()
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
    logger.info(f"Executable {name} not found. Attempting to download PEST++ binaries...")
    try:
        download_pestpp_binaries(bin_dir)
    except Exception as e:
        logger.error(f"Failed to download PEST++ binaries: {e}")
        # Continue to check if it exists anyway (e.g. partial success)

    if local_path.exists():
        return str(local_path.absolute())
    
    return name

def download_pestpp_binaries(target_dir: Path):
    """Download and extract PEST++ binaries."""
    target_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine URL based on platform
    system = platform.system()
    if system == "Windows":
        filename = f"pestpp-{PESTPP_VERSION}-win.zip"
        zip_path_internal = "bin"
    elif system == "Linux":
        filename = f"pestpp-{PESTPP_VERSION}-linux.zip"
        zip_path_internal = "bin"
    elif system == "Darwin":
        filename = f"pestpp-{PESTPP_VERSION}-mac.zip"
        zip_path_internal = "bin"
    else:
        raise OSError(f"Unsupported platform for auto-download: {system}")

    url = f"{PESTPP_BASE_URL}/{PESTPP_VERSION}/{filename}"
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

    @classmethod
    def from_adjustable_parameter(cls, p: AdjustableParameter) -> "PestParameter":
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
        )

    def format(self) -> str:
        return f"{self.name} {self.trans} {self.change_limit} {self.value} {self.lower_bound} {self.upper_bound} {self.group} {self.scale} {self.offset}"


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
        self.tpl_file: str = "case.tpl"
        self.ins_file: str = "case.ins"
        self.model_input_file: str = "params.csv"
        self.model_output_file: str = "results.dat"
        self.pestpp_options: Dict[str, Any] = {}

    def add_parameter(self, param: PestParameter):
        self.parameters.append(param)
        self.control_data.npar = len(self.parameters)

    def add_observation(self, obs: PestObservation):
        self.observations.append(obs)
        self.control_data.nobs = len(self.observations)

    def write(self, filepath: Path):
        """Write the complete control file to disk."""
        # Sync dimensions
        self.control_data.npar = len(self.parameters)
        self.control_data.nobs = len(self.observations)

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

            f.write("* model command line\n")
            f.write(f"{self.model_command_line}\n")

            f.write("* model input/output\n")
            f.write(f"{self.tpl_file} {self.model_input_file}\n")
            f.write(f"{self.ins_file} {self.model_output_file}\n")

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
    ):
        self.problem = problem
        self.engine = engine
        self.work_dir = Path(work_dir)
        self.case_name = case_name
        self.max_nfev = max_nfev
        self.n_workers = n_workers
        self.pestpp_options = pestpp_options or {}

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

        if self.engine == "pestpp-glm":
            pst.pestpp_options.setdefault("lambdas", "1.0")

        pst.control_data.noptmax = self.max_nfev

        for p in params:
            pst.add_parameter(PestParameter.from_adjustable_parameter(p))
        for o in obs:
            pst.add_observation(PestObservation.from_observation(o))

        pst.write_tpl(self.work_dir / pst.tpl_file)
        pst.write_ins(self.work_dir / pst.ins_file)
        pst.write(self.work_dir / f"{self.case_name}.pst")

        return pst

    def run(self) -> Dict[str, Any]:
        """Execute PEST++ binary and parse the outputs."""
        pst = self.prepare()
        pst_file = f"{self.case_name}.pst"

        # Find Executable
        exe_path = get_executable_path(self.engine)

        logger.info(f"Starting {self.engine}...")
        cmd = [exe_path, pst_file]

        try:
            subprocess.run(cmd, cwd=str(self.work_dir), check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"PEST++ failed: {e}")
            return {"success": False, "phi": -1}
        except FileNotFoundError:
            logger.error(f"Executable {exe_path} not found.")
            return {"success": False, "phi": -1}

        # Parse Results
        result = {"success": True, "phi": 0.0}

        # Try reading phi.csv (IES)
        phi_file = self.work_dir / f"{self.case_name}.phi.actual.csv"
        if phi_file.exists():
            df = pd.read_csv(phi_file)
            result["phi"] = df.iloc[-1].mean()

        # Try reading .par (GLM best parameters)
        par_file = self.work_dir / f"{self.case_name}.par"
        if par_file.exists():
            try:
                optimal_params = {}
                with open(par_file, "r") as f:
                    next(f)  # skip header
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            optimal_params[parts[0]] = float(parts[1])
                result["optimal_parameters"] = optimal_params
            except Exception as e:
                logger.warning(f"Failed to read par file: {e}")

        return result


def run_pestpp(
    problem,
    engine: str = "pestpp-ies",
    work_dir: str = "pest_workspace",
    case_name: str = "case1",
    max_nfev: int = 100,
    n_workers: int = 1,
    pestpp_options: Optional[Dict[str, Any]] = None,
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
    )
    return runner.run()
