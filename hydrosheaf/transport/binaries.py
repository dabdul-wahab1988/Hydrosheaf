"""
MODFLOW/MT3DMS Binary Manager.
Handles auto-download of FloPy-compatible executables.
"""

import os
import shutil
import platform
import zipfile
import urllib.request
import stat
from pathlib import Path

from ..log import get_logger

logger = get_logger("transport.binaries")

# Release from https://github.com/MODFLOW-ORG/executables
MF_release_TAG = "23.0"
MF_BASE_URL = f"https://github.com/MODFLOW-ORG/executables/releases/download/{MF_release_TAG}"

def get_bin_dir() -> Path:
    # Use the shared bin directory in the project root
    # Assuming this file is in hydrosheaf/transport/
    root = Path(__file__).resolve().parent.parent.parent
    return root / "bin"

def get_platform_asset_name() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    if system == "windows":
        return "win64.zip"
    elif system == "linux":
        return "linux.zip"
    elif system == "darwin":
        return "mac.zip"
    else:
        raise OSError(f"Unsupported platform: {system}")

def get_executable_name(base_name: str) -> str:
    system = platform.system().lower()
    if system == "windows":
        return f"{base_name}.exe"
    return base_name

def download_modflow_binaries(target_dir: Path):
    """Download and extract MODFLOW executables."""
    target_dir.mkdir(parents=True, exist_ok=True)
    
    asset_name = get_platform_asset_name()
    url = f"{MF_BASE_URL}/{asset_name}"
    temp_zip = target_dir / "temp_mf.zip"
    
    logger.info(f"Downloading MODFLOW binaries from {url}...")
    try:
        urllib.request.urlretrieve(url, temp_zip)
    except Exception as e:
        logger.error(f"Download failed: {e}")
        return

    logger.info("Extracting MODFLOW binaries...")
    try:
        with zipfile.ZipFile(temp_zip, 'r') as z:
            # We specifically want mf2005 and mt3dms (and maybe mfnwt, mf6 later)
            # The zip structure is typically flat or has a folder.
            # Let's verify by checking names.
            for name in z.namelist():
                filename = Path(name).name
                if not filename: continue
                
                # Check for executables we care about
                # mf2005, mt3dms
                base = filename.lower()
                if base.startswith("mf2005") or base.startswith("mt3dms"):
                     if base.endswith(".exe") or "." not in base:
                        # Extract
                        source = z.open(name)
                        target_path = target_dir / filename
                        with open(target_path, "wb") as t:
                            shutil.copyfileobj(source, t)
                        
                        if platform.system() != "Windows":
                            st = os.stat(target_path)
                            os.chmod(target_path, st.st_mode | stat.S_IEXEC)
                            
    except Exception as e:
        logger.error(f"Extraction failed: {e}")
    finally:
        if temp_zip.exists():
            try:
                os.remove(temp_zip)
            except Exception as e:
                logger.warning(f"Could not remove temp zip: {e}")

def get_executable_path(name: str) -> str:
    """
    Get path to a binary, downloading if missing.
    name: e.g. 'mf2005', 'mt3dms'
    """
    bin_dir = get_bin_dir()
    exe_name = get_executable_name(name)
    local_path = bin_dir / exe_name
    
    if local_path.exists():
        return str(local_path.absolute())
        
    # Check PATH
    path_exe = shutil.which(exe_name) or shutil.which(name)
    if path_exe:
        return path_exe
        
    # Download if not found
    logger.info(f"Binary {exe_name} not found. Attempting download...")
    download_modflow_binaries(bin_dir)
    
    if local_path.exists():
        return str(local_path.absolute())
        
    # Fallback to name (will likely fail in flopy but provides clear error)
    return name
