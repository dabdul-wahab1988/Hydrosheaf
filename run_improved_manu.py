import sys
from pathlib import Path

# Add project root to path to ensure hydrosheaf is importable
sys.path.insert(0, str(Path.cwd()))

try:
    from hydrosheaf.workflows.manu import run_manu
    run_manu("manu.xlsx")
except ImportError as e:
    print(f"Import Error: {e}")
    # Fallback if package structure issue
    sys.path.append("hydrosheaf")
    from workflows.manu import run_manu
    run_manu("manu.xlsx")
except Exception as e:
    print(f"Execution Error: {e}")
