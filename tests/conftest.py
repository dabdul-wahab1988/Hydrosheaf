import os
import sys

# Fix for PyMC/MKL crash on Windows (Heap Corruption 0xc0000374)
# Must be set before numpy/pymc are imported
if sys.platform == "win32":
    os.environ.setdefault("MKL_THREADING_LAYER", "GNU")
