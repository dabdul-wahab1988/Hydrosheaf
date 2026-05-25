"""Run the real M4 topology-validation benchmark.

This entrypoint intentionally delegates to the Phase 2B public-archive
validation pipeline. The previous controlled toy benchmark has been removed
from the main M4 path so `independent_graph_vs_modpath.csv` is always generated
from the Savage MODPATH reference graph and Hydrosheaf's independent graph
inference scenarios.
"""
from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from phase2b_independent_validation import main


if __name__ == "__main__":
    main()
