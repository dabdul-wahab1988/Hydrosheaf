# Reproducibility Guide: Hydrosheaf PhD Benchmarks

This document provides instructions for regenerating the benchmarks, figures, and tables presented in the PhD thesis "Methodological Framework for Integrated Hydrogeochemical and Flow Network Analysis".

## 1. Environment Setup

To ensure reproducibility and avoid binary dependency issues (especially with BLAS/LAPACK on Linux/WSL), use the included `uv` manager.

### Using `uv` (Recommended)

```bash
# Clean existing virtual environment if having C-extension issues
rm -rf .venv
uv venv
source .venv/bin/activate

# Install the package and dependencies
uv pip install -e .[all] pytest
```

## 2. Core Science Verification

Before running benchmarks, verify the core package functionality. From the activated virtual environment, run:

```bash
# Run unit tests to ensure the physics and inference engines are correct
.venv/bin/pytest tests/ -v
```

## 3. Manuscript Benchmarks (M2-M5)

Each manuscript has a dedicated benchmark package. You can run them individually or use the master script.

### Master Reproduction Script

The master script runs the benchmarks in sequence, providing a single entry point to reproduce the thesis evidence.

```bash
python run_all_benchmarks.py
```

### Individual Benchmarks

#### M2: Framework & Multi-Tracer Validation
- **Goal:** Validate joint LPM fitting and tracer agreement.
- **Command:** `python M2/m2_benchmark/scripts/run_m2_benchmark.py`
- **Key Output:** `M2/m2_benchmark/docs/m2_results_summary.md`

#### M3: Nuclear Tracer Age Inference
- **Goal:** Test graph-regularized age inference against single-node estimates.
- **Command:** `python M3/m3_age_benchmark/scripts/run_m3_age_benchmark.py`
- **Key Output:** `M3/m3_age_benchmark/results/graph_regularization_scenarios.csv`

#### M4: Topology & MODPATH Validation
- **Goal:** Compare reduced-order graph topology with MODPATH reference advection.
- **Command:** `python M4/m4_topology_benchmark/scripts/run_m4_topology_benchmark.py`
- **Key Output:** `M4/m4_topology_benchmark/results/independent_graph_vs_modpath.csv`

#### M5: Inverse Hydrogeochemical Reaction
- **Goal:** Screen sparse linear inverse models with PHREEQC thermodynamics.
- **Command:** `python M5/m5_inverse_reaction_benchmark/scripts/run_m5_inverse_reaction_benchmark.py`
- **Key Output:** `M5/m5_inverse_reaction_benchmark/results/thermodynamic_bound_violations.csv`

## 4. Data Inventory

- `data/synthetic/`: Ground-truth datasets for benchmarks.
- `data/NorthenGhana/`: Field hydrochemistry demonstration data.
- `M2/m2_benchmark/external/`: External archives (USGS, MODPATH) for tier-2 validation.

## 5. Known Limitations & Guardrails

- Hydrosheaf is a **research methodological framework**, not a production-grade forward simulator.
- For M5, inverse results are "screened" by PHREEQC, not "solved" by a fully coupled nonlinear solver.
- Figure generation requires `matplotlib` and `seaborn`.

---
**Last Updated:** May 2026
**Contact:** Dickson Abdul-Wahab
