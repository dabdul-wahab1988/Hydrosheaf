# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Hydrosheaf is a Python framework for inverse geochemical modeling in groundwater hydrogeochemistry. It determines optimal transport processes (evaporation, mixing) and geochemical reactions that explain observed chemical evolution along flow paths in aquifer networks.

The project has two main components:
1. **Core Python Package** (`hydrosheaf/`) - The scientific modeling engine
2. **Web Application** (`web/`) - FastAPI backend + React frontend for interactive analysis

## Common Commands

### Core Package

```bash
# Install the core package
pip install .

# Install with plotting support
pip install .[plot]

# Install with PHREEQC integration
pip install .[phreeqc]

# Run CLI
hydrosheaf --help

# Run all tests
python -m pytest tests/

# Run a single test file
python -m pytest tests/test_edge_fit.py

# Run a specific test
python -m pytest tests/test_edge_fit.py::EdgeFitTests::test_edge_fit_evap_and_halite
```

### Web Backend (FastAPI)

```bash
cd web/backend

# Install dependencies
pip install -r requirements.txt
pip install ../../   # Install Hydrosheaf core (required!)

# Start development server
uvicorn app.main:app --reload --port 8000

# Run integration tests
python test_integration.py
```

### Web Frontend (React + Vite)

```bash
cd web/frontend

# Install dependencies
npm install

# Start development server
npm run dev

# Build for production
npm run build

# Lint
npm run lint
```

## Architecture

### Core Package Structure (`hydrosheaf/`)

The core modeling engine is organized by functionality:

- **`config.py`** - Central `Config` dataclass with all model parameters (~100 configurable options)
- **`api.py`** - High-level pipeline functions: `fit_network_pipeline()`, `fit_network_with_priors()`
- **`inference/`** - Core fitting algorithms
  - `edge_fit.py` - Single edge transport+reaction fitting (`fit_edge()`, `EdgeResult`)
  - `network_fit.py` - Network-level orchestration (`fit_network()`, `infer_edges()`)
- **`models/`** - Geochemical models
  - `transport.py` - Evaporation/mixing models
  - `reactions.py` - LASSO reaction fitting with mineral dictionary
  - `sheaf.py` - Sheaf-theoretic residual computation
- **`graph/`** - Flow network construction
  - `build.py` - `infer_edges_probabilistic()` for Bayesian flow inference
  - `types.py` - `Edge` dataclass
- **`graph3d/`** - 3D layered aquifer support
- **`phreeqc/`** - Thermodynamic constraint integration via PHREEQC
- **`temporal/`** - Time-series analysis and residence time estimation
- **`uncertainty/`** - Bootstrap and Bayesian uncertainty quantification
- **`reactive_transport/`** - Forward validation with kinetic rate laws
- **`vadose/`** - Vadose zone physics priors

### Key Data Flow

1. **Input**: Water samples with ion concentrations (Ca, Mg, Na, HCO3, Cl, SO4, NO3, F, Fe, PO4)
2. **Edge Inference**: `infer_edges_probabilistic()` builds flow network from coordinates/heads
3. **Transport Fitting**: For each edge, evaluates evaporation (`gamma`) vs mixing (`f`) models
4. **Reaction Fitting**: LASSO coordinate descent finds sparse mineral reactions explaining residuals
5. **Output**: `EdgeResult` objects with transport parameters, reaction extents, and uncertainties

### Web Application Structure

```
web/
├── backend/
│   ├── app/
│   │   ├── main.py              # FastAPI app entry point
│   │   ├── hydrosheaf_adapter.py # Data transformation layer (mg/L → mmol/L)
│   │   └── routers/
│   │       ├── analysis.py      # /api/analysis/* - Hydrosheaf integration
│   │       ├── network.py       # /api/network/* - Flow inference
│   │       ├── samples.py       # /api/samples/* - Data management
│   │       └── projects.py      # /api/projects/* - Project management
│   └── requirements.txt
└── frontend/
    └── src/                     # React + Vite application
```

The `hydrosheaf_adapter.py` module handles conversion between:
- Frontend format (lowercase fields, mg/L units)
- Hydrosheaf format (capitalized fields, mmol/L units)

## Key Concepts

### Config Object
All model behavior is controlled through the `Config` dataclass. Key parameters:
- `lambda_l1` - LASSO regularization strength (sparsity)
- `phreeqc_enabled` - Enable thermodynamic constraints
- `isotope_enabled` - Enable δ18O/δ2H analysis
- `active_minerals` - List of minerals to include in reaction dictionary

### EdgeResult
The primary output object containing:
- `transport_model` - "evap" or "mix"
- `gamma` - Evaporation factor (≥1.0)
- `f` - Mixing fraction (0-1)
- `z_extents` - Molar reaction extents for each mineral
- `anomaly_norm` - Residual unexplained by model

### Ion Order Convention
The default ion order is: `["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]`

All concentration vectors must follow this 10-element order. The `Config.ion_order` field can customize this.

### Units
- Internal calculations use **mmol/L**
- Web frontend accepts **mg/L** (converted by adapter)
- Use `hydrosheaf.data.units.mgL_to_mmolL()` for manual conversion

## Testing

Tests use Python's `unittest` framework. Most tests are self-contained and don't require external dependencies.

```bash
# Run specific test categories
python -m pytest tests/test_phreeqc*.py  # PHREEQC integration (requires PHREEQC)
python -m pytest tests/test_accuracy*.py # Accuracy/regression tests
python -m pytest tests/test_doc_examples.py # Validate documentation examples
```

## API Endpoints Summary

- `POST /api/analysis/run` - Start Hydrosheaf analysis
- `GET /api/analysis/results/{job_id}` - Get analysis results
- `POST /api/network/create` - Create flow network
- `POST /api/network/{id}/infer-flow` - Bayesian flow inference
- `POST /api/samples/upload` - Upload sample data
- `GET /api/health` - Health check
