# Hydrosheaf

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18157915.svg)](https://doi.org/10.5281/zenodo.18157915)

**Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry: A Mathematical Framework for Inverse Geochemical Modeling**

## Overview

Hydrosheaf is a Python-based framework for solving inverse problems in groundwater hydrogeochemistry. It determines optimal transport processes (evaporation, mixing) and geochemical reactions that explain observed chemical evolution along flow paths in aquifer networks.

The framework integrates:

- **Weighted Least Squares Optimization** for transport model selection.
- **Sparse LASSO Regression** with coordinate descent for parsimonious reaction fitting.
- **Thermodynamic Constraints** utilizing PHREEQC for saturation index calculations.
- **Isotope Hydrogeology** for independent process validation (Deuterium Excess, LMWL).
- **Graph Theory** for probabilistic flow network inference.

## Features

- **Transport Modeling**: Distinguishes between evaporation (Rayleigh distillation-like) and end-member mixing.
- **Reaction Path Modeling**: Identifies mineral dissolution/precipitation, redox reactions (denitrification), and ion exchange.
- **Sparsity**: Uses L1 regularization to find the simplest chemical explanation for observed data.
- **Physical Consistency**: Enforces thermodynamic bounds (e.g., minerals cannot precipitate from undersaturated solutions).
- **Network Inference**: Infers flow direction probabilities from hydraulic head data with uncertainty.
- **Nitrate Source Discrimination**: Distinguishes between manure and fertilizer sources using a hybrid approach: Dual Isotope Bayesian Mixing ($\delta^{15}\text{N}, \delta^{18}\text{O}$) prioritized over CoDA-based hydrochemical statistics.
- **Reactive Transport**: Validates inverse results against kinetic rate laws (Arrhenius) and Damköhler numbers.
- **3D Flow Networks**: Analyzes layered aquifer systems with vertical anisotropy and topographic Bayesian priors.
- **Temporal Dynamics**: Resolves time-variant signals with cross-correlation residence times and seasonal decomposition.
- **Uncertainty Quantification**: Provides rigorous confidence intervals via Bayesian MCMC (NUTS) and Bias-Corrected Bootstrap (BCa).
- **Web Application**: Full-stack web interface with FastAPI backend and React frontend that supports **Real-time Analysis** monitoring via WebSockets.
- **Verified Documentation**: All mathematical examples in the technical reference are computationally verified by the test suite.

## Installation

### Core Package (CLI)

```bash
git clone https://github.com/dabdul-wahab1988/Hydrosheaf.git
cd Hydrosheaf
pip install .
```

### External Dependencies

**For PEST++ Calibration (Auto-handled):**
Hydrosheaf now automatically manages external dependencies for advanced calibration.
When you run a calibration routine that requires PEST++, MODFLOW, or MT3DMS, the necessary binaries will be **automatically downloaded** and configured for your operating system (Windows, Linux, or macOS).

No manual compilation or download is required.

**Note:**
1. **Core inverse modeling** works entirely within Python and needs no external tools.
2. **Advanced calibration** (using PEST++) triggers the auto-download on first use.

### Web Application

The web application provides a modern browser-based interface for Hydrosheaf. It consists of a FastAPI backend and a React frontend.

#### Prerequisites

- Python >= 3.8
- Node.js >= 18.x
- npm >= 9.x

#### Backend Setup

```bash
# Navigate to the backend directory
cd web/backend

# Create a virtual environment (recommended)
python -m venv venv

# Activate the virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Start the backend server
uvicorn app.main:app --reload --port 8000
```

The API will be available at `http://localhost:8000` with interactive documentation at `http://localhost:8000/api/docs`.

#### Frontend Setup

```bash
# Navigate to the frontend directory
cd web/frontend

# Install dependencies
npm install

# Start the development server
npm run dev
```

The frontend will be available at `http://localhost:5173` and expects the backend at `http://localhost:8000` (override with `VITE_API_BASE`).

#### Development Roadmap

For details on the current integration status, known issues, and future plans for the web interface, please refer to the [Integration Roadmap](web/INTEGRATION_ROADMAP.md).

#### Running Both Services

For development, run the backend and frontend in separate terminal windows:

**Terminal 1 (Backend):**
```bash
cd web/backend
uvicorn app.main:app --reload --port 8000
```

**Terminal 2 (Frontend):**
```bash
cd web/frontend
npm run dev
```

Open your browser and navigate to `http://localhost:5173` to access the Hydrosheaf web application.

## Quick Start

### ✅ Works Out-of-the-Box

When you do `pip install .`, you get **full functionality** for:
- Core inverse geochemical modeling (transport + reactions)
- PHREEQC thermodynamic constraints
- Isotope analysis and forensics
- Nitrate source discrimination (Bayesian MCMC)
- Network inference and topology refinement
- 3D flow networks
- Temporal dynamics and residence time estimation
- Uncertainty quantification (Bootstrap, MCMC)
- **Advanced Calibration (PEST++)**: Binaries are auto-downloaded when needed.

**This is sufficient for 95% of use cases!**

### Python API Usage

The recommended entry point is `fit_network_pipeline()` from `hydrosheaf.api`:

```python
from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import Config

# Load your data
samples = {...}  # List of sample dictionaries
edges = [...]    # List of edge definitions
config = Config()

# Run the pipeline
results, diagnostics = fit_network_pipeline(
    samples=samples,
    edges=edges,
    config=config,
    auto_disable_missing=True
)

# Access results
for edge_result in results:
    print(f"Edge {edge_result.u} → {edge_result.v}:")
    print(f"  Transport model: {edge_result.transport_model}")
    print(f"  Reactions: {edge_result.z_labels}")
    print(f"  Reaction extents: {edge_result.z_extents}")
```

### ⚙️ Optional: Manual Compilation

**Only needed if you want to modify PEST++ source code:**

#### System Requirements for Compilation
- Windows 10+ with Visual Studio 2022 Build Tools
- CMake and Ninja
- ~500 MB disk space

#### Build PEST++
```bash
# Compile from source (NOT necessary for standard usage)
compile_pestpp.bat
```

### Configuration

Hydrosheaf uses a central `Config` dataclass. Key settings:

```python
from hydrosheaf.config import Config

config = Config(
    # Chemistry
    ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
    charge_balance_limit=0.1,
    
    # Inference
    lambda_l1=0.01,  # LASSO penalty
    transport_models_enabled=["evap", "mix"],
    
    # Thermodynamics
    phreeqc_enabled=True,
    si_threshold_tau=0.2,
    
    # Isotopes & Nitrate
    isotope_enabled=True,
    nitrate_source_enabled=True,
    nitrate_source_isotope_mcmc_enabled=True,
    
    # Temporal
    residence_time_method="bayesian_lag",
    
    # 3D Flow
    network_3d_enabled=False,
    vertical_anisotropy=0.1,
)
```

For detailed configuration reference, see [HYDROSHEAF_USER_MANUAL.md](HYDROSHEAF_USER_MANUAL.md#16-configuration-reference).

### Data Input Format

Samples should be provided as a list of dictionaries with keys:

```python
samples = [
    {
        "site_id": "W001",
        "sample_id": "W001_2023-01",
        "Ca": 80.0,      # mg/L (converted internally to mmol/L)
        "Mg": 20.0,
        "Na": 15.0,
        "K": 2.0,
        "HCO3": 250.0,
        "Cl": 30.0,
        "SO4": 50.0,
        "NO3": 5.0,
        "F": 1.0,
        "Fe": 0.1,
        "PO4": 0.05,
        "pH": 7.5,
        "d18o": -8.5,    # Optional: δ18O (‰)
        "d2h": -60.0,    # Optional: δ2H (‰)
        "latitude": 40.5,
        "longitude": -74.0,
        "elevation": 50.0,
        "head": 45.0,    # Optional: hydraulic head (m)
    },
    # ... more samples
]
```

Edges should specify connectivity:

```python
edges = [
    {"u": "W001", "v": "W002", "distance_km": 2.5, "delta_h": 1.5},
    {"u": "W002", "v": "W003", "distance_km": 3.0, "delta_h": 2.0},
]
```



## Documentation

For comprehensive reference, see:

- **[HYDROSHEAF_USER_MANUAL.md](HYDROSHEAF_USER_MANUAL.md)**: Complete technical reference for all modules and functions (v0.1.0)
- **[DEVELOPMENT.md](DEVELOPMENT.md)**: For developers and contributors (building from source, running tests, contributing)
- **[Technical Document (PDF)](docs/papers/hydrosheaf_technical_document.pdf)**: Mathematical theory and proofs
- **[User Guide](docs/USER_GUIDE.md)**: Extended usage instructions and workflows
- **[Mathematical Reference](docs/math.md)**: Compact math notes aligned with the codebase
- **[Technical Reference](docs/TECHNICAL_REFERENCE.md)**: Code architecture and module details
- **[PHREEQC Integration](docs/phreeqc.md)**: Setting up thermodynamic constraints
- **[Examples](docs/examples.md)**: Walkthroughs of common scenarios

## Troubleshooting & Common Questions

### Q: Do I need to compile anything to use Hydrosheaf?
**A: No!** Simply run `pip install .` and everything works. Compilation is only optional for PEST++ advanced calibration.

### Q: What if I can't compile PEST++?
**A: That's fine!** Core hydrosheaf inverse modeling works perfectly without PEST++. Download pre-compiled PEST++ binaries from [USGS/pestpp releases](https://github.com/usgs/pestpp/releases) if needed.

### Q: Will it work on Linux/Mac?
**A: Yes!** The Python package is cross-platform. The `bin/` directory contains Windows executables, but:
- Core hydrosheaf works on Windows/Linux/Mac
- Download Linux/Mac binaries separately if needed
- Most users don't need the `bin/` files

### Q: Can I use Hydrosheaf without PHREEQC?
**A: Yes!** While the PHREEQC library is installed by default, you can disable thermodynamic constraints in your configuration if you don't need them.

### Q: I got an error about missing dependencies
**A:** Ensure you have installed the package using `pip install .`. All necessary dependencies are now included by default.

### Q: How do I run the web interface?
**A: Follow the Web Application section above.** It requires Node.js, but Python setup is the same.

---

## Authors

**Dickson Abdul-Wahab**,
**Ebenezer Aquisman Asare**,
**Abdul Rashid Dickson**

## License

[MIT License](LICENSE) (See LICENSE file for details)
