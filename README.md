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

## Installation

```bash
git clone https://github.com/dabdul-wahab1988/Hydrosheaf.git
cd Hydrosheaf
pip install .
```

To install with plotting dependencies:

```bash
pip install .[plot]
```

## Requirements

- Python >= 3.8
- PHREEQC (optional, for thermodynamic constraints)

## Usage

Hydrosheaf provides a command-line interface (CLI) for running the modeling pipeline.

```bash
hydrosheaf --help
```

### Basic Workflow

1. Prepare your input data (chemical concentrations, hydraulic heads, coordinates).
2. Configure the model parameters (weights, penalties, active minerals).
3. Run the inference pipeline.

## Documentation

**[Technical Document (PDF)](hydrosheaf_technical_document.pdf)**: Comprehensive mathematical theory and proofs.

Detailed guides are available in the `docs/` directory:

- **[User Guide](docs/USER_GUIDE.md)**: Extended usage instructions and workflows.
- **[Mathematical Reference](docs/math.md)**: Compact math notes aligned with the codebase implementation.
- **[Technical Reference](docs/TECHNICAL_REFERENCE.md)**: Deep dive into the code architecture and modules.
- **[PHREEQC Integration](docs/phreeqc.md)**: Setting up thermodynamic constraints.
- **[Examples](docs/examples.md)**: Walkthroughs of common modeling scenarios.

## Authors

**Dickson Abdul-Wahab**
**Ebenezer Aquisman Asare**

## License

[MIT License](LICENSE) (See LICENSE file for details)
