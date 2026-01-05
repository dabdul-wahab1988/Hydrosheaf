# Hydrosheaf

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

For a detailed mathematical description of the algorithms, optimization problems, and theoretical proofs, please refer to the [Technical Document](hydrosheaf_technical_document.pdf).

## Author

**Dickson Abdul-Wahab**

## License

[MIT License](LICENSE) (See LICENSE file for details)
