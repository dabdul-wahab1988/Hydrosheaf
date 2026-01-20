# Hydrosheaf Functional Specification Document (FSD)

**Version:** 1.0.0
**Status:** Active Development
**Target Audience:** Scientific Supervisors, Hydrogeologists, and Research Engineers.

---

## 1. Executive Summary

**Hydrosheaf** is an advanced Python-based framework for **Inverse Hydrogeochemical Modeling**. It solves the "Sheaf" problem: given a set of water quality samples scattered in space and time, what is the optimal network of flow paths and geochemical reactions that explains the observed data?

Unlike traditional forward models (e.g., PHREEQC, MODFLOW) that require exhaustive parameterization upfront, Hydrosheaf works backwards from the data. It uses **sparse optimization (L1 regularization)** and **Bayesian inference** to reconstruct the simplest, most scientifically probable explanation for groundwater evolution.

### Key Capabilities
1.  **Network Inference**: Automatically identifies flow connections between wells based on hydraulic heads and chemical similarity.
2.  **Process Identification**: Quantifies mixing ratios, evaporation extents, and mineral dissolution/precipitation rates for every connection.
3.  **Source Discrimination**: Distinguishes nitrate sources (Manure vs. Fertilizer) using a multi-evidence Bayesian engine.
4.  **Self-Calibration**: Includes a built-in PEST++ engine to automatically optimize physical parameters (dispersivity, reaction rates) against field data.

---

## 2. Scientific Workflow

The Hydrosheaf pipeline follows a structured scientific logic:

1.  **Data Ingestion**: Loads hydrochemistry (major ions, isotopes) and spatial data (coordinates, screen depths).
2.  **Topology Construction**:
    *   Builds a 3D graph of potential connections based on aquifer layers and spatial proximity.
    *   Infers missing hydraulic heads using a Hierarchical Bayesian Model constrained by topography.
3.  **Inverse Inference (The Core)**:
    *   For every potential connection $u \to v$, solves the inverse problem:
        $$v = \text{Mixing}(u, \text{Endmember}) + \text{Reactions} + \text{Evaporation}$$
    *   Uses L1-regularization to enforce sparsity (Occam's Razor: favor fewer reactions unless data demands them).
4.  **Temporal Analysis**:
    *   If time-series data exists, uses Cross-Correlation and Deconvolution to estimate groundwater residence times ($ \tau $).
5.  **Validation**:
    *   Runs forward kinetic simulations (PHREEQC) to verify that the inferred reactions are thermodynamically plausible.
6.  **Reporting**: Outputs detailed JSON results, geochemical interpretation reports, and visualization plots.

---

## 3. Core Architecture & Modules

### 3.1 Command Line Interface (`hydrosheaf.cli`)
*   **Purpose**: The primary entry point for executing the pipeline.
*   **Key Commands**:
    *   `hydrosheaf`: Runs the full inverse modeling pipeline.
    *   `hydrosheaf-cal`: Runs the calibration subsystem.
*   **Inputs**: CSV samples, Edge definitions (optional), YAML configuration.

### 3.2 High-Level API (`hydrosheaf.api`)
*   **Purpose**: Orchestrates the complex interaction between sub-modules.
*   **Key Component**: `fit_network_pipeline` manages the flow of data from ingestion $\to$ topology $\to$ inference $\to$ validation.

### 3.3 Configuration (`hydrosheaf.config`)
*   **Purpose**: Centralized control of scientific parameters.
*   **Features**: Strongly typed configuration for ~100 parameters (e.g., `lambda_l1` penalty strength, `dispersivity_m`). Supports dynamic updates from calibration results.

### 3.4 Logging & Transparency (`hydrosheaf.log`)
*   **Purpose**: Ensures scientific reproducibility.
*   **Feature**: Dual-channel logging records a high-level summary to the console and a granular **Audit Trail** to a file, documenting every decision (e.g., "Edge A->B rejected due to thermodynamic infeasibility").

---

## 4. Inverse Modeling Engine

### 4.1 Inference Logic (`hydrosheaf.inference`)
*   **Purpose**: The mathematical heart of the system.
*   **`edge_fit.py`**: Solves the inverse problem for a single edge. It evaluates multiple competing hypotheses (e.g., "Is this change due to Evaporation or Mixing?") and selects the best one based on the Akaike Information Criterion (AIC) and residual error.
*   **`network_fit.py`**: Applies `edge_fit` across the entire graph topology.

### 4.2 Geochemical Models (`hydrosheaf.models`)
*   **`reactions.py`**: Implements the sparse linear solver for mineral mass transfer ($\Delta = S \cdot \xi$).
*   **`mixing.py`**: Handles affine transformations for conservative mixing and Rayleigh fractionation (evaporation).
*   **`gibbs.py`**: Implements Gibbs diagrams to classify waters (Rock/Rain/Evaporation dominance) as a prior constraint.
*   **`ec_tds.py`**: Machine learning (Linear Regression) to fill missing TDS data from Electrical Conductivity.

### 4.3 Compositional Data Analysis (`hydrosheaf.coda_sbp`)
*   **Purpose**: robust statistics for chemical ratios.
*   **Innovation**: Uses **Sequential Binary Partitions (SBP)** to transform ion concentrations into orthogonal **isometric log-ratio (ilr)** coordinates. This allows standard statistical methods (Z-scores, PCA) to be applied to compositional data without generating spurious correlations.

---

## 5. Calibration & Optimization (`hydrosheaf.calibration`)

*   **Purpose**: Automates the search for unknown physical parameters.
*   **Engine**: **PESTGLM**, a native Python implementation of the industry-standard PEST++ algorithm.
*   **Capabilities**:
    *   **Parallelization**: Distributed Jacobian calculation across multiple CPU cores.
    *   **Regularization**: Tikhonov regularization to prevent overfitting.
    *   **Adapters**: Plug-and-play interfaces to calibrate different model components (Vadose Zone, Transport, Kinetics).

---

## 6. Physics & Transport Modules

### 6.1 Reactive Transport (`hydrosheaf.reactive_transport`)
*   **Purpose**: Forward validation.
*   **Method**: Generates PHREEQC input files (`KINETICS` blocks) based on the inverse model's results. If the forward simulation matches the observed data, the inverse result is validated.

### 6.2 Saturated Transport (`hydrosheaf.transport`)
*   **Purpose**: Simulates physical flow in the aquifer.
*   **Integration**: Wraps **FloPy** (MODFLOW/MT3DMS) to simulate 1D advection-dispersion along flow paths.

### 6.3 Vadose Zone (`hydrosheaf.vadose`)
*   **Purpose**: Simulates recharge from the soil surface.
*   **Method**: Solves the 1D Richards Equation to estimate recharge timing and magnitude.
*   **`calibrate.py`**: Special routine to fit soil hydraulic properties ($K_s$, $\alpha$) to moisture observations.

### 6.4 Physics Priors (`hydrosheaf.physics`)
*   **Purpose**: Constrains the graph using physical laws.
*   **`modpath.py`**: Imports particle tracking results from external MODFLOW models to define valid hydraulic connections.

---

## 7. Specialized Scientific Modules

### 7.1 Nitrate Source Tracking (`hydrosheaf.nitrate_source_v2`)
*   **Purpose**: Distinguishes between Manure and Fertilizer nitrate sources.
*   **Algorithm**: **Bayesian Evidence Fusion**. It combines three independent lines of evidence into a single probability score:
    1.  **Chemical Ratios**: CoDA-based Z-scores of NO3/Cl, NO3/K.
    2.  **Dual Isotopes**: MCMC mixing model of $\delta^{15}N$ and $\delta^{18}O$.
    3.  **Process Evidence**: Denitrification extent derived from the inverse graph.

### 7.2 Head Inference (`hydrosheaf.graph.head_inference`)
*   **Purpose**: Estimates groundwater levels in unmonitored wells.
*   **Method**: **Hierarchical Bayesian Model**. It uses surface topography as a "prior" for hydraulic head, updating the estimate based on nearby measurements and depth-to-water data.

### 7.3 Temporal Analysis (`hydrosheaf.temporal`)
*   **Purpose**: Estimates groundwater age.
*   **Methods**: Cross-correlation of tracer time-series and Deconvolution of Transit Time Distributions (TTDs).

---

## 8. Data & Outputs

### 8.1 Data Layer (`hydrosheaf.data`)
*   Handles data ingestion, unit conversion (mg/L $\leftrightarrow$ mmol/L), and Quality Control (Charge Balance Error checks).

### 8.2 Outputs (`hydrosheaf.outputs`)
*   **`interpret.py`**: Generates text-based "Geochemical Reports" for decision makers.
*   **`science_plots.py`**: Generates publication-ready figures (Piper diagrams, mixing plots, uncertainty ridges).
*   **Exports**: Full results available as JSON/CSV for integration with GIS or other tools.
