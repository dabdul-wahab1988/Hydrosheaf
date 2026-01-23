# Hydrosheaf Reference Manual

**Advanced Groundwater Inverse Modeling Framework**

**Version:** 0.1.0  
**Date:** January 2026  
**License:** MIT

---

## Table of Contents

1.  [Overview](#1-overview)
2.  [High-Level Pipeline (`api`)](#2-high-level-pipeline)
3.  [Module 1: Inference Engine (`inference`)](#3-module-inference)
4.  [Module 2: Network Topology & Sheaf Theory (`graph`, `sheaf`)](#4-module-topology)
5.  [Module 3: 3D Topology (`graph3d`)](#5-module-3d-topology)
6.  [Module 4: Vadose Zone Physics (`vadose`)](#6-module-vadose)
7.  [Module 5: Reaction Thermodynamics (`phreeqc`, `models`)](#7-module-reactions)
8.  [Module 6: Calibration & Uncertainty (`calibration`, `uncertainty`)](#8-module-calibration)
9.  [Module 7: Nitrate Forensics (`nitrate_source_v2`)](#9-module-nitrate)
10. [Module 8: Temporal Dynamics (`temporal`)](#10-module-temporal)
11. [Module 9: Reactive Transport Validation (`reactive_transport`)](#11-module-validation)
12. [Module 10: Data Integrity & Outputs (`data`, `outputs`)](#12-module-data--outputs)
13. [Module 11: Physics Priors & Forward Transport (`physics`, `transport`)](#13-module-physics--transport)
14. [Module 12: Hyperparameter Tuning (`tuning`)](#14-module-hyperparameter-tuning)
15. [Configuration Reference (`config`)](#15-configuration-reference)

---

## 1. Overview

Hydrosheaf is a Python library for inverse hydrogeochemical modeling. It infers flow paths, reaction networks, and transport parameters from sparse field data. This manual provides a comprehensive reference for the software's architecture, documenting every major module and function available to the researcher.

---

## 2. High-Level Pipeline

**Path:** `hydrosheaf.api`  
**Purpose:** Streamlined entry points for running complex modeling workflows.

### 2.1 Network Pipeline
**Function:** `fit_network_pipeline(samples, edges, config, *, auto_disable_missing=True, physics_priors=None, physics_priors_mode="override", temporal_nodes=None, temporal_hydraulic_params=None, phreeqc_results=None)`

**Returns:** `Tuple[List[EdgeResult], Dict[str, object]]`

**Description:** The recommended high-level interface. It orchestrates the entire process:
1.  **Validation:** Automatically validates inputs or disables missing modules via `auto_disable_missing_modules`.
2.  **Physics Priors:** Applies optional `PhysicsPrior` objects (e.g., from MODPATH) to the graph topology.
3.  **Temporal Fitting:** If `temporal_nodes` are provided, it fits residence times using `fit_temporal_edges`.
4.  **Inverse Modeling:** Calls `fit_network` to solve transport and reactions for every edge.
5.  **Coupling:** Attaches temporal summaries to the final `EdgeResult` objects.

### 2.2 Network Fitting with Priors
**Function:** `fit_network_with_priors(samples, edges, config, *, auto_disable_missing=True, physics_priors=None, physics_priors_mode="override", phreeqc_results=None, residence_time_overrides=None)`

**Returns:** `Tuple[List[EdgeResult], List[Edge]]`

**Description:** Alternative interface that applies physics priors without temporal fitting.

### 2.3 Utility Functions
**Function:** `validate_required_inputs(samples, config)`  
**Description:** Raises `ValueError` if required inputs for enabled modules are missing.

**Function:** `auto_disable_missing_modules(samples, config) -> Config`  
**Description:** Returns a modified config with feature flags disabled if required inputs are missing.

**Function:** `build_vadose_priors(profile, forcing, links, *, config=None, water_table_depth_m=None)`  
**Returns:** `Tuple[List[PhysicsPrior], Dict[str, Dict[str, object]]]`  
**Description:** Run vadose zone physics and return `PhysicsPrior` objects for use in network fitting.

**Function:** `fit_temporal_edges(temporal_nodes, edges, config, *, hydraulic_params_by_edge=None)`  
**Returns:** `Tuple[List[TemporalEdgeResult], Dict[str, float]]`  
**Description:** Fit temporal edges and return residence time overrides for use in `fit_network`.

---

## 3. Module: Inference

**Path:** `hydrosheaf.inference`  
**Purpose:** Core inverse solver that fits transport and reaction models to observed data.

### 3.1 Network Fitting
**Function:** `fit_network(samples, edges, config, ...)`  
**Description:** Iterates over the provided graph topology and solves the inverse problem for every edge.
*   **Algorithm:** 
    1.  Validates input samples against schema.
    2.  Builds thermodynamic bounds (using PHREEQC) for each edge.
    3.  Calls `fit_edge` for each connection.
    4.  Aggregates results into a coherent network solution.
    5.  Optionally runs uncertainty propagation (Bootstrap/Monte Carlo).

### 3.2 Edge Fitting
**Function:** `fit_edge(x_u, x_v, config, ...)`  
**Path:** `hydrosheaf.inference.edge_fit`  
**Returns:** `EdgeResult`  
**Description:** Solves the inverse problem for a single pair of nodes ($u \to v$).
*   **Mathematical Model:** Minimizes $\|\mathbf{x}_v - \text{Transport}(\mathbf{x}_u) - \mathbf{S}\mathbf{z}\|_{\mathbf{W}}^2 + \lambda \|\mathbf{z}\|_1$.
*   **Key Logic:**
    *   **Transport Separation:** Uses `conservative_weights` (default: Cl at 1.0, others at 0.01) to fit mixing/evaporation parameters first, preventing reaction bias.
    *   **Reaction Fitting:** Uses LASSO coordinate descent to fit mineral mass transfers ($\mathbf{z}$) subject to thermodynamic bounds.
    *   **Candidate Selection:** Evaluates multiple transport hypotheses (evaporation, mixing) and selects the best based on weighted residual scoring.

### 3.3 EdgeResult Dataclass
**Class:** `EdgeResult`  
**Key Fields:**
- `edge_id`, `u`, `v`: Edge identifier and node endpoints
- `transport_model`: Selected transport model ("evap" or "mix")
- `gamma`: Evaporation parameter (0-1, where 1 = no evaporation)
- `f`: Mixing fraction (0-1, where 0 = pure upstream, 1 = 100% endmember)
- `z_extents`, `z_labels`: Reaction extents and mineral/redox names
- `transport_residual_norm`, `objective_score`: Goodness-of-fit metrics
- `l1_norm`: L1 penalty (sparsity measure)
- `si_u`, `si_v`: Saturation indices at upstream/downstream nodes
- Temporal fields: `temporal_residence_time_days`, `temporal_reaction_extents_mean`, etc.
- Uncertainty fields: `gamma_std`, `gamma_ci_low`, `gamma_ci_high`, `extents_std`, `extents_ci_low`, `extents_ci_high`
- Nitrate source fields: `nitrate_source_p_manure`, `nitrate_source_logit`, `nitrate_source_evidence`

---

## 4. Module: Topology

**Path:** `hydrosheaf.graph`, `hydrosheaf.sheaf`  
**Purpose:** Infers the flow network structure from sparse head data and refines it using geochemical consistency.

### 4.1 Graph Building
**Function:** `build_edges(edges)`  
**Path:** `hydrosheaf.graph.build`  
**Returns:** `List[Edge]`  
**Description:** Normalizes edge input format to standardized `Edge` dataclass.

**Function:** `infer_edges_from_coordinates(node_ids, lat, lon, config, ...)`  
**Description:** Simple nearest-neighbor edge inference based on Euclidean distance.

**Function:** `infer_edges_probabilistic(node_ids, lat, lon, head_meas, elevation, config, ...)`  
**Description:** Constructs a flow graph where edge existence is probabilistic.
*   **Logic:** Uses a Linear-Gaussian Bayesian model to estimate hydraulic head distributions at every node. Computes $P(\text{flow } u \to v) = P(h_u > h_v)$ by integrating the joint posterior.
*   **Key Parameters (from Config):**
    - `edge_p_min`: Minimum probability threshold for edge inclusion (default: 0.75)
    - `edge_radius_km`: Maximum distance for candidate edges (default: 5.0)
    - `edge_max_neighbors`: Max neighbors per node (default: 3)
    - `edge_head_inference`: Method for head estimation ("heuristic", "bayesian", or "bayesian_mcmc")

### 4.2 Sheaf Refinement
**Function:** `refine_edges_with_sheaf(edges, samples, config, ...)`  
**Path:** `hydrosheaf.sheaf.topology_refine`  
**Returns:** `List[Edge]`  
**Description:** Filters the probabilistic graph using isotope and chloride consistency.
*   **Theory:** Uses isotope and chloride gradients to identify physically inconsistent connections.
*   **Algorithm:**
    *   For each edge, computes isotope and chloride cost functions.
    *   Calculates local and global "sheaf scores" indicating consistency.
    *   Returns edges ranked by consistency.
*   **Key Parameters (from Config):**
    - `sheaf_isotope_enabled`, `sheaf_cl_enabled`: Toggle isotope/Cl consistency checks
    - `sheaf_iso_sigma_d18o`, `sheaf_iso_sigma_d2h`: Measurement uncertainty for isotopes
    - `sheaf_weight_*`: Relative weights for different consistency terms

### 4.3 Head Estimation
**Class:** `HeadPosterior`  
**Path:** `hydrosheaf.graph.head_inference`  
**Description:** Holds posterior distribution of hydraulic head across nodes.
*   **Fields:**
    - `node_ids`: List of node identifiers
    - `head_mean`: Mean head estimate (vector)
    - `head_cov`: Covariance matrix of head uncertainty
    - `mu_dtw_mean`, `mu_dtw_sigma`: Mean depth-to-water and its uncertainty

**Function:** `bayesian_head_estimate(nodes, config, ...)`  
**Description:** Computes joint posterior head distribution using hierarchical Bayesian model.

---

## 5. Module: 3D Topology

**Path:** `hydrosheaf.graph3d`  
**Purpose:** Advanced network inference for multi-aquifer systems, accounting for vertical gradients and aquitard leakage.

### 5.1 3D Edge Inference
**Function:** `infer_edges_3d_probabilistic(nodes, config, ...)`  
**Description:** Extends the 2D probabilistic inference to 3D space.
*   **Logic:** Computes a composite probability $P_{uv} = P_{head} \times P_{dist} \times P_{layer} \times P_{screen}$.
*   **$P_{dist}$:** Anisotropic distance decay (vertical distances are penalized by `vertical_anisotropy`).
*   **$P_{layer}$:** Probability of crossing aquitards (leakage probability).

---

## 6. Module: Vadose

**Path:** `hydrosheaf.vadose`  
**Purpose:** Models vertical recharge and contaminant transport through the unsaturated zone.

### 6.1 Richards Equation Solver
**Function:** `run_richards_column(profile, forcing, config, ...)`  
**Description:** Solves the 1D Richards Equation for variably saturated flow using Van Genuchten-Mualem parameterization.

### 6.2 Dual-Domain Transport
**Function:** `mixture_ttd_from_series(...)`  
**Description:** Generates the Travel Time Distribution (TTD) kernel using a Dual-Domain model (fast/slow flow paths).

### 6.3 Calibration
**Function:** `calibrate_ks_and_kc(...)`  
**Description:** Calibrates $K_s$ (Conductivity) and $K_c$ (Evapotranspiration) against Soil Moisture ($\theta$) and Tracer data.

---

## 7. Module: Reactions

**Path:** `hydrosheaf.models`, `hydrosheaf.phreeqc`  
**Purpose:** Handles chemical logic and thermodynamic constraints.

### 7.1 Thermodynamic Bounds
**Function:** `build_edge_bounds(...)`  
**Description:** Calls PHREEQC to calculate Saturation Indices (SI) to constrain mineral dissolution/precipitation.

### 7.2 Isotope Utilities
**Path:** `hydrosheaf.isotopes`  
**Description:** Utilities for $\delta^{18}$O and $\delta^{2}$H analysis, including LMWL fitting (`fit_lmwl`) and Deuterium Excess computation.

---

## 8. Module: Calibration & Uncertainty

**Path:** `hydrosheaf.calibration`, `hydrosheaf.uncertainty`  
**Purpose:** Parameter estimation and uncertainty quantification.

### 8.1 PEST++ Integration
**Path:** `hydrosheaf.calibration`  
**Description:** Wraps USGS PEST++ for high-dimensional parameter estimation across multiple edges. 
*   **Key Functions:**
    - `pestpp_glm_calibration(...)`: Generalized Linear Model calibration via pestpp-glm
    - `pestpp_ies_calibration(...)`: Iterative Ensemble Smoother calibration via pestpp-ies
*   **Typical Workflow:**
    1. Generate PEST control file from edge results and observations
    2. Run pestpp executable (if available locally)
    3. Parse resulting parameter ensemble
    4. Update edge parameters with optimized values

### 8.2 Uncertainty Quantification
**Function:** `bootstrap_edge_fit(x_u, x_v, reaction_matrix, reaction_labels, config, ...)`  
**Path:** `hydrosheaf.uncertainty.bootstrap`  
**Returns:** `UncertaintyResult`  
**Description:** Bias-Corrected and Accelerated (BCa) bootstrap to estimate confidence intervals for transport and reaction parameters.

**Function:** `bayesian_edge_fit(x_u, x_v, reaction_matrix, reaction_labels, config, ...)`  
**Path:** `hydrosheaf.uncertainty.bayesian`  
**Returns:** `UncertaintyResult`  
**Description:** Bayesian MCMC using HMC/NUTS (via PyMC5) to estimate posterior distributions.
*   **Parameters:**
    - `n_samples`: Number of MCMC draws per chain (default: 5000)
    - `n_chains`: Number of independent chains (default: 4)
    - `target_accept`: Target acceptance rate for HMC (default: 0.95)
*   **Requires:** `pip install pymc>=5.0 arviz>=0.15`

**Function:** `monte_carlo_propagate(edge_results, config, ...)`  
**Path:** `hydrosheaf.uncertainty.propagation`  
**Description:** Propagates uncertainty from individual edges through a network ensemble.

---

## 9. Module: Nitrate Forensics

**Path:** `hydrosheaf.nitrate_source_v2`  
**Purpose:** Identification of nitrate sources using isotopes and compositional chemistry.

### 9.1 Nitrate Source Discrimination
**Function:** `infer_node_posteriors(samples, nodes, config, ...)`  
**Path:** `hydrosheaf.nitrate_source_v2`  
**Returns:** `Dict[str, NitrateSourceResult]`  
**Description:** Estimates nitrate source fractions per node.
*   **Algorithm:**
    1. **Isotope MCMC (Primary):** If $\delta^{15}$N and $\delta^{18}$O available, runs Bayesian mixing model (PyMC) to estimate manure vs. synthetic fertilizer fractions.
    2. **Compositional Data Analysis (Secondary):** Uses CoDA-based isometric log-ratio (ILR) transforms via `hydrosheaf.coda_sbp` to analyze hydrochemical ratios (K/NO3, Cl/NO3, etc.).
    3. **Evaporation Gating:** Uses Deuterium Excess to detect evaporation and gate unreliable signals.
    4. **Edge Coupling:** Incorporates denitrification extents inferred from the flow network to adjust source estimates.

### 9.2 NitrateSourceResult
**Dataclass:** `NitrateSourceResult`  
**Key Fields:**
- `p_manure`: Probability or fraction of manure source (0-1)
- `p_fertilizer`: Probability of synthetic fertilizer (computed as 1 - p_manure)
- `logit_score`: Raw logit-scale evidence (-inf to +inf)
- `top_evidence`: List of primary evidence lines (e.g., "MCMC mixing estimate")
- `gating_flags`: List of flags indicating data quality issues
- `ilr_valid`: Boolean indicating valid CoDA transformation
- `source_fractions`: Dictionary of source fractions (MCMC output)
- `ci_lower`, `ci_upper`: Credible intervals per source

### 9.3 Compositional Data Analysis
**Function:** `ilr_from_sbp(component_fractions, sbp_matrix)`  
**Path:** `hydrosheaf.coda_sbp`  
**Description:** Isometric log-ratio transformation for compositional analysis.

---

## 10. Module: Temporal Dynamics

**Path:** `hydrosheaf.temporal`  
**Purpose:** Estimates groundwater age and dynamics from time-series tracers.

### 10.1 Temporal Edge Fitting
**Function:** `fit_temporal_edge(node_u, node_v, config, edge_id="", hydraulic_params=None)`  
**Path:** `hydrosheaf.temporal.temporal_edge_fit`  
**Returns:** `TemporalEdgeResult`  
**Description:** Fits transport and reaction models across time-series data points.
*   **Algorithm:**
    1. Estimate residence time $\tau$ via multiple methods
    2. Align upstream and downstream time series using $\tau$ (with interpolation)
    3. For each time point, fit transport (evaporation/mixing) and reactions
    4. Aggregate to time-averaged parameters with uncertainty estimates
*   **Key Parameters:**
    - `residence_time_method`: Method for estimating $\tau$ ("cross_correlation", "bayesian_lag", "ttd", "recharge_piston")
    - `temporal_window_days`: Time window for analysis
    - `ttd_max_lag_days`: Maximum lag for TTD deconvolution

### 10.2 Residence Time Estimation Methods
**Function:** `estimate_residence_time_with_details(node_u, node_v, config, hydraulic_params=None)`  
**Path:** `hydrosheaf.temporal.residence_time`  
**Returns:** `Dict[str, object]` containing `residence_time_days`, `method`, `uncertainty`, etc.
**Methods:**
- **cross_correlation:** Lag-1 cross-correlation of time series (fastest, least accurate)
- **bayesian_lag:** Physics-informed Bayesian model with hydraulic constraints (preferred)
- **ttd:** Travel Time Distribution deconvolution (requires multiple tracers)
- **recharge_piston:** Piston-flow assumption for volcanic/karst aquifers

### 10.3 TemporalEdgeResult
**Dataclass:** `TemporalEdgeResult`  
**Key Fields:**
- `edge_id`, `u`, `v`: Edge identifier and endpoints
- `residence_time_days`: Estimated mean residence time
- `residence_time_method`: Method used for estimation
- `residence_time_uncertainty`: Standard deviation or CI
- `transport_model`: Selected transport model across time
- `gamma_mean`, `gamma_std`: Mean and std of evaporation parameter
- `f_mean`, `f_std`: Mean and std of mixing fraction
- `reaction_extents_mean`, `reaction_extents_std`: Time-averaged reaction extents
- `timestamps`: List of time points included in fit
- `residence_time_flags`: QC flags (e.g., "insufficient_data", "low_correlation")
- `residence_time_details`: Dictionary with method-specific diagnostic information

---

## 11. Module: Reactive Transport Validation

**Path:** `hydrosheaf.reactive_transport`  
**Purpose:** Forward verification of inverse modeling results using kinetic rate laws.

### 11.1 Forward Validation
**Function:** `validate_network_forward(edge_results, config, ...)`  
**Path:** `hydrosheaf.reactive_transport.validation`  
**Returns:** `Dict[str, object]` with validation metrics  
**Description:** Runs forward kinetic simulation using inferred reaction extents to verify physical consistency.
*   **Approach:**
    1. For each edge, extract inferred transport and reaction parameters
    2. Set up reactive transport system with Arrhenius kinetics
    3. Solve for concentration profiles over residence time
    4. Compare forward modeled outlet concentrations to observed values
    5. Calculate fit metrics (RMSE, NSE, Damköhler numbers)
*   **Key Parameters:**
    - `rt_simulator`: Backend ("phreeqc" or "mt3dms")
    - `rt_default_rate_constant`: Default kinetic rate constant (s⁻¹)

### 11.2 Kinetic Rate Laws
**Module:** `hydrosheaf.reactive_transport.rate_laws`  
**Functions:** Arrhenius temperature-dependent rate law implementations

### 11.3 Validation Metrics
**Module:** `hydrosheaf.reactive_transport.metrics`  
**Metrics:** Nash-Sutcliffe Efficiency (NSE), RMSE, Damköhler number, reaction progress ratio

---

## 12. Module: Data Integrity & Outputs

**Path:** `hydrosheaf.data`, `hydrosheaf.outputs`  
**Purpose:** Input validation, quality control, and results formatting.

### 12.1 Schema & Validation
**Function:** `normalize_sample(sample, config)`  
**Path:** `hydrosheaf.data.schema`  
**Description:** Normalizes a single sample dict to standard format, handling units and missing values.

**Function:** `qc_flags(sample, config)`  
**Path:** `hydrosheaf.data.qc`  
**Returns:** `List[str]`  
**Description:** Computes QC flags for a sample including:
- Charge Balance Error (CABE) 
- Detection limit warnings
- Physical range violations (negative concentrations, impossible pH, etc.)

**Function:** `normalize_samples(samples, config)`  
**Description:** Normalizes a sequence of samples.

### 12.2 Data Units
**Module:** `hydrosheaf.data.units`  
**Key Functions:**
- `mgL_to_mmolL(conc_mgL, ion_name, config)`: Convert mg/L to mmol/L
- `meqL_to_mmolL(meq, ion_name)`: Convert meq/L to mmol/L
- Hydrosheaf internally uses **mmol/L** for all concentration computations; web interfaces convert to mg/L for user display

### 12.3 Scientific Visualization
**Path:** `hydrosheaf.outputs.science_plots`  
**Functions:**
*   `plot_ttd_kernel(edge_results, max_tau=365.0, path=None)`: Visualizes transit time distributions (aquifer memory).
*   `plot_breakthrough(edge_results, path=None)`: Observed vs. Modeled contaminant breakthrough curves.
*   `plot_posterior_ridges(edge_results, path=None)`: Bayesian posterior distributions for transport and reaction parameters.
*   `plot_reactive_transport_validation(edge_results, path=None)`: Forward model fit metrics (RMSE/NSE per edge).

### 12.4 Network Visualization
**Path:** `hydrosheaf.outputs.plots_3d`  
**Function:** `plot_network_3d(edge_results, nodes, ..., path=None)`  
**Description:** 3D interactive visualization of flow network with color-coded properties (residence time, reaction extent, etc.).

### 12.5 Results Export
**Path:** `hydrosheaf.outputs.export`  
**Functions:**
- `export_edge_results_csv(edge_results, path)`: Export to CSV format
- `export_edge_results_json(edge_results, path)`: Export to JSON format
- `temporal_edge_results_to_rows(temporal_results)`: Convert temporal results to tabular format

### 12.6 Interpretation Reports
**Path:** `hydrosheaf.outputs.interpret`  
**Function:** `print_interpretation_report(edge_results, config, ...)`  
**Description:** Generates human-readable summary of flow paths, reactions, and uncertainties.

---

## 13. Module: Physics Priors & Forward Transport

**Path:** `hydrosheaf.physics`, `hydrosheaf.transport`  
**Purpose:** Integration with external physics engines (MODFLOW/MODPATH) for constraining the inverse model.

### 13.1 Physics Priors
**Class:** `PhysicsPrior`  
**Path:** `hydrosheaf.physics.priors`  
**Fields:**
- `u`, `v`: Source and target node IDs
- `weight`: Confidence in prior (0-1)
- `tau_mean_days`, `tau_std_days`: Residence time estimate and uncertainty
- `p10_days`, `p90_days`: Percentile bounds

**Function:** `priors_from_modpath_endpoints(path, nodes, config, ...)`  
**Description:** Reads MODPATH particle endpoint files to extract residence times and flow paths, creating `PhysicsPrior` objects for network inference.

**Function:** `apply_physics_priors(edges, physics_priors, mode="override")`  
**Description:** Applies physics priors to edge weights or topology.
*   **Modes:**
    - `"override"`: Replace edge attributes with physics priors
    - `"weight"`: Soft constraint, weighted into optimization

### 13.2 Forward Transport (1D)
**Path:** `hydrosheaf.transport`  
**Functions:**
- `build_1d_transport_model(config, ...)`: Constructs 1D advection-dispersion-reaction system
- `analytical_1d_transport(c0, q, D, L, t, reaction_rate)`: Ogata-Banks solution for linear reactive transport
*   **Description:** Used to verify inverse modeling results by forward simulation along inferred flow paths.

### 13.3 Vadose Zone Priors
**Module:** `hydrosheaf.vadose`  
**Function:** `build_vadose_edge_priors(profile, forcing, links, config=None, ...)`  
**Path:** `hydrosheaf.vadose.run`  
**Returns:** `Tuple[List[VadosePrior], Dict]`  
**Description:** Solves Richards Equation for unsaturated flow through vadose column, producing residence time and mixing/transport priors for saturated zone network fitting.
*   **Key Components:**
    - Richards solver (Van Genuchten-Mualem)
    - Dual-domain mass transfer
    - TTD kernel generation from soil moisture and tracer data

---

## 14. Module: Hyperparameter Tuning

**Path:** `hydrosheaf.tuning`  
**Purpose:** Automated selection of LASSO sparsity penalties ($\lambda$) using Stability Selection.

### 14.1 Reaction Network Tuning
**Function:** `stability_selection_lambda(x_u, x_v, reaction_matrix, config, ...)`  
**Path:** `hydrosheaf.tuning.reaction_tuning`  
**Returns:** `Dict[str, object]` with optimal $\lambda$ and robustness metrics  
**Description:** Uses subsampling and LASSO regression to identify robust reaction networks.
*   **Algorithm:**
    1. Subsample the data (e.g., 80% of samples)
    2. Fit LASSO over a grid of $\lambda$ values
    3. Track which reactions are selected as $\lambda$ varies
    4. Calculate selection frequency for each reaction
    5. Return $\lambda$ that yields high-frequency (stable) reactions
*   **Key Parameters:**
    - `n_subsamples`: Number of subsampling iterations
    - `subsample_fraction`: Fraction of data per subsample
    - `lambda_grid`: Grid of $\lambda$ values to evaluate

---

## 15. Module: 3D Topology

**Path:** `hydrosheaf.graph3d`  
**Purpose:** Advanced network inference for multi-aquifer systems with vertical gradients.

### 15.1 3D Edge Inference
**Function:** `build_network_3d(nodes_3d, config, ...)`  
**Path:** `hydrosheaf.graph3d.build_3d`  
**Returns:** `List[Edge3D]`  
**Description:** Constructs probabilistic flow network in 3D with vertical anisotropy.
*   **Composite Probability:** $P_{uv} = P_{\text{head}} \times P_{\text{distance}} \times P_{\text{layer}} \times P_{\text{screen}}$
*   **Key Features:**
    - Anisotropic distance decay (vertical distances penalized by `vertical_anisotropy`)
    - Aquitard leakage probability (`layer_enabled`)
    - Screen/perforation matching (`screen_match_enabled`)

### 15.2 3D Edge Type
**Class:** `Edge3D`  
**Fields:**
- Base fields from `Edge` (u, v, distance_km, delta_h)
- `distance_3d_m`: True 3D distance
- `horizontal_distance_m`, `vertical_distance_m`: Components
- `type_3d`: Classification (e.g., "horizontal", "vertical", "cross-layer")
- `layer_u`, `layer_v`: Aquifer layer assignments
- `prob_head`, `prob_distance`, `prob_layer`: Component probabilities

### 15.3 3D Constraints
**Module:** `hydrosheaf.graph3d.constraints`  
**Description:** Thermodynamic and physical constraints for layered systems (aquitard leakage, confining layer permeability).

---

## 16. Configuration Reference

**File:** `hydrosheaf/config.py`  
**Class:** `Config` (dataclass)  
**Version:** 0.1.0

Key configuration categories:

| Category | Key Parameters | Defaults |
| :--- | :--- | :--- |
| **Chemistry** | `ion_order` | DEFAULT_ION_ORDER (Ca, Mg, Na, HCO3, Cl, SO4, NO3, F, Fe, PO4) |
| | `unit`, `unit_mode` | "mmol/L", "mmol_L" |
| | `weights` | [1.0] * 10 |
| | `conservative_weights` | [0.01]*9 + [1.0] (Cl weight = 1.0) |
| | `charge_balance_limit` | 0.1 (meq/L) |
| | `detection_limit_policy` | "half" |
| **Inference** | `lambda_l1` | 0.0 (set manually or via tuning) |
| | `reaction_max_iter` | 300 |
| | `reaction_tol` | 1e-6 |
| | `transport_models_enabled` | ["evap", "mix"] |
| | `allow_signed_reactions` | False |
| **Thermodynamics** | `phreeqc_enabled` | True |
| | `phreeqc_mode` | "phreeqpython" |
| | `phreeqc_database` | Built-in phreeqc.dat |
| | `si_threshold_tau` | 0.2 |
| | `constraints_hard` | True |
| **Topology** | `edge_head_inference` | "heuristic" |
| | `edge_p_min` | 0.75 |
| | `edge_radius_km` | 5.0 |
| | `edge_max_neighbors` | 3 |
| **Isotopes** | `isotope_enabled` | False |
| | `isotope_d18o_key`, `isotope_d2h_key` | "d18o", "d2h" |
| | `lmwl_defined`, `lmwl_a`, `lmwl_b` | False, 8.0, 10.0 |
| **Nitrate** | `nitrate_source_enabled` | False |
| | `nitrate_source_isotope_mcmc_enabled` | True |
| | `min_mg_L` | 1.0 |
| **Temporal** | `residence_time_method` | "bayesian_lag" |
| | `residence_time_hydraulic_k` | 0.0 (m/day) |
| | `residence_time_porosity` | 0.0 |
| | `ttd_max_lag_days` | 365 |
| **Validation** | `reactive_transport_validation` | False |
| | `rt_simulator` | "phreeqc" |
| | `rt_default_rate_constant` | 1e-9 (s⁻¹) |
| **3D Flow** | `network_3d_enabled` | False |
| | `vertical_anisotropy` | 0.1 |
| | `layer_enabled` | True |

### 16.1 Config Creation & Usage
```python
from hydrosheaf.config import Config

# Default config
config = Config()

# Custom config
config = Config(
    ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
    lambda_l1=0.01,
    phreeqc_enabled=True,
    isotope_enabled=True,
    nitrate_source_enabled=True,
)
```

---

## 17. CLI Tools

**Entry Points:**
- `hydrosheaf`: Main CLI for network fitting, visualization, and export
- `hydrosheaf-cal`: PEST++ calibration interface

**Key Commands:**
```bash
# Fit network from CSV
hydrosheaf fit --samples samples.csv --edges edges.csv --config config.yaml

# Export results
hydrosheaf export --results results.json --format csv

# Plot results
hydrosheaf plot --results results.json --type ttd
```

---

## 18. Spatial Units & Conventions

**Internal Units:**
- Concentrations: **mmol/L** (all computations)
- Distances: **km** (edges), **m** (3D distances)
- Time: **days** (residence time, temporal analysis)
- Temperature: **°C**

**Web Frontend Units (converts automatically):**
- Concentrations: **mg/L** (displayed to users)
- Distances: **km** (same)

**Ion Order Convention:**
- Hydrosheaf uses a fixed 10-ion order: **Ca, Mg, Na, HCO3, Cl, SO4, NO3, F, Fe, PO4**
- Customizable via `Config.ion_order`
- Ensure consistency across input data

---

## 19. Web Application Architecture

**Backend:** FastAPI (Python)
- Location: `web/backend/`
- Run: `uvicorn app.main:app --reload --port 8000`
- Dependencies: `requirements.txt` (includes hydrosheaf core)

**Frontend:** React + Vite (TypeScript/JavaScript)
- Location: `web/frontend/`
- Run: `npm run dev`
- Expected backend: `http://localhost:8000`
- Converts mg/L ↔ mmol/L via adapters

**WebSocket Support:** Real-time analysis monitoring (defined but implementation details TBD)

---

## 20. Dependencies & Installation

### Core Package
```bash
pip install .
```

### With Plotting
```bash
pip install .[plot]
```
Adds: seaborn, pyvista, vtk

### With PHREEQC
```bash
pip install .[phreeqc]
```
Adds: phreeqpy

### Everything
```bash
pip install .[all]
```
Adds: flopy, seaborn, phreeqpy, pyvista, vtk

### Optional: Bayesian MCMC
For uncertainty quantification with MCMC:
```bash
pip install pymc>=5.0 arviz>=0.15
```

---

**End of Reference Manual**
