# Hydrosheaf Calibration Guide

This guide provides a comprehensive reference for configuring, running, and troubleshooting model calibration in **Hydrosheaf**.

---

## 1. Quick Start

Calibration in Hydrosheaf is run via the console command `hydrosheaf-cal`:

```bash
# Generate a starter configuration file of type 'composite'
hydrosheaf-cal --write-template composite > calibration_config.yaml

# Run a dry-run validation to check inputs, parameters, and observations
hydrosheaf-cal calibration_config.yaml --dry-run

# Run the actual optimization
hydrosheaf-cal calibration_config.yaml
```

---

## 2. Configuration Schema Reference

The calibration settings are defined in a standard YAML configuration file.

### Top-Level Keys

| Key | Type | Description |
|---|---|---|
| `calibration.type` | String | Problem type: `kinetic`, `transport`, `vadose`, `nitrate`, `age_temporal`, `topology`, or `composite`. |
| `calibration.settings` | Dict | Engine settings (optimizer, iterations, parallel workers, loss, etc.). |
| `calibration.model` | Dict | Adapter-specific settings (used when running a single problem type). |
| `calibration.sub_models` | List | A list of sub-model configuration dictionaries (used for `composite` calibration). |
| `calibration.parameters` | List | Global or generic adjustable parameter definitions. |

### Settings Parameters (`calibration.settings`)

| Key | Type | Default | Description |
|---|---|---|---|
| `engine` | String | `"internal"` | Optimization engine: `"internal"` (Python GLM solver) or `"pestpp"` / `"pestpp-ies"`. |
| `max_iterations` | Integer | `50` | Maximum number of model evaluations (`max_nfev`). |
| `n_workers` | Integer | `1` | Number of parallel worker processes / threads. |
| `output_dir` | String | `"calibration_results"` | Directory where results and `results.json` will be saved. |
| `loss` | String | `"linear"` | Robust loss function for `internal` engine: `linear`, `huber`, `soft_l1`, `cauchy`. |

---

## 3. Which Calibration Adapter Should I Use?

| Scenario | Calibration Type | Key Adjustable Parameters | Observations |
|---|---|---|---|
| Fitting mineral reaction rates/surface areas | `kinetic` | Mineral rates (`k`), surface areas (`A`), extents (`xi`) | Mineral solute concentrations (mg/L or mmol/L) |
| Fitting solute travel time & aquifer properties | `transport` | Dispersivity, decay, effective porosity | Concentration time-series at observation points |
| Estimating soil parameters, recharge & ET crop coefficient | `vadose` | Conductivity (`Ks`), VG `alpha`/`n`, crop coefficient `Kc`, preferential flow | Water content (`theta`) or suction pressure (`psi`) at depths |
| Source attribution and mixing weights | `nitrate` | Priors, evidence weights, QC thresholds | Manure/Fertilizer post-attribution probabilities |
| Estimating travel time distributions (TTD) | `age_temporal` | Mean travel time (`tau`), TTD CV (`ttd_cv`), decay | Solute concentration/groundwater age time-series |
| Building directed network flow topology | `topology` | Edge probabilities/presence | Edge presence, hydraulic heads, layers |

---

## 4. Adapter Configuration Examples

### A. Kinetic Calibration
```yaml
calibration:
  type: kinetic
  settings:
    engine: internal
    max_iterations: 30
    n_workers: 4
  model:
    minerals: ["calcite", "dolomite"]
    fit_parameters: ["calcite:k:global", "dolomite:A:per_edge"]
    observations_file: "data/kinetic_observations.csv"
```

### B. Vadose Zone Calibration
```yaml
calibration:
  type: vadose
  settings:
    engine: internal
  model:
    profile_file: "configs/vadose_profile.json"
    forcing_file: "data/forcing.csv"
    observations_file: "data/vadose_observations.csv"
    layers_to_fit: [0]
    fit_parameters: ["ks_L0", "alpha_L0", "n_L0", "kc", "preferential_flow_fraction"]
```

### C. Flow Network Topology Search (Outer-Loop)
```yaml
calibration:
  type: topology
  settings:
    engine: internal
  model:
    outer_loop: true
    samples_file: "data/network_samples.csv"
    edges_file: "data/candidate_edges.csv"
    observations_file: "data/topology_observations.csv"
    max_iterations: 10
    max_neighbors: 2
    head_key: "hydraulic_head"
```

---

## 5. Output Interpretation

The optimization run outputs a summary directly to the console and saves a detailed `results.json` file inside the configured `output_dir`.

### Key Metrics
1. **Objective Function Value ($\Phi$)**: Sum of weighted squared residuals. A lower value indicates a better fit.
2. **Covariance Diagnostics**:
   - `inverse`: Matrix inverted successfully (well-conditioned parameter space).
   - `svd_pseudoinverse`: Jacobian was rank-deficient; SVD pseudoinverse was used to compute parameter uncertainties.
3. **AIC / BIC**: Akaike and Bayesian Information Criteria. Useful for comparing models with different numbers of parameters (especially for topology search).

---

## 6. Troubleshooting & Common Failures

### 1. "NameError / ModuleNotFoundError" during run
- **Cause**: Missing dependencies or incorrect environment active.
- **Solution**: Run the CLI using the package manager runner (e.g. `uv run hydrosheaf-cal config.yaml`).

### 2. SVD Pseudoinverse flag active
- **Cause**: Some parameters are non-identifiable or highly correlated.
- **Solution**: Reduce the number of parameters to fit, apply parameter tying, or add prior constraints with small standard deviations (`prior_sigma`).

### 3. Outer-loop topology search selecting 0 edges
- **Cause**: Constraints (like acyclicity or head direction) are too strict for the input candidate edge set.
- **Solution**: Review the hydraulic heads in your samples file and ensure flow direction from higher head to lower head is physically possible for candidate edges.
