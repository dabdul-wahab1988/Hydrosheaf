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
| `engine` | String | `"internal"` | Optimization engine: `"internal"` (PESTGLM), `"pestpp-glm"`, `"pestpp-ies"`, `"pestpp-sen"`, `"pestpp-swp"`, `"pestpp-mou"`, `"pestpp-opt"`, `"pestpp-da"`. |
| `max_iterations` | Integer | `50` | Maximum number of model evaluations (`max_nfev`, `noptmax` in PEST++). |
| `n_workers` | Integer | `1` | Number of parallel workers. For PEST++, activates manager/agent PANTHER architecture when > 1. |
| `output_dir` | String | `"calibration_results"` | Directory where results and `results.json` will be saved. |
| `loss` | String | `"linear"` | Robust loss function for `internal` engine: `linear`, `huber`, `soft_l1`, `cauchy`. |
| `work_dir` | String | `"pest_workspace"` | Workspace for PEST++ input/output files. |
| `pestpp_options` | Dict | `{}` | Generic `++` options passed to all PEST++ engines. |
| `ies` | Dict | `{}` | PESTPP-IES specific options (see Section 5). |
| `sen` | Dict | `{}` | PESTPP-SEN specific options. |
| `swp` | Dict | `{}` | PESTPP-SWP specific options. |
| `opt` | Dict | `{}` | PESTPP-OPT / PESTPP-MOU specific options. |
| `da` | Dict | `{}` | PESTPP-DA specific options. |

### Parameter Definitions (`calibration.parameters`)

| Key | Type | Required | Description |
|---|---|---|---|
| `name` | String | Yes | Parameter name. Must match adapter expectations or be unique. |
| `initial` | Float | Yes | Starting value. |
| `bounds` | `[lower, upper]` | Yes | Lower and upper bounds. |
| `log` | Boolean | No | Log10-transform the parameter. Default `false`. |
| `prior_mean` | Float | No | Prior mean for regularization. |
| `prior_sigma` | Float | No | Prior standard deviation for regularization. |
| `fixed` | Boolean | No | If `true`, the parameter is held constant during calibration. |
| `tied_to` | String | No | If set, this parameter is tied to another parameter's value. |

Example with `fixed` and `tied_to`:

```yaml
calibration:
  parameters:
    - name: dispersivity
      initial: 1.0
      bounds: [0.1, 10.0]
      log: true
    - name: decay
      initial: 0.0
      bounds: [0.0, 1.0]
      fixed: true               # do not adjust during calibration
    - name: velocity
      initial: 0.1
      bounds: [0.01, 1.0]
      tied_to: dispersivity     # tied to another parameter's value
```

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
    fit_parameters: ["calcite:k:global", "dolomite:k:global"]
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

## 5. PEST++ Engine Specifics

When `engine` is set to a PEST++ tool, Hydrosheaf automatically generates `.pst`, `.tpl`, `.ins`, ensemble files, and invokes the native binary.

### PESTPP-GLM (Parameter Estimation + FOSM)

Suitable for: single-point optimization with sensitivity/uncertainty analysis.

```yaml
calibration:
  settings:
    engine: pestpp-glm
    max_iterations: 20
    pestpp_options:
      lambdas: "1.0,10.0"
      uncertainty: true
```

### PESTPP-IES (Iterative Ensemble Smoother)

Suitable for: ensemble-based history matching, uncertainty reduction.

```yaml
calibration:
  settings:
    engine: pestpp-ies
    max_iterations: 10
    pestpp_options:
      max_run_fail: 3
      panther_agent_restart_on_error: true
    ies:
      ies_num_reals: 200
      ies_num_threads: 8
      ies_lambda_mults: "0.1,1.0,10.0"
      ies_subset_size: 20
      ies_accept_phi_fac: 1.01
      par_sigma_range: 4.0
      # Optional: auto-generate covariance from bounds
      parcov: true
      # Optional: localization matrix (auto-generates identity)
      ies_localizer: true
      # Optional: restart from previous run
      # ies_restart_parameter_ensemble: "restart.0.par.csv"
      # ies_restart_observation_ensemble: "restart.0.obs.csv"
      # Optional: forecast observations to summarize
      forecasts: "gw_age,mean_conc"
```

**Key IES Options:**

| Option | Type | Default | Description |
|---|---|---|---|
| `ies_num_reals` | int | 100 | Number of realizations in the ensemble. |
| `ies_num_threads` | int | -1 | Number of threads for matrix ops. `-1` = all cores. |
| `ies_lambda_mults` | str | `"0.1,1.0,10.0"` | Marquardt lambda multipliers to test. |
| `ies_subset_size` | int | 4 | Size of realization subset. |
| `ies_accept_phi_fac` | float | 1.05 | Phi acceptance factor. Lower = stricter. |
| `par_sigma_range` | float | 4.0 | Number of standard deviations spanned by bounds. |
| `parcov` | bool/str | — | `true` to auto-generate from bounds, or path to `.csv`. |
| `obscov` | bool/str | — | `true` to auto-generate from weights, or path to `.csv`. |
| `ies_localizer` | bool/str | — | `true` for identity localizer, or path to `.csv`. |

### PESTPP-SEN (Sensitivity Analysis)

```yaml
calibration:
  settings:
    engine: pestpp-sen
    sen:
      gsa_method: morris           # or sobol
      gsa_morris_r: 100
      gsa_morris_delta: 0.1
      # Sobol alternatives:
      # gsa_sobol_samples: 1000
      # gsa_sobol_par_dist: unif     # or norm
```

Hydrosheaf still accepts older aliases such as `sen_method`,
`sen_num_samples`, and `sen_morris_delta`, but writes the PEST++ control file
with the real `gsa_*` option names.

### PESTPP-SWP (Parameter Sweep)

```yaml
calibration:
  settings:
    engine: pestpp-swp
    swp:
      sweep_n_runs: 100
      sweep_sampler: latin_hypercube  # or grid, random
      sweep_forgive: true
```

### PESTPP-MOU (Multi-Objective Optimization)

```yaml
calibration:
  settings:
    engine: pestpp-mou
    max_iterations: 50
    opt:
      mou_population_size: 100
      mou_generator: de               # or pso
      mou_objectives: "cost,rmse"
      opt_constraint_groups: "less_than_limit,greater_than_target"
      opt_risk: 0.95
```

PESTPP-MOU uses `max_iterations`/`NOPTMAX` for generations. Objectives and
constraints must be observations or prior-information equations whose groups
start with `less_`, `less_than`, `l_`, `greater_`, `greater_than`, or `g_` so
PEST++ can infer the minimization/maximization direction.

### PESTPP-OPT (Optimization Under Uncertainty)

```yaml
calibration:
  settings:
    engine: pestpp-opt
    opt:
      opt_risk: 0.95
      opt_dec_var_groups: "pargp"
      opt_constraint_groups: "less_than_limit"
      opt_direction: min
```

### PESTPP-DA (Data Assimilation)

```yaml
calibration:
  settings:
    engine: pestpp-da
    max_iterations: 5
    da:
      da_hotstart_cycle: 0
      # da_stop_cycle: 5
      # da_noptmax_schedule: "da_noptmax_schedule.dat"
      # da_parameter_ensemble: "initial.par.csv"
```

PESTPP-DA cycle count is not controlled by a `da_num_cycles` option. Cycles are
defined by cycle-aware control-file sections and optional DA cycle tables; the
number of ensemble-smoother iterations per cycle comes from `max_iterations`
or `da_noptmax_schedule`.

---

## 6. Parallel Execution with `n_workers`

For **PEST++** engines, setting `n_workers > 1` activates the PANTHER manager/agent architecture:

- A **manager** process is launched: `pestpp-<engine> case.pst /h :<port>`
- `n_workers` **agent** processes connect: `pestpp-<engine> case.pst /h localhost:<port>`
- Agents execute model runs in parallel and report back to the manager.

The runner automatically finds a free ephemeral port, handles agent lifecycle, and cleans up processes on completion or failure (Windows-safe `terminate` → `wait(5s)` → `kill` cascade).

On Windows, the PANTHER manager can occasionally keep running after writing
valid outputs. Set the Hydrosheaf-only runtime option `panther_timeout_secs` to
make the runner terminate the manager after the timeout and parse outputs if
the expected files were already written. This option is not written to the
PEST++ control file.

### Example

```yaml
calibration:
  settings:
    engine: pestpp-ies
    n_workers: 4
    max_iterations: 20
    pestpp_options:
      panther_timeout_secs: 600
```

On a 4-core machine, this launches 1 manager + 4 agents, giving 4 concurrent model evaluations.

---

## 7. Pipeline Coupling (`--load-calibration`)

After calibration, you can load optimized parameters back into the core `Config`:

```python
from hydrosheaf.config import Config
config = Config()
config.load_from_calibration_json("calibration_results/results.json")
```

**Supported mappings:**
- **Direct matches:** If the JSON parameter name matches a `Config` attribute, it is cast and set directly.
- **Aliases:** Common names like `dispersivity` → `dispersivity_m`, `decay` → `denitrification_k_1_day`, `porosity` → `aquifer_porosity`, etc.
- **Ensemble support:** If the JSON contains `posterior_parameters` (IES) but no `optimal_parameters`, the method automatically computes posterior means and maps those.
- **Kinetic parameters:** `log_k_mineral` names are recognized and mapped into `mineral_rate_constants`.

---

## 8. Output Interpretation

The optimization run outputs a summary directly to the console and saves a detailed `results.json` file inside the configured `output_dir`.

### Key Metrics
1. **Objective Function Value ($\Phi$)**: Sum of weighted squared residuals. A lower value indicates a better fit.
2. **Covariance Diagnostics**:
   - `inverse`: Matrix inverted successfully (well-conditioned parameter space).
   - `svd_pseudoinverse`: Jacobian was rank-deficient; SVD pseudoinverse was used to compute parameter uncertainties.
3. **AIC / BIC**: Akaike and Bayesian Information Criteria. Useful for comparing models with different numbers of parameters.

### IES-Specific Outputs
- `prior_parameters`: Ensemble of prior parameter realizations.
- `posterior_parameters`: Ensemble of posterior parameter realizations.
- `phi_history`: Mean Phi per iteration.
- `posterior_forecast_summaries`: Mean, std, min, max, median for each forecast observation.
- `parcov_path`, `obscov_path`, `localizer_path`: Paths to covariance/localizer files used.

---

## 9. Troubleshooting & Common Failures

### 1. "NameError / ModuleNotFoundError" during run
- **Cause**: Missing dependencies or incorrect environment active.
- **Solution**: Run the CLI using the package manager runner (e.g. `uv run hydrosheaf-cal config.yaml`).

### 2. SVD Pseudoinverse flag active
- **Cause**: Some parameters are non-identifiable or highly correlated.
- **Solution**: Reduce the number of parameters to fit, apply parameter tying, or add prior constraints with small standard deviations (`prior_sigma`).

### 3. Outer-loop topology search selecting 0 edges
- **Cause**: Constraints (like acyclicity or head direction) are too strict for the input candidate edge set.
- **Solution**: Review the hydraulic heads in your samples file and ensure flow direction from higher head to lower head is physically possible for candidate edges.

### 4. PEST++ agents not connecting (Windows)
- **Cause**: Windows firewall or port conflict.
- **Solution**: Ensure the ephemeral port is not blocked. The runner automatically scans for a free port, but if you see `ConnectionRefused`, try setting `n_workers: 1` to run serially, or manually specify a known-open port via `pestpp_options.panther_master_port`.

### 5. IES ensemble convergence issues
- **Cause**: `ies_num_reals` too low, or `ies_accept_phi_fac` too high.
- **Solution**: Increase `ies_num_reals` to 200–500. Lower `ies_accept_phi_fac` to 1.01. Add `parcov: true` to provide a structured prior covariance.

### 6. Missing `parcov` / `obscov` files
- **Cause**: User-provided path does not exist in `work_dir`.
- **Solution**: Use absolute paths, or set `parcov: true` / `obscov: true` to auto-generate diagonal matrices from current parameter bounds and observation weights.
