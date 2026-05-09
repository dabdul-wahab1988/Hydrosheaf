# Hydrosheaf User Guide

## 1. Quick Start

### Basic Run

Run the model on your prepared dataset using default minerals:

```powershell
python -m hydrosheaf.cli --samples data_prepared.csv --output results.json --infer-edges --phreeqc-enabled
```

### Advanced Run (Expert Mode)

Enable isotopes, auto-calibration, and interpretation report:

```powershell
python -m hydrosheaf.cli \
  --samples data_prepared.csv \
  --output results_final.json \
  --infer-edges \
  --phreeqc-enabled \
  --isotope-enabled \
  --auto-lmwl \
  --iso-consistency-weight 2.0 \
  --interpret
```

---

## 2. CLI Options Reference

### Data Inputs

* `--samples`: Path to CSV file containing ion concentrations (Site ID, Ca, Mg, Na... etc).
* `--unit`: Input units. Options: `mmol/L` (default), `mg/L`, `meq/L`.

### Geochemical Settings

* `--minerals`: Comma-separated list of allowed minerals.
  * *Example*: `"calcite,dolomite,gypsum,halite,albite,pyrite_oxidation_aerobic"`
* `--list-minerals`: Print all valid minerals in the library and exit.
* `--phreeqc-enabled`: Check saturation indices to constrain dissolution/precipitation.

### Isotope Settings

* `--isotope-enabled`: Turn on isotope physics.
* `--auto-lmwl`: Automatically fit the Local Meteoric Water Line from your data.
* `--iso-consistency-weight`: Weight for the penalty connecting isotope enrichment to Chloride increase. Set to `2.0` for strict physics.

### Nitrate Source Discrimination

* `--nitrate-source-enabled`: Enable the module to distinguish between manure and fertilizer sources of nitrate.
* `--nitrate-source-min-conc`: Minimum nitrate concentration (mg/L) required to attempt source discrimination (default: 10 mg/L). Samples below this threshold are classified as background and skipped.
* Uses geochemical ratios (NO3/Cl, NO3/K) and robust statistics to assign a probability of manure contamination.

### Temporal Dynamics

* `--temporal-data <file.csv>`: Time-series chemistry CSV (timestamp + ions).
* `--temporal-enabled`: Turn on temporal edge fitting and residence time estimation.
* `--residence-method`: `cross_correlation`, `gradient`, `bayesian_lag`, `ttd`, or `tracer_decay`.
* `--residence-tracer`: Conservative tracer for time-lag methods (e.g., `Cl`, `18O`, `2H`, `auto`).
* `--temporal-output <file.json>`: Optional file with temporal summaries.

### Uncertainty Quantification

* `--uncertainty`: `none`, `bootstrap`, `bayesian`, or `monte_carlo`.
* `--bootstrap-n`, `--bootstrap-ci`: Resampling settings for bootstrap CIs.
* `--bayesian-samples`, `--bayesian-chains`, `--bayesian-accept`: MCMC settings (PyMC required).
* `--monte-carlo-n`, `--input-uncertainty`: Monte Carlo settings.

### Vadose Zone & Physics Priors

* `--vadose-enabled`: Run a 1D Richards column and generate travel-time priors.
* `--vadose-forcing`, `--vadose-profile`, `--vadose-links`: Forcing, profile, and link files.
* `--vadose-no3-loading`: Optional nitrate loading CSV for breakthrough curves.
* `--physics-priors`: Apply external physics priors CSV/JSON.
* `--modpath-endpoints`: Build priors from MODPATH endpoints (requires FloPy).

### Reactive Transport Validation

* `--validate-forward`: Run kinetic forward validation of inverse results.
* `--rt-simulator`: `phreeqc_kinetic` or `mt3dms`.
* `--rt-time-steps`, `--rt-residence-time`: Control kinetic runtime.

### 3D Flow Networks

* `--3d`: Enable 3D edge inference.
* `--z-key`: Column name for vertical coordinate (e.g., `screen_depth`).
* `--layer-file`: Layer definition YAML for multi-aquifer systems.
* `--anisotropy`: Vertical anisotropy factor αᵥ.

### Interpretation & Output

* `--interpret`: Print a text-based geochemical report after the run finishes.
* `--report-only <file.json>`: Generate a report from an existing results file without re-running the model.

### Plotting

* `--plot-ilr`: Generate an Isometric Log-Ratio (ILR) plot for hydrochemical facies analysis (4-panel diagram).
* `--plot-gibbs`: Generate a Gibbs diagram (TDS vs Na/(Na+Ca) and Cl/(Cl+HCO3)) to identify precipitation, rock, or evaporation dominance.
* `--plot-ttd`: Plot Transit Time Distribution (TTD) kernels (Gamma/Dispersion) for temporal edges.
* `--plot-breakthrough`: Plot input vs. output breakthrough curves (convolution) to visualize signal propagation.
* `--plot-uncertainty`: Plot Bayesian posterior ridge distributions for reaction extents and transport parameters.
* `--plot-rt-validation`: Plot modeled vs. observed concentrations (RMSE/NSE) for reactive transport validation.
* `--plot-output <path>`: Prefix/path for plot files (default: derived from output filename).


---

## 3. Interpreting the Report

The interpretation report provides high-level insights:

### Transport Physics

* **Evap-Dominant**: Water is concentrating via evaporation.
* **Mix-Dominant**: Water is mixing with a saline end-member (e.g., seawater).

### Isotope Consistency

* **Enrichment Slope**: ~3.5-6.0 indicates strong evaporation. ~8.0 indicates rain/recharge.
* **Consistency Penalty**: High values warn that chemical salinity is increasing faster than physical evaporation can explain (Process: **Rock Interaction**).

### Redox State

* **Oxic**: High Nitrate, Low Iron. Pyrite burns to Sulfate.
* **Reducing**: Low Nitrate, High Iron. Pyrite is stable; Aerobic oxidation is disabled.

### Nitrate Sources

* **Manure-Dominated (`p_manure > 0.5`)**: High NO3/Cl and NO3/K ratios, often coupled with high denitrification signals.
* **Fertilizer-Dominated (`p_manure < 0.5`)**: Very high NO3/Cl ratios, low K, lack of organic signals.
* **Evidence**: The system provides the top reasoning for its classification (e.g., `NO3/K_high_manure`, `denitrif_strong`).
