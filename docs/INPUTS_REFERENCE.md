# Hydrosheaf Inputs Reference

This guide lists the inputs each module expects so users can quickly prepare data.
All concentrations are assumed to be in mmol/L unless noted otherwise.

## Global Concepts

- **Sample identifiers**: `site_id` is used to connect samples to edges.
  `sample_id` is accepted and often required for file I/O, but `site_id` is the
  canonical node key during inference.
- **Ion order**: Uses `config.ion_order` (default: Ca, Mg, Na, HCO3, Cl, SO4, NO3, F, Fe, PO4).
- **Units**: Default is mmol/L. The CLI can convert from mg/L or meq/L. Web requests are in mg/L and convert via adapters.

- **Missing data**:
  - If `missing_policy="skip"`, missing ions skip the sample.
  - If `missing_policy="impute_zero"`, missing ions are treated as 0.
  - Detection limits like "<0.01" use `detection_limit_policy`.

## Core Network Fit (fit_network / fit_network_pipeline)

Minimum inputs per sample (for core chemistry fit):
- `site_id`
- Every ion in `config.ion_order`

Optional but common:
- `pH`, `EC`, `TDS`, `temp_c`, `K`, `lat`, `lon`, `elevation`


Edge inputs:
- `(u, v)` tuples or edges with `edge_id`, `u`, `v`.
- Optional edge attrs for probabilistic edges: `distance_km`, `delta_h`, `edge_confidence`.

## Edge Inference (infer_edges / probabilistic graph)

For probabilistic inference:
- `site_id`
- `lat`, `lon`
- `head_meas` (or `edge_head_key`)
- Optional: `dtw`, `elevation`, `aquifer_unit`, `screen_depth`, `well_depth`

For 3D inference:
- `screen_depth` (or `z_mASL` depending on `config.z_coordinate_key`)
- Optional: `aquifer_unit`, `well_depth`, `layer` fields from layer YAML


## PHREEQC (thermodynamic constraints)

Requires (per sample):
- `pH`
- Major ions used in PHREEQC input: Ca, Mg, Na, K, Cl, SO4, NO3, F, HCO3
- Optional: `temp_c`

Config flags:
- `config.phreeqc_enabled=True`
- `config.phreeqc_mode="phreeqpython"` or `"subprocess"`
- `config.phreeqc_database` (optional path)

If `pH` is missing or PHREEQC is unavailable, constraints are skipped for those samples.

## Isotopes (delta18O / delta2H)

Requires (per sample):
- `18O` and `2H` (configurable via `config.isotope_d18o_key` and `config.isotope_d2h_key`)

Config flags:
- `config.isotope_enabled=True`
- LMWL parameters: `config.lmwl_a`, `config.lmwl_b`, `config.lmwl_defined=True`

If isotope pairs are missing, isotope penalties are skipped.

## Nitrate Source (v2)

Minimum:
- `NO3` (for threshold check)

Common hydrochemistry (improves results):
- `Cl`, `K`, `PO4`, `Fe`, `HCO3`
- `d_excess` (if available)

Optional isotopes (dual-isotope mixing):
- `d15N`, `d18O_NO3` (configurable via `nitrate_isotope_n15_col` / `nitrate_isotope_o18_col`)
- Requires PyMC/ArviZ for Bayesian mixing; otherwise analytical probabilities are returned.


Config flag:
- `config.nitrate_source_enabled=True`

If NO3 is below threshold, inference is skipped and a reason is provided.

## Temporal Module

Time-series CSV required fields:
- `timestamp` (default format: YYYY-MM-DD)
- `site_id` or `node_id` or `sample_id` (the loader tries these in order)
- All ions in `config.ion_order`

Optional:
- `temp_c` or `temperature`, `pH`, `EC`, isotopes (`18O`, `2H`, `d18O`, `d2H`)

Config flags:
- `config.temporal_enabled=True` (CLI) or use `fit_temporal_edges` (API)
- `config.residence_time_method` and `config.residence_time_tracer`

If the time series cannot be aligned across edges, temporal fits are skipped for those edges.

## Vadose Module (unsaturated zone)

Inputs:
- Forcing CSV: `timestamp`, `P_mm_day`, `ET0_mm_day` (optional: `I_mm_day`)
- Profile JSON/YAML:
  - `profile_id`, `depth_m`
  - `layers`: per-layer parameters such as `theta_r`, `theta_s`, `alpha`, `n`, `ks`
  - Optionally use `texture` to auto-assign parameters
- Links CSV: `u`, `v`, `u_depth_m`, `p_uv`
- Optional water table CSV: `timestamp`, `water_table_depth_m`

Outputs:
- Physics priors with travel-time statistics, suitable for `fit_network_with_priors`.
- Optional nitrate breakthrough curves when `vadose-no3-loading` is supplied.


## Physics Priors (external)

CSV/JSON fields:
- `u`, `v`
- Optional: `p_uv`, `tt_mean_days`, `tt_std_days`, `tt_p10_days`, `tt_p90_days`

Apply with:
- `fit_network_with_priors(..., physics_priors=..., physics_priors_mode="override|merge|only")`

## Uncertainty Quantification

Config flags:
- `config.uncertainty_method` in `none`, `bootstrap`, `bayesian`, `monte_carlo`

No extra input columns required, but uncertainty settings control runtime. Bayesian mode requires PyMC.


## Reactive Transport Validation

Config flag:
- `config.reactive_transport_validation=True`

Uses the same chemistry inputs as core fit; requires PHREEQC kinetic support if enabled. `rt_residence_time` defaults to temporal results if available.


## Convenience API Pipeline

If you use `fit_network_pipeline(...)` with `auto_disable_missing=True` (default):
- PHREEQC, isotope penalties, and nitrate source are automatically disabled if the dataset does not contain the required inputs.
- Temporal results are only computed if `temporal_nodes` is provided.
- Vadose priors are only computed if you call `build_vadose_priors(...)`.
