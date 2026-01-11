# Hydrosheaf Input Interconnection Map (Primary vs Derived)

Definitions
- Primary inputs: raw data supplied by the user (CSV/JSON rows, config, and optional files).
- Derived inputs: values computed by a module and then used as inputs by downstream modules.

Global flow (high level): raw samples -> normalization -> ion vectors -> edge inference (optional) -> edge fitting (transport + reactions) -> optional modules (PHREEQC, isotopes, nitrate source, temporal, vadose, UQ, reactive transport).

Normalization (`hydrosheaf/data/schema.py`)
- Primary inputs: raw sample rows (strings or numbers), `config.detection_limit_policy`, `config.missing_policy`, `config.ion_order`
- Derived outputs: numeric `normalized` sample dict and ion vector `x` in model order
- Downstream use: core fit (`fit_edge`, `fit_network`), EC/TDS penalty, PHREEQC runner, QC

Graph/edge inference (`hydrosheaf/graph/build.py`, `hydrosheaf/graph3d/build_3d.py`)
- Primary inputs: `site_id`, `lat`, `lon`, `head_meas`, optional `dtw`, `elevation`, `screen_depth`, `well_depth`
- Derived outputs: edges with `edge_id`, distances/gradients, `edge_confidence`
- Downstream use: `fit_network` uses edges and attrs for per-edge fitting

Reaction dictionary (`hydrosheaf/models/reactions.py`)
- Primary inputs: `config.active_minerals`, `config.exchange_enabled`, mineral library
- Derived outputs: reaction matrix, labels, mineral mask
- Downstream use: reaction fit in `fit_edge`, PHREEQC constraints mapping

Transport fit (`hydrosheaf/models/transport.py`)
- Primary inputs: `x_u`, `x_v`, endmembers, weights
- Derived outputs: `gamma` or `f` and transport residual
- Downstream use: reaction fit uses residual and transport selection

Reaction fit (`hydrosheaf/models/reactions.py` via `fit_edge`)
- Primary inputs: transport residual, reaction matrix, bounds, lambda
- Derived outputs: reaction extents `z_extents`, reaction residual
- Downstream use: edge summaries, nitrate source, reactive transport validation

PHREEQC runner (`hydrosheaf/phreeqc/runner.py`)
- Primary inputs: `pH`, `temp_c`, Ca/Mg/Na/K/Cl/SO4/NO3/F/HCO3
- Derived outputs: ionic strength, saturation indices, charge balance, `phreeqc_ok`
- Downstream use: `hydrosheaf/phreeqc/constraints.py` converts SI to bounds for `fit_edge`

Isotope penalties (`hydrosheaf/isotopes.py`, used in `fit_edge`)
- Primary inputs: `18O`, `2H` from normalized `obs_u` and `obs_v`
- Derived outputs: d-excess, evaporation index, isotope penalty
- Downstream use: adds penalty to objective in `fit_edge`

Gibbs metrics (`hydrosheaf/models/gibbs.py`)
- Primary inputs: Na/Ca/Cl/HCO3, optional TDS
- Derived outputs: Gibbs ratios and evaporation penalty
- Downstream use: adds penalty to objective when enabled

EC/TDS penalty (`hydrosheaf/models/ec_tds.py`)
- Primary inputs: modeled ion vector, observed `EC` and `TDS` (if present)
- Derived outputs: EC/TDS penalty
- Downstream use: adds penalty to objective in `fit_edge`

Nitrate source v2 (`hydrosheaf/nitrate_source_v2.py`)
- Primary inputs: raw node chemistry (NO3, Cl, K, PO4, Fe, HCO3), optional `d_excess`, optional nitrate isotopes
- Derived outputs: posterior probabilities and evidence per node
- Downstream use: attached to edges in `fit_network` for reporting

Temporal module (`hydrosheaf/temporal/time_series.py`, `hydrosheaf/temporal/temporal_edge_fit.py`)
- Primary inputs: time-series CSV with `timestamp`, node id, and ions
- Derived outputs: residence time per edge, time-averaged gamma/f/extents
- Downstream use: residence time overrides in `fit_network_pipeline` and attached summaries

Vadose module (`hydrosheaf/vadose/run.py`)
- Primary inputs: forcing CSV (P, ET0), vadose profile JSON/YAML, links CSV
- Derived outputs: travel-time priors (mean/std/quantiles)
- Downstream use: converted to `PhysicsPrior` for `fit_network_with_priors`

Physics priors (`hydrosheaf/physics/priors.py`)
- Primary inputs: CSV/JSON with `u`, `v`, `p_uv`, `tt_*`
- Derived outputs: edge attrs (`physics_tau_*`, `edge_confidence`)
- Downstream use: edges updated before `fit_network`

Uncertainty (`hydrosheaf/uncertainty/*`)
- Primary inputs: edge fit inputs/outputs + config settings
- Derived outputs: confidence intervals for gamma/f/extents
- Downstream use: attached to `EdgeResult`

Reactive transport validation (`hydrosheaf/reactive_transport/validation.py`)
- Primary inputs: inverse results (extents), residence time, sample chemistry
- Derived outputs: RMSE/NSE and thermo consistency flags
- Downstream use: attached to `EdgeResult`

Auto-disable linkage (`hydrosheaf/api.py`)
- If required primary inputs are missing, PHREEQC, isotopes, and nitrate source are disabled automatically to avoid invoking modules without data.
