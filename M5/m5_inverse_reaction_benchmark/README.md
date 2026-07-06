# M5 Inverse Reaction Benchmark

This folder contains the complete M5 hydrogeochemical inverse-modelling
validation benchmark. It supports the defensible Hydrosheaf framing:

> a sparse linear inverse reaction model screened and stress-tested using
> PHREEQC thermodynamic, Hydrosheaf-Core sparse evidence-gate, calibrated
> support-selection, and forward-validation diagnostics.

It does not frame Hydrosheaf as a fully coupled nonlinear PHREEQC inverse
solver.

The locked benchmark includes a conventional PHREEQC `INVERSE_MODELING`
baseline. Strict 5% PHREEQC inverse models were feasible for 45.8% of the
clean synthetic pairs; a documented relaxed 20% uncertainty fallback increased
feasibility to 99.6% but produced high model multiplicity and lower
first-minimal equivalence-class recovery than the Hydrosheaf guarded model.

The current benchmark separates three Hydrosheaf-relevant layers:
`thermo_elastic_net` is the thermodynamic sparse baseline,
`hydrosheaf_guarded` is the primary calibrated guarded performance model, and
`hydrosheaf_core` is the Hydrosheaf-Core evidence-gated comparator/audit layer
that adds Gibbs, CAI/exchange, redox, EC/TDS, saturation-index, and optional
Ghana SiO2/Sr/water isotope evidence without requiring rare tracers.

The workflow also includes a formal data-tier experiment. `core` uses sparse
major-ion evidence only; `plus_lite` adds controlled synthetic SiO2, Sr, and
water-isotope diagnostics; `enhanced` adds controlled synthetic Br, DO, DOC,
sulphate-isotope, and nitrate-isotope diagnostics. These optional diagnostics
are a measurement-design experiment, not field-measured reaction truth.

Run the complete analysis, table, database, and figure workflow from the
repository root:

```powershell
.\.venv\Scripts\python.exe M5\m5_inverse_reaction_benchmark\scripts\run_m5_all.py
```

For display-only refreshes after the synthetic CSV files exist:

```powershell
.\.venv\Scripts\python.exe M5\m5_inverse_reaction_benchmark\scripts\run_m5_all.py --reuse-synthetic
```

For model-comparison refreshes that reuse existing PHREEQC ground truth but
regenerate the sparse-inversion outputs, tables, and figures:

```powershell
.\.venv\Scripts\python.exe M5\m5_inverse_reaction_benchmark\scripts\run_m5_all.py --reuse-phreeqc
```

The full run requires USGS PHREEQC 3.7.3 or environment variables
`PHREEQC_EXE` and `PHREEQC_DATABASE`. Outputs are written to `results/`,
`tables/`, `figures/`, `phreeqc_inputs/`, and `docs/`.

The workflow now writes a complete DuckDB evidence store at
`results/m5_results.duckdb` plus `results/m5_results_database_catalog.csv`.
The intended manuscript-grade display layer is the R Nature-style figure set in
`figures/r_publication/` and `figures/r_publication/supplementary/`, saved as
PNG, TIFF, and PDF. Python database-backed figures in `figures/publication/`
remain a reproducible fallback and quality-control display layer.

To regenerate only the R publication figures after the database/CSV outputs
exist, run from `M5/m5_inverse_reaction_benchmark`:

```powershell
& "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" r_figures\plot_m5_publication_figures.R
```

Use `--skip-r-figures` with `run_m5_all.py` only on machines without R.
