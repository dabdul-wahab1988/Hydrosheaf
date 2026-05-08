# M2 Hydrosheaf Benchmark Package

This directory is the M2 benchmark package for generating and auditing the
tables and figures specified in the M2 outline, table guide, and synthetic-data
guide.

## Rebuild

From the repository root:

```powershell
python M2\m2_benchmark\scripts\run_m2_benchmark.py
```

The default run uses the locked `config/ground_truth.yaml` setting of 100 Monte
Carlo realisations.

For a quick smoke test:

```powershell
python M2\m2_benchmark\scripts\run_m2_benchmark.py --realisations 3
```

## Output Structure

- `config/ground_truth.yaml` - locked benchmark design.
- `data/` - generated ground-truth node/edge tables and synthetic realisations.
- `results/` - quantitative recovery and validation CSVs.
- `tables/` - main-text and supplementary table CSVs.
- `figures/` - main-text and supplementary figure PNGs.
- `external/` - external validation workspace for public age data, MODPATH,
  PHREEQC, and the Northern Ghana data-limited field-hydrochemistry
  demonstration.
- `docs/m2_results_summary.md` - concise result summary for manuscript drafting.
- `docs/external_validation_plan.md` - required external validation tasks and
  source mapping.
- `MANIFEST.md` - generated file inventory.

## Main Required Outputs

- Table 1: `tables/table1_module_architecture.csv`
- Table 2: `tables/table2_input_fields.csv`
- Table 3: `tables/table3_residence_time_options.csv`
- Table 4: `tables/table4_validation_design_and_results.csv`
- Table 5: `tables/table5_method_comparison.csv`
- Figure 1: `figures/figure1_architecture_workflow.png`
- Figure 2: `figures/figure2_process_network.png`
- Figure 3: `figures/figure3_synthetic_benchmark_recovery.png`
- Figure 4: `figures/figure4_residence_time_network_update.png`
- Figure 5: `figures/figure5_sensitivity_uncertainty.png`

## Repository Policy For Raw Inputs

Bulky public/raw inputs are not committed to Git. This includes downloaded USGS
archives, extracted public data tables, MODPATH endpoint/pathline files, DGMETA
workbooks, and superseded local pilot input files. The committed package keeps
the scripts, source manifests, README files, generated result tables, figures,
and interpretation guardrails.

If a source file is needed for a full rebuild, use the relevant runner script
or the source DOI/URL recorded in the result manifest for that tier.

## External Validation Status

The committed result artefacts include the synthetic benchmark and external
validation tiers:

- public tracer-age agreement from USGS DOI `10.5066/P9W7T0DN`;
- MODFLOW/MODPATH topology comparison from USGS DOI `10.5066/F7J102FK`;
- live PHREEQC forward validation from PHREEQC v3 examples/databases, DOI
  `10.3133/tm6A43`;
- a data-limited field-hydrochemistry demonstration using the corrected
  Northern Ghana workbook at `data/NorthenGhana/NorthernGhana.xlsx`.

These are summarized in `tables/table4_validation_design_and_results.csv` and
`docs/m2_results_summary.md`. Retain the guardrails in those files when writing
the manuscript.

## E1 Public Age Validation

Run the E1 validation from the repository root:

```powershell
python M2\m2_benchmark\scripts\run_e1_usgs_age_validation.py
```

If the USGS ScienceBase endpoint is unavailable but the seven release tables have
already been downloaded into `external/usgs_age/input`, use:

```powershell
python M2\m2_benchmark\scripts\run_e1_usgs_age_validation.py --no-download
```

The runner writes a failure report instead of Figure 4A when the actual source
tables cannot be read.
