# O4 manuscript package

Objective 4 (Chapter 1) benchmark-synthesis manuscript for **Computers &
Geosciences**. Sibling paper to `O3` (Objective 3) and companion to `M2.3`
(the HydroSheaf framework paper). O4 harmonizes three already-locked,
independently designed component benchmarks -- `M6` (field-transfer and
robustness under data scarcity), `M7` (conditional integration and
non-identifiability, M7.3-M7.6), `M8` (calibration is not identification:
parameter recovery, uncertainty calibration, optimal measurement design) --
under one taxonomy of whether internally-generated confidence signals track
externally-verifiable truth under stress. It performs no PHREEQC, MODFLOW/
MODPATH, PEST-GLM, or Bayesian active-learning re-runs; see `DECISIONS.md`
D2-D3.

## Working title

**Fit quality is not identifiability: a harmonized robustness, integration-
value and calibration audit of groundwater inference under data limitation**

## Layout

```
O4/
  DECISIONS.md              locked decisions and claim ledger (incl. staleness verification, D3)
  Outline.md                positioning, argument spine, word budget
  proposal.normalized.json  objectives, inputs, risks
  analysis/
    python/                 computation authority; reads M6/M7/M8 locked
                             result files only, writes read-only exports
      _common.py
      derive_taxonomy.py
      derive_robustness_gradient.py
      derive_integration_value.py
      derive_calibration_gap.py
      derive_benchmark_scale.py
    r/                      figure authority; consumes exports only
      _theme.R
      fig01..fig06
      make_all_figures.R
  manuscript/
    Manuscript-O4.md         assembled manuscript
    supplementary/           Supplementary Methods and Information
    artifacts/
      data/                  tidy CSV exports (the evidence interface)
      figures/                PDF and PNG at 600 dpi
      TAB-*.md                generated tables
    review/                  internal adversarial reviewer report
    LITERATURE.bib
```

## Reproduction

From the repository root:

```
.venv/Scripts/python.exe O4/analysis/python/derive_taxonomy.py
.venv/Scripts/python.exe O4/analysis/python/derive_robustness_gradient.py
.venv/Scripts/python.exe O4/analysis/python/derive_integration_value.py
.venv/Scripts/python.exe O4/analysis/python/derive_calibration_gap.py
.venv/Scripts/python.exe O4/analysis/python/derive_benchmark_scale.py
"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" O4/analysis/r/make_all_figures.R
```

Python reads only files already committed or generated under
`M6/m6_field_transfer_benchmark/results/`,
`M7/m7_nonuniqueness_benchmark/results/` and
`M8/m8_calibration_benchmark/results/` (plus `M8`'s own registered manuscript
tables) and writes CSV exports under `manuscript/artifacts/data/`. R consumes
those exports and recomputes no reported statistic.

## Staleness verification

Before any number was written into the manuscript, the import graph of every
`M6`/`M7`/`M8` script producing a cited result was traced against three
commits made after those results were locked (`2bd5db0`, `8718d66`,
`2d4b8af`, 2026-08-01/02), and three headline numbers were independently
recomputed from raw locked CSVs rather than taken on trust from prose
summaries. Full detail, including the one number chain found stale and
excluded (`M7`'s public-pipeline system-acceptance run), is in
`DECISIONS.md` D3.
