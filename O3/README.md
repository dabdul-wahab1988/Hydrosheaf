# O3 manuscript package

Objective 3 (Chapter 1) benchmark-synthesis manuscript for **Computers &
Geosciences**. Companion paper to `M2.3` (the HydroSheaf framework paper).
O3 harmonizes the three already-locked, independently designed component
benchmarks -- `M3` (age/residence time), `M4` (topology), `M5` (reaction) --
under one common identifiability taxonomy. It performs no PHREEQC, MODFLOW/
MODPATH or TracerLPM re-runs; see `DECISIONS.md` D2-D3.

## Working title

**HydroSheaf benchmarks: a harmonized identifiability comparison of directed
connectivity, residence-time, and reaction inference in data-limited
aquifers**

## Layout

```
O3/
  DECISIONS.md              locked decisions and claim ledger
  Outline.md                positioning, argument spine, word budget
  proposal.normalized.json  objectives, inputs, risks
  analysis/
    python/                 computation authority; reads M3/M4/M5 locked
                             result files only, writes read-only exports
      derive_taxonomy.py
      derive_headline_metrics.py
      derive_calibration_gap.py
      derive_benchmark_scale.py
      derive_field_transfer.py
    r/                      figure authority; consumes exports only
      _theme.R
      fig01..fig06
      make_all_figures.R
  manuscript/
    Manuscript-O3.md         assembled manuscript
    sections/                authoritative section sources
    supplementary/           Supplementary Methods
    artifacts/
      data/                  tidy CSV exports (the evidence interface)
      figures/                PDF and PNG at 600 dpi
      TAB-*.md                generated tables
    review/                  reviewer report and resolution ledger
    LITERATURE.bib
```

## Reproduction

From the repository root:

```
.venv/Scripts/python.exe O3/analysis/python/derive_taxonomy.py
.venv/Scripts/python.exe O3/analysis/python/derive_headline_metrics.py
.venv/Scripts/python.exe O3/analysis/python/derive_calibration_gap.py
.venv/Scripts/python.exe O3/analysis/python/derive_benchmark_scale.py
.venv/Scripts/python.exe O3/analysis/python/derive_field_transfer.py
"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" O3/analysis/r/make_all_figures.R
```

Python reads only files already committed under `M3/m3_age_benchmark/`,
`M4/m4_topology_benchmark/` and `M5/m5_inverse_reaction_benchmark/` and
writes CSV exports under `manuscript/artifacts/data/`. R consumes those
exports and recomputes no reported statistic.
