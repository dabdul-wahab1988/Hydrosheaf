# M6 — Hydrosheaf Field-Transfer & Robustness Benchmark

Field transfer of the identifiability-aware Hydrosheaf reaction-inference workflow to
three real Ghanaian aquifer datasets, under data scarcity. M6 asks whether the M5
workflow degrades **gracefully, honestly and reproducibly** as evidence is removed —
it does **not** claim field validation of true reactions, flow paths or ages.

See `Outline.md` (Q1 Nature-Portfolio outline) and `docs/m6_locked_analysis_plan.md`.

## What is reused (frozen, not re-fit)
- M5 inverse-reaction primitives (`M5/.../scripts/m5_common.py`): FISTA sparse inverse,
  reaction dictionary, identifiability diagnostics, equivalence classes, SI bounds.
- The **frozen M5 MRS calibration** (`M5/.../results/mrs_calibration_model.json`) — the
  *transferred* identifiability classifier. This is the scientific core of a transfer study.

## Datasets (real)
All three datasets are sourced strictly from `data/FieldData/`
(`NorthenGhana/NorthernGhana.xlsx`, `Talensi_MiningArea/talensi.csv`,
`LowerAnayari/manu.csv`). An earlier revision additionally read a different,
antecedent study's own derived workbook (`Aquifers_Dataset_Mendeley.xlsx` — aquifer/
geology/lithology metadata, workbook saturation indices, a provided graph-edge sheet,
inferred process labels for the same Northern Ghana boreholes); that workbook is not
this project's field data and has been removed entirely (see `DECISIONS.md`).

| Dataset | n | Native tier | Role |
|---|---:|---|---|
| Northern Ghana (canonical raw workbook) | 320 (160 wells × wet/dry) | Tier 4 (M6 chemistry ladder) | Primary field transfer + robustness |
| Talensi (mining area) | 63 | Tier 1 | External sparse transfer |
| Lower Anayari (manu) | 41 | Tier 2 | External sparse transfer |

`Tier 4` is the maximum tier of the M6 reaction-diagnostic ladder: major ions + F +
isotopes + Sr/SiO2 + Hydrosheaf-computed PHREEQC saturation indices. It does not
imply that environmental age tracers, repeated heads, screen intervals, independent
flow-path truth, or any aquifer/geology/lithology metadata are available. No
independent aquifer-type classification exists for any of the three datasets;
stratified reporting uses administrative region/district for Northern Ghana instead.
See
[`docs/objective6_data_limited_synthesis.md`](docs/objective6_data_limited_synthesis.md)
for the revised Objective 6 interpretation.

## Reproduce
```powershell
# Complete deterministic Q1 workflow: analysis, tables, main, Extended Data,
# supplementary figures, and publication cartography.
.venv\Scripts\python.exe scripts\run_m6_q1.py

# Faster modes when locked results already exist.
.venv\Scripts\python.exe scripts\run_m6_q1.py --analysis-only
.venv\Scripts\python.exe scripts\run_m6_q1.py --figures-only
```

## Layout
```
Outline.md                    Q1 manuscript outline (main + supplementary plan)
docs/                         locked analysis plan, results summary, audits, manifest
scripts/
  m6_common.py                harmonisation, tiers, transport correction, transferred classifier, edges
  synthetic_reaction_truth_model.py
                              independent extended model for synthetic validation only
  run_m6_field_transfer.py    six experiments -> results CSVs
  run_m6_robustness_diagnostics.py
                              gate, transport, CBE and discrimination diagnostics
  run_m6_null_sensitivity.py  supplementary competing-no-flow sensitivity
  run_m6_chemistry_robustness.py
                              supporting regularisation/PSI results used by M2
  run_m6_chemistry_stress_tests.py
                              supporting LOO and bias checks used by M2
  export_m6_figure_data.py    map coords, tier ladder, hydrochem context
  make_m6_tables.py           5 main + 8 supplementary tables (CSV + MD)
  run_m6_q1.py                authoritative one-command Q1 orchestrator
r_figures/
  theme_m6.R                  shared Nature-style theme + palettes
  plot_m6_publication_figures.R      core main figures + External Data 1
  plot_m6_validation_figures.R       Extended Data 2–3
  plot_m6_supplementary_figures.R    11 supplementary figures
results/                      all derived CSVs
tables/                       main + supplementary tables
figures/r_publication/        six main figures: PDF, 600-dpi PNG, 300-dpi TIFF
figures/extended_data/        three Extended Data figures in the same formats
```

The two `run_m6_chemistry_*` scripts and their three `results/m6_*` chemistry
robustness CSVs were consolidated from the former separate M6 scaffold. They
preserve M2 cross-milestone reproducibility but are deliberately excluded from
`run_m6_q1.py` and from the M6 Q1 submission assets.

## Experiments
1. Dataset readiness / missingness (evidence tiers).
2. Northern Ghana seasonal × region transfer (Tier 4).
3. Tier ablation (Tier 4 → Tier 0).
4. Edge / flow-path uncertainty (3 Hydrosheaf-generated edge-set strategies).
5. External sparse transfer (Talensi, Lower Anayari).
6. Limitation map + conservative-vs-best-fit reporting.

## Honesty boundary
No dataset provides independent reaction, flow-path or age truth. An earlier revision
of this workflow additionally read a different, antecedent study's own inferred labels
(`Dominant_Process`, `Aquifer_Evolution_Label`, `Graph_Edges`) for the Northern Ghana
boreholes; those fields, and the workbook that supplied them, have been removed from
this study entirely rather than retained as a concordance reference (`DECISIONS.md`).

The supplementary null-model sensitivity treats scores above 0.5 and 0.8 as
screening flags for competing and dominant no-flow explanations. They are not
calibrated field probabilities. High flag rates further restrict interpretation
of inferred edges as physical flow paths.

Figure 1 is generated with the R spatial stack (`sf` + `ggspatial`). The
repository contains the unsimplified 2021 geoBoundaries ADM1 Shapefile for all
16 Ghana regions; its aligned national outline is derived with `sf::st_union()`.
Both the boundaries and WGS84 well coordinates are transformed to WGS84 / UTM
zone 30N (EPSG:32630) before plotting. Boundary provenance and licensing are in
`data/ghana_geoboundaries_2021/SOURCE.md`. All displayed well coordinates remain
masked. Talensi longitudes are converted from the source's positive degrees-west
convention to signed WGS84 longitudes. Northern Ghana masked positions falling
outside the higher-resolution border are constrained to a 5 km interior buffer
for display only. The map is a country-scale sampling locator, not a source for
recovering exact well positions or inferring local flow paths.
