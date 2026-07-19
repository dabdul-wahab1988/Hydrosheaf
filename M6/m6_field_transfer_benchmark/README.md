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
| Dataset | n | Native tier | Role |
|---|---:|---|---|
| Northern Ghana (Mendeley workbook) | 320 (160 wells × wet/dry) | Tier 4 (M6 chemistry/metadata ladder) | Primary field transfer + robustness |
| Talensi (mining area) | 63 | Tier 1 | External sparse transfer |
| Lower Anayari (manu) | 41 | Tier 2 | External sparse transfer |

`Tier 4` is the maximum tier of the M6 reaction-diagnostic ladder. It does not
imply that environmental age tracers, repeated heads, screen intervals or
independent flow-path truth are available. See
[`docs/objective6_data_limited_synthesis.md`](docs/objective6_data_limited_synthesis.md)
for the revised Objective 6 interpretation.

## Reproduce
```bash
# analysis + auxiliary data + tables (deterministic, seed=1234)
python scripts/run_m6_all.py
# Nature-style figures (R 4.3+)
Rscript r_figures/plot_m6_publication_figures.R
Rscript r_figures/plot_m6_supplementary_figures.R
```

## Layout
```
Outline.md                    Q1 manuscript outline (main + supplementary plan)
docs/                         locked analysis plan, results summary, audits, manifest
scripts/
  m6_common.py                harmonisation, tiers, transport correction, transferred classifier, edges
  run_m6_field_transfer.py    six experiments -> results CSVs
  export_m6_figure_data.py    map coords, tier ladder, hydrochem context
  make_m6_tables.py           4 main + 8 supplementary tables (CSV + MD)
  run_m6_all.py               one-command orchestrator
r_figures/
  theme_m6.R                  shared Nature-style theme + palettes
  plot_m6_publication_figures.R      6 main figures
  plot_m6_supplementary_figures.R    10 supplementary figures
results/                      all derived CSVs
tables/                       main + supplementary tables
figures/r_publication/        PNG/PDF/TIF at 300 dpi
```

## Experiments
1. Dataset readiness / missingness (evidence tiers).
2. Northern Ghana seasonal × aquifer transfer (Tier 4).
3. Tier ablation (Tier 4 → Tier 0).
4. Edge / flow-path uncertainty (4 edge-set strategies).
5. External sparse transfer (Talensi, Lower Anayari).
6. Limitation map + conservative-vs-best-fit reporting.

## Honesty boundary
Mendeley `Dominant_Process`, `Aquifer_Evolution_Label`, `Graph_Edges` and heuristic
residence times are **inferred references, never validation targets**. Concordance with
them is reported as agreement with prior labels, not accuracy against ground truth.
