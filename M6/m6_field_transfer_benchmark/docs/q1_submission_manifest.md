# M6 Q1 submission manifest

This is the authoritative display-asset inventory for the M6 manuscript.
Generate it with:

```powershell
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\run_m6_q1.py
```

Every figure is retained as vector PDF, 600-dpi PNG and LZW-compressed 300-dpi
TIFF. PDF is the editable/vector master, TIFF is the production upload, and PNG
is the review/preprint copy.

## Main figures

1. `figures/r_publication/figure1_dataset_tier_design`
2. `figures/r_publication/figure2_workflow`
3. `figures/r_publication/figure3_ng_stability`
4. `figures/r_publication/figure4_tier_ablation`
5. `figures/r_publication/figure5_field_prequential`
6. `figures/r_publication/figure6_limitation_map`

## Extended Data figures

1. `figures/extended_data/figureED1_external_transfer`
2. `figures/extended_data/figureED2_synthetic_validation`
3. `figures/extended_data/figureED3_circularity_sensitivity`

## Supplementary figures

`figures/r_publication/supplementary/figureS1_hydrochem_context` through
`figureS11_null_sensitivity`.

## Main tables

`tables/table1_dataset_readiness` through `table5_field_prequential`, each as
CSV and Markdown.

## Supplementary tables

`tables/tableS2_qc_missingness` through `tableS9_null_sensitivity`, each as CSV
and Markdown.

## Required claim boundaries

- Tier 4 means the maximum M6 chemistry/metadata tier, not complete
  age-head-screen field evidence.
- Ghana coordinates are masked and support country-scale sampling context only.
- Seasonal hold-forward is within-campaign chemistry prediction, not
  independent-basin validation.
- Null-model scores are screening flags for competing no-flow explanations,
  not calibrated field probabilities.
- Field labels, residence times and candidate edges are inferred references,
  not independent reaction, age or topology truth.

## Excluded from the submission package

- Python bytecode caches and other generated runtime caches.
- SVG duplicates; PDF is the maintained vector format.
- Superseded names `figure5_external_transfer`, `figure7_*`, `figure8_*` and
  `figure9_field_prequential`.
- The removed compatibility entry point `make_m6_dataset_figure.py`.
- Cross-milestone `run_m6_chemistry_*` scripts and their three chemistry
  robustness CSVs; these are retained only to reproduce legacy M2 assets.
