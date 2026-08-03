# M3 final accuracy audit

Accuracy lock: 2026-07-28  
Branch: `codex/m3-correctness`  
Starting commit: `47fbf7ca385b3210547673519ac2b30f910a80d9`

## Verdict and scope

The current M3 code, frozen analysis outputs, publication tables and figures,
main manuscript, supplement, highlights, cover letter, and graphical abstract
are internally synchronized to the same corrected evidence package.

This is an accuracy statement about reproducibility, internal consistency,
reference-agreement calculations, and honest claim boundaries. It is **not** a
claim that published USGS LPM ages are error-free true ages, that Hydrosheaf has
been independently field-validated, or that graph regularisation is generally
beneficial.

## Authoritative claim boundary

Current claims must be taken from:

- `DECISIONS.md`;
- `docs/m3_results_summary.md`;
- `results/m3_manuscript_analysis_manifest.json`;
- `results/m3_design_matrix_manifest.json`;
- `results/m3_real_usgs_graph_benchmark_manifest.json`;
- `results/m3_cv_benchmark_*_manifest.json`;
- `tables/Manuscript_Ready/`; and
- `figures/Manuscript_Ready/`.

Older QA and planning notes are retained for provenance but are explicitly
labelled historical or superseded. They are not eligible sources for current
manuscript numbers.

## Locked numerical results

### Published-configuration agreement

The canonical dataset contains 1,272 harmonised rows. The reference-free
identifiability/reportability guard retained 329 strict reported-configuration
fits.

| Result | Reportable N | Median absolute log10 discrepancy | Log10 RMSE | Log10 R2 | Within factor 2 |
|---|---:|---:|---:|---:|---:|
| Strict reported configuration | 329 | 0.027932 | 0.276882 | 0.937147 | 0.875380 |
| Reported age-fraction sensitivity | 289 | 0.021613 | 0.196357 | 0.969780 | 0.916955 |
| Hydrosheaf-selected configuration | 309 | 0.130446 | 0.609790 | 0.763909 | 0.673139 |

The age-fraction result is a shared-provenance emulation sensitivity. The
fractions and reference ages come from the same USGS LPM release, so it is not
independent validation. The Hydrosheaf-selection and strict results have
different reportable supports; their aggregate metrics do not establish a
general ranking. Paired effects are confined to the 40-row common support in
Manuscript Table 6.

### Real-USGS graph benchmark

The graph contains 329 reportable nodes and 1,245 edge rows. The weak
hydraulic-proxy candidate changed log10 RMSE by -0.001349 but reduced
within-factor-two agreement from 0.875380 to 0.869300. It therefore failed the
predeclared robust-improvement rule. A wrong-direction weak control increased
RMSE by 0.100802 and a randomised weak control increased it by 0.654969.

No candidate graph supports a general reference-agreement improvement claim.

### Leakage-guarded tracer withholding

Every target tracer was withheld from node features before prediction.

| Tracer | Eligible rows | Reportable rows | Baseline RMSE | Hydraulic RMSE | Depth RMSE | Random RMSE | Eligible graph edges |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3H | 794 | 121 | 20.818104 | 20.818055 | 20.841494 | 21.566597 | 4 / 14 / 19 |
| SF6 | 262 | 75 | 2.843988 | 2.857993 | 2.928528 | 2.928427 | 5 / 12 / 14 |
| 14C | 1,103 | 169 | 26.830072 | 26.829250 | 26.902373 | 27.819135 | 11 / 28 / 30 |
| CFC-11 | 28 | 4 | non-estimable graph effect | — | — | — | 0 |
| CFC-12 | 16 | 6 | non-estimable graph effect | — | — | — | 0 |

The candidate graphs show no meaningful target-withheld predictive gain. The
CFC graph effects are non-estimable because no reportable CFC fits formed an
eligible graph edge.

## Withdrawn or non-estimable results

- The hierarchical old-water prior is withdrawn. Its audit performance was
  log10 RMSE 1.310 and R2 0.004, and its pooled priors were estimated from the
  same release rather than independent evidence.
- The gas-correction audit found zero supported corrected-versus-raw pairs.
  The effect is non-estimable; Supplementary Table S4 and Supplementary Figure
  S1 were removed.
- The MRVA release lacks a like-for-like reported-LPM table and is not used as
  cross-aquifer graph-age replication.
- Synthetic network examples demonstrate controlled behaviour and ambiguity;
  they are not field validation.
- The short 300-draw Bayesian network smoke test recovers the expected age
  ordering but is not convergence-qualified. The API retains the strict
  convergence flag (R-hat <= 1.01, bulk and tail ESS >= 400, zero divergences),
  and the manuscript does not use this smoke run as a performance claim.

## Publication-artifact integrity

| Artifact | Pages | QA result |
|---|---:|---|
| `M3_geochemistry/Manucript_3.docx` | 31 | Rendered and visually inspected; no clipping, blank pages, or red revision text |
| `M3_geochemistry/Supplementary_Information_M3.docx` | 7 | Rendered and visually inspected; no clipping, blank pages, or red revision text |
| `M3_geochemistry/Highlights-AGeochemistry.docx` | 1 | Rendered and visually inspected |
| `M3_geochemistry/Cover Letter-AGeochemistry.docx` | 2 | Rendered and visually inspected |
| `M3_geochemistry/Graphical Abstract.png` | — | Visually inspected |

All seven figures embedded in the manuscript and all three figures embedded in
the supplement are byte-for-byte matches to the current regenerated files in
`figures/Manuscript_Ready/`. References to true ages are expressly negated or
confined to controlled synthetic examples. The only remaining mention of
Supplementary Table S4 states that it is not reported because the gas effect is
non-estimable.

## Verification record

- `scripts/check_manuscript_artifacts.py`: M3 and M4 tracked tables reproduce
  from the frozen results.
- Selected M3 computational tests: 72 passed.
- Bayesian network structural/diagnostic smoke test: 1 passed.
- Repository-wide suite: 652 passed, 2 skipped, 4 failed, 13 errors.
- Every remaining failure/error is in `tests/synthetic_data_tests` and is caused
  by absent legacy fixtures under `data/synthetic/`, beginning with
  `water_chem_full.csv`. These are not used by the M3 analysis or manuscripts,
  but the repository as a whole must not be described as fully green until the
  fixtures are restored or the legacy test group is retired.
- `git diff --check`: no whitespace errors.

## SHA-256 provenance

```text
2fc28234ef9d5d6d4b24d533114b6f1b1f8abac6bb2aff62a5f3f96db8ed50bb  M3/M3_geochemistry/Manucript_3.docx
ef0e5e56b81e53d220de7ada5194574d5ef6f234ea54f48e045840b8a5aba955  M3/M3_geochemistry/Supplementary_Information_M3.docx
fc2be7bcc07bd5e3ba69ac3524541b8f32c892fd64a1f5187ef920777c6c3494  M3/M3_geochemistry/Highlights-AGeochemistry.docx
28278790d80dadfbf0c6a1070eaa37791dc34df0595202ea477d40ef260e49bb  M3/M3_geochemistry/Cover Letter-AGeochemistry.docx
50925d358ecd1ff2f6704089b48f52560a6955237bef0cedbad8ec30f21feab7  M3/M3_geochemistry/Graphical Abstract.png
847840ebed1c843eb564eea5b9d9cc1ba7385f6fdcb46756364653923116b5f3  M3/m3_age_benchmark/results/m3_manuscript_analysis_manifest.json
3c97d1f2d92c6d598bc2cf4bf5942636d33d4bfd484beee48e768b1b688edb1f  M3/m3_age_benchmark/results/m3_design_matrix_manifest.json
fc44bad8e49b9275956fa02811dbc91685747588fec95ae3f14055f65b1f2daa  M3/m3_age_benchmark/results/m3_real_usgs_graph_benchmark_manifest.json
f1e1075386471408a42be69e1bdfb236fab982ae87739d653521426c33d7ab0a  M3/m3_age_benchmark/tables/Manuscript_Ready/Manuscript_Table2_Design_Matrix_Performance.csv
cead42007d1b8c78c142b6eef3d2b4a424169ec4467e17d4907ec03c7300815e  M3/m3_age_benchmark/tables/Manuscript_Ready/Manuscript_Table4_Real_USGS_Graph_Benchmark.csv
a4c90f0030085c03dbc9fea76905c5743f0c678792f73663592a4b2a706d0c66  M3/m3_age_benchmark/tables/Manuscript_Ready/Manuscript_Table6_Statistical_Significance.csv
7fe1c3d30bf5c65e1b704d5e80fd91f38d267fb9fcd881f6f9dc001fac356065  M3/m3_age_benchmark/figures/Manuscript_Ready/Manuscript_Fig3_USGS_Benchmark_Parity.png
29897d3348bbef0e5c397c97095345431b23594d969d9f93b71c27f782ad81c0  M3/m3_age_benchmark/figures/Manuscript_Ready/Manuscript_Fig4_Real_USGS_Graph_Benchmark.png
168ae9dfdc2bc91eecb160a855a1d28d82eb05bd43de345ac5042cf025f7b999  M3/m3_age_benchmark/figures/Manuscript_Ready/Manuscript_Fig6_Cross_Validation_Results.png
```

These hashes describe the working-tree accuracy lock. Commit provenance must be
updated after the files are committed.
