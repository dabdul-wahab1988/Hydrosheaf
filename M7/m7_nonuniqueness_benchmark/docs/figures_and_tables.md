# M7.3-M7.5 manuscript figures and tables

## Main figures

### Figure 1 | Benchmark architecture and claim boundary

The independent synthetic-truth branch separates official flow/pathline
generation, nonlinear chemistry and tracer generation, blinded HydroSheaf
inference, and locked negative-control scoring. The Ghana branch distinguishes
supportable component diagnostics from field interpretations that are
non-identifiable under the available evidence.

### Figure 2 | Evidence integration is conditional

Panel a compares native hydraulic (H), age (A), chemistry (C), pairwise and
fully integrated topology-ranking performance. Panel b shows case-block
incremental PR-AUC contrasts. Panels c and d separate the misspecification
stress test into certainty and discrimination shifts. Their aligned 2 × 2
layout makes false confidence explicit: permuted evidence can reduce entropy
while degrading PR-AUC.

### Figure 3 | Correct topology improves age inference

Panels a and b show changes in age MAE and 95% interval width for correct,
partial and reversed topology assumptions. Panel c reports importance-sampling
ESS fractions, including the reversed-graph failures. Panel d shows coverage
changes. Reversed tritium-only accuracy contrasts are not interpreted because
8/12 cases failed the ESS rule.

### Figure 4 | Reaction recovery remains process dependent

The 2 × 2 panels compare core and enhanced evidence for modal-family recovery,
probability assigned to the true family, support entropy and effective support
size. Carbonate weathering and precipitation remain unrecovered despite low
support entropy, demonstrating that stability can be confidently wrong.

### Figure 5 | M7-only Ghana evidence and claim boundary

Panels a and b map observed evidence to defensible claims. Panel c reports the
M7 field-audit sample (160 wells, 140 complete seasonal pairs and 320 seasonal
observations), and panel d reports the truth-free seasonal hold-forward test.
No M6 result is included.

### Figure 6 | Incremental contribution of affine sheaf structure

The balanced 2 x 2 display compares overall discrimination and calibration,
scenario-specific PR-AUC, and paired case-block contrasts among the edge-local
weighted graph, identity graph Laplacian, native affine sheaf, and permuted-map
control. Error bars are simultaneous intervals controlling all 120 published
contrasts. The intended printed width is 7.08 in and the minimum label size is
8 pt.

### Figure 7 | Robust and hybrid sheaf estimator diagnostic

The balanced 2 x 2 display reports fresh locked-test PR-AUC and log loss,
selected local-first/global-fallback contrasts against the edge-local graph,
and mechanism contrasts separating LOO robustification and native-map
semantics. Error bars are simultaneous intervals controlling all 560 published
contrasts. The intended printed width is 7.08 in and the minimum label size is
8 pt.

## Supplementary figures

### Figure S1 | Locked synthetic model domain

The representative confirmatory realization shows the MODFLOW/MODPATH truth
network, observation nodes and hydraulic heads in synthetic model-space
coordinates. It is not a geographic map and must not be placed on a Ghana
basemap.

## Main tables

1. Seven-audit design and claim map.
2. Compact primary M7.3 decision table.
3. Competence-matched M7.4 sheaf-versus-graph locked-test means.
4. Fresh-seed M7.5 estimator locked-test means (selected local weight 1.0).

## Supplementary tables

- Table S1: every evidence panel under every native/adverse condition.
- Table S2: all case-block evidence contrasts.
- Table S3: per-case topology-to-age sensitivity and ESS.
- Table S4: edgewise reaction support and entropy.
- Table S5: the predeclared conflict diagnostic, including its zero sensitivity.
- Table S6: multiplicity-correction robustness check.
- Table S7: complete sheaf-versus-graph paired contrast family.
- Table S8: strict public-pipeline system-acceptance conditions and contrasts.
- Table S9: complete M7.5 robust/hybrid paired contrast family.
- Table S10: all 120 M7.4 contrasts with family-wise simultaneous intervals.
- Table S11: all 560 M7.5 contrasts with family-wise simultaneous intervals.
- Table S12: post-review precision and future-replication planning audit.
- Table S13: public-pipeline selection rule and confusion-count audit.

## Reproduction

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_3_publication_assets.py
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_sheaf_vs_graph_assets.py
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_robust_hybrid_assets.py
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\assemble_m7_manuscript.py
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\assemble_m7_supplement.py
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' M7\m7_nonuniqueness_benchmark\scripts\make_supporting_model_domain_map.R
```

Each figure is exported as vector PDF, 600-dpi PNG and LZW-compressed 300-dpi
TIFF. The source manifest records the locked input files for every figure.
