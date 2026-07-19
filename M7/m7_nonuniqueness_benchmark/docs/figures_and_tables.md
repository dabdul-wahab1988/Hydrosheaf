# M7.3 manuscript figures and tables

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
incremental PR-AUC contrasts. Panel c demonstrates false confidence: permuted
evidence reduces entropy while degrading PR-AUC.

### Figure 3 | Correct topology improves age inference

Panels a and b show changes in age MAE and 95% interval width for correct,
partial and reversed topology assumptions. Panel c reports importance-sampling
ESS fractions, including the reversed-graph failures. Panel d shows coverage
changes. Reversed tritium-only accuracy contrasts are not interpreted because
8/12 cases failed the ESS rule.

### Figure 4 | Reaction recovery remains process dependent

Core and enhanced panels are compared for modal-family accuracy, probability
on the true family and effective supported families. Carbonate weathering and
precipitation remain unrecovered despite low support entropy, demonstrating
that stability can be confidently wrong.

### Figure 5 | Ghana data support component diagnostics, not complete field truth

Panels a and b map observed evidence to defensible claims. Panel c shows the M6
evidence-tier ablation and panel d shows the truth-free M7.2 seasonal
hold-forward test. The figure operationalizes revised Objective 6 without
equating chemistry prediction with residence-time or topology validation.

## Main tables

1. Locked benchmark design and computational scale.
2. Native locked-test performance for all seven evidence panels.
3. Case-block evidence contrasts, including adverse controls.
4. Topology-conditioned age contrasts.
5. Reaction-family recovery and non-uniqueness.
6. Ghana data scope and claim boundary.

## Supplementary tables

- Table S1: every evidence panel under every native/adverse condition.
- Table S2: all case-block evidence contrasts.
- Table S3: per-case topology-to-age sensitivity and ESS.
- Table S4: edgewise reaction support and entropy.
- Table S5: the predeclared conflict diagnostic, including its zero sensitivity.

## Reproduction

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_3_publication_assets.py
```

Each figure is exported as vector PDF, 600-dpi PNG and LZW-compressed 300-dpi
TIFF. The source manifest records the locked input files for every figure.
