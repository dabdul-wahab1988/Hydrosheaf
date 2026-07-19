# M6 manuscript-ready asset plan

## Recommended six-figure main narrative

1. **Dataset and evidence ceiling:** existing Figure 1, explicitly describing
   Northern Ghana coordinates as masked and Tier 4 as the maximum M6
   chemistry/metadata tier.
2. **Field-transfer workflow:** existing Figure 2.
3. **Northern Ghana reaction-class stability:** existing Figure 3.
4. **Diagnostic ablation:** existing Figure 4, the principal M6 result.
5. **Within-campaign seasonal hold-forward:** new Figure 9. Graph ridge is
   superior to persistence but not distinguishable from the expanding
   mean-delta baseline.
6. **Limitation map:** existing Figure 6, showing conservative equivalence-class
   reporting and external-data non-identifiability.

Existing Figure 5 (external transfer) is suitable as Extended Data or may
replace Figure 3 if the target journal prioritizes generalization. Existing
Figures 7 and 8 should remain Extended Data because they answer validation and
circularity objections rather than advancing the main field narrative.

## Main tables

1. Dataset readiness and claim strength.
2. Variable availability.
3. Northern Ghana evidence-tier ablation.
4. External sparse-data transfer.
5. Northern Ghana truth-free seasonal hold-forward.

The `Tier 4` label must always be expanded as “maximum M6
chemistry/metadata tier”. It is not a complete field age–head–screen evidence
tier.

## Reproduction

```powershell
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\make_m6_dataset_figure.py
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\make_m6_tables.py
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\make_objective6_prequential_figure.py
```

The Python Figure 1 replaces the blank coordinate field with a sourced Ghana
outline and makes the coordinate mask explicit. The existing R figure suite
remains authoritative for M6 Figures 2–8. Python Figure 9 uses the same
two-column dimensions, accessible palette and vector/raster export conventions.
