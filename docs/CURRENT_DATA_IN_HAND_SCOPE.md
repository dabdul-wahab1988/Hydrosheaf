# Current data-in-hand execution scope

**Scope lock:** 2026-09-08  
**Purpose:** keep the active HydroSheaf work centred on the data and
generators already present in this repository.

The active analyses use three complementary sources:

1. the four separate packages under `data/FieldData/`:
   `LowerAnayari`, `NorthenGhana`, `NorthernGhanaNew`, and
   `Talensi_MiningArea`;
2. the supplied Aiken County release under `data/AikenCounty`, read through
   the calibrated MODFLOW/MODPATH reference adapter; and
3. the independent controlled-synthetic generators used for topology,
   age/transport, and reaction component tests.

## Analyses supported by the available material

- canonical-unit and missingness audits for each field package;
- geochemical-ratio diagnostics, with the ratio covariance/missingness kept
  explicit and no duplicate counting of the source ions;
- mapped-geology context priors, optionally joined through a user-supplied
  `geology_key` dictionary and retained as metadata rather than screened
  lithology;
- Aiken model-conditioned endpoint and pathline transport summaries,
  direction/censoring diagnostics, and CFC apparent-age interval
  concordance;
- controlled direct-versus-indirect age Bayes-factor scoring with sealed
  synthetic truth, complete reachability denominators, adverse controls, and
  case-block uncertainty intervals; and
- truth-free transfer/screening summaries from the observed field packages.

The field packages are loaded separately. A combined analysis requires an
explicit harmonisation declaration and a new run identifier; no workbook is
silently merged with another workbook that has a different schema or sampling
frame.

## Reporting rule

The active report should lead with the results that can be computed from these
assets: data coverage, ratios, mapped-geology sensitivity, transport
concordance, age identifiability strata, abstention rates, and controlled
synthetic directness performance. Machine-readable manifests retain the source
kind and scoring flags so that observed, calibrated-model, and synthetic rows
cannot be pooled accidentally. The current locked M2/M2.3/M7.3 artifacts are
not overwritten by these analyses.

The reproducible entry point for the age/transport work is:

```powershell
& '.venv\Scripts\python.exe' `
  M7\m7_nonuniqueness_benchmark\scripts\run_age_adjacency_two_tier.py `
  --aiken-root data\AikenCounty `
  --output .codex_work\runs\RUN-AGE-ADJACENCY-TWO-TIER-YYYYMMDD-01
```

The output contains separate synthetic and Aiken child panels and a manifest
that records the non-pooling rule.
