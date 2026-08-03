# M3 benchmark decisions

## 2026-07-28: reported-configuration emulation and scalar-age reportability

Status: adopted for the regenerated manuscript analysis.

The USGS TracerLPM Table 4 release publishes tracer-specific input scale factors,
the parameters selected for optimization, and the initial or fixed LPM
configuration. The benchmark now uses those fields in the
`tracerlpm_strict_parity` scenario. For carbon-14, the published scale is applied
as the initial activity (`A0 = 100 * scale`); for the other supported tracers it
scales the modeled concentration. Only the published initial/fixed
configuration and optimization declaration are used. Reported final ages and
final fitted parameters are never supplied to the optimizer.

This scenario is therefore described as **USGS reported-configuration parity
emulation**, not independent validation. The age-fraction scenario additionally
uses USGS-reported age fractions and is labeled a reported-output-constrained
sensitivity analysis.

Raw scalar ages are retained in `est_age_multi_raw`. A scalar estimate is copied
to `est_age_multi` only when all reference-free reportability criteria pass:

- the optimizer converged;
- a requested reported-model emulation used an exactly supported LPM rather
  than a fallback model;
- the number of fitted tracer observations exceeds the effective number of free
  parameters;
- standardized tracer-space RMSE is at most 2.0; and
- grid candidates within four objective units of the best grid candidate span
  no more than 0.5 log10 age units (approximately a factor of 3).

Failure reasons are retained in `fit_identifiability_reason`. Reference ages are
used only after fitting to evaluate reportable estimates and do not participate
in the guard.

## 2026-07-28: hierarchical old-water prior withdrawn

Status: withdrawn from the active design matrix; retained as a disabled audit
scenario invokable with `--include-withdrawn`.

The corrected full rerun of `tracerlpm_parity_hier_oldwater` produced log10 RMSE
1.310 and log10 R2 0.004. More importantly, its empirical pooling priors were
estimated from the same evaluation release, without an external or held-out
calibration split. The scenario is therefore not eligible for manuscript tables,
figures, pairwise comparisons, or primary claims.

A redesign may re-enter the active matrix only after its priors are learned from
an external or explicitly held-out calibration set, its hyperparameters and
acceptance criteria are declared before evaluation, and it improves a locked
test set without reducing identifiability coverage or introducing catastrophic
young-water branches.

## Result lock

The manuscript pipeline is the result-locking command:

```powershell
.venv\Scripts\python.exe M3\m3_age_benchmark\scripts\run_m3_manuscript_analysis.py --full --age-steps 90
```

The resulting `m3_manuscript_analysis_manifest.json` records file sizes and
SHA-256 hashes for the canonical results, QA reports, tables, and figures.

## 2026-07-28: graph and tracer-withholding leakage controls

Status: adopted for the final M3 regeneration.

The main graph benchmark now consumes only reportable
`tracerlpm_strict_parity` rows from the canonical design-matrix result. Earlier
screened-scenario graph outputs are historical diagnostics and are not eligible
for manuscript claims. Because the USGS ages are model-derived reference values,
graph results are described as agreement diagnostics rather than field-validated
predictive accuracy.

Tracer-withholding cross-validation now fits every age entering a graph with the
target tracer withheld and then applies the same reference-free scalar-age
reportability guard as the main benchmark. Full-fit neighbouring ages are no
longer allowed into the graph because they can carry the withheld tracer back
into the prediction through regularisation.

The legacy calibrated-emulation output is withdrawn from manuscript Table 5.
It predates the reportability guard and used a reference-age-derived age class as
a predictor. The legacy national-plus-MRVA graph table is also withdrawn because
MRVA lacks the reported-LPM table required for a like-for-like strict-emulation
replication.
