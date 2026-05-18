# M3 TracerLPM-Parity Implementation Plan

This document is a step-by-step implementation handoff for improving Hydrosheaf M3 public-age benchmarking from screening-level agreement toward a defensible TracerLPM-parity benchmark.

The goal is not to make a prettier Fig. 5. The goal is to make Hydrosheaf estimate the same benchmark quantity as the USGS reported-age table, then clearly separate that parity mode from Hydrosheaf's own model-selection mode.

## Non-Negotiable Rules

1. Do not mix parity mode and Hydrosheaf-selection mode in the same metric.
2. Do not use `hydrosheaf_selection_corrected` for the main Fig. 5 public-age validation until paired diagnostics prove it improves agreement.
3. Do not claim independent truth. `Rpt_TotAge_yrs` is a USGS model-reported age, not a direct observation.
4. Every change must be tested on a small subset before any full 1,272-row run.
5. Keep generated output filenames explicit so old and new metrics are not accidentally compared.
6. Do not change manuscript claims until the new benchmark summary is produced and reviewed.

## Current Problem Summary

The current M3 public benchmark is limited because:

- `LPM_Name` and `LPM_TracersMod` are available in the USGS source table, but Hydrosheaf-selection mode can ignore the reported model family.
- `Rpt_TotAge_yrs` may include unsaturated-zone travel time, while Hydrosheaf currently mostly estimates saturated-zone tracer age.
- `get_localized_input()` currently returns a generic tritium input and does not use latitude, longitude, study unit, or recharge year.
- Atmospheric gas histories in `hydrosheaf/nuclear/multi_tracer.py` are compact screening histories.
- SF6/CFC contamination is mostly handled by soft reliability weights, not by a censored or contaminated likelihood.
- 14C and 4He corrections are still too simple for old groundwater.
- Age fractions in the USGS table are not used as constraints.
- Current BMA averages scalar ages and uses plain AIC, which is unstable when there are only 2-5 tracers.
- Binary-mixture model refinement is incomplete because the local optimizer assumes a `mean_age_years` parameter.

## Key Source Files

Edit these files first:

- `hydrosheaf/nuclear/joint_lpm.py`
- `hydrosheaf/nuclear/multi_tracer.py`
- `hydrosheaf/nuclear/tracer_inputs.py`
- `hydrosheaf/nuclear/old_groundwater.py`
- `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py`
- `M3/m3_age_benchmark/scripts/run_m3_design_matrix.py`
- `M3/m3_age_benchmark/configs/design_matrix.yaml`
- `M3/m3_age_benchmark/scripts/make_publication_figures.py`
- `M3/m3_age_benchmark/scripts/make_publication_tables.py`

Add or update tests here:

- `tests/test_joint_lpm.py`
- `tests/test_multi_tracer_age.py`
- `tests/test_tracer_inputs.py`
- `tests/test_old_groundwater.py`
- `tests/test_m3_usgs_benchmark.py`
- `tests/test_m3_real_graph_benchmark.py`

## Success Metrics

The primary benchmark metric is log-age parity against `Rpt_TotAge_yrs`.

For any paired comparison, report:

- `N`
- log10 `R²`
- median absolute log10 error
- log10 RMSE
- within factor 2
- within factor 10
- number of gained factor-2 rows
- number of lost factor-2 rows
- improved fraction

Minimum target before manuscript upgrade:

- TracerLPM-parity mode must improve or preserve the reported-model baseline on a paired test.
- Full run should not degrade median absolute log10 error relative to the current reported-model parity mode.
- Hydrosheaf-selection/BMA must be reported separately unless it improves on paired rows.

## Required Commands

Run small checks first:

```powershell
python -m py_compile hydrosheaf\nuclear\joint_lpm.py hydrosheaf\nuclear\multi_tracer.py M3\m3_age_benchmark\scripts\run_m3_usgs_benchmark.py
python -m pytest tests\test_joint_lpm.py tests\test_multi_tracer_age.py tests\test_tracer_inputs.py tests\test_old_groundwater.py tests\test_m3_usgs_benchmark.py -q
python M3\m3_age_benchmark\scripts\run_m3_design_matrix.py --max-rows 80 --age-steps 90 --scenario parity_reported_corrected --output M3\m3_age_benchmark\results\m3_parity_smoke.csv
```

Only after the smoke tests pass, run a full parity benchmark:

```powershell
python M3\m3_age_benchmark\scripts\run_m3_design_matrix.py --full --age-steps 90 --scenario tracerlpm_strict_parity --output M3\m3_age_benchmark\results\m3_tracerlpm_strict_parity_full.csv
```

Do not run full Hydrosheaf-selection/BMA until checkpointing is added.

## Phase 1: Strict TracerLPM-Parity Mode

### Objective

Make a benchmark mode that follows the USGS reported model family and reported tracer set as closely as Hydrosheaf supports.

This mode should drive Fig. 5.

### Implementation Steps

1. Add a new scenario to `M3/m3_age_benchmark/configs/design_matrix.yaml`.

Use this scenario:

```yaml
  - scenario_id: tracerlpm_strict_parity
    label: Strict TracerLPM parity
    purpose: Match USGS reported LPM family, reported tracer set, corrected inputs, and reported total-age convention as closely as Hydrosheaf supports.
    gas_correction_mode: usgs_dgm
    c14_correction_mode: selected
    he4_mode: calibrated
    tracer_set: reported
    lpm_strategy: reported
    age_target_mode: reported_total
    uztt_mode: add_reported
```

2. In `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py`, preserve these columns from `Table_2_Ages.txt`:

- `LPM_Name`
- `LPM_TracersMod`
- `Rpt_TotAge_yrs`
- `Rept_TotAge_Err_yrs`
- `Rpt_UZtt_yrs`
- `Rpt_UZtt_Err_yrs`
- `Rpt_ChiSquare`
- `Rpt_Probability`
- `FracAnthropocene`
- `FracHolocene`
- `FracPleistocene`
- `AgeCat`
- `Comments`

3. Add a helper:

```python
def _reported_tracer_tokens(value: Any) -> set[str]:
    ...
```

It must map the USGS reported tracer text into Hydrosheaf keys:

- `3H` -> `tritium_TU`
- `3He(trit)` or `3H/3He` -> `he3_trit_TU` plus `tritium_TU`
- `SF6` -> `sf6_pptv`
- `CFC-11` or `CFC11` -> `cfc11_pptv`
- `CFC-12` or `CFC12` -> `cfc12_pptv`
- `CFC-113` or `CFC113` -> `cfc113_pptv`
- `14C` -> `c14_pmc`
- `4He` -> `he4_ccpg`
- `39Ar` -> `ar39_pmc`
- `85Kr` -> `kr85_pptv`

4. Add a helper:

```python
def _apply_reported_tracer_mask(obs: Mapping[str, Any], lpm_tracers_mod: Any) -> dict[str, Any]:
    ...
```

It must set non-reported tracer values and sigmas to `None`.

Example: if `LPM_TracersMod` is `3H, 3He(trit), 14C`, then SF6, CFCs, 4He, 39Ar, and 85Kr must be removed.

5. In `_fit_prepared_benchmark_row()`, apply `_apply_reported_tracer_mask()` when:

```python
factors.get("tracer_set") == "reported"
```

6. Make sure `_supported_reported_model()` maps USGS reported models to Hydrosheaf model names correctly.

Minimum mapping:

- `PFM` -> `PFM`
- `EM` -> `EM`
- `DM` -> `DM`
- `GA` -> `GA`
- `EPM` -> `EPM`
- `PEM` -> `PEM`
- `BMM-DM-DM` -> `BMM-DM-DM`
- `BMM-PEM-PEM` -> `BMM-PEM-PEM`

If unsupported, fall back to `GA`, but record a diagnostic column:

- `reported_model_supported`
- `reported_model_fallback_reason`

### Tests

Add tests to `tests/test_m3_usgs_benchmark.py`:

1. `test_reported_tracer_mask_keeps_only_usgs_reported_tracers`
2. `test_strict_parity_uses_reported_lpm_name`
3. `test_unsupported_reported_model_records_fallback`

### Acceptance Criteria

- A row with `LPM_TracersMod = "3H, 3He(trit), 14C"` must not fit SF6 or CFCs.
- A row with `LPM_Name = "DM"` must fit `DM`, not model-selected `PFM` or `EM`.
- `tracerlpm_strict_parity` must produce 80 finite rows in a smoke test.

## Phase 2: Unsaturated-Zone Travel-Time Handling

### Objective

Compare the same age quantity as the USGS reported age.

If `Rpt_TotAge_yrs` includes unsaturated-zone travel time, Hydrosheaf must either:

- add reported UZ travel time to its saturated tracer age, or
- compare against a saturated-only reference age.

For Fig. 5 parity, use `add_reported` first because it is closest to the USGS reported total-age convention.

### Implementation Steps

1. Add helper in `run_m3_usgs_benchmark.py`:

```python
def _reported_uztt_years(row: Mapping[str, Any]) -> float:
    ...
```

Return `0.0` when missing, non-finite, or negative.

2. Add helper:

```python
def _apply_age_target_mode(est_saturated_age: float, row: Mapping[str, Any], factors: Mapping[str, Any]) -> tuple[float, dict[str, Any]]:
    ...
```

Supported modes:

- `reported_total`: compare estimated total age to `Rpt_TotAge_yrs`
- `saturated_only`: compare saturated estimate to `Rpt_TotAge_yrs - Rpt_UZtt_yrs`

Supported UZ modes:

- `add_reported`: `est_total = est_saturated_age + Rpt_UZtt_yrs`
- `ignore`: `est_total = est_saturated_age`

3. Add output columns:

- `est_age_saturated_years`
- `est_age_total_years`
- `reported_uztt_years`
- `age_target_mode`
- `uztt_mode`
- `ref_age_saturated_proxy`

4. In metrics, use `est_age_total_years` as `est_age_multi` only for `reported_total/add_reported`.

### Tests

Add tests:

1. Missing `Rpt_UZtt_yrs` returns `0.0`.
2. `est_saturated=20`, `Rpt_UZtt_yrs=5`, `uztt_mode=add_reported` gives `est_total=25`.
3. `saturated_only` reference subtracts UZ travel time but never below `0.1`.

### Acceptance Criteria

- Output CSV contains both saturated and total estimated ages.
- Summary clearly states whether Fig. 5 uses total or saturated age.

## Phase 3: Site-Specific Input Histories

### Objective

Replace compact screening atmospheric histories with site/year-aware input histories.

This is necessary because tritium, CFC, SF6, and recharge temperature vary by hemisphere, latitude, altitude, recharge year, and dissolved-gas correction.

### Implementation Steps

1. Create or extend `hydrosheaf/nuclear/tracer_inputs.py` with:

```python
@dataclass(frozen=True)
class SiteInputContext:
    site_id: str
    sample_year: float
    latitude: float | None = None
    longitude: float | None = None
    study_unit: str = ""
    aquifer_group: str = ""
    recharge_temperature_c: float | None = None
    elevation_m: float | None = None
```

2. Add:

```python
def build_site_tracer_histories(context: SiteInputContext) -> dict[str, InputHistory]:
    ...
```

First implementation can be simple but must be structured:

- choose Northern Hemisphere or Southern Hemisphere by latitude
- return tritium history
- return SF6, CFC11, CFC12, CFC113, 85Kr histories
- record `history_source` metadata

3. Replace `get_localized_input()` in `run_m3_usgs_benchmark.py` with a real call to `build_site_tracer_histories()`.

4. Pass histories into `fit_lpm_model()` and `fit_lpm_models()`.

5. Add output columns:

- `input_history_mode`
- `input_history_region`
- `input_history_source`

### Tests

Add tests:

1. Northern latitude returns Northern Hemisphere histories.
2. Southern latitude returns Southern Hemisphere histories or a recorded fallback.
3. `fit_benchmark_row()` passes a non-empty `histories` mapping to `fit_lpm_model()`.

### Acceptance Criteria

- No benchmark run silently uses generic histories without recording it.
- Fig. 5 caption can state which input-history mode was used.

## Phase 4: Censored and Contamination Likelihoods for Gases

### Objective

Handle contaminated SF6/CFC values as censored or contaminated observations rather than ordinary Gaussian observations.

Current soft downweighting can still let contaminated gas values pull the fit toward unrealistic ages.

### Implementation Steps

1. Add a new field to `TracerFitObservation` in `hydrosheaf/nuclear/joint_lpm.py`:

```python
likelihood: str = "gaussian"
```

Supported values:

- `gaussian`
- `upper_censored`
- `lower_censored`
- `contaminated_mixture`

2. Add function:

```python
def _standardized_observation_loss(obs: TracerFitObservation, predicted: float) -> float:
    ...
```

Rules:

- `gaussian`: existing squared standardized residual
- `upper_censored`: no penalty when predicted <= observed; penalty only when predicted > observed
- `lower_censored`: no penalty when predicted >= observed; penalty only when predicted < observed
- `contaminated_mixture`: use a robust loss, for example Student-t or Huber loss

3. In `build_lpm_tracer_observations()`, mark SF6/CFC as censored/contaminated when diagnostics say:

- above historical maximum
- impossible modern age relative to independent proxy
- strong disagreement with tritium/3He age

4. Add output columns:

- `gas_likelihood_mode`
- `censored_gas_count`
- `contaminated_gas_count`

### Tests

Add tests in `tests/test_joint_lpm.py`:

1. Upper-censored observation gives zero loss when prediction is below observed.
2. Upper-censored observation penalizes prediction above observed.
3. Contaminated robust loss is smaller than Gaussian loss for a very large residual.

### Acceptance Criteria

- Supersaturated SF6 no longer drives the fitted age to a wrong young/old solution.
- Benchmark rows record how gas observations were treated.

## Phase 5: Hierarchical 14C and 4He Corrections

### Objective

Improve old-groundwater performance by using study-unit/aquifer-level priors for:

- initial 14C activity `A0`
- dead-carbon correction
- 4He background
- 4He accumulation rate

### Implementation Steps

1. In `hydrosheaf/nuclear/old_groundwater.py`, add:

```python
@dataclass(frozen=True)
class OldGroundwaterPrior:
    study_unit: str
    aquifer_group: str
    a0_pmc_mean: float
    a0_pmc_sigma: float
    he4_background_mean: float
    he4_background_sigma: float
    he4_rate_mean: float
    he4_rate_sigma: float
    n_support: int
```

2. Build priors from available USGS rows:

```python
def build_old_groundwater_priors(df: pd.DataFrame) -> dict[str, OldGroundwaterPrior]:
    ...
```

Use group keys:

1. `StudyUnit`
2. `AqGroup`
3. global fallback

3. Apply priors in `prepare_c14_observation()` and `apply_he4_uncertainty_mode()`.

4. Add scenario:

```yaml
  - scenario_id: tracerlpm_parity_hier_oldwater
    label: Strict parity plus hierarchical old-water priors
    gas_correction_mode: usgs_dgm
    c14_correction_mode: hierarchical
    he4_mode: hierarchical
    tracer_set: reported
    lpm_strategy: reported
    age_target_mode: reported_total
    uztt_mode: add_reported
```

5. Add output columns:

- `oldwater_prior_scope`
- `oldwater_prior_n_support`
- `a0_prior_mean`
- `a0_prior_sigma`
- `he4_background_prior_mean`
- `he4_rate_prior_mean`
- `he4_rate_prior_sigma`

### Tests

Add tests:

1. Study-unit prior is preferred over aquifer fallback.
2. Aquifer fallback is preferred over global fallback.
3. Missing old-water support records fallback scope.

### Acceptance Criteria

- Old-water rows show lower median error without damaging modern rows.
- Improvements must be reported by age class.

## Phase 6: Age Fractions as Fitting Targets

### Objective

Use USGS age fractions as additional constraints on the transit-time distribution.

The source table includes:

- `FracAnthropocene`
- `FracHolocene`
- `FracPleistocene`

These constrain TTD shape, not only scalar mean age.

### Implementation Steps

1. In `hydrosheaf/nuclear/joint_lpm.py`, add:

```python
def age_fraction_predictions(model: str, parameters: Mapping[str, float], max_age_years: float = 50000.0) -> dict[str, float]:
    ...
```

Use transit-time PDF:

- Anthropocene: `tau <= sample_year - 1950` approximately, or use a fixed 70-year cutoff for first implementation.
- Holocene: `tau <= 11700`
- Pleistocene: `tau > 11700`

For a simple first version:

- `anthropocene_fraction = integral tau <= 70`
- `holocene_fraction = integral 70 < tau <= 11700`
- `pleistocene_fraction = integral tau > 11700`

2. Add optional observation fields:

- `frac_anthropocene`
- `frac_holocene`
- `frac_pleistocene`

3. Add loss terms to `fit_lpm_model()` only when scenario factor `use_age_fractions=True`.

Suggested standard deviation:

- `sigma_fraction = 0.10` initially

4. Add output columns:

- `pred_frac_anthropocene`
- `pred_frac_holocene`
- `pred_frac_pleistocene`
- `age_fraction_loss`

5. Add scenario:

```yaml
  - scenario_id: tracerlpm_parity_agefractions
    label: Strict parity plus age-fraction constraints
    gas_correction_mode: usgs_dgm
    c14_correction_mode: selected
    he4_mode: calibrated
    tracer_set: reported
    lpm_strategy: reported
    age_target_mode: reported_total
    uztt_mode: add_reported
    use_age_fractions: true
```

### Tests

Add tests:

1. A young PFM predicts high Anthropocene fraction.
2. An old PFM predicts high Pleistocene fraction.
3. Fraction predictions sum to approximately 1.
4. Age-fraction loss changes objective when enabled.

### Acceptance Criteria

- Shape-sensitive rows improve without overfitting scalar age.
- Report performance by `AgeCat`.

## Phase 7: Replace Plain BMA with Gated AICc/LOO BMA

### Objective

Prevent BMA from damaging benchmark agreement when model ranking is weak or model ages are divergent.

### Implementation Steps

1. In `JointLpmFit`, add fields:

- `aicc`
- `loo_score`
- `effective_n_params`

2. Replace current plain AIC weighting with AICc:

```python
def _aicc(aic: float, n: int, k: int) -> float:
    if n <= k + 1:
        return float("inf")
    return aic + (2 * k * (k + 1)) / (n - k - 1)
```

3. Add:

```python
def compute_gated_bma_age(fits: list[JointLpmFit], *, max_delta_aicc: float = 4.0, max_log_age_span: float = 0.7) -> dict[str, Any]:
    ...
```

Rules:

- Ignore non-converged fits.
- Use only fits within `delta_aicc <= 4`.
- If contributing model ages span more than `0.7 log10 years`, do not average; use top model.
- If fewer than two models pass gates, use top model.
- Record whether BMA was used or skipped.

4. Add output columns:

- `bma_used`
- `bma_skip_reason`
- `bma_n_models`
- `bma_log_age_span`
- `top_model_age_years`
- `bma_age_years`

5. In `run_m3_usgs_benchmark.py`, use gated BMA only for Hydrosheaf-selection scenarios, not strict parity.

### Tests

Add tests:

1. Divergent model ages trigger top-model fallback.
2. Close model ages allow BMA.
3. AICc is infinite when `n <= k + 1`.
4. Selection mode records `bma_used`.

### Acceptance Criteria

- Selection/BMA no longer severely degrades paired benchmark rows.
- BMA is diagnostic unless it improves paired metrics.

## Phase 8: Proper Binary-Mixture Model Refinement

### Objective

Fix local continuous optimization for BMM models.

Current refinement assumes `mean_age_years`; BMM models use:

- `mean_age_1_years`
- `mean_age_2_years`
- `binary_fraction`
- possible `dispersion`, `dispersion_2`
- possible `capture_fraction`, `capture_fraction_2`

### Implementation Steps

1. In `joint_lpm.py`, refactor local refinement into:

```python
def _pack_parameters(model: str, params: Mapping[str, float]) -> tuple[np.ndarray, list[str], list[tuple[float, float]]]:
    ...

def _unpack_parameters(model: str, x: np.ndarray, keys: list[str], bounds: list[tuple[float, float]]) -> dict[str, float]:
    ...
```

2. For non-BMM models:

- pack `mean_age_years`
- pack secondary parameters

3. For BMM models:

- pack `log10(mean_age_1_years)`
- pack `log10(mean_age_2_years)`
- pack `logit(binary_fraction)`
- pack dispersion/capture parameters if present

4. Enforce constraints:

- `mean_age_1_years >= 0.01`
- `mean_age_2_years > mean_age_1_years`
- `0.001 <= binary_fraction <= 0.999`
- dispersion bounds: `1e-4` to `2.0`
- capture bounds: `0.05` to `1.0`

5. Use `scipy.optimize.minimize(..., method="Nelder-Mead")` first.

6. If Nelder-Mead fails, fallback to `L-BFGS-B` on bounded transformed parameters.

7. Remove noisy warning:

```text
Refinement failed for BMM-DM-DM: 'mean_age_years'
```

Replace with a structured diagnostic in the returned fit:

- `refinement_attempted`
- `refinement_success`
- `refinement_message`

### Tests

Add tests:

1. BMM synthetic data refines without `mean_age_years`.
2. BMM refinement preserves `mean_age_2 > mean_age_1`.
3. BMM refinement objective is not worse than grid objective.

### Acceptance Criteria

- No BMM refinement warnings in an 80-row smoke test.
- BMM rows have finite refined parameters.

## Phase 9: USGS-Calibrated Parity Layer

### Objective

Create a legitimate benchmark-emulation layer that improves R² against USGS reported ages while being clearly labelled as calibrated emulation, not independent validation.

### Implementation Steps

1. Add script:

```text
M3/m3_age_benchmark/scripts/run_m3_usgs_calibrated_emulator.py
```

2. Inputs:

- pointwise strict-parity output CSV
- source USGS metadata

3. Features:

- `log10(est_age_total_years)`
- `multi_model`
- `n_tracers`
- `young_gas_proxy_coherence`
- `old_groundwater_case`
- `age_class`
- `StudyUnit`
- `AqGroup`
- `Depth_m`
- `Rpt_Probability`
- `Rpt_ChiSquare`
- `FracAnthropocene`
- `FracHolocene`
- `FracPleistocene`

4. Target:

```text
log10(Rpt_TotAge_yrs)
```

5. Validation:

Use leave-study-unit-out cross-validation.

Do not use random row split as primary evidence.

6. Models:

Start simple:

- ridge regression
- isotonic calibration by age class
- optional gradient boosting only as secondary

7. Outputs:

- calibrated prediction
- raw prediction
- residual
- fold id
- study unit held out
- calibration model

8. Required label:

```text
USGS-calibrated benchmark emulation
```

### Tests

Add tests:

1. Leave-study-unit-out split never puts the same `StudyUnit` in train and test.
2. Calibration script writes fold predictions.
3. Calibration does not train on missing target rows.

### Acceptance Criteria

- Calibrated layer improves held-out log-age R² over raw strict-parity mode.
- Manuscript text clearly labels it as calibrated emulation.

## Phase 10: Figure and Table Updates

### Objective

Make Fig. 5 and M3 tables report the correct mode and avoid overclaiming.

### Implementation Steps

1. Update `M3/m3_age_benchmark/scripts/make_publication_figures.py` to prefer:

1. `m3_tracerlpm_strict_parity_full.csv`
2. `m3_tracerlpm_parity_hier_oldwater_full.csv`
3. `m3_phase4_screened_full_results.csv`

2. Update M2 Fig. 5 loader in:

```text
M2/m2_benchmark/scripts/make_publication_figures.py
M2/m2_benchmark/scripts/make_supplementary_figures.py
```

3. Update table scripts to include:

- strict parity metrics
- Hydrosheaf-selection metrics separately
- calibrated-emulation metrics separately, if available

4. Update docs:

- `M3/m3_age_benchmark/docs/m3_results_summary.md`
- `M2/m2_benchmark/docs/02_figures.md`
- `M2/m2_benchmark/docs/m2_results_summary.md`

### Acceptance Criteria

- Fig. 5 title or caption states the exact mode.
- Tables do not mix strict parity, selection/BMA, and calibrated emulation.

## Recommended Implementation Order

Do not implement all phases at once.

Implement in this order:

1. Strict TracerLPM-parity mode.
2. Unsaturated-zone travel-time handling.
3. Reported tracer mask and metadata diagnostics.
4. BMM refinement fix.
5. Censored gas likelihoods.
6. Hierarchical 14C/4He priors.
7. Age-fraction constraints.
8. Gated AICc/LOO BMA.
9. USGS-calibrated parity layer.
10. Figure/table updates.

## First Pull Request Scope

The first PR should include only:

- strict parity scenario
- reported tracer mask
- UZ travel-time handling
- tests
- 80-row smoke result

Do not include BMA, hierarchical priors, or calibrated emulation in the first PR.

## First PR Acceptance Checklist

- [ ] `python -m py_compile ...` passes.
- [ ] `python -m pytest tests\test_m3_usgs_benchmark.py tests\test_joint_lpm.py -q` passes.
- [ ] `m3_parity_smoke.csv` has 80 rows.
- [ ] `tracerlpm_strict_parity` has finite estimates for most rows.
- [ ] Output has `est_age_saturated_years`, `est_age_total_years`, `reported_uztt_years`.
- [ ] Output records `reported_model_supported`.
- [ ] Output records the tracer set actually used.
- [ ] No manuscript claim is changed.

## Full Benchmark Checklist

After first PR:

```powershell
python M3\m3_age_benchmark\scripts\run_m3_design_matrix.py --full --age-steps 90 --scenario tracerlpm_strict_parity --output M3\m3_age_benchmark\results\m3_tracerlpm_strict_parity_full.csv
python M3\m3_age_benchmark\scripts\make_publication_tables.py
python M3\m3_age_benchmark\scripts\make_publication_figures.py
python M2\m2_benchmark\scripts\make_publication_figures.py
python M2\m2_benchmark\scripts\make_supplementary_figures.py
```

Then report:

- Full-run strict parity metrics.
- Full-run screened benchmark metrics.
- Pairwise deltas where site IDs overlap.
- Age-class metrics.
- Tracer-count metrics.
- Study-unit held-out metrics if calibration layer is implemented.

## Expected R² Improvement Route

The most likely legitimate R² improvement path is:

1. Reported tracer mask prevents extra tracers from pulling away from USGS reported fit.
2. Reported LPM model preserves the USGS target convention.
3. UZ travel time aligns estimated age with reported total age.
4. Site-specific histories reduce systematic young-tracer bias.
5. Censored gas likelihood prevents contamination from dominating.
6. Hierarchical old-water priors reduce 14C/4He old-age bias.
7. Age fractions constrain TTD shape.
8. Calibrated emulation improves R² against USGS reported ages, but must be labelled as calibrated emulation.

## Claim Discipline

Use these labels:

- `TracerLPM-parity benchmark`: Hydrosheaf follows USGS reported model/tracer conventions.
- `Hydrosheaf-selection benchmark`: Hydrosheaf chooses model structure itself.
- `USGS-calibrated benchmark emulation`: Hydrosheaf output is calibrated to match USGS reported ages under held-out validation.

Avoid these claims unless proven by paired full-run metrics:

- "equivalent to TracerLPM"
- "better than TracerLPM"
- "validated against true age"
- "guaranteed age accuracy"

Acceptable wording after strict parity improves:

> Under a strict reported-model parity configuration, Hydrosheaf reproduces USGS model-reported groundwater ages with quantified log-age agreement. Hydrosheaf-selection and calibrated-emulation modes are reported separately because they answer different validation questions.

