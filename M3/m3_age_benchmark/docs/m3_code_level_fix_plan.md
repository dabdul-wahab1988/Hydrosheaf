# Hydrosheaf: Code-Level Fix Plan

Date: 2026-05-15

## The Root Cause of Limited Improvement
My re-examination of the codebase reveals that the "limited improvement" is directly caused by how the objective function calculates residuals in `hydrosheaf/nuclear/joint_lpm.py` combined with the binary masking in `run_m3_usgs_benchmark.py`.

In `fit_lpm_model`, the objective is a simple L2 sum of standardized residuals:
```python
standardized = (obs.value - pred) / obs.sigma
objective += standardized * standardized
```
Because there is no concept of "tracer reliability" or "weight", a slightly contaminated SF6 sample (that wasn't caught by the binary mask) will generate a massive standardized residual. Since it's squared, the optimizer is forced to pull the entire age distribution toward the contaminated tracer just to minimize the explosive penalty.

## The Code-Level Fixes Required

### 1. Introduce `weight` into `TracerFitObservation`
**File:** `hydrosheaf/nuclear/joint_lpm.py`

Modify the `TracerFitObservation` dataclass to support soft weighting:
```python
@dataclass(frozen=True)
class TracerFitObservation:
    tracer: str
    value: float
    sigma: float
    units: str
    note: str = ""
    weight: float = 1.0  # NEW
```

### 2. Update the Objective Function
**File:** `hydrosheaf/nuclear/joint_lpm.py`

In `fit_lpm_model` (around line 730), update the inner loop to apply the weight:
```python
        for obs in fit_observations:
            pred = predicted.get(obs.tracer)
            if pred is None or not math.isfinite(pred):
                missing_prediction = True
                break
            standardized = (obs.value - pred) / obs.sigma
            objective += (standardized * standardized) * getattr(obs, "weight", 1.0) # UPDATED
            residuals.append(...)
```

### 3. Read Weights from the Observations Dictionary
**File:** `hydrosheaf/nuclear/joint_lpm.py`

In `build_lpm_tracer_observations`, read `{tracer}_weight` keys for all tracers. For example, for gases:
```python
    for tracer, value_key, sigma_key in gas_specs:
        value = _finite_float(observations.get(value_key))
        if value is None or value < 0:
            continue
        weight = _finite_float(observations.get(f"{tracer.lower()}_weight"))
        if weight is None or weight < 0:
            weight = 1.0
            
        out.append(
            TracerFitObservation(
                tracer,
                value,
                _sigma(observations, sigma_key, value, relative=0.10, floor=0.05),
                "pptv",
                note="expects atmospheric-equivalent corrected gas input",
                weight=weight, # NEW
            )
        )
```

### 4. Replace Binary Masking with Soft Weighting
**File:** `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py`

Rewrite `_apply_young_gas_masking` to output weights instead of dropping the tracer entirely. 
```python
        # INSTEAD OF:
        # if value > max_hist * 1.02: out[value_key] = None
        
        # USE SOFT WEIGHTING:
        tracer_weight = 1.0
        if value > max_hist * 1.02:
            tracer_weight = 0.05  # Strongly down-weight supersaturated gas instead of dropping
            tracer_reasons.append("above_historical_max_downweighted")
            
        if (math.isfinite(reference_age) and ...):
            tracer_weight *= 0.1  # Down-weight due to proxy conflict
            tracer_reasons.append("excessively_modern_vs_independent_proxy_downweighted")
            
        out[f"{tracer.lower()}_weight"] = tracer_weight
```

### 5. Internal Dissolved Gas Optimization
Currently, Hydrosheaf uses the USGS DGM output statically. To fix this fundamentally:
1. In `hydrosheaf/nuclear/dissolved_gas.py`, expose `fit_recharge_conditions_from_nobles(observations)`.
2. In `_fit_prepared_benchmark_row` (in the benchmark script), call this *before* fitting the LPM.
3. If noble gases (Ne, Ar) exist, calculate the optimal Recharge Temp and Excess Air, then calculate the Atmospheric Equivalent of SF6 and CFCs dynamically.
4. Pass these *internal* atmospheric equivalents into the `fit_lpm_model` instead of the static `dgm_sf6` from the USGS tables.

## Summary
By changing the objective function to support `weight`, the LPM optimizer will stop allowing contaminated or questionable tracers to hijack the multi-tracer fit. This is the exact code change needed to convert Hydrosheaf from a brittle tool into a robust, probabilistic engine.