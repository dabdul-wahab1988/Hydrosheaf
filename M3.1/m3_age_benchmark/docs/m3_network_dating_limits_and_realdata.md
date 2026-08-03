# Network-Enhanced Dating — Bayesian Limits and Real-Data Findings

> **Historical diagnostic — superseded 2026-07-28.** The controlled-twin
> experiments remain useful implementation history, but they are not the
> current real-USGS graph benchmark or a field-validation result.

Two investigations requested to give the network-dating result full weight.

## 1. Can the Bayesian MCMC (`infer_network_ages_bayesian`) date old water? No.
Attempted fixes on the controlled twin (true ages to ~360 yr; real tritium model):
- tighten the velocity prior (σ 50 → 6 → 1.5 m/yr);
- reparameterise age as `recharge_age + cumulative_flow_distance / velocity` so flow
  geometry dates old water.

**Outcome: not fixable from tritium alone.** In every configuration the sampler drove
velocity far outside its prior (to 340–4900 m/yr), collapsing old ages to < 30 yr
(within-factor-2 ≈ 0). Root cause is structural, not a bug: tritium is *dead* beyond
~60–70 yr, so the likelihood is flat there and the joint velocity/input-scale/age
posterior is **degenerate** — the model can trade a huge velocity + small input-scale for
young ages with no likelihood penalty. No reparameterisation can manufacture information
the tracer does not carry.

**Conclusion / action taken.** The framework model was restored to its published form with
an explicit identifiability caveat added in `network_aging.py`. Bayesian network dating is
reliable in the **tracer-informative / ambiguity-resolution** regime; accurate old-water
dating requires an independent velocity (hydraulic) constraint or a longer-lived tracer
(¹⁴C). The robust network-dating method for the tritium bomb-peak regime is the
**ordering-based alias resolution** demonstrated in `run_m3_network_dating_demo.py`
(within-factor-2 0.63 → 0.84 on ambiguous nodes).

## 2. Real bomb-peak tracer data (`m3_cv_benchmark_3H.csv`, 794 USGS sites)
**Fundamental limitation: real groundwater sites have no independent *true* age**, so a
true-accuracy demonstration of ambiguity resolution is impossible on real data (this is the
same reason the whole field is hard). What the real data *can* show:
- **Prevalence & consequence:** graph priors changed the inferred ³H age by > 2 yr on
  **5% (42/794)** of sites, including **12 sites flipped young (< 30 yr) → old (> 50 yr)**
  — real bomb-peak alias resolutions.
- **Validation is inconclusive on real data:** among the 42 changed sites, the withheld-
  tracer prediction error was lower under the graph prior on only 9/42 (mean 1.99 vs 2.47),
  i.e. no decisive real-data confirmation — because there is no reference age and the
  tracer-prediction proxy is noisy.

**Honest position.** The accuracy benefit of network alias-resolution is demonstrated on the
controlled twin (known truth); on real data it is shown to be *prevalent and consequential*
(5% of sites, 12 flips) but *not independently validatable* without matched multi-tracer or
hydraulic data. Carrying it to full weight is a **data-availability step** (paired ³H + a
monotonic tracer such as CFC/SF6, or ³H/³He, on the same wells), not an analysis step.
