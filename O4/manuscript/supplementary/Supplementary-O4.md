# Supplementary Information — Fit quality is not identifiability: a harmonized robustness, integration-value and calibration audit of groundwater inference under data limitation

This Supplementary Information carries the full per-component stress-design
detail, equations, algorithms, extended tables and extended within-component
results that the main text (Section 3) can only summarise. Every number
below traces to the same locked `M6`/`M7`/`M8` result files cited in the main
text (Section 3.1, `O4/DECISIONS.md` D2-D3); no PHREEQC, MODFLOW/MODPATH,
PEST-GLM, or Bayesian active-learning re-run was performed to produce it.

## S1. Extended data provenance

### S1.1 M6: the five-level evidence-tier ladder

The Northern Ghana workbook (`data/FieldData/NorthenGhana/NorthernGhana.xlsx`)
provides 160 boreholes across four administrative regions (Northern, North
East, Upper East, Upper West), each sampled once, reported on Dry and
reconstructed-Wet worksheet panels; this paper, consistent with `M6`'s own
disclosure, uses the Dry panel as the primary measured chemistry and reports
no seasonal-change finding. The tier ladder is defined as:

- **T0 (majors):** calcium, magnesium, sodium, potassium, bicarbonate,
  chloride, sulphate, nitrate.
- **T1 (+isotopes):** T0 plus stable water isotopes (δ18O, δ2H).
- **T2 (+fluoride):** T1 plus fluoride.
- **T3 (+Sr/SiO2):** T2 plus strontium and dissolved silica.
- **T4 (full):** T3 plus Hydrosheaf-computed PHREEQC saturation indices
  (calcite, dolomite, gypsum, fluorite and related phases).

Each of the 160 wells is scored independently at each of the five tiers
(800 well-tier rows total), producing a resolution class, a mechanism-
resolution score (MRS), a bootstrap support-stability estimate, and a
dominant Cl-corrected process label. Talensi (63 samples, native Tier 1) and
Lower Anayari (41 samples, native Tier 2) enter only at their measured tier;
no tier ablation is performed on the external datasets, since their
measurement panels are fixed by what each survey originally collected.

### S1.2 M7: four independent generators

`M7.3`'s generator produces fresh MODFLOW 6/MODPATH 7 cases with paired
hydraulic, tracer (tritium and argon-39) and reaction-family evidence; six
development cases (seeds 5201-5206) and twelve untouched locked-test cases
(5301-5312) were used, with 50,000 age-importance particles per case/tracer
regime and 64 reaction bootstraps per case. `M7.4`/`M7.5` use a
code-independent scalar-section generator (`independent_sheaf_graph_
generator.py`) producing 64 held-out cases across four scenario strata
(native, identity-limit, heterogeneous-affine, incompatible-cycles,
noisy-missing) for `M7.4`, and a further 64 development plus 128 locked-test
cases for `M7.5`'s two-stage development/confirmatory design. `M7.6` uses a
third-generation generator (`independent_modflow_generator_v2.py`) informed
only by declared synthetic nuisance parameters, with six development cases
(5401-5406) and twelve locked-test cases (5501-5512), to test one specific
selective-mechanism hypothesis about the age layer's own infeasibility
pattern (Section S4.3). None of these four generators imports HydroSheaf
inference code; each was verified independently before any locked-test case
was scored.

### S1.3 M8: three controlled-synthetic truth mechanisms

The matched-model transport experiment uses an analytical one-dimensional
advection-dispersion forward model with declared synthetic truth
(dispersivity 2.0 m, first-order decay 0.005 d⁻¹, domain 10 m, velocity
0.1 m d⁻¹, unit source concentration, observation error
σᵢ = max(0.01, 0.03|Cᵢ|)). The independent-model experiment regenerates truth
with an implicit finite-volume/upwind numerical solver implemented only in
the benchmark script (240 grid cells, accepted only if its maximum deviation
from a 480-cell reference run is within 0.25 declared observation standard
deviations); the *same* production analytical adapter is then calibrated
against this independent truth, so any divergence between matched-model and
independent-model results is attributable to model-form mismatch, not to a
different calibration engine. The kinetic experiment drives the production
PHREEQC kinetic adapter against declared truth (calcite rate constant
k = 1×10⁻¹⁰ mol m⁻² s⁻¹, reactive surface area A = 0.1 m² L⁻¹, concentration
error 0.005 mmol L⁻¹). The frontier active-learning experiment uses a fourth,
independent MODFLOW 6/MODPATH 7 heterogeneous-aquifer generator with
nonlinear synthetic geochemistry, 256 topology particles per case, and three
model-discrepancy scenarios (nominal, separation-stress, noise-stress)
combined with a robust acquisition weight of 0.75 on the worst scenario.

## S2. Extended methods: equations

### S2.1 Mechanism-resolution score and resolution class (M6)

`M6`'s mechanism-resolution score (MRS) is a continuous 0-100 index computed
from the evidence-gated sparse-inversion fit at each well-tier, combining
normalised reconstruction residual, bootstrap support stability across
64 resamples, and thermodynamic-bound consistency. A well-tier is classified
`non_identifiable` when no candidate reaction family clears the evidence
gate at a fixed confidence threshold, `partially_identifiable` when one or
more families clear the gate but an equivalence class of two or more
families remains indistinguishable, and `identifiable` when a single family
is uniquely supported (this class does not occur in the Northern Ghana panel
at any tier; Table S2). Full source: `M5/m5_inverse_reaction_benchmark/
scripts/m5_common.py::fit_inverse`, reused unmodified by `M6`.

### S2.2 Posterior edge entropy (M7)

For an inferred directed edge *e* with posterior existence probability
*p*(*e*), mean edge entropy over a case's candidate edge set *E* is

H = −(1/|E|) Σ_{e∈E} [p(e) log p(e) + (1−p(e)) log(1−p(e))]

reported in nats. A reduction in H indicates the posterior has moved away
from 0.5 (maximum uncertainty) toward 0 or 1 for more candidate edges,
irrespective of whether that movement is toward the edge's true label.

### S2.3 Ranking and calibration metrics (M7, M8 frontier)

PR-AUC is the area under the precision-recall curve for the binary
true/false edge-existence label, computed by trapezoidal integration over
the full threshold range [@davisgoadrich2006prcurves]. Brier score is the
mean squared difference between predicted probability and the binary
outcome [@brier1950verification]. Log-loss is the mean negative
log-likelihood of the true label under the predicted probability, clipped at
a small epsilon to avoid unbounded loss for exact 0/1 predictions.

### S2.4 Whitened Fisher information and linearised coverage (M8)

Sensitivities are central finite differences of the forward model with
respect to log10(parameter) (step 1×10⁻⁴ log10 units), divided by the
declared Gaussian observation-error standard deviation, giving a
noise-whitened Jacobian **J**. The Fisher information matrix is
**F** = **J**ᵀ**J**; its condition number, minimum eigenvalue, and numerical
rank are reported directly. Linearised 95% marginal parameter uncertainty is
±1.96 × the square root of the corresponding diagonal of **F**⁻¹ (or its
SVD pseudoinverse when **F** is singular), in log10-parameter space.
Empirical coverage at a given design is the fraction of replicates for which
the known true log10-parameter value falls inside this nominal 95% interval;
because the interval is linearised and truth is known exactly in every
experiment, this is a direct, unconfounded test of whether the local
Gaussian approximation used for uncertainty reporting is honest at that
design, not an estimate of the interval's field performance.

### S2.5 Bootstrap intervals

`M6` reports bootstrap support stability from 64 within-well resamples
(no case-level resampling; the field-transfer point estimates in Table 2 of
the main text carry no analogous interval and are reported as point values,
consistent with `M6`'s own convention). `M7` reports 95% intervals from
10,000 case-block paired bootstrap resamples (`M7.3`) or 10,000 paired
case-block bootstrap replicates (`M7.4`/`M7.5`), resampling whole independent
cases rather than individual edges, so that within-case correlation does not
inflate apparent precision. `M8` reports 2,000-resample percentile bootstrap
intervals for the fixed-design sweep and 5,000-resample paired bootstrap
intervals for the optimal-design, independent-model and frontier
active-learning contrasts. No interval in this paper is computed by O4's own
harmonization layer; all are read directly from each component's own locked
result file (`O4/DECISIONS.md` D2).

## S3. Full taxonomy table (Table S1)

| Component | Experiment | Stress axis | Internal signal | External signal | Neg./structural control | Synthetic ground truth | Source |
|---|---|---|---|---|:---:|:---:|---|
| M6 | Tier ablation T4→T0 | data limitation | mean MRS | verified identifiability class | yes | yes | `results/m6_tier_ablation.csv`; `results/m6_field_gate_structural.csv` |
| M6 | Edge-set perturbation | model-form (connectivity) | per-edge MRS | network process-composition TV distance | no | no | `results/m6_edge_sensitivity.csv` |
| M6 | External sparse transfer | data limitation + domain shift | mean MRS | % non-identifiable | no | no | `results/m6_external_summary.csv` |
| M6 | Synthetic validation | model-form | exact-mineral F1 | true dominant-process recovery | no | yes | `results/m6_synthetic_recovery_by_tier.csv` |
| M7 | Native evidence-panel integration | evidence combination | entropy reduction | PR-AUC/Brier/log-loss | yes | yes | `results/m7_3_locked/evidence_case_bootstrap_contrasts.csv` |
| M7 | Adverse controls | evidence misspecification | entropy reduction | PR-AUC/Brier/log-loss | self | yes | `results/m7_3_locked/evidence_case_bootstrap_contrasts.csv` |
| M7 | M7.4 sheaf vs graph | model-form | log-loss/ECE | PR-AUC/selected-F1 | yes | yes | `results/RUN-M7-SHEAF-VS-GRAPH-20260729-01/paired_bootstrap_contrasts.csv` |
| M7 | M7.5 robust hybrid | evidence misspecification | Brier/log-loss | PR-AUC | yes | yes | `results/RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01/locked_test/paired_bootstrap_contrasts.csv` |
| M8 | Fixed 16-design sweep | data limitation | success rate; coverage | recovery error | no | yes | `manuscript/artifacts/m8_transport_parameter_summary.csv` |
| M8 | Optimal design, independent truth | model-form | coverage | recovery error | yes | yes | `results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv` |
| M8 | Kinetic structural confound | data limitation | convergence; objective value | Fisher rank | yes | yes | `manuscript/artifacts/m8_kinetics_structural_summary.csv` |
| M8 | Frontier active learning | evidence-value | entropy reduction/cost | Brier/PR-AUC | yes | yes | `provenance/runs/RUN-M8-FRONTIER-AL-20260728-01/strategy_summary.csv` |

## S4. Extended within-component detail

### S4.1 M6: edge-set perturbation and external-transfer detail

Three Hydrosheaf-generated candidate edge sets (chemistry-similarity kNN,
geographic-nearest, random-perturbed) were compared against each other; the
network-level process composition shifts materially with the assumed edge
set (total-variation distance versus the chemistry-kNN reference:
geographic-nearest 0.12, random-perturbed 0.05), but per-edge identifiability
class is edge-invariant (99.3-100% partially identifiable across all three
sets; mean MRS 72.6-73.3), so connectivity assumptions in this layer affect
network-level attribution, not point-reaction identifiability — the same
internal-signal-invariance-under-a-varying-external-condition pattern
documented for the tier ladder, at smaller magnitude. External transfer at
matched tiers: Talensi (Tier 1) shows 37.2% non-identifiable at mean MRS
68.8; the matched Northern Ghana Tier-1 reference shows 54.2% non-
identifiable at a *higher* mean MRS of 70.9. Lower Anayari (Tier 2) shows
95.3% non-identifiable at mean MRS 69.9; the matched Northern Ghana Tier-2
reference shows 53.3% non-identifiable at mean MRS 73.1. In both external
comparisons, the dataset with the *higher* internal signal (mean MRS) is the
Northern Ghana reference, which also has the *lower* external non-
identifiable fraction in the Tier-1 case but is directionally uninformative
in the Tier-2 case, where the gap in mean MRS (73.1 vs. 69.9, a 3.2-point
difference) does not scale with the gap in non-identifiable fraction (53.3%
vs. 95.3%, a 42-point difference).

### S4.2 M7.5 and M7.6: two further instances

`M7.5`'s locked 128-case test found the selected local-first/global-fallback
hybrid improved PR-AUC over the edge-local baseline (+0.0200, 95% CI +0.0073
to +0.0324) while its Brier-score interval crossed zero (−0.00151, −0.00419
to +0.00105) and mean log-loss was slightly worse (+0.00333, −0.00341 to
+0.01009); the prespecified overall-superiority gate, requiring all three
primary outcomes to improve, failed even though the single most commonly
reported outcome (PR-AUC) did improve. `M7.6`'s locked mechanism-diagnostic
run found that a declared synthetic shared nuisance increased full
nuclear-panel infeasibility (+0.2882, 95% CI +0.2118 to +0.3646) — the
internal, model-side signal moved in the direction consistent with the
hypothesis under test — but the predeclared CFC-12 specificity control moved
at least as strongly as the CFC-11 contrast it was meant to isolate (+0.7396
versus +0.7188), so the selective-mechanism claim the internal signal seemed
to support was not confirmed once the external, predeclared specificity
check was applied.

### S4.3 M8: frontier active-learning detail beyond the main text

Across 24 untouched cases, the robust information-and-decision policy was
actionable in every case (rate 1.0) with mean candidate recall 0.9861. All
120 selected actions across all 24 cases were chemistry panels; age-tracer
and directed-connectivity actions were available under the declared cost
structure but never selected, so the policy's multimodal-selection capacity
is unexercised by this run and is reported as a capability, not an empirical
contribution. Against the strong legacy uncertainty-chemistry policy, the
robust policy was noninferior within a frozen ±0.01 margin on both Brier
score (difference +0.00148, 95% CI 0.00000 to 0.00582) and entropy-reduction-
per-cost (−0.00422, −0.00955 to 0.00000, a narrow pass with 0.00045 nats of
margin remaining) — a result this paper reports as a passed decision gate,
not as evidence of superiority over the legacy heuristic.

## S5. Full benchmark-scale table (Table S4, unabridged)

| Component | Design unit | Primary scale | Replicate/bootstrap count | External reference | Field/archive transfer |
|---|---|---|---|---|---|
| M6 robustness | wells × evidence tiers | 160 × 5 = 800 rows | not applicable (deterministic, seed 1234) | independent CBE QC; synthetic extended model | Northern Ghana (160), Talensi (63), Lower Anayari (41) |
| M7 identifiability (M7.3) | independent MODFLOW6/MODPATH7 cases | 6 dev + 12 locked-test | 50,000 particles/case-tracer; 64 reaction bootstraps/case; 10,000 case-block bootstrap | fresh independent generator | Northern Ghana workbook, readiness audit only |
| M7 identifiability (M7.4) | independent scalar-section cases | 64 held-out, 4 strata | 10,000 paired case-block bootstrap | independent sheaf/graph generator | none |
| M7 identifiability (M7.5) | independent scalar-section cases | 64 dev + 128 locked-test | 10,000 paired case-block bootstrap | independent sheaf/graph generator | none |
| M7 identifiability (M7.6) | independent nuisance-mechanism cases | 6 dev + 12 locked-test | 10,000 paired case-block bootstrap | independent generator, declared nuisance | none |
| M8 calibration (matched-model) | calibrations | 4,000 fixed + 4,500 optimal-design = 8,500 | 250 paired replicates/design; 2,000-resample bootstrap | known synthetic truth | none |
| M8 calibration (independent-model) | calibrations | 1,000 locked-test | 80 dev replicates; 5,000 paired bootstrap | independent finite-volume/upwind solver | none |
| M8 calibration (kinetic) | designs | 3 (single-time, multi-time, +surface area) | not applicable (deterministic Fisher analysis) | known synthetic truth | none |
| M8 calibration (frontier AL) | independent MODFLOW6/MODPATH7 cases | 24 locked-test, ≤5 actions each | 5,000 paired case-bootstrap | independent generator, nonlinear geochemistry | none |

## S6. Extended staleness-verification detail

The main text (Section 3.1) and `O4/DECISIONS.md` D3 summarise a check
traced at the level of individual Python import statements. In full: the
three post-lock commits (`2bd5db0`, `8718d66`, `2d4b8af`) changed
`hydrosheaf/inference/network_fit.py` (the `infer_edges` function gained a
new `method="null_aware_sheaf"` branch and two new optional parameters, with
the pre-existing `max_neighbors` default behaviour explicitly preserved by a
`1 if max_neighbors is None else max_neighbors` guard; the body of
`fit_network`, which `M7.3` imports, is textually unchanged in the diff),
`hydrosheaf/sheaf/topology_refine.py` (the `refine_edges_with_sheaf`
function gained a new candidate-annotation and posterior-selection path,
entirely gated behind `config.topology_posterior_enabled`, default `False`),
`hydrosheaf/models/reactions.py` (the `build_reaction_dictionary` function's
concentration-based phase-pruning gate changed from treating a missing
SO4/Cl/F/NO3 measurement as zero, which could silently trigger pruning, to
skipping the pruning check when the measurement is genuinely absent),
`hydrosheaf/config.py` (thirty added lines, all new fields with declared
defaults and validation, no existing field's default changed),
`hydrosheaf/api.py`, `hydrosheaf/calibration/benchmark.py`,
`hydrosheaf/calibration/validation_workflow.py`, and the `hydrosheaf/
validation/` package. Repository-wide grep restricted to
`M6/.../scripts/*.py`, `M7/.../scripts/*.py` and `M8/.../scripts/*.py`
confirmed that `models/redox.py`, `calibration/benchmark.py`,
`calibration/validation_workflow.py`, and every `hydrosheaf/validation/`
module are imported by none of the three components' benchmark scripts;
that `models/reactions.py`'s changed pruning logic is imported only by
`M5`'s synthetic benchmark (fully-specified, no missing ions, so inert) and
by two `M6` scripts (`run_m6_chemistry_robustness.py`,
`run_m6_chemistry_stress_tests.py`) explicitly excluded from `M6`'s own Q1
submission pipeline per `M6/README.md`; and that the one script found to
import a changed module for a cited result (none — the sole hit,
`M7/.../run_m7_system_acceptance.py` importing `hydrosheaf.api`, produces a
result excluded from every claim in this paper). Three headline numbers were
independently recomputed from raw locked CSVs as a cross-check against
possible transcription drift between a component's own prose summary and its
underlying data: `M6`'s tier-ablation table (`derive_robustness_gradient.py`,
reading `results/m6_tier_ablation.csv` directly) reproduced
`docs/m6_results_summary.md` exactly; `M7.3`'s native and adverse-control
contrasts (`derive_integration_value.py`, reading
`evidence_case_bootstrap_contrasts.csv` directly) reproduced
`docs/m7_3_results.md` exactly; `M7.4`'s sheaf-versus-graph contrasts
reproduced `M7/README.md`'s summary exactly. No discrepancy was found.

## S7. Reproducibility appendix

All five Python harmonization scripts and six R figure scripts referenced in
this Supplementary Information and in the main text are archived under
`O4/analysis/` and run end to end from the repository root:

```
.venv/Scripts/python.exe O4/analysis/python/derive_taxonomy.py
.venv/Scripts/python.exe O4/analysis/python/derive_robustness_gradient.py
.venv/Scripts/python.exe O4/analysis/python/derive_integration_value.py
.venv/Scripts/python.exe O4/analysis/python/derive_calibration_gap.py
.venv/Scripts/python.exe O4/analysis/python/derive_benchmark_scale.py
"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" O4/analysis/r/make_all_figures.R
```

Each Python script prints the row count and destination path of every CSV it
writes; each R figure script prints the pixel dimensions of every figure it
saves. No script accepts a random seed or configuration argument, because
none performs new inference; every number is either passed through from a
named source file unchanged or computed by a disclosed, auditable arithmetic
operation (row selection, relabelling, or unit conversion) on already-locked
values.

## Supplementary references

Citation keys follow the main-text bibliography (`O4/manuscript/LITERATURE.bib`).
