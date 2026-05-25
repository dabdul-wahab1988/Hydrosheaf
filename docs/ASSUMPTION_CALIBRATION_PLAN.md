# Assumption Calibration Plan

## Document Status

| Field | Value |
|---|---|
| Status | Revised draft after assumption review |
| Date | 2026-05-24 |
| Target | Hydrosheaf v0.6.x or next minor release |
| Purpose | Convert Hydrosheaf topology inference from implicit assumption scoring into explicit, falsifiable, calibrated edge evidence |

## Review Conclusion

Hydrosheaf's major assumption flaw is not a single parameter. It is the way topology evidence is interpreted:

> Hydrosheaf currently turns hydrogeologic compatibility into connectivity evidence without first testing whether the same observations are expected under non-connectivity.

This creates a false-positive risk. Similar chemistry, similar isotope signatures, short distance, common layer membership, or monotonic age can all occur without direct flow between two wells. They may reflect shared recharge, common lithology, regional gradients, dispersion, sampling bias, long screens, pumping effects, or the same end-member mixture.

The current code exposes this issue in three places:

- `hydrosheaf/graph3d/build_3d.py` multiplies `P_head * P_dist * P_layer * P_screen`, which treats correlated variables as separable evidence.
- `hydrosheaf/sheaf/topology_refine.py` penalizes isotope, chloride, age, and global-section mismatch, but it has no explicit no-flow/null model. Missing evidence can also contribute zero cost, which may look neutral rather than uncertain.
- `hydrosheaf/nuclear/network_aging.py` constrains scalar mean ages to increase downstream, which can fail in mixed, fractured, bypass, or long-screen systems where the observed age is a distribution.

The practical advice is: do not claim that Hydrosheaf "proves connectivity" from sparse chemistry or graph scores. Claim that it proposes candidate edges whose assumptions are scored, falsified, or validated against independent evidence.

## Revised Scientific Position

Hydrosheaf should be positioned as an assumption-audited topology inference framework:

1. Candidate edges are hypotheses, not facts.
2. Edge probabilities are not posterior probabilities unless calibrated against independent labels.
3. Chemical and isotope similarity must be compared with a no-flow explanation before supporting connectivity.
4. MODPATH-assisted priors and MODPATH validation must remain separated. A benchmark used to train or tune an edge model cannot also be reported as independent validation.
5. Sparse-data outputs should include evidence class, uncertainty, and failure reason rather than only a selected edge list.

## Priority Fixes

### 1. Evidence Ladder

Add edge-level classes before adding more model complexity.

| Class | Meaning |
|---|---|
| `FALSIFIED` | Rejected by null model, benchmark, tracer contradiction, or impossible head/age relation |
| `AMBIGUOUS` | Insufficient independent evidence |
| `PRIOR_ASSISTED` | Supported mainly by imported model or MODPATH prior |
| `PROBABLE` | Multiple independent evidence streams support the edge and no strong null explanation dominates |
| `VALIDATED` | Supported by an independent benchmark, tracer test, or field-confirmed topology |

This should live in `hydrosheaf/validation/claims.py` or a nearby validation module and be written into every final edge's attributes.

### 2. Null Models for Non-Connectivity

Implement explicit null penalties before treating chemistry or isotopes as edge support.

Minimum null explanations:

- common lithology without direct flow;
- shared recharge or end-member mixture;
- regional evaporation/isotope trend;
- spatial autocorrelation or dispersion-like similarity;
- common anthropogenic source, especially nitrate and chloride.

The first implementation should be simple and testable: compute `P(observed similarity | no direct edge)` and downgrade or reject the edge when the no-flow explanation is more plausible than the flow explanation.

### 3. Calibrated Joint Edge Model

Replace the independent product rule with a calibrated joint model only after labels are available.

Candidate features:

- hydraulic-head difference and uncertainty;
- 3D distance and anisotropy;
- layer separation and screen overlap;
- chemistry/isotope residuals after null correction;
- residence-time compatibility;
- geological barriers or lithology tags;
- pumping or seasonal-head flags.

Use logistic regression, probit, or another transparent probabilistic classifier first. Reserve more complex Bayesian/copula models until benchmarks show the simple model is insufficient.

Validation must use held-out labels by site, archive, or scenario. Random edge-level splits are not enough if neighboring edges share the same hydrogeologic structure.

### 4. Distribution-Aware Age Constraints

Replace scalar monotonic-age penalties with distribution compatibility.

Required behavior:

- store each well's transit-time distribution or sampled posterior, not only mean age;
- test whether the downstream distribution could be produced by upstream water plus a physically plausible travel-time distribution;
- treat reversals as evidence of mixing, bypass, wrong direction, or ambiguous topology rather than automatic failure.

This is especially important for fractured, karst, crystalline, long-screen, and heavily pumped aquifers.

### 5. Global Flow Consistency Audit

Add a graph-level consistency check after local edge scoring.

The first useful version does not need a full mathematical Hodge implementation. It should detect:

- high-confidence cycles inconsistent with head gradients;
- edges that create most of the non-gradient residual;
- fragile topology where small evidence changes alter selected paths.

If a Hodge decomposition is added later, report it as a topology audit metric, not as proof that the graph is physically true.

### 6. Dynamic Head Treatment

Well head should not be treated as aquifer head unless measurement conditions are known.

Add optional uncertainty inflation or correction for:

- pumping status and pumping rate;
- seasonal water-level variation;
- datum/elevation uncertainty;
- nearby production wells;
- stale or mixed-date measurements.

This can be implemented as a head-uncertainty model before building more elaborate dynamic-head priors.

## Implementation Order

| Phase | Work | Deliverable | Status |
|---|---|---|---|
| 0 | Evidence semantics and guardrails | Edge evidence classes; validation wording; no-flow null design | ✅ Complete (2026-05-24) |
| 1 | Null-model filtering | Chemistry/isotope/Cl null penalties; synthetic negative controls | ✅ Complete (2026-05-24) |
| 2 | Joint edge calibration | Trained and held-out calibrated edge model; Brier/ECE/AUROC metrics | Not started |
| 3 | Age and head correction | Distributional age compatibility; head uncertainty inflation | Not started |
| 4 | Global topology audit | Cycle/non-gradient diagnostics; class-wise benchmark reports | Not started |
| 5 | Active learning | Rank next measurements by expected reduction in edge uncertainty | Not started |

The first release should stop after Phase 1 or Phase 2 if validation does not show a clear reduction in false positives.

## Validation Requirements

| Risk | Required Test |
|---|---|
| Chemistry false positives | Synthetic no-flow pairs with same lithology, similar recharge, and similar chemistry |
| Correlated predictors | Compare product rule vs. joint model on held-out MODPATH labels |
| Benchmark leakage | Separate prior-assisted runs from independent validation runs |
| Age-model overconstraint | Test scalar age vs. distributional age on mixed and bypass-flow cases |
| Global inconsistency | Report cycle or non-gradient residual before and after pruning |
| Sparse evidence overconfidence | Report edge classes and uncertainty, not only precision/recall |

Core metrics should include precision, recall, F1, false-positive count, Brier score, calibration error, and class-wise precision. AUROC alone is insufficient because Hydrosheaf users need calibrated confidence, not only ranking skill.

## Documentation Changes

Update these files with the revised claim language:

- `docs/hydrosheaf_model_assumptions.md`: add no-flow null assumptions, evidence classes, and dynamic-head caveats.
- `docs/INPUTS_REFERENCE.md`: document optional fields for lithology, pumping status, screen interval, sample date, uncertainty, and age distributions.
- `M4/m4_topology_benchmark/docs/`: separate independent validation from prior-assisted or trained modes.
- `M3/m3_age_benchmark/docs/`: distinguish scalar apparent age, mean residence time, and full transit-time distribution.

Remove or verify paper-specific claims and exact expected performance gains before manuscript use. The previous draft included precise false-positive reductions and citation-specific claims that should not remain unless they are independently checked.

## Recommended Manuscript Claim

Use this conservative claim:

> Hydrosheaf infers candidate groundwater-flow topology under sparse data by scoring hydraulic, hydrochemical, isotopic, age, spatial, and geological evidence. The revised framework treats every edge as an assumption-bearing hypothesis: competing no-flow explanations are tested through null models, correlated predictors are calibrated with held-out reference labels, and final edges are reported by evidence class rather than as deterministic connectivity.

Avoid this claim:

> Hydrosheaf proves groundwater connectivity from chemistry or sheaf consistency.

## Immediate Advice

The next technical step should be the evidence ladder plus null-model filter, not the full seven-module expansion. That is the smallest change that directly addresses the major flaw. Without it, adding Hodge decomposition, active learning, or more Bayesian machinery may make the system look more rigorous while preserving the same overclaim: compatibility is being treated as causality.
