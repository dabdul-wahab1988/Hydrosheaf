**Abstract**

**Purpose.** The problem is that groundwater evidence integration is often
assumed to reduce interpretive non-uniqueness, although lower uncertainty does
not prove better inference, and graph-based connectivity prediction does not
isolate a sheaf layer's contribution. We tested both propositions.

**Design and methods.** Seven linked audits were conducted. The process-based
integration and public-pipeline audits used an independent MODFLOW 6/MODPATH 7
flow, tracer and nonlinear-chemistry generator; the competence-matched
representation benchmark and follow-up estimator diagnostic used a separate
scalar affine graph generator. Locked tests comprised 12, 6, 64 and 128 cases,
respectively. A descriptive Northern Ghana audit defined the field-data claim
boundary. The full representation-benchmark and estimator-diagnostic contrast
families were evaluated with
10,000-resample max-z simultaneous intervals.

**Primary representation result.** The competence-matched representation
benchmark passed identity-limit nesting
and the native affine sheaf outperformed the identity graph on PR-AUC (+0.0854,
simultaneous 95% CI [0.0539, 0.1169]) and the permuted-map control (+0.0909
[0.0571, 0.1246]). Planted-cycle conflict localisation also improved over the
edge-local graph (+0.0689 [0.0318, 0.1061]). However, overall PR-AUC versus
edge-local was +0.0097 [-0.0149, 0.0343], and neither Brier score nor log loss
improved. The representation benchmark failed the prespecified complete gate.

**Follow-up result.** In the estimator diagnostic, development selected a local
weight of 1.0, so
the estimator was local-first/global-fallback. Against edge-local, differences
were +0.0200 [-0.0055, 0.0454] for PR-AUC, -0.00151 [-0.00691, 0.00389] for
Brier score and +0.00333 [-0.0102, 0.0168] for log loss. The estimator diagnostic failed the
prespecified complete gate. Native maps nevertheless beat the permuted-map
control on all three outcomes, and conflict localisation survived
correction.

**Scope and significance.** In the process-based integration benchmark,
chemistry improved the topology-ranking
outcomes, correct topology reduced age MAE by 0.062--0.164 years, adverse
controls reduced uncertainty while worsening skill, and carbonate reactions
were not recovered under either tested indicator panel. These controlled-
synthetic benchmarks support a conditional representation claim: affine maps
encode non-identity relations and global compatibility supplies a conflict
diagnostic and missing-endpoint fallback. They do not establish general
predictive superiority, field validity, or performance for temporal,
three-dimensional, vadose-zone, vector-stalk or active-learning capabilities.
