### Estimator clarification and reviewer-requested sensitivities

PR-AUC was computed on all candidate edges pooled within each test case and
then contrasted by resampling the 12 independent cases as blocks. This preserves
the edge-ranking estimand while preventing edges within a generated aquifer
from being treated as independent replicates. The reported interval is the
percentile interval of the paired case-block contrast; direct pooled-edge
contrasts are supplied only as descriptive checks. Posterior graph sampling by
generative flow networks provides a relevant modern alternative
[@deleu2022daggflownet], but this benchmark intentionally evaluates a fixed,
auditable candidate-edge scoring task rather than full graph-posterior learning.

The effective-sample-size threshold was varied from 200 to 400 and 1,000. All
conditions remained stable in every case except reversed topology: under
tritium only, the stable-case fraction changed from 0.583 to 0.333 and 0.083;
under two tracers it changed from 1.000 to 1.000 and 0.833. Thus, the qualitative
diagnosis of reversed-topology incompatibility is stable and becomes more
conservative as the threshold rises.

Logistic L2 regularisation was varied over C = 0.1, 0.3, 1, 3, and 10. Native
HAC PR-AUC changed narrowly from 0.477 to 0.486; Brier score improved from
0.096 to 0.083 and log loss from 0.305 to 0.247. Standardised H/A/C
coefficients are archived for every C in
`results/m7_3_locked/review_sensitivity.csv`, so the stability claim is not
based on performance values alone.

The probability-span conflict threshold was varied from 0.10 to 0.50. At 0.10
it flagged 39.8% of native edges and 59.5%-61.2% under adverse conditions, but
its error rate remained condition-dependent; thresholds of 0.30 or higher
flagged no edges. This confirms that the originally predeclared univariate
diagnostic is not a reliable standalone safeguard. The 3% reaction-noise level
and five-year topology order scale remain predeclared fixed settings, not
validated universal constants, and are retained as explicit process-based-generator
limitations.

A reviewer-requested case-block bootstrap placed the enhanced-minus-core modal
family accuracy difference at 0.0556 (95% CI 0.0278-0.0833; 12 cases). Under
tritium-only evidence, reversed topology also produced a falsely narrow mean
interval relative to complete topology (-3.394 years; 95% CI -5.777 to
-1.069), while failing the ESS diagnostic in most cases. Values in prose are
rounded only after analysis; CSV outputs retain full precision.

Recent hydrochemical metamodelling has used machine-learning ensembles to
estimate groundwater-age distributions [@tschritter2023age]. That work predicts
ages from chemistry, whereas the present benchmark tests whether adding an
evidence stream improves or harms a separately scored topology, age, or
reaction inference under adverse controls.
