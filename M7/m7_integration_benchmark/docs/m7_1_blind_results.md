# M7.1 Blind Replicated Results

## Source-frozen run

The definitive run used Git commit
`0a9c16777cf630f84c86628fe2af3f945766f49f`, a clean pre-run source tree,
and 171 hashed runtime inputs. The Git commit and every runtime-input hash were
unchanged at the end of the 898.84-second run.

- 20 development, 100 disjoint test, and 8 heavy-audit aquifers completed.
- Candidate-universe recall was 1.000 on the test set.
- The truth-poisoning leakage test passed.
- No heavy-audit case raised an execution exception.
- Module provenance recorded 128 `infer_edges`, 264
  `refine_edges_with_sheaf`, and 136 `fit_network` calls, plus 8 calls each to
  the topology posterior, PHREEQC, and Bayesian network aging.

## Held-out Tier-A performance

Values are means across 100 aquifers. F1 parentheses are 95% aquifer-bootstrap
intervals from 2,000 resamples.

| Method | F1 | MCC | PR-AUC | Brier | ECE10 | SHD |
|---|---:|---:|---:|---:|---:|---:|
| Hydraulic-spatial | 0.659 (0.649–0.670) | 0.579 | 0.627 | 0.112 | 0.112 | 7.89 |
| Age-only | 0.314 (0.311–0.317) | -0.001 | 0.219 | 0.165 | 0.091 | 42.12 |
| Sheaf-multievidence | 0.355 (0.347–0.364) | 0.149 | 0.341 | 0.159 | 0.105 | 28.19 |
| Chemistry | 0.462 (0.445–0.480) | 0.318 | 0.422 | 0.145 | 0.104 | 14.85 |
| Equal-weight joint | 0.567 (0.552–0.583) | 0.461 | 0.527 | 0.134 | 0.107 | 10.94 |
| Development-trained logistic joint | 0.692 (0.680–0.704) | 0.619 | 0.700 | 0.098 | 0.111 | 6.61 |

For the logistic joint score, paired method-minus-hydraulic differences were
F1 +0.0332 (95% CI +0.0211 to +0.0456), MCC +0.0401 (+0.0243 to +0.0550),
PR-AUC +0.0721 (+0.0609 to +0.0830), and Brier -0.0133 (-0.0151 to
-0.0116; lower is better). The ECE10 difference was -0.0018 (-0.0072 to
+0.0033), so improved calibration error is not established.

Equal weighting was harmful relative to the hydraulic baseline: paired F1 was
-0.0922 (-0.1102 to -0.0754), MCC -0.1187 (-0.1422 to -0.0962), and
PR-AUC -0.1004 (-0.1152 to -0.0858).

The development-only logistic coefficients were +6.296 for hydraulic, +1.562
for chemistry, and -5.215 for age-only, with intercept -2.269. Thus the positive
fusion result is primarily a learned hydraulic/chemistry combination that
downweights the present age score; it is not evidence that the current age
component adds positive topology information. The age-only F1 remained poor in
both the tracer-informative subgroup (0.314) and the deliberately old,
tracer-uninformative subgroup (0.313).

## Held-out chemistry result

Strictly held-out Mg/SO4/Fe reaction-delta RMSE was 0.0487 mmol/L. Only 58.8% of
candidate-covered true processes belonged to a reaction family represented in
the fitted dictionary, and supported-family accuracy was 38.7%. Reaction-family
identification is therefore not validated. An expanded process dictionary and
an identifiability analysis are required, especially for sulfate reduction,
unmodelled mixing, and observationally equivalent reaction combinations.

## Tier-B execution and diagnostic audit

All eight heavy cases executed the intended sequential path:
PHREEQC-constrained refitting, posterior-MAP topology into Bayesian aging, and
posterior-age feedback into a second age-only sheaf pass. Execution is not the
same as convergence or material influence.

| Audit gate | Valid cases | Exact interpretation |
|---|---:|---|
| Sequential path executed | 8/8 | All required module calls and hand-offs occurred. |
| Posterior MAP structure valid | 8/8 | Each MAP graph met edge-count, DAG, weak-connectivity, out-degree, and root-reachability constraints. |
| Topology posterior converged | 0/8 | Confirmatory posterior probability claims are withheld. |
| PHREEQC execution/constrained fit valid | 8/8 | All samples succeeded; 413 constrained edge fits were performed. |
| Bayesian-age posterior converged | 0/8 | Confirmatory age MAE and interval coverage are withheld. |

Topology edge-count R-hat ranged from 1.002 to 1.037 and edge-count ESS from
125.4 to 166.1. Worst-edge R-hat ranged from 1.017 to 1.240 and minimum edge ESS
from 82.2 to 225.4. These fail the predeclared R-hat ≤ 1.01 and ESS ≥ 400 gates
despite acceptance rates of 0.222–0.382.

PHREEQC sample success was 100%, but the mean absolute change in the inverse-fit
objective was zero in seven cases and approximately `5.7e-14` in one. The
thermodynamic bounds were therefore non-binding in these synthetic cases.

Bayesian-age worst R-hat ranged from 1.0 to 2.3, minimum bulk ESS from 5.2 to
160, and the eight runs contained 660 divergences in total. Consequently, the
publication-valid age MAE and coverage fields are null. The exploratory,
non-confirmatory averages were 21.43 years MAE and 0.70 interval coverage. Age
feedback changed zero selected edges in all eight cases (Jaccard 1.0), so this
run demonstrates data flow but not a beneficial operational update.

## Defensible conclusion

Under this fixed synthetic generator, a development-trained logistic stacker
improved selected held-out topology metrics over hydraulic-spatial inference.
The result does not establish equal-weight integration, a positive contribution
from the current age score, unique reaction recovery, converged posterior
probabilities, useful Bayesian-age updating, field transfer, or a validated
digital twin.

The next model-development priorities are:

1. redesign tracer-age likelihoods for multimodal/pre-bomb ambiguity and verify
   age evidence adds out-of-sample information;
2. improve topology-posterior mixing before reporting marginal probabilities;
3. reparameterize Bayesian network aging to eliminate divergences and raise
   effective sample sizes;
4. broaden and test the reaction/process dictionary, with explicit
   non-identifiability reporting;
5. create PHREEQC stress cases whose thermodynamic constraints are active; and
6. perform external simulator and field prequential validation before using
   operational digital-twin language.
