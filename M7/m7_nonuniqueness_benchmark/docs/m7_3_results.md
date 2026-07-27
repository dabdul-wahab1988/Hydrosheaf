# M7.3 locked result — conditional integration and non-identifiability

## Confirmatory status

The protocol and executable analysis were committed as `d336e87` before any
declared 5301-series test case was generated. The locked run used:

- six development cases (seeds 5201–5206);
- twelve untouched MODFLOW 6/MODPATH 7 test cases (5301–5312);
- 50,000 age-importance particles per case and tracer regime;
- 64 reaction bootstraps per case; and
- 10,000 case-block bootstrap replicates.

Candidate generation retained every true test edge. All local age analyses met
their declared convergence rule. Synthetic truth is model-conditioned and is
not evidence about a Ghanaian aquifer.

## Decision summary

| Question | Locked result | Decision |
|---|---|---|
| Does age add topology-ranking value beyond hydraulics and chemistry? | HAC minus HC PR-AUC = **−0.0060**, 95% CI [−0.0122, −0.0011]. Entropy fell by only 0.00062; Brier and log-loss changes included zero. | **No positive incremental age value.** The small entropy reduction does not compensate for poorer ranking. |
| Does chemistry add topology-ranking value? | HAC minus HA PR-AUC = **+0.447**, Brier = −0.0196 and log loss = −0.0791; all 95% intervals excluded zero. | **Strong positive contribution in this generator.** |
| Do hydraulics add topology-ranking value? | HAC minus AC PR-AUC = **+0.0091**, Brier = −0.00100 and log loss = −0.00593; all 95% intervals excluded zero. | **Small but calibrated positive contribution.** |
| Does correct topology improve age estimation? | Correct graph minus no graph MAE = **−0.0619 y** with 3H+39Ar and **−0.1642 y** with 3H only; both intervals excluded zero. Intervals narrowed by 0.252 y and 0.912 y, respectively, without a coverage decrease. | **Yes, conditionally, when topology is correct.** |
| Can a wrong graph harm age inference? | For 3H+39Ar, reversed minus correct MAE = **+0.282 y** and coverage = −0.0278; both intervals excluded zero. In 3H-only runs, reversed-graph importance sampling failed the ESS rule in 8/12 cases. | **Yes.** The ESS collapse is an incompatibility warning; unstable reversed-graph age summaries are not inferential results. |
| Does the enhanced chemistry panel resolve reactions? | Overall modal-family recovery rose descriptively from 0.556 to 0.611. Carbonate weathering and precipitation recovery remained **0/36 edges** under both panels. | **Only partial family recovery; carbonate remains non-identifiable/misclassified.** |
| Does lower uncertainty always mean better integration? | Permuted streams reduced entropy while worsening PR-AUC, Brier, log loss and overconfident error. | **No.** Misspecification can create false confidence. |

## Evidence-panel interpretation

On native test data, chemistry was the strongest individual edge-ranking
stream (PR-AUC 0.459), hydraulics was weaker (0.176), and age alone was below
the useful ranking level in this candidate set (0.111). The best pair was
hydraulics plus chemistry (HC, PR-AUC 0.485). Adding age produced HAC PR-AUC
0.480 rather than a gain.

The adverse controls are an important positive result for the framework's
guardrails. Permuted age reduced mean edge entropy by 0.0207 but reduced
PR-AUC by 0.0754 and increased Brier score by 0.00341. Permuted hydraulics
reduced entropy by 0.0482 while reducing PR-AUC by 0.0686 and increasing log
loss by 0.0745. Joint misspecification reduced entropy most strongly
(−0.0706) while reducing PR-AUC by 0.139 and raising overconfident error by
0.0387. Thus uncertainty reduction is not accepted without predictive and
calibration checks.

The predeclared univariate probability-span conflict flag identified no edges.
That threshold was insensitive even though the case-level negative-transfer
metrics detected misspecification clearly. It must be reported as a failed
conflict heuristic, not silently replaced after seeing the test outcomes.

## Topology-to-age interpretation

Correct and partial true topology assumptions were numerically stable in all
cases. Correct topology produced small, calibrated improvements when the local
tracer likelihood was already informative and larger improvements when only
tritium was available. This is the defensible positive integration finding:
topology can reduce age non-uniqueness when it contributes correct information
that is not already supplied by the tracer panel.

Reversed topology is not merely an alternative estimate. With two tracers it
worsened MAE and coverage. With tritium alone it caused severe importance-weight
degeneracy in most cases. The latter is evidence of model–data incompatibility,
not evidence that the reversed graph gives narrower or better age estimates.

## Reaction interpretation

The enhanced panel improved aggregate recovery only descriptively; no
case-block interval was predeclared for the tier contrast. Denitrification,
sulfate reduction and silicate weathering were recovered reliably, while iron
reduction remained mixed. Carbonate weathering and precipitation received zero
probability as the true family on all relevant edges.

Some carbonate rows had very low support entropy because the model was
consistently confident in the wrong family. This is a central
non-identifiability lesson: bootstrap stability or low entropy cannot establish
mechanistic truth when the reaction dictionary or observation model is
misspecified.

## Ghana scope decision

The audited canonical Northern Ghana workbook (`data/FieldData/NorthenGhana/
NorthernGhana.xlsx`; an earlier revision audited a different, antecedent
study's own derived workbook instead, since removed — DECISIONS.md) has major
chemistry, Sr, SiO2, stable water isotopes and one static-water-level field
per well. It has no independent aquifer-type/geology/land-use classification,
no environmental age-tracer panel, screen intervals, repeated head series,
exact public coordinates, intra-season sampling dates, processed graph edges,
independently verified connectivity or reaction truth.

Consequently, the Ghanaian data support component diagnostics, seasonal
chemistry hold-forward under a disclosed arbitrary well-revelation order,
alternative-edge sensitivity, reaction-family plausibility/equivalence
classes and explicit non-identifiability mapping. They do not support field
residence-time validation, exact directed connectivity, screen-resolved flow,
unique reaction mechanisms or a fully observed field
digital twin.

The field evidence is synthesized in
[`objective6_data_limited_synthesis.md`](../../../M6/m6_field_transfer_benchmark/docs/objective6_data_limited_synthesis.md).

## Honest scientific conclusion

M7.3 supports a conditional-integration paper, not a universal-integration
claim. Chemistry and hydraulics provide complementary topology information in
the independent generator. Correct topology improves groundwater-age inference,
especially under a weak tracer panel. Age does not improve topology ranking in
the reverse direction, misspecified streams create false confidence, and
carbonate mechanisms remain unresolved.

This pattern is aligned with the revised Objective 6: the contribution is a
framework that distinguishes supportable interpretations from
non-identifiable ones rather than forcing a complete field reconstruction from
data that cannot identify it.
