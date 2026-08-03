# M3 accuracy-locked results summary

Locked on 2026-07-28. Published USGS LPM ages are model-derived reference
products, not error-free true ages. All reference-agreement metrics below use
only fits passing the reference-free identifiability/reportability guard.

## Reference-workflow agreement

- Strict reported-configuration emulation: 329 reportable fits from 1,272 rows;
  median absolute log10 discrepancy 0.027932; log10 RMSE 0.276882; log10 R2
  0.937147; within-factor-2 agreement 0.875380.
- Reported-output-constrained age-fraction sensitivity: 289 reportable fits;
  median discrepancy 0.021613; RMSE 0.196357; R2 0.969780; within-factor-2
  agreement 0.916955.
- The age-fraction result is not independent validation because reported
  fractions and reference ages share USGS LPM provenance.
- Hydrosheaf model selection: 309 reportable fits; median discrepancy 0.130446;
  RMSE 0.609790; R2 0.763909; within-factor-2 agreement 0.673139. Unequal support
  prevents a general ranking against the emulation modes.

## Graph-age benchmark

- Support: 329 strict reportable nodes.
- Smallest candidate RMSE change: -0.001349 log10 for weak hydraulic-proxy
  regularisation, while within-factor-2 agreement declined from 0.875380 to
  0.869301. It fails the joint robust-improvement rule.
- Weak wrong-direction and randomised controls increased RMSE by 0.100802 and
  0.654969 log10, respectively. This supports topology falsification, not
  confirmation of candidate flow paths.
- No MRVA graph-age replication is claimed because the required reported-LPM
  reference table is unavailable.

## Leakage-guarded target-tracer withholding

Every graph-node age was refitted without the target tracer before graph
regularisation.

- 3H: 121 reportable of 794 eligible rows. RMSE was 20.818104 TU at baseline,
  20.818055 hydraulic, 20.841494 depth, and 21.566597 randomised.
- SF6: 75 of 262. RMSE was 2.843988 pptv at baseline, 2.857993 hydraulic,
  2.928528 depth, and 2.928427 randomised.
- 14C: 169 of 1,103. RMSE was 26.830072 pmC at baseline, 26.829250 hydraulic,
  26.902373 depth, and 27.819135 randomised.
- CFC-11 and CFC-12 retained 4 and 6 reportable fits and zero eligible graph
  edges. Their graph effects are non-estimable, not equivalent or neutral.

No candidate graph provides a meaningful target-withheld predictive gain.

## Withdrawals and bounded claim

- The hierarchical old-water prior is withdrawn from active outputs because its
  same-release pooled priors are non-independent and its audit result was
  catastrophic (RMSE 1.310; R2 0.004).
- The gas-correction comparison is withdrawn because zero supported paired raw
  and corrected atmospheric-equivalent values remained.
- The controlled-synthetic ambiguity demonstration supports only a capability
  claim for known simulated flow ordering. It is not field validation.

Hydrosheaf M3 is a controlled public-data diagnostic and topology-falsification
benchmark. It does not establish true groundwater ages, verified groundwater
flow paths, universal graph benefit, or complete-framework validation.
