# M7 — Conditional integration and non-uniqueness

This is the single maintained M7 code package for the Q1 manuscript. The
superseded original M7/M7.1 implementation has been removed. The independent
generator and confirmatory validation formerly developed as M7.2 are retained
here only because they are direct dependencies of the final benchmark and its
supporting evidence.

M7.3 tests the thesis claim that integration is useful when independent
evidence is complementary, but can be redundant or harmful when an evidence
stream is weak or misspecified. It does not require every additional stream to
increase point accuracy.

The benchmark has seven linked audits:

1. Fresh independent MODFLOW 6/MODPATH 7 cases compare hydraulic, age and
   chemistry evidence alone, in pairs and together.
2. Correct, partial and reversed flow graphs are imposed on the same local
   tracer likelihoods to quantify how connectivity assumptions change age
   uncertainty and accuracy.
3. Reaction-family recovery is bootstrapped under core and enhanced chemistry
   panels, with special reporting for carbonate processes.
4. The Northern Ghana workbook is audited as a data-limited application:
   supported component diagnostics are separated from unavailable field truth.
5. The strict public pipeline is exercised on fresh independent cases, with
   mandatory nuclear-age, sheaf-refinement, and network-fit stage assertions.
6. A prospectively locked M7.4 comparator isolates what affine restriction
   maps and a global section add beyond edge-local and identity-Laplacian
   weighted graphs.
7. A separately locked M7.5 diagnostic tests local-first/global-fallback
   hybrids, leave-one-edge-out false-edge downweighting, and separately
   cross-fitted versus shared calibration on fresh cases.

The final locked protocol is in
[docs/m7_3_protocol.md](docs/m7_3_protocol.md). Results must be written only
after the protocol and runner have been committed.

The locked interpretation and decision table are in
[docs/m7_3_results.md](docs/m7_3_results.md).

## Package layout

- `scripts/run_m7_3_nonuniqueness.py`: final conditional-integration benchmark.
- `scripts/independent_modflow_generator.py`: truth generator used by the final
  benchmark; it intentionally imports no HydroSheaf code.
- `scripts/strong_inference.py`: blind HydroSheaf inference adapter.
- `scripts/run_supporting_validation.py`: confirmatory supporting experiment,
  retained for reproducibility rather than as a separate M7 version.
- `results/m7_3_locked`: final locked M7 results.
- `results/supporting_validation`: confirmatory evidence used by M7 and M6.
- `scripts/run_m7_sheaf_vs_graph.py`: competence-matched M7.4 comparator.
- `scripts/independent_sheaf_graph_generator.py`: code-independent M7.4
  scalar-section generator.
- `results/RUN-M7-SHEAF-VS-GRAPH-20260729-01`: locked M7.4 results and claim
  decision.
- `scripts/run_m7_robust_hybrid_sheaf.py`: two-stage M7.5 development and
  one-time confirmatory runner.
- `results/RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01`: frozen development
  settings and the single 128-case locked-test result.

The supporting experiment's
[protocol](docs/supporting_validation_protocol.md),
[fresh-seed amendment](docs/supporting_validation_amendment.md), and
[results](docs/supporting_validation_results.md) remain as the audit trail for
the negative incremental-age finding.

## Replay the final benchmark

Place `mf6.exe` and `mp7.exe` in `.codex_work/modflow-bin`, then run:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_m7_3_nonuniqueness.py
```

A non-claim-bearing smoke run is available with `--quick`.

## Replay the supporting validation

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_supporting_validation.py --confirmatory
```

To regenerate Supplementary Figure S1 for representative locked realization
4101:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' M7\m7_nonuniqueness_benchmark\scripts\make_supporting_model_domain_map.R
```

This figure uses synthetic model-space metres. It is not georeferenced and
must not be presented on a Ghana or other geographic basemap.

## Public-pipeline system acceptance

`RUN-M7-SYSTEM-20260728-01` exercises the public `fit_network_pipeline` on
six fresh cases from the independent MODFLOW 6/MODPATH 7 generator. The run
uses strict nuclear-age, sheaf-refinement, and network-fit stage assertions and
compares hydraulic-only, local-age, full-sheaf, and age-permuted conditions.

The system gate passed: every requested stage completed, mean candidate recall
was 0.9815, all metrics were finite, and generator independence was verified.
The full global sheaf improved PR-AUC over local age alone by 0.0586 (95% paired
bootstrap interval 0.0386 to 0.0777), but it had lower PR-AUC than the hydraulic
baseline and did not beat the age-permuted adverse control. The result therefore
supports verified controlled-synthetic system execution and a conditional
global-section increment, not overall topology superiority.

The locked protocol is
[`M7_SYSTEM_ACCEPTANCE_PROTOCOL.md`](M7_SYSTEM_ACCEPTANCE_PROTOCOL.md), and the
runner is `scripts/run_m7_system_acceptance.py`.

## What the sheaf contributes beyond a weighted graph

M7.4 compares an edge-local weighted graph, an ordinary identity-coupled graph
Laplacian, the HydroSheaf affine section solver, and a within-case permuted-map
control under equal inputs, calibration form, optimisation routine, iteration
count, threshold budget, and locked cases. The execution/equivalence gate
passed on 64 held-out cases. The affine sheaf improved PR-AUC over the identity
Laplacian by 0.0854 (95% CI 0.0666 to 0.1050) and Brier score by -0.0193
(-0.0235 to -0.0152), exactly equalled it in the identity limit, and beat the
permuted-map control. It did not beat the stronger edge-local weighted graph
overall in both primary outcomes (PR-AUC +0.0097, CI -0.0054 to +0.0248;
Brier +0.00053, CI -0.00330 to +0.00443).

The allowed conclusion is conditional: non-identity restriction maps and the
global section improve representation and localisation of planted incompatible
cycles, but the tested sheaf is not generally superior to a strong local graph
score. See [`M7_SHEAF_VS_GRAPH_PROTOCOL.md`](M7_SHEAF_VS_GRAPH_PROTOCOL.md),
[`docs/m7_sheaf_vs_graph_results.md`](docs/m7_sheaf_vs_graph_results.md), and
[`M7_SHEAF_VS_GRAPH_CLAIM_DECISION.md`](M7_SHEAF_VS_GRAPH_CLAIM_DECISION.md).

M7.5 then tested whether the M7.4 heterogeneous-affine calibration loss arose
from discarded local evidence, candidate self-influence, or calibration. On
128 fresh locked cases, the selected local-first/global-fallback hybrid
improved PR-AUC over edge-local by +0.0200 (95% CI +0.0073 to +0.0324), but
the Brier interval crossed zero and mean log loss was slightly worse, so the
strict overall-superiority gate failed. Restoring local evidence significantly
improved the original global estimator, whereas leave-one-edge-out
robustification worsened Brier score and log loss. Conditional PR-AUC gains
were supported for incompatible cycles and noisy/missing observations, and
native maps beat the permuted control on all three primary outcomes. See
[`M7_ROBUST_HYBRID_SHEAF_PROTOCOL.md`](M7_ROBUST_HYBRID_SHEAF_PROTOCOL.md),
[`docs/m7_robust_hybrid_results.md`](docs/m7_robust_hybrid_results.md), and
[`M7_ROBUST_HYBRID_CLAIM_DECISION.md`](M7_ROBUST_HYBRID_CLAIM_DECISION.md).

Replay the locked comparator with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_m7_sheaf_vs_graph.py --overwrite
```

Regenerate Figure 6, main Table 3, the complete unadjusted M7.4 contrast table
(Table S7), the family-wise table (Table S10), and the public-pipeline
acceptance table (Table S8) from immutable result
tables only with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_sheaf_vs_graph_assets.py
```

The confirmatory M7.5 test must not be replayed after viewing its result.
Regenerate only Figure 7, main Table 4 and Tables S9/S11 from the immutable
locked-test tables with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\make_m7_robust_hybrid_assets.py
```
