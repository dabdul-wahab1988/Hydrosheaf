# M3 identified-TTD development rerun: full-result QA and claim decision

**Run date:** 2026-07-31  
**Protocol:** `m3-identified-ttd-development-20260730`  
**Status:** development; implementation-only evidence  
**Graph mode:** disabled (local inference only)

## Decision

The rerun was scientifically necessary, but it does not replace or alter the
previous scalar M3 outputs. It asks a harder and more relevant question:

> After withholding one measured tracer, what broad age fractions are actually
> determined by the remaining tracer constraints, and can the resulting
> tracer-consistent set accommodate the withheld observation?

The answer materially narrows the M3 claim. HydroSheaf now demonstrates the
ability to compute identified sets, expose non-identifiability, and abstain.
The field benchmark does **not** show that broad groundwater-age fractions are
generally identified from the available tracer combinations.

## Reproducible run identity

The full national dataset contained 1,272 sites. Each of six tracers was
withheld in turn, producing the expected 7,632 site-tracer folds.

```powershell
.\.venv\Scripts\python.exe M3\m3_age_benchmark\scripts\run_m3_identified_ttd_benchmark.py `
  --sources national `
  --withhold-tracer all `
  --output M3\m3_age_benchmark\results\m3_identified_ttd_development_full_20260731.csv
```

Artifacts:

- `results/m3_identified_ttd_development_full_20260731.csv`
- `results/m3_identified_ttd_development_full_20260731_manifest.json`
- protocol SHA-256: `31ee761d52582082a983d2688e9e88e091cddab837ae3aeb2f59c12907e5b179`
- runner SHA-256: `53ea46e0f7a5f5afd9b86f49bf23cfd89f504bb1e581573c94e9987513a372c9`
- result SHA-256: `42b9fb4ae541e6270f0e23fcaecbc9a374a35110bee2564b4f4b7503db14a892`
- recorded Git commit: `53beb46034d5230c1a061341a5cf2175d9af858e` — **incorrect; see correction below**
- freeze commit: `2e73d51` (retroactive WP0 freeze, 2026-07-31)

The three recorded hashes match the current files. The result was written to a
new path, so the previous M3 benchmark outputs were not overwritten.

### Provenance correction (2026-07-31)

The recorded Git commit above is wrong and must not be used to reproduce this
run. None of the code that produced these results existed in `53beb46`:
`ttd_identified.py`, `ttd_graph.py`, `ttd_design.py`, the runner, and the
protocol were all untracked working-tree files at run time, and
`joint_lpm.py::tracer_response_kernel` — which `ttd_identified.py` imports — is
absent from that commit entirely.

WP0 of the implementation plan ("Before coding: commit or otherwise isolate the
current intended M3 and nuclear changes") was skipped. The stage was developed
and executed in an uncommitted tree, so at the time this QA document was
written the run was not reproducible from version control.

Commit `2e73d51` retroactively freezes that working tree and is the correct
provenance reference. It is an honest freeze, not a reconstruction: it also
contains the separate 2026-07-28 M3-correctness changes to `joint_lpm.py`,
because no intermediate commit separates the two and the split cannot be
verified without rerunning the 17,808-row correctness benchmark. Anyone
reproducing this run should be aware that the frozen `joint_lpm.py` carries both
changes.

Focused tests at freeze time: 36 passed across `test_ttd_identified.py`,
`test_ttd_graph.py`, `test_ttd_design.py`, and
`test_m3_identified_ttd_benchmark.py`.

## Integrity and leakage audit

- 7,632 rows are present, with 7,632 unique site-withheld-tracer pairs.
- All six intended withheld tracers are represented.
- `reference_fields_used` is false for every row. Published ages and published
  age fractions were not used to construct the feasible set.
- `conditioned_on_held_out_observation` is false wherever recorded. The held-out
  observation was used only after fitting for predictive assessment.
- No row has an unrecorded `evaluation_error`.
- Five otherwise feasible folds could not form a supported held-out likelihood
  and were preserved as explicit held-out abstentions rather than invented
  predictions.
- Graph information was disabled, so this run cannot support a graph-gain claim.
  To be precise about *how* it was disabled: `graph_mode: disabled_local_only`
  is a literal string written to the protocol and output columns, not a runtime
  switch. `ttd_graph.py` has no callers anywhere in the M3 pipeline, and no code
  path exists that could enable graph conditioning. The guarantee that no graph
  information entered this run is therefore stronger than a configuration flag —
  the capability is simply not wired in. The corollary is that WP3 remains
  unstarted in practice, despite the module existing.

The first full-run attempt also exposed a software edge case: an unsupported
held-out contaminated-mixture likelihood returned an abstention without a
prediction interval. The runner was repaired to preserve that state as an
audited abstention, and unexpected row failures are now converted to explicit
row-level abstentions instead of terminating the benchmark.

## Identifiability result

Of 7,632 total folds, 3,501 had enough eligible observations to enter the local
assessment. Their outcomes were:

| Outcome among eligible folds | Count | Percent |
|---|---:|---:|
| All three age fractions identified at the declared width gates | 269 | 7.68% |
| At least one, but not all, requested fractions identified | 391 | 11.17% |
| Feasible set, but no requested fraction met its width gate | 1,794 | 51.24% |
| Inconsistent constraints / solver could not produce a feasible set | 975 | 27.85% |
| Insufficient conditioning tracers | 72 | 2.06% |
| **Total eligible** | **3,501** | **100.00%** |

Thus, 2,454 of 3,501 eligible folds (70.09%) had a non-empty local
tracer-consistent set. Only 269 of those 2,454 feasible folds (10.96%)
identified all three requested fractions.

The Holocene fraction was the principal bottleneck:

| Age functional | Identified folds among 2,454 feasible sets | Median bound width | 90th-percentile width |
|---|---:|---:|---:|
| Anthropocene fraction | 477 | 0.745 | approximately 1.00 |
| Holocene fraction | 269 | 0.952 | approximately 1.00 |
| Pleistocene fraction | 453 | 0.629 | approximately 1.00 |

Widths near 1 span nearly the full possible fraction range and are therefore
not decision-informative.

## Held-out tracer assessment

There were 2,449 supported held-out prediction intervals. The measured value's
uncertainty interval was compatible with the feasible prediction range in
1,819 cases: **74.28%** (95% Wilson interval **72.51% to 75.97%**).

This is an interval-overlap compatibility rate, not conventional nominal 95%
coverage. It should not be described as a calibrated 95% predictive interval.

| Withheld tracer | Eligible folds | Feasible local sets | All fractions identified | Supported predictions | Compatible | Compatibility | Median prediction width / observation sigma |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3H | 1,033 | 816 | 77 | 816 | 639 | 78.31% | 237.8 |
| SF6 | 949 | 749 | 150 | 745 | 487 | 65.37% | 2.21 |
| CFC11 | 108 | 49 | 7 | 48 | 21 | 43.75% | 67.2 |
| CFC12 | 107 | 48 | 11 | 48 | 22 | 45.83% | 93.9 |
| CFC113 | 98 | 36 | 8 | 36 | 20 | 55.56% | 78.2 |
| 14C | 1,206 | 756 | 16 | 756 | 630 | 83.33% | 21.7 |
| **Total** | **3,501** | **2,454** | **269** | **2,449** | **1,819** | **74.28%** | — |

The very large uncertainty-normalized widths for several tracers show that
compatibility was often obtained with weak prediction bounds. Consequently,
74.28% compatibility should not be read as strong predictive performance.
The manifest's pooled native-unit width statistics are intentionally not
interpreted because they combine TU, pptv, and pmc and have no physical meaning.

An important diagnostic is that narrow age-fraction bounds did not guarantee
held-out predictive compatibility. Across all tracers, folds classified as
fully `IDENTIFIED` had 131 compatible predictions among 268 supported
predictions (48.88%). This aggregate is strongly influenced by SF6 and the CFCs:
identified SF6 folds had 53/150 compatibility, whereas identified 3H and 14C
folds had 58/77 and 14/16, respectively. The result points to tracer-specific
model discrepancy, input-history or correction uncertainty, and/or a mismatch
between the age-functional width gate and predictive reliability. It is not
evidence that a linear-program solver can make an incorrect set uniquely true.

## What changed relative to the previous M3 result?

The older strict scalar benchmark remains exactly what it was: an emulation
comparison against model-derived reference products. Its reported strict
configuration had 329 estimates, median absolute log10 discrepancy 0.027932,
log10 RMSE 0.276882, R-squared 0.937147, and 87.538% within a factor of two.
Those metrics are not recalculated or overwritten by this run.

The scientific interpretation changes as follows:

| Question | Previous scalar M3 evidence | Identified-TTD rerun |
|---|---|---|
| Can HydroSheaf reproduce model-derived reference summaries? | Strong agreement in the strict reported configuration | Not retested; previous result remains |
| Are field age fractions uniquely determined by the available tracers? | Not established by reference emulation | Usually no under the declared development model and width gates |
| Is unsupported precision made visible? | Limited scalar reportability guards | Yes: identified, partially identified, abstain, infeasible, and ineligible states are explicit |
| Is there independent predictive evidence? | Sparse scalar cross-validation | Held-out tracer set-compatibility is evaluated, but is only 74.28% and often uses broad bounds |
| Does graph structure improve inference? | No robust general gain in the prior M3 audit | Not assessed; graph mode was disabled |

The rerun therefore does not reverse the code-emulation result. It prevents that
result from being mistaken for proof of field identifiability.

## Claim decision

Defensible now:

1. HydroSheaf can form tracer-consistent TTD identified sets, bound broad age
   functionals, validate against an untouched tracer, and explicitly abstain.
2. Under this development protocol, only 7.68% of eligible field folds identify
   all three requested broad age fractions at the declared width gates.
3. The current forward assumptions and uncertainty intervals show substantial
   model-data tension: 27.85% of eligible folds do not admit a feasible local
   TTD, while supported held-out compatibility is 74.28% and is often broad.
4. The previous high scalar agreement is implementation/emulation evidence, not
   proof that field TTDs or age fractions are uniquely recoverable.

Not defensible yet:

1. that the age distribution or all broad age fractions are generally
   identifiable from these field data;
2. that an `IDENTIFIED` age-fraction label guarantees calibrated held-out
   prediction;
3. that the graph layer improves identifiability or prediction; or
4. that this development run is confirmatory evidence.

## Recommended next analysis

1. Diagnose the 975 inconsistent folds by tracer combination and by aquifer,
   separating input-history, gas-correction, observation-error, and forward-model
   discrepancy hypotheses.
2. Recalibrate uncertainty models using a predeclared procedure, then rerun this
   same held-out design without changing the decision rules after seeing results.
3. Examine whether functional width gates predict held-out compatibility; revise
   the gate only in a new development protocol, not retrospectively in this run.
4. Freeze a confirmatory protocol and environment after the model-discrepancy
   audit passes.
5. Introduce graph constraints only with independently sourced edge provenance
   and compare them against this locked local baseline.

This is the valuable outcome of the rerun: HydroSheaf becomes scientifically
stronger by reporting what the tracers do **not** determine, rather than forcing
a precise-looking age estimate where the evidence cannot support one.
