# M7.6 claim decision — auxiliary M3 mechanism diagnostic

**Run:** `RUN-M7-6-M3-MECHANISM-20260731-01`  
**Runner commit:** `b06523eef09104d3ac6a4a3d74bfc26576029a37`  
**Protocol:** `e2473186169fa4cc4fb5d699dac243677ce6867a` plus the frozen
pre-test amendment `m7_6_execution_amendment_20260731.md`  
**Seeds:** development 5401–5406; locked test 5501–5512  
**Bootstrap:** 10,000 paired case-block resamples

## Execution gate

The run passed the binding environmental-isotope age control. Adding E to N
for T2 produced an exact MAE difference of 0.0000 years with a 95% interval of
[0.0000, 0.0000]. All locked-test predictions were written before the truth
files were opened for scoring.

## M3 mechanism results

Severe shared nuisance increased full nuclear-panel infeasibility relative to
none by **+0.2882** (95% CI **+0.2118 to +0.3646**). This demonstrates that the
declared synthetic nuisance can create additional infeasible age constraints.
It did not, however, reproduce the M3 signature in a selective way.

Under the severe condition, reducing minus non-reducing pair-infeasibility was:

| Diagnostic pair family | Difference | 95% CI | Decision |
|---|---:|---:|---|
| CFC-11-containing pairs | +0.7188 | [+0.6667, +0.7500] | positive association |
| CFC-12-containing pairs (specificity control) | +0.7396 | [+0.7240, +0.7500] | specificity fails |
| ³H--³H/³He pair | +0.3229 | [+0.1354, +0.4896] | positive association |

Because the CFC-12 specificity control moved at least as strongly as the
CFC-11 contrast, the predeclared selective-degradation gate failed. The
result therefore does **not** support a controlled-synthetic CFC-11-specific
explanation under this generator.

## Decision and claim boundary

**Execution decision:** PASS.  
**Selective CFC-11 mechanism decision:** NOT SUPPORTED.  
**Shared-nuisance finding:** supported only as a controlled-synthetic increase
in full-panel infeasibility, not as an M3-like family-specific explanation.

The authorised statement is:

> In the tested controlled-synthetic generator, the declared shared nuisance
> increased full nuclear-panel infeasibility, but the redox-stratified CFC
> pattern failed the predeclared CFC-12 specificity control; the experiment
> therefore did not identify a selective CFC-11 mechanism.

The following statements remain prohibited:

- the tested mechanism caused the USGS or Ghana infeasibility pattern;
- M7.6 is field validation of age, topology, or tracer corrections;
- any general superiority claim for HydroSheaf or multi-stream fusion;
- revision of any locked M7.3–M7.5 number.

M7.6 is retained as an auxiliary diagnostic audit. The unresolved USGS cause
remains unresolved and requires co-located field measurements or a different
validated error model.
