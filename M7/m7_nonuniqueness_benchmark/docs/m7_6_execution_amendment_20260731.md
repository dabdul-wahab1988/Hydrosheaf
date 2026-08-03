# M7.6 execution amendment — M3 mechanism diagnostic

**Status:** pre-test execution amendment; frozen before any 5501-series case is
generated.

This amendment narrows the purpose of M7.6 to the unresolved M3 mechanism
question. It does not alter M7.3 results and does not authorise field-
validation or USGS-cause claims.

## Locked implementation choices

- Development seeds: 5401–5406.
- Locked-test seeds: 5501–5512.
- Nuisance levels: `none = 0.0`, `mild = 0.5`, `severe = 1.0`, where the scale
  multiplies the generator's declared recharge-temperature and 14C-A0 error
  standard deviations.
- All three nuisance levels are paired conditions of each seed. A seed is the
  resampling unit; the nuisance level is never treated as an additional
  independent case.
- The four stream codes are H (hydraulic), C (hydrochemistry), N (nuclear),
  and E (environmental isotope). All 15 non-empty subsets are evaluated.
- T2's E feature mask is empty because the generated stable isotopes are fixed
  at recharge and transported conservatively. The binding gate is the exact
  native `E_added_to_N` T2 age contrast; any non-zero mean or interval fails the
  run before an M3 interpretation is accepted.
- T3's synthetic source fraction is the recharge x-position divided by the
  maximum recharge x-position in that case. This is a generator-defined target,
  not a field endmember estimate.

## M3-compatible feasibility diagnostic

The claim-bearing diagnostic is a local convex-mixture transit-time test over
240 ages from 0.25 to 100 years. The reference responses are the independent
v2 analytic 3H, 39Ar, 14C, CFC-11/12/113, SF6, 4He and 3H/3He formulations. A
panel is feasible when a non-negative age-mixture with unit mass satisfies all
observed tracer intervals at `k = 1.96`. Every CFC pair and the 3H--3H/3He pair
are tested separately, as well as the full nuclear panel.

Redox labels are assigned from the withheld generator process history only
after truth-blind locked-test predictions are complete:

- `reducing`: sulfate-reduction or iron-reduction has occurred;
- `suboxic`: denitrification has occurred without the reducing processes;
- `nonreducing`: neither condition has occurred.

The primary mechanism contrasts are:

1. severe minus none full-panel infeasibility;
2. reducing minus nonreducing CFC-11-containing pair infeasibility at severe;
3. the CFC-12-containing pair specificity control;
4. the analogous 3H versus 3H/3He contrast.

The M3 interpretation is licensed only as a controlled-synthetic possibility:
the tested mechanism can or cannot generate an M3-like signature under this
generator. It never establishes that the mechanism caused the USGS pattern.

## Reporting gate

The M7.6 claim decision is written before any manuscript edit. Only after the
decision is complete may the M7 manuscript add an auxiliary diagnostic note;
M7.3–M7.5 locked numbers and their claim boundaries remain unchanged.
