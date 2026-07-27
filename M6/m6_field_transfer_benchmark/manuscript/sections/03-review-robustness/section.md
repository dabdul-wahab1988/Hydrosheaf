### Reviewer-requested robustness analyses

Recent Ghanaian studies use machine learning to predict groundwater-quality
outcomes in the Nabogo Basin and crystalline aquifer terrain
[@apogba2024nabogo; @abu2025ghana]. They complement, but do not duplicate,
this study: prediction of a quality index does not test whether a particular
reaction family remains identifiable when a diagnostic tracer is withheld.

The tier ladder is descended by removing fluoride from Tier 2 to Tier 1 and
then isotopes from Tier 1 to Tier 0; strontium and silica are removed earlier,
from Tier 3 to Tier 2. The combined evidence-gate sensitivity used low,
predeclared, and high corroboration profiles. Across these profiles, the Tier
2 non-identifiable fraction ranged only from 0.513 to 0.531 (predeclared
0.519), whereas Tier 3 ranged from 0.000 to 0.025 (predeclared 0.006). Thus,
the location and direction of the Tier 3-to-2 collapse are robust to these
reasonable threshold changes, although its exact magnitude remains conditional
on the operational gate.

The near-prohibitive penalty applied to reactions requiring unmeasured ions was
also varied from 20 to 10 and 50 in 40 fixed well pairs at Tiers 2 and 4. Mean
support Jaccard similarity was 0.990 for 10 versus 20 and 1.000 for 50 versus
20; mean maximum absolute extent change was 0.0091 and 0.0000 mmol/L,
respectively. The penalty therefore does not explain the principal tier
contrast over the tested range.

Bootstrap counts were originally set for runtime rather than precision. A
200-resample audit of 20 fixed Tier-4 pairs estimated the mean standard error
of stability at 0.0123 for B = 8, 0.0050 for B = 48, and 0.0025 for B = 200.
Comparisons using B = 8 should consequently be read as coarse diagnostics, not
high-precision estimates. Full unit-level and aggregate values are archived in
`results/m6_review_sensitivity.csv`.

The synthetic validation's exact-mineral F1 applies only to the fixed reaction
dictionary and simulated archetypes; family- and class-level recovery are the
transfer targets. Figure and table notes distinguish within-well wet-to-dry
pairs from chemistry-kNN candidate edges. The latter are unobserved screening
constructions, so a viable no-flow explanation remains a competing account of
every external-transfer contrast. Talensi's result cannot be attributed
causally to mining or any single process because charge-balance failure and
missing tracers make it an exploratory screening case only. Likewise, the
external datasets' smaller candidate pools can alter nearest-neighbour graph
structure independently of tracer coverage; prospective matched-size sampling
is needed to separate this sample-size confound from the tier effect.
