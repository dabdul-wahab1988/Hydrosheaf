# Results

## Benchmark overview and numerical validity

The live-PHREEQC factorial benchmark generated 240 scenarios across four
hydrochemical archetypes (carbonate, crystalline/silicate, evaporitic, and
mixed), each with known active phases, directions, and extents, yielding
21,600 factorial inverse fits across seven compared methods, three
analytical-noise levels, and five ion panels [[TAB:TAB-1]]. Agreement
between the PHREEQC-simulated endpoint chemistry and the concentration
change implied by the sampled stoichiometry and extents was high: the
maximum PHREEQC-to-stoichiometric generation residual across all retained
scenarios was 6.56 x 10^-4 mmol/L, several orders of magnitude below the
smallest applied analytical-noise level, confirming that the generator
itself introduced negligible error relative to the noise conditions under
test [[FIG:FIG-1]]. In well-conditioned, noise-free cases, exact recovery
of the true reaction set was confirmed before any stress test was applied.
Convergence was high but not universal: the guarded and thermodynamic
elastic-net models converged in 82.5% of fits, lasso and elastic net in
76.7-77.9%, and bounded least squares in only 35.4%, reflecting the
comparative difficulty of the unregularised, unbounded optimisation problem
[[TABREF:TAB-1]]. Non-convergent fits were retained and flagged rather than
discarded, so all reported recovery statistics reflect the full factorial
design rather than a convergence-filtered subset, and every summary
statistic below traces to an explicit fit-level record carrying its own
convergence flag, runtime, and identifiability diagnostics [[TABREF:TAB-1]].
Across all 3,600 bounded-least-squares fits, converged fits (n = 1,592) had
higher mean phase F1 (0.438 versus 0.378), higher class F1 (0.474 versus
0.440), and lower false-discovery rate (0.673 versus 0.745) than the 2,008
flagged non-converged fits. Retaining the latter is conservative, but the
convergence flag is substantively informative.

## Good concentration fits can conceal incorrect pathways

Across the benchmark, concentration-reconstruction accuracy and phase-support
accuracy diverged substantially. Among the lowest concentration-residual
quartile of fits, 55.0% still had exact-phase F1 below 0.80, directly
supporting the central claim that a low residual does not imply a correctly
recovered mechanism [[FIG:FIG-2]]. At 3% analytical noise with the full
11-ion panel, the training-calibrated Hydrosheaf guarded model achieved a
mean exact-phase F1 of 0.563 [95% CI 0.525-0.602] and a mean equivalence-class
F1 of 0.609 [0.569-0.648], a materially higher equivalence-class than
exact-phase score that is itself evidence of structural, not merely
numerical, ambiguity [[TABREF:TAB-1]]. The guarded model's false-discovery
rate of 0.361 [0.319-0.402] was lower than the legacy thermodynamic
elastic-net model's 0.465 [0.428-0.503] and the unconstrained lasso and
elastic-net models' 0.505 and 0.516 respectively, showing that
thermodynamic guarding and penalty scaling reduced spurious reaction
selection without eliminating it. Reaction-equivalence classes recovered a
consistently higher fraction of scenarios than exact-phase matching across
every comparator, and the gap between the two was widest for the three
exactly collinear class pairs identified in the dictionary's null-space
structure, confirming that part of the
observed exact-phase shortfall reflects information-limited ambiguity rather
than solver failure [[FIG:FIG-3]]. Bounded least squares, the only
comparator without a sparsity penalty, illustrated the opposite failure
mode. Its mean phase F1 of 0.435 [0.413-0.457] and false-discovery rate of
0.697 [0.679-0.715] were the worst of any comparator, despite the lowest
mean reconstruction RMSE of all seven methods (0.008 mmol/L)
[[TABREF:TAB-1]]. An unpenalised least-squares fit can therefore
reconstruct concentrations essentially exactly while still selecting a
materially wrong reaction set. Lasso and elastic net narrowed this gap
relative to bounded least squares but still returned false-discovery rates
above 0.50, and only the guarded and Hydrosheaf-Core models, which combined
sparsity with thermodynamic and evidence-based screening, brought
false-discovery rate below 0.40 while maintaining comparable or better
phase and equivalence-class F1. The Hydrosheaf-Core evidence gate itself had
a mean reaction-evidence score of 0.622 across the synthetic benchmark, but
only 23.0% of reaction rows it flagged as high evidence corresponded to a
truly active reaction in the known ground truth, so a high evidence score
from this gate should be read as a moderate prior in favour of a reaction,
not as a reliable positive indicator on its own.

## Effects of sparsity and stoichiometric identifiability

Under the full ion panel the reaction dictionary had rank 8 and nullity 4,
so at least four independent reaction directions could be reallocated among
themselves without changing any predicted concentration, regardless of
noise level, solver choice, or sample size. Three exact reaction-equivalence
classes accounted for this null space: the two oppositely signed cation
exchange pairs (calcium-sodium and magnesium-sodium exchange, each exact
algebraic negatives of one another) and, less expectedly, an anorthite/
calcite pair, in which silicate weathering and carbonate dissolution proved
exactly indistinguishable under the measured ion panel because both
reactions were written to load identically onto calcium and bicarbonate,
the only ions either reaction affects that the panel measures
[[FIGREF:FIG-3]]. This structural ceiling, established independently of any
specific inverse fit, sets an upper bound on exact-phase recoverability
under the full panel that no comparator method exceeded. Conditioning on
additional diagnostic evidence separated the two exchange pairs
substantially (mean evidence-lifted resolution index, ELRI, around 0.3) but
barely separated anorthite from calcite (mean ELRI below 0.1); correctly
identifying the active member of
the anorthite/calcite pair when truly active succeeded in only 49% of
core-tier cases, rising to 70% with plus-lite diagnostics and falling back
to 65% with the enhanced tier, a non-monotonic pattern [[FIGREF:FIG-3]].
Regularisation paths showed stable, well-separated support in
well-conditioned archetypes and abrupt, penalty-sensitive support changes in
archetypes with higher mutual coherence, consistent with this rank
deficiency rather than with noise alone (Supplementary Methods S4;
Supplementary Figure S6). Reducing $\lambda_1$ increased recall at the cost
of false discoveries; increasing $\lambda_2$ stabilised support across
bootstrap resamples without resolving the null-space ambiguity, improving
stability in correlated but non-degenerate subsets while leaving the three
collinear pairs unresolved by construction. Bootstrap support-stability
analysis independently confirmed this: well-separated reactions were
selected with high, penalty-insensitive frequency across resamples, while
the three collinear pairs shifted sharply with small changes in
$\lambda_1$. The Mechanism Resolution Score (MRS), calibrated on
the carbonate, crystalline, and evaporitic archetypes, achieved 48.9%
four-class classification accuracy on the untouched mixed archetype
(Supplementary Methods S4; Supplementary Figure S7). This is materially
better than the 25% expected under uniform random assignment but only 11.3
percentage points above the 37.6% majority-class baseline. This indicates
that structural and stability diagnostics carry some transferable
information about recovery reliability, but 48.9% is a moderate, not a
strong, transfer result and is reported without inflation: the MRS narrows
expectations about recovery reliability but does not by itself certify a
specific interpretation as correct.

## Missing ions and minimum analytical information

Ion-ablation results showed that removing individual diagnostic ions
produced abrupt, reaction-specific changes in phase support even when the
residual on the retained ions remained low, consistent with the
introduction's prediction that conventional missing-ion tests can appear
favourable for the wrong reason (Supplementary Figure S3). Held-out-ion
prediction, which tested whether the fitted reaction predicted a withheld
ion's concentration change rather than merely omitting it from the
objective, distinguished reaction combinations that residual comparison on
the retained ions alone could not. Retrospective next-best-measurement
selection based on expected ambiguity reduction was checked against every
non-selected candidate within each of 40 scenarios. The top-ranked ion
produced mean realised class-F1 gain 0.099 versus 0.013 for the other 400
candidates, and mean support change 0.429 versus 0.151. Because this was a
retrospective, same-benchmark comparison, it supports prioritisation for
prospective testing rather than claiming proven monitoring benefit
(Supplementary Figure S5).

The data-tier
experiment [[FIG:FIG-4]] showed a consistent, modest improvement in mean class F1 moving
from the core panel (0.606) through plus-lite (0.610) to the enhanced panel
(0.614), with false-discovery rate falling from 0.383 to 0.364 to 0.361
[[TABREF:TAB-1]]. ELRI values rose more
sharply, from a mean of 0.225 (core) to 0.243 (plus-lite) and a
conditionally-preferred-or-resolved fraction rising from 55.1% to 62.8% for
both plus-lite and enhanced tiers, indicating that additional diagnostics
narrowed evidential support within ambiguous equivalence classes more than
they changed the underlying mass-balance recovery rate [[FIGREF:FIG-4]]. The
divergence between the modest gain in class F1 (0.606 to 0.614, a difference
of 0.008 across the full core-to-enhanced range) and the larger gain in
conditionally-preferred-or-resolved fraction (55.1% to 62.8%) is itself
informative: it shows that most of the practical value of the additional
diagnostics tested here lay in helping an analyst choose among already
identified equivalence-class members, rather than in expanding the set of
scenarios for which the correct class was found at all. Realised
next-best-measurement value tracked expected value closely, but diverged for
reactions whose diagnostic ion carried elevated analytical noise, so
expected-information-gain rankings should be read as a priority ordering
rather than a guaranteed reduction in ambiguity.

## Value and limits of thermodynamic screening

Comparing unconstrained and PHREEQC-bounded fits showed that thermodynamic
bounds reliably eliminated fitted reaction directions inconsistent with the
local saturation state, removing a class of false positives that no purely
statistical penalty could distinguish from a feasible alternative
[[FIG:FIG-5]]. This screening effect was real but partial: in the predeclared
3%-noise, full-11-ion subset, 174 of 240 thermodynamically constrained fits
(72.5%; 95% Wilson interval 66.5%-77.8%) remained below an equivalence-class
F1 of 0.80, showing that feasible, non-unique solutions persisted after
infeasible directions were removed. The conventional PHREEQC
`INVERSE_MODELING` baseline illustrated the same limit from a different
angle. Under strict 5% analytical-uncertainty bounds it was feasible for
45.8% of scenarios; relaxing the bound to 20% raised feasibility to 99.6%,
but at the cost of generating a mean of 185.8 feasible models and 6.8
minimal models per scenario [[TABREF:TAB-1]]. Its first-minimal-model
equivalence-class F1 of 0.512 was below every guarded comparator's
equivalence-class F1, showing that reporting one arbitrarily selected
feasible model, even from a fully thermodynamically consistent solver,
understated the true model multiplicity that the same solver's own output
revealed. Aggregating across the strict and relaxed uncertainty settings,
the guarded model's mean phase F1 of 0.563 and equivalence-class F1 of 0.609
both exceeded the conventional baseline's first-minimal-model values of
0.510 and 0.512 respectively, while the conventional baseline's own
minimal-model union achieved a higher equivalence-class F1 of 0.539, and an
oracle-selected minimal model, an idealised upper bound obtained only by
already knowing which feasible model was correct and therefore not
achievable by any practical selection rule, reached 0.586. This ordering
shows that the conventional solver's thermodynamic completeness was not
itself the limiting factor on its theoretical ceiling; the practical
limiting factor was which one of its many internally consistent outputs a
residual-based rule reported as the answer, and no residual-based rule
closed the gap to the oracle ceiling in this benchmark [[TABREF:TAB-1]].

## Northern Ghana chemistry-only demonstration

The Northern Ghana chemistry-only demonstration used all 160 wet-to-dry
borehole pairs across four administrative regions (Northern, North East,
Upper East, Upper West); no independent aquifer-type or lithology
classification is available for these boreholes (Methods). Median
Hydrosheaf-Core evidence score was 0.690 and median total-dissolved-solids
consistency score was 0.942, indicating generally coherent, chemically
plausible fits [[FIG:FIG-6]]. Under the transferred MRS calibration, all 160
pairs were classified as partially identifiable, none as fully identifiable,
supporting conservative
class-level rather than single-phase reporting for this dataset. Fluoride
carried the highest retrospective measurement-value score among tested
ions in the Northern Ghana transfer analysis, consistent with its
diagnostic specificity for fluorite/apatite-related reactions that the
major-ion panel alone constrains only weakly. The external field-transfer
extension, repeating the ELRI audit on three independently sourced
chemistry datasets, computed saturation indices independently by live
PHREEQC for all three datasets (Talensi and Lower Anayari measure complete
pH, temperature and major-ion panels, the only inputs this computation
needs, so it was not withheld from them as if it required strontium or
silica), and returned a median ELRI of 0.072 for `NorthernGhana.xlsx` (160
edges) and Lower Anayari (49 edges), and 0.147 for Talensi (85 edges),
consistently low values that reinforce the chemistry-only, non-validating
status of these comparisons [[TABREF:TAB-1]]. Comparing the identifiability-aware output against a
conventional single best-fit interpretation on the same field pairs showed
materially fewer confidently stated exact mechanisms and correspondingly
more explicit equivalence-class and alternative-pathway statements. Each
field inference unit was one wet-to-dry pair, so the design cannot estimate
separate wet-season and dry-season mechanism supports. Region-stratified pair
diagnostics are reported, but no seasonal mechanism contrast is claimed from
these paired residuals; repeated within-season observations or independent
flow-path evidence would be needed for that test. All three external
field-transfer datasets returned ELRI
values well below the enhanced-tier synthetic benchmark values, indicating
that real-world diagnostic coverage narrowed equivalence-class ambiguity
considerably less than the controlled synthetic diagnostics did, so the
Ghana and external field results demonstrate transfer under realistic data
sparsity rather than an upper bound on what the framework can achieve.

## Computational performance

Median runtime per fit ranged from 16.8 ms (Hydrosheaf-Core) to 28.0 ms
(bounded least squares) across the 21,600 factorial fits, with the guarded
and thermodynamic elastic-net models both at 17.2 ms [[TABREF:TAB-1]]. This
runtime scale, combined with linear scaling in the number of ions and
reactions used by the coordinate-descent solver, supports the framework's
feasibility for screening groundwater networks substantially larger than
the 160-pair Northern Ghana demonstration without a qualitative change in
per-edge computational cost. The conventional PHREEQC inverse baseline,
by contrast, required an explicit combinatorial search over feasible phase
subsets to report its full model set, consistent with its much higher
reported model multiplicity (185.8 feasible models per scenario on average)
relative to the single regularised solution returned by every guarded
comparator. Runtime differences among the six sparse-regression comparators
were small in absolute terms (17-28 ms median per fit) but proportionally
large, with bounded least squares taking approximately 1.6 times longer per
fit than the guarded and thermodynamic elastic-net models despite converging
on only about 35% of fits, indicating additional unproductive iterations
rather than a more thorough search. No comparator showed runtime scaling
that would prevent extending the benchmark, or an operational deployment, to
substantially more than the 160 field edges evaluated here.
