## 3. Methods

### 3.1 Computation and figure authority; no re-simulation

This paper performs no PHREEQC, MODFLOW/MODPATH, or TracerLPM re-run. Every
number reported below was read from a result file that the age, topology, or
reaction benchmark had already written and locked, under each benchmark's own predeclared
design, before this comparison began; residual staleness risk from code
changes made after those files were locked was checked directly against the
specific modules each benchmark exercises and is disclosed in Limitations
(Section 5.5) rather than assumed away. A Python harmonization layer reads
these files and performs only disclosed, auditable arithmetic on
already-published values -- selecting a row, summing a column, grouping and
taking a median, or, in one place, pooling published per-reaction confusion
counts into a precision/recall pair -- and writes the results as tidy CSV
exports. R consumes those exports and owns every figure; no figure recomputes
a reported statistic. Where an operation could plausibly disagree with a
value the source component already reports as its own headline number, both
values are retained and the discrepancy is written to an evidence-discrepancy
record rather than silently resolved.

### 3.2 Each benchmark's own design

**Age/residence time.** Candidate age fits are re-derived from the
USGS national release using only the published initial or fixed
lumped-parameter-model configuration (Section 2.1). A fit is admitted to the
reportable set only if the optimizer converged, the requested reported-model
emulation used an exactly supported lumped-parameter model rather than a
fallback, the number of fitted tracer observations exceeded the number of
free parameters, standardized tracer-space error was bounded, and competing
grid candidates near the optimum did not span an excessive range of
candidate ages; reference ages are used only after fitting, to score
reportable estimates, never to select or constrain them. A second, held-out
cross-validated comparison partitions sites into folds, fits a ridge
calibration on a training fold, and scores both the uncalibrated estimate and
the calibrated estimate against the ages in the held-out fold; only the
uncalibrated held-out result is treated as independent (Section 3.4). A
graph-regularization design compares single-node age estimates against
multiple graph-regularized candidates spanning a range of regularization
strengths, scored against a joint robust-improvement rule, and against
randomized negative-control graphs, to test whether network information
helps or harms inference rather than assuming it helps; Results Section 4.4
reports the smallest candidate improvement found across that range and the
randomized-control comparison, not a full strength-by-strength breakdown.

**Topology.** Independent graph-inference scenarios (head-gradient,
depth-tiered, and hydrostratigraphic edge construction) infer directed edges
from elevation and depth proxies with no access to the MODPATH reference,
and are scored against the Savage archive's 174 reference edges by standard
edge-classification precision, recall and F1. Four control families bound
this independent result: a spatial-proximity-only control tests whether
geometry alone recovers topology; a randomized-edge control and a
wrong-direction control test sensitivity to graph identity and orientation;
a misrouted-sink control reassigns each source to the wrong reference
receptor at constant graph size, testing whether the benchmark rewards
correct receptor attribution rather than merely graph shape; and a
sink-aware structural baseline, informed by the three-receptor set alone with
no hydraulic information, establishes the floor that any inference scenario
must be judged against. A separate, explicitly non-independent mode enters
MODPATH edges into the Hydrosheaf graph as prior information (override,
merge, and only-prior variants) and is reported on its own terms, never as
independent validation.

**Reaction.** Six comparator inversion methods -- bounded least
squares, LASSO [@tibshirani1996lasso], elastic net [@zou2005elasticnet], a
thermodynamically bounded elastic net, a calibrated Hydrosheaf-guarded
model, and an evidence-gated Hydrosheaf-Core comparator -- are fitted to
every one of the 21,600 factorial scenario/noise/panel combinations and
scored against the known PHREEQC ground truth on phase-support precision,
recall, F1, false-discovery rate, extent error, and reconstruction error. A
conventional PHREEQC `INVERSE_MODELING` baseline is run in parallel at a
strict 5% uncertainty tolerance and a relaxed 20% fallback, reporting model
multiplicity (the number of thermodynamically feasible models per scenario)
alongside its own phase-recovery accuracy. Support selection uses a
predeclared regularization criterion rather than the minimum-residual model
alone, following stability-selection practice for sparse regression
[@meinshausen2010stability]. A calibrated Mechanism Resolution Score is
trained on three of the four hydrochemical archetypes and evaluated,
untouched, on the fourth, so that its reported accuracy is a genuine transfer
test rather than a within-sample fit. A data-tier ablation compares three
nested measurement panels -- a core major-ion panel, a "plus-lite" panel
adding optional SiO2, strontium, and water-isotope diagnostics, and an
"enhanced" panel further adding bromide, dissolved oxygen, dissolved organic
carbon, and sulphate/nitrate isotope diagnostics -- to quantify how much
additional analytical coverage narrows the recovery gap, and a retrospective
next-best-measurement analysis, in the spirit of Bayesian experimental design
[@chaloner1995design], ranks candidate ions by their expected reduction in
reaction-class ambiguity.

### 3.3 A common evidence taxonomy

Every retained row of the three benchmarks was classified into one of seven
categories, applied without altering any component's own claim: independent
(no access to the reference used for scoring); non-independent (the fitting
procedure shares provenance with the reference, for example the age
benchmark's age-fraction scenario, which uses USGS-reported age fractions); informed
control (uses partial reference-derived structure deliberately, to establish
a floor rather than to claim skill); negative control (designed to fail);
sensitivity diagnostic (varies one design factor, such as node density, and
is not comparable one-for-one with single-run scenario metrics); prior
informed (reference information enters as an input, not as a score); and
conventional baseline (an established, non-HydroSheaf comparator method).
This taxonomy is applied to every row retained in the comparison and
plotted against each component's own headline metric in Section 4.1 and
Section 4.4.

### 3.4 Capture-type versus correctness-type metrics, and their limits

The three benchmarks do not share a statistic. Topology is a literal edge
classification, so precision and recall are exact. Reaction reports a
macro-averaged phase F1 and false-discovery rate per scenario as its own
published headline result; because no macro-averaged recall is published
alongside it, this paper additionally pools the published per-reaction
confusion counts (true and false positives and negatives, already reported
for every method and reaction in the reaction benchmark's own supplementary table) into a
micro-averaged recall and precision, and reports both the macro-averaged
precision and the micro-averaged pair, disclosing that the two conventions
differ numerically (Section 4.2). Age has no classification truth at all; its
two reported axes are the reportability rate -- the fraction of attempted
fits that clear the reference-free reportability guard, an abstention design
rather than a defect -- and the within-factor-2 agreement of the fits that
are reported, scored on an independent held-out cross-validation split. These
three pairs are presented on each component's own native scale in parallel
panels (Section 4.2), not merged onto one shared numeric axis, and the analogy
between them is stated as a disclosed methodological choice: a "capture-type"
axis (how much of the target space yields a claim) against a
"correctness-type" axis (how trustworthy that claim is), not a claim that the
three statistics are mathematically identical.

### 3.5 Reproducibility

The full harmonization layer -- five Python scripts and six R figure scripts
-- is archived under `O3/analysis/`, runs end to end from the repository
root with fixed, published inputs, and writes nothing back to any age,
topology, or reaction benchmark result file.
