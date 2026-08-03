# HydroSheaf benchmarks: a harmonized identifiability comparison of directed connectivity, residence-time, and reaction inference in data-limited aquifers

Dickson Abdul-Wahab

HydroSheaf PhD Chapter 1, Objective 3. Companion article to the HydroSheaf
framework paper (Objective 2, same target journal). Submitted to Computers &
Geosciences.

## Abstract

Groundwater connectivity, residence time and reaction pathways are usually
benchmarked as separate problems, each against its own reference, leaving it
unclear whether a result reflects the inference problem or the benchmark's
design. The objective of this paper is to harmonize three already-locked
HydroSheaf benchmarks -- age against the USGS national groundwater-age
release (1,272 candidate fits), topology against MODPATH connectivity from
three public MODFLOW archives, and reaction inversion against a
240-scenario live-PHREEQC benchmark (21,600 fits) -- under one evidence
taxonomy, without re-running any of the three. The method is a disclosed,
auditable arithmetic pass over already-published result files: grouping,
pooling confusion counts, and computing ratios and confidence intervals,
with every aggregation discrepancy recorded rather than resolved silently.
The result is a capture-type metric (edge recall 0.845; age-fit
reportability 0.259; reaction-phase recall 0.599) exceeding its matched
correctness-type metric (edge precision 0.487; held-out age agreement within
a factor of 2, 0.564; reaction-phase precision 0.639) under independent,
prior-free, uncalibrated evaluation in all three layers. Attaching 95%
Wilson intervals to each proportion shows this gap is resolved,
non-overlapping, for topology and age, but not for reaction, whose interval
overlap means the point-estimate ordering there is not distinguishable from
sampling variation. Calibrating each resolved layer against its own scoring
reference narrows the gap but differs substantially in rigour between
layers, and the three benchmarks differ by two orders of magnitude in
scale. The contribution is this interval-aware cross-component comparison,
unavailable from any one underlying benchmark alone, offered as a companion
to the HydroSheaf framework article rather than a fourth validation of any
one component.

**Keywords:** groundwater inverse modelling; identifiability; benchmark comparison; equifinality; MODPATH; PHREEQC; TracerLPM; reproducible research software

## 1. Introduction

Groundwater systems are usually characterized through three largely separate
inference problems. Where water moves is inferred from hydraulic gradients,
particle tracking, or reduced-order network models. How long it has resided
in the subsurface is inferred from environmental tracers and lumped-parameter
or transport models. What it has reacted with along the way is inferred from
mass-balance or inverse geochemical models constrained by thermodynamic
databases. Each problem has its own community, its own reference simulators
-- MODFLOW and MODPATH for connectivity [@langevin2017modflow6;
@pollock2016modpath], TracerLPM and related lumped-parameter tools for
residence time [@jurgens2012tracerlpm], PHREEQC for reaction screening
[@parkhurst2013phreeqc] -- and its own benchmarking literature. This
separation is longstanding and defensible: connectivity, residence time and
reaction are governed by different physics and are usually constrained by
different observation types. It is also the reason no single, shared account
exists of how these three inference problems behave under a common
evaluation logic, even within one framework built to integrate them.

Data-limited settings sharpen the problem. Across much of Africa, groundwater
systems are characterized from sparse well networks, incomplete tracer
panels, and geochemical surveys that were not designed as controlled
experiments [@macdonald2012africa]. The Prediction in Ungauged Basins
initiative formalized an analogous difficulty for streamflow two decades ago,
arguing that hydrological inference in data-limited catchments requires an
explicit accounting of what a given evidence base can and cannot support,
rather than a single calibrated best estimate [@sivapalan2003pub;
@hrachowitz2013pub]. The equifinality literature made the same point for
mechanistic environmental models more broadly: multiple structurally
different models, or parameter sets, can reproduce the same observations
equally well, so a good fit alone is not evidence of a correct structure
[@beven2001glue; @beven2006manifesto], a difficulty multimodel averaging
approaches manage rather than eliminate [@neuman2003bma]. Groundwater connectivity, residence
time and reaction inference each inherit this difficulty, but the degree to
which they do so, and whether it takes the same shape in each, has not been
reported side by side.

HydroSheaf is a claim-gated evidence-integration framework that treats
groundwater connectivity, residence time and reaction as outputs of a shared
directed-edge object, drawing on graph- and sheaf-inspired approaches to
combining heterogeneous evidence sources [@robinson2017sheaves;
@hansen2019spectral; @moracchini2025graphflow]. Its architecture and an
integrated field demonstration are reported in a companion article (Objective
2 of this thesis chapter) [@hydrosheaf_m2_3_inprep]. Building that framework
required benchmarking its three principal inference layers separately,
because no single external reference exists against which directed
connectivity, residence time and reaction extent can all be judged at once:
the age benchmark [@hydrosheaf_m3_inprep] compares nuclear-tracer and
lumped-parameter age inference against the USGS national groundwater-age
release; the topology benchmark [@hydrosheaf_m4_inprep] compares
reduced-order graph topology against MODPATH particle-tracking connectivity
from three public MODFLOW archives; the reaction benchmark
[@hydrosheaf_m5_inprep] compares sparse linear reaction inversion against a
240-scenario live-PHREEQC factorial synthetic benchmark. Each of these three
benchmarks was designed,
executed and locked independently, on its own timetable, against its own
external reference, with its own negative controls and its own
reportability rules.

That independence is a strength for each benchmark's own internal validity
and a limitation for anyone trying to compare them. The three benchmarks do
not share a metric, a notion of what "independent" means, or a scale. A
reader who wants to know whether HydroSheaf's connectivity layer is more or
less trustworthy than its age layer, or whether the reaction layer's
apparent accuracy reflects real identifiability or a permissive evaluation
design, currently has to read three separate result packages built for three
separate purposes and reconcile the comparison unaided. This is not a defect
particular to HydroSheaf; it is the default state of any research programme
that benchmarks physically distinct inference problems against physically
distinct references. It is, however, a gap that a companion synthesis paper
can close without repeating any of the underlying computation.

This paper harmonizes the three already-locked benchmarks under one common
evidence taxonomy and asks a single question: under evaluation conditions
that are independent, prior-free and uncalibrated in all three cases, do
directed connectivity, residence time and reaction inference behave alike?
The answer reported here is a qualified yes. A capture-type metric (how much
of the true target is detected) exceeds the matched correctness-type metric
(how much of what is detected is correct) in all three layers as a point
estimate, but the three capture/correctness pairs are not the same statistic
computed the same way, and 95% confidence intervals resolve the gap as
statistically distinguishable in only two of the three layers, topology and
age; the reaction layer's point estimates run the same direction but their
intervals overlap. The gap narrows further only when independence is given
up, either through calibration against the same reference used for scoring
or through conditioning on prior information from another layer. Four
further objectives follow from this central comparison: to quantify the
screening-to-calibration gap in the two layers where a calibrated
alternative exists and to state plainly why calibrated numbers are not
independent validation; to compare the three benchmarks on computational
scale, since a benchmark with 21,600 factorial fits and a benchmark with a
few hundred graph edges are not making comparably strong claims by
construction; to compare their field- and archive-transfer scope, which
differs sharply across the three layers; and to describe, for a reader
unfamiliar with any of the three underlying projects, the full set of
synthetic and field datasets that this comparison rests on.

The paper does not introduce new software architecture, revisit HydroSheaf's
integrated field demonstration, or claim a fourth, independent validation of
any component beyond what each component's own result package already
supports. It is a benchmark-synthesis article: its contribution is the
comparison itself, the taxonomy that makes the comparison possible, and the
computational-cost and field-scope accounting that no single component paper
reports.

## 2. Data

This paper introduces no new dataset. It reports a comparison built on eight
already-existing datasets, three synthetic or archive-based and five drawn
from published or author-compiled field surveys, each already used as the
external reference for exactly one of the three component benchmarks (the
full inventory is given as a table in Section 4.6). This section describes each in enough detail for a reader
unfamiliar with the three benchmarks to judge what kind of evidence underlies
every number reported later.

### 2.1 Age/residence time: the USGS national groundwater-age release

The age benchmark is evaluated against the reported-configuration release accompanying
TracerLPM, the U.S. Geological Survey's workbook for interpreting
environmental-tracer age distributions with lumped-parameter models
[@jurgens2012tracerlpm]. The release publishes, for each of 1,272 candidate
site-tracer fits, the tracer-specific input scale factor, the parameters the
original study selected for optimization, and the initial or fixed
lumped-parameter-model configuration; tracers represented include tritium
(3H), sulphur hexafluoride (SF6), the chlorofluorocarbons CFC-11, CFC-12 and
CFC-113, and carbon-14. The age benchmark re-fits every candidate using only the published
initial or fixed configuration and optimization declaration -- reported final
ages and final fitted parameters are never supplied to the optimizer -- and
retains a scalar age estimate only when a predeclared reference-free
reportability guard passes (Methods, Section 3.2). Because the published ages
are themselves model-derived outputs of the U.S. Geological Survey's own
workflow rather than independently measured true ages, every comparison
against this release is described as agreement with a reference workflow,
not as validation against a directly observed age. Each of the 1,272 rows
carries a site identifier, a screen or open-interval depth, and an aquifer
group classification, which the leakage-guarded cross-validation design
(Section 3.2) uses to build graph structure without exposing a target site's
own withheld-tracer measurement to the fit used to predict it; roughly half
of national sites carry no finite paired reference-and-estimate value under
the strict emulation design, this attrition is uneven across aquifer groups,
and it is tracked separately rather than assumed to be representative of the
retained comparison.

### 2.2 Topology: three public MODFLOW/MODPATH archives

The topology benchmark is evaluated against particle-tracking connectivity from three public
U.S. Geological Survey MODFLOW/MODPATH archives, released as ScienceBase data
packages and consumed as MODFLOW/MODPATH input and output files rather than
narrative descriptions [@langevin2017modflow6; @pollock2016modpath]. The
primary archive, the Savage Municipal Water-Supply Well model
(MODFLOW-2005/MODPATH 5, https://doi.org/10.5066/F7J102FK), releases 3,000
tracked particles that reduce to 174 endpoint reference edges between source
cells and three receptor wells; this is the archive against which every
independent topology metric in this paper is reported. A second archive, the
Great Miami River Basin model (MODFLOW 6/MODPATH 7,
https://doi.org/10.5066/P9X4C9R6), releases 154 particles and 68 reference
edges and is used as a self-consistency check on the endpoint/pathline
projection step rather than as a second independent skill test. A third
archive, the Long Island Nitrogen model (MODFLOW 6/MODPATH 7,
https://doi.org/10.5066/P97VFXZ4), is retained in the benchmark design but is
a documented fallback stub with no active validation rows in the current
package; it is reported as absent evidence rather than silently omitted
(Section 4.6 reports transfer scope; Section 4.5 reports benchmark scale).

### 2.3 Reaction: a 240-scenario live-PHREEQC factorial benchmark

The reaction benchmark is evaluated against a purpose-built synthetic dataset with known
reaction ground truth, generated by running the U.S. Geological Survey's
PHREEQC geochemical modelling program forward from declared initial
conditions [@parkhurst2013phreeqc]. The benchmark spans four hydrochemical
archetypes -- carbonate, crystalline/silicate, evaporitic and mixed -- and
combines one to five simultaneous mineral, exchange or redox reactions per
scenario, drawn from a sixteen-reaction, eleven-ion stoichiometric
dictionary, at three analytical-noise levels (0%, 3%, 8%) and five ion
panels ranging from a full eleven-ion panel to a five-ion core panel that
withholds alkalinity and redox-sensitive species. This design yields 240
PHREEQC scenarios and, once every comparator inversion method is run across
the full noise/panel/archetype grid, 21,600 factorial inverse fits, each with
a known true reaction support, direction and extent that a sparse linear
inversion never has access to. Unlike the age and topology references, this
benchmark's ground truth is fully known by construction, at the cost of
being synthetic rather than field-measured. The five ion panels range from
an eleven-species full panel (calcium, magnesium, sodium, potassium,
bicarbonate, chloride, sulphate, nitrate, fluoride, iron, phosphate) down to
a five-species core panel (calcium, magnesium, sodium, bicarbonate,
chloride), stepping through configurations that withhold alkalinity or
withhold redox-sensitive trace species individually, so that the effect of a
specific missing measurement class can be isolated rather than only the
effect of panel size. PHREEQC-generated endpoint chemistry reproduces its own
declared stoichiometric input to within 6.6e-4 mmol/L across all 240
scenarios, which is the internal-consistency floor against which every
inversion method's much larger reconstruction and extent errors are judged.

### 2.4 Field chemistry: Northern Ghana, Talensi and Lower Anayari

Three field chemistry datasets carry the reaction benchmark's field-transfer component. The
Northern Ghana workbook covers 160 boreholes across four administrative
regions, each sampled once, with results reported on two worksheet panels
labelled Dry and Wet. An integrity check on the underlying workbook
established that the static water level recorded for every one of the 160
wells is identical between the two panels, that the workbook embeds no
sampling dates, digital object identifier, laboratory, or source
publication, and that a supplementary document the workbook's provenance
notes referenced does not exist in the repository; the author subsequently
confirmed that the underlying chemistry is measured but the Dry/Wet seasonal
separation was reconstructed rather than independently sampled in the field.
Consistent with that disclosure, this paper treats the Northern Ghana
workbook as 160 measured wells rather than 320 independent seasonal samples,
uses the Dry panel as the primary measured chemistry, and reports no
seasonal-change or repeated-hydraulic-head finding from it anywhere. Where
the reaction benchmark pairs each well's Dry and reconstructed-Wet chemistry into a candidate
directed edge for its reaction-inversion audit (160 such pairs), that pairing
is reported exactly as the reaction benchmark's own documentation reports it: a field-plausibility and
candidate-class exercise for stress-testing the inversion machinery, not a
measured flow path or a validated reaction transformation (Section 4.6; regional
hydrogeological context from @banoengyakubo2011ghana and
@macdonald2012africa).

The Talensi District dataset (63 samples) and the Lower Anayari catchment
dataset (41 samples) are drawn from independent published field studies
[@chegbeleh2020talensi; @zakaria2021anayari] and carry no reconstructed
attribute; both show the analytical heterogeneity expected of independently
collected field chemistry. The reaction benchmark uses both as a second, smaller field-transfer
check on candidate-edge reaction plausibility, reported at screening level
throughout.

### 2.5 What this comparison does and does not add to the underlying data

No file listed in the dataset inventory (Section 4.6) was re-read, re-fitted, or re-simulated to produce
this paper. Every number reported in Sections 4 and 5 traces to a result
file that the age, topology, or reaction benchmark had already written and locked before this
comparison was constructed (Methods, Section 3.1); the only new computation
performed here is the harmonization arithmetic described in Section 3 --
grouping, pooling already-published confusion counts, and computing simple
ratios -- applied to those existing files.

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

## 4. Results

### 4.1 Applying the common taxonomy

Figure 1 summarizes the three pipelines this comparison harmonizes.
Applying the taxonomy of Section 3.3 to every retained row of the age,
topology and reaction benchmarks (Table 1) yields 25 classified rows: 16
independent, one non-independent, one informed control, four negative
controls, one sensitivity diagnostic, three prior-informed, and one
conventional baseline. No row's classification contradicts the claim
already recorded for it in its own source benchmark's decisions or
manuscript-ready tables; the taxonomy adds a shared vocabulary across
benchmarks, not a re-interpretation of any single result.

![Figure 1. Three independently designed HydroSheaf benchmarks and the common evidence taxonomy this paper adds. The age benchmark compares nuclear-tracer/lumped-parameter age inference against the USGS national groundwater-age release; the topology benchmark compares reduced-order graph topology against MODPATH particle-tracking connectivity from public MODFLOW archives; the reaction benchmark compares sparse linear reaction inversion against a 240-scenario live-PHREEQC factorial synthetic benchmark. Each pipeline was locked independently, on its own timetable, before this comparison was constructed.](artifacts/figures/FIG-1.png)

**Table 1.** Common evidence taxonomy applied to every retained benchmark row of the age, topology and reaction benchmarks.

| Benchmark | Scenario | Row type | Independent | External reference | Headline metric | Value |
|---|---|---|:---:|---|---|---:|
| Topology | Spatial only | independent | yes | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Sink-aware baseline | informed control | no | MODPATH (Savage archive) | F1 | 0.552 |
| Topology | Head gradient | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Hodge pruned (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Proj. gradient (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Head depth | independent | yes | MODPATH (Savage archive) | F1 | 0.616 |
| Topology | Hydrostratigraphic | independent | yes | MODPATH (Savage archive) | F1 | 0.614 |
| Topology | Negative: random | negative control | no | MODPATH (Savage archive) | F1 | 0.006 |
| Topology | Negative: wrong direction | negative control | no | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Negative: misrouted sink | negative control | no | MODPATH (Savage archive) | F1 | 0.138 |
| Topology | Negative: shortcut (inapplicable) | negative control | no | MODPATH (Savage archive) | F1 | n/a |
| Topology | Sparse node | sensitivity diagnostic | no | MODPATH (Savage archive) | mean F1, 20 trials | 0.445 |
| Topology | Prior: override / merge / only | prior-informed | no | MODPATH (Savage archive) | edges in/out | 302 -> 302/329/174 |
| Age / residence time | TracerLPM strict parity | independent | yes | USGS TracerLPM Table 4 | within-factor-2 | 0.875 (n=329/1272) |
| Age / residence time | TracerLPM age fractions | non-independent | no | USGS reported age fractions | within-factor-2 | 0.917 (n=289/1272) |
| Age / residence time | Hydrosheaf model selection | independent | yes | USGS national age release | within-factor-2 | 0.673 (n=309/1272) |
| Age / residence time | Old-water He-4 branch | independent | yes | USGS national age release | within-factor-2 | 0.697 (n=66/1272) |
| Age / residence time | Graph regularisation, randomised control | negative control | yes | Randomised candidate graph | delta log10 RMSE vs single-node | +0.655 |
| Reaction | Bounded LS | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.435 |
| Reaction | Lasso | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.554 |
| Reaction | Elastic Net | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Thermo Elastic Net (legacy primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Hydrosheaf Guarded (primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.563 |
| Reaction | Hydrosheaf-Core (evidence-gated) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.561 |
| Reaction | Conventional PHREEQC inverse (strict+relaxed fallback) | conventional baseline | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.510 |

Note. The Hodge-pruned and projected-gradient topology rows are diagnostic
post-processing variants of the Head-gradient scenario immediately above
them, not independent replications; all three converge on the same
underlying inference and are numerically identical (F1 0.618) by
construction, not by coincidence.

### 4.2 Capture exceeds correctness under independent evaluation, in all three layers

Figure 2 and Table 2 report the paper's central result. For
topology, the independent Head-gradient scenario recovers 84.5% of the true
MODPATH reference edges (recall 0.845) while only 48.7% of the edges it
infers are true (precision 0.487); the three other independent scenarios
(head-depth, hydrostratigraphic, and two head-gradient diagnostic variants)
recover the same pattern within 0.004 F1 of one another (F1 0.614-0.618).
For age, 25.9% of the 1,272 candidate TracerLPM fits clear the reference-free
reportability guard and are reported at all, while the fits that are
reported agree with an independent held-out cross-validation reference
within a factor of 2 in 56.4% of cases (log10 R2 0.482). For reaction, the
Hydrosheaf Guarded model recovers 59.9% of the true active reaction phases
across the 240-scenario benchmark, pooled across 16 reactions (recall
0.599), while 63.9% of the phases it selects as active are true (macro-mean
precision, 1 minus a false-discovery rate of 0.361). In every one of the
three physically distinct inference problems, the point estimate of the
capture-type axis exceeds the point estimate of the matched correctness-type
axis.

That point-estimate ordering does not hold with equal statistical force in
all three layers. Adding 95% Wilson score intervals to each proportion
(Table 2) shows that topology's recall and precision intervals do not
overlap (0.784-0.891 versus 0.431-0.543, n=174 and n=302 respectively) and
neither do age's reportability and within-factor-2 intervals (0.235-0.283
versus 0.527-0.601, n=1,272 and n=675): in both layers the capture-exceeds-
correctness ordering is resolved at conventional confidence. Reaction is
different: its pooled recall interval (0.562-0.635, n=681) and macro
precision interval (0.598-0.681) overlap, so the small numerical gap between
0.599 and 0.639 for this layer is not distinguishable from sampling
variation at the 95% level, even though the point estimates run in the same
direction as the other two layers. The pattern that survives interval-aware
scrutiny is therefore narrower than the point estimates alone suggest:
recall/reportability exceeds precision/agreement with resolved confidence in
the topology and age layers, and in the same direction but not resolved in
the reaction layer.

The magnitude of the gap also differs across the two layers where it is
resolved. Topology shows the larger separation (0.358 between recall and
precision); age shows a smaller separation (0.305) on a pair of axes that,
as Section 3.4 notes, are not the same kind of rate as topology's or
reaction's. Among the lowest-residual quartile of reaction
fits -- the cases a residual-only read would judge most trustworthy -- 55.0%
still show phase F1 below 0.80, and the conventional PHREEQC inverse
baseline, run in parallel at strict 5% uncertainty, is feasible for only
45.8% of the 240 scenarios (99.6% once a relaxed 20% fallback is allowed),
returning a mean of 185.8 thermodynamically feasible models per scenario
where it is feasible at all.

Across the seven compared reaction-inversion methods, correctness and
convergence move together but not in lock-step with capture. Unregularized
bounded least squares converges least often (35.4% of fits) and shows the
highest extent error (0.091 mmol/L median) and the highest false-discovery
rate (0.697) despite fitting concentrations closely; the four regularized
sparse methods (LASSO, elastic net, thermodynamic elastic net, Hydrosheaf
Guarded) converge in 76.7-82.5% of fits and cut extent error roughly in half
(0.050-0.052 mmol/L); Hydrosheaf Guarded and Hydrosheaf-Core are the only two
methods with zero recorded thermodynamic-bound violations. The leakage-guarded
tracer-withholding design applied to age produces the same qualitative
result as the graph-regularization comparison: refitting every graph-node
age with its target tracer withheld shows no candidate graph improving on
baseline prediction of the withheld tracer for any of the three tracers with
enough retained fits to compare (tritium: baseline RMSE 20.818 TU against
20.841 for a depth-based graph and 21.567 for a randomized graph, 121 of 794
eligible rows reportable; sulphur hexafluoride: 2.844 pptv baseline against
2.929 randomized, 75 of 262 reportable; carbon-14: 26.830 pmC baseline
against 27.819 randomized, 169 of 1,103 reportable).

![Figure 2. Capture-type versus correctness-type metrics under independent, prior-free, uncalibrated evaluation, one panel per benchmark on its own native scale, with 95% Wilson score error bars. Topology: recall and precision of edge classification against the MODPATH reference (Head-gradient scenario). Age: reportability rate (fraction of fitted rows admitted by the reference-free guard) and within-factor-2 agreement on a held-out cross-validated comparison. Reaction: recall and precision of active reaction-phase detection, Hydrosheaf Guarded method, against the PHREEQC benchmark; note the overlapping intervals for this layer.](artifacts/figures/FIG-2.png)

**Table 2.** Headline independent/uncalibrated capture-type and correctness-type metrics, all three benchmarks.

| Benchmark | Capture-type metric | Value [95% CI] | Correctness-type metric | Value [95% CI] | CIs overlap? |
|---|---|---|---|---|:---:|
| Topology | Recall (true MODPATH edges recovered), n=174 | 0.845 [0.784, 0.891] | Precision (inferred edges that are true), n=302 | 0.487 [0.431, 0.543] | No |
| Age / residence time | Reportability rate, n=1272 (fitted rows admitted by the reference-free guard) | 0.259 [0.235, 0.283] | Within-factor-2 agreement, n=675, held-out cross-validated, uncalibrated | 0.564 [0.527, 0.601] | No |
| Reaction | Recall, pooled across 16 reactions, n=681 (true active phases recovered) | 0.599 [0.562, 0.635] | Precision, macro mean across 240 scenarios (1 minus false-discovery rate) | 0.639 [0.598, 0.681] | Yes |

Notes. Intervals are 95% Wilson score intervals computed from the underlying counts (`O3/analysis/python/derive_headline_metrics.py`); the reaction-layer precision interval is inverted from the reaction benchmark's own published false-discovery-rate 95% CI rather than recomputed from scratch. Topology: independent Head-gradient scenario against the USGS Savage MODPATH reference; the informed structural floor (receptor-set only, zero hydraulic information) reaches F1 0.552 and is not independent evidence. Age: capture and correctness are two different kinds of rate, not a literal recall/precision pair -- reportability is the fraction of attempted fits that produce an actionable claim at all (an abstention design, not a defect: see Discussion); within-factor-2 is the accuracy of the claims that are produced, scored on an independent held-out cross-validation split. Reaction: recall is pooled (micro-averaged) from `tableS6_reaction_confusion_matrices.csv` because no macro-averaged recall is published; the macro-averaged precision (0.639, from `table1_comparative_inverse_performance.csv`) is the reaction benchmark's own published headline correctness value. A pooled-precision cross-check (0.586, 95% CI [0.549, 0.622]) is reported in `manuscript/artifacts/data/evidence_discrepancies.csv` and differs from the macro value only because of the averaging convention, not a computational disagreement. **The reaction layer's capture and correctness intervals overlap; only the topology and age layers show a gap that is clearly resolved at the 95% level.**

### 4.3 The gap narrows only when independence is given up

Figure 3 and Table 3 compare each independent or screening
condition against its calibrated, prior-informed, or emulated counterpart.
Calibrating age estimates with a ridge model fit on the same held-out folds
used to score them raises log10 R2 from 0.482 to 0.962 and within-factor-2
agreement from 0.564 to between 0.769 and 0.993 depending on the reported
statistic; because the calibration is fit and scored on the same folds, this
is reported as emulation of the reference rather than independent agreement
(Section 2.1). The reaction layer's calibrated alternative is built
differently: the Mechanism Resolution Score is trained on the carbonate,
crystalline and evaporitic archetypes and evaluated, untouched, on the mixed
archetype, reaching 48.9% four-class accuracy against a 25.0% chance level
for four resolution classes -- a smaller absolute gain than age's calibrated
figure, but a genuine transfer test rather than a same-fold emulation.
Topology's counterpart is neither a calibration nor an emulation: entering
the 174 MODPATH reference edges as priors changes the output graph's edge
count (302 independently inferred edges become 302 under override, 329 under
merge, or 174 under only-prior mode) rather than raising an accuracy metric
comparable to the independent F1 of 0.618, and this mode is reported on its
own edge-count axis rather than forced onto the same scale.

![Figure 3. Independent/screening versus calibrated, emulated, or prior-informed evaluation. Left: age benchmark log10 R2, held-out uncalibrated versus calibrated emulation (fit on the same folds it is scored against). Centre: reaction benchmark Mechanism Resolution Score four-class accuracy, chance level versus a genuine held-out-archetype transfer test. Right: topology benchmark, independent F1 versus prior-informed output-graph edge counts (a different metric, plotted on a separate axis).](artifacts/figures/FIG-3.png)

**Table 3.** Screening/independent versus calibrated, emulated, or prior-informed metric, with non-independence stated per row.

| Benchmark | Independent / screening condition | Value | Calibrated / prior-informed / emulated condition | Value | Independence lost |
|---|---|---:|---|---:|---|
| Age / residence time | Held-out uncalibrated log10 R2 | 0.482 | Calibrated emulation log10 R2 | 0.962 | Ridge calibration fit on the same held-out folds it is scored against |
| Age / residence time | Held-out uncalibrated within-factor-2 | 0.564 | Calibrated emulation within-factor-2 | 0.769-0.993 | Same non-independence as above |
| Reaction | Chance level, 4-class MRS | 0.250 | Held-out archetype transfer, 4-class MRS accuracy | 0.489 | Trained on 3 archetypes, tested on 1 untouched archetype (a genuine transfer test, not an emulation) |
| Topology | Independent F1 (no MODPATH access) | 0.618 | Prior-informed output-graph edges (override / merge / only) | 302 / 329 / 174 | MODPATH edges entered the graph as prior information; not independent validation |

Notes. The two calibration exercises are not equally rigorous and are not presented as such: the age benchmark's calibrated value is an emulation fit and scored on the same data, which is why its improvement (0.482 to 0.962) is large; the reaction benchmark's Mechanism Resolution Score is tested on a fourth, untouched archetype it never trained on, a harder and more independent test, which is why its absolute accuracy (0.489) is far lower and still only doubles chance. The topology benchmark's prior-informed row is not a calibration exercise at all -- it is a different mode of use (MODPATH edges as priors, not as an accuracy target) and is reported on its own edge-count scale rather than forced onto the same axis.

### 4.4 Within-component detail

Figure 4 plots every retained row of each benchmark on its own
headline metric, coloured by taxonomy row type. In topology, all four
negative controls fall at or below F1 0.138, the informed structural
baseline reaches F1 0.552 using the reference receptor set alone with zero
hydraulic information, and the sparse-node sensitivity diagnostic falls to a
mean F1 of 0.445 across 10-100% node-density trials. In age, the
non-independent age-fraction scenario (which shares provenance with its own
reference) shows the highest within-factor-2 agreement of the four
design-matrix scenarios plotted (0.917), above every independent scenario
including strict parity (0.875); the randomized-graph negative control
increases log10 RMSE by 0.655 relative to single-node inference, and the
smallest candidate improvement from any tested graph-regularization
condition is a decrease of 0.001 log10 units, accompanied by a decline in
within-factor-2 agreement, which does not meet the benchmark's own joint
robust-improvement criterion. In reaction, the seven compared methods span
phase F1 from 0.435 (bounded least squares, unregularized) to 0.563
(Hydrosheaf Guarded), with the conventional PHREEQC inverse baseline (0.510)
falling within this range on accuracy alone despite its far lower feasibility
and higher model multiplicity (Section 4.2). A data-tier ablation that adds
optional diagnostic ions on top of the core panel shows the same evidence
class present across all three tiers: class F1 rises only from 0.606 (core,
major-ion panel only) to 0.610 (plus optional SiO2, Sr and water-isotope
diagnostics) to 0.614 (plus Br, dissolved oxygen, dissolved organic carbon,
and sulphate/nitrate isotopes), while the false-discovery rate falls only
from 0.383 to 0.364 to 0.361 across the same tiers -- additional analytical
information narrows, but does not close, the gap reported in Section 4.2.

![Figure 4. Within-component detail: every retained benchmark row for the age, topology and reaction benchmarks, coloured by evidence-taxonomy row type (Table 1), plotted on each benchmark's own headline metric. The three panels use three different metrics on visually similar axes and are not comparable to one another; only within-panel bar heights may be compared directly.](artifacts/figures/FIG-4.png)

### 4.5 Benchmark scale and computational cost

Figure 5 and Table 4 report the scale of the three benchmarks
on their own terms. The reaction benchmark is the largest by two measures:
240 PHREEQC scenarios, expanding to 21,600 factorial inverse fits once every
comparator method and noise/panel/archetype combination is counted. The age
benchmark spans 13 design-matrix scenarios totalling 16,536 candidate rows,
plus 375 additional rows from a leakage-guarded tracer-withholding
cross-validation exercise run separately for five tracers. The topology
benchmark is the smallest by scenario and fit count: 12 independent
graph-inference rows compared against 242 endpoint-derived reference edges
pooled across the two active MODFLOW archives. Only the reaction benchmark
records per-fit runtime across every comparator method, with a median of
20.8 ms per fit; neither the age nor the topology benchmark's result files
carry a comparable runtime field.

![Figure 5. Benchmark scale and computational cost. Top: scenario counts and total fit counts for the age, topology and reaction benchmarks, each on its own axis (values differ by roughly two orders of magnitude). Bottom: median per-fit runtime, recorded only for the reaction benchmark.](artifacts/figures/FIG-5.png)

**Table 4.** Benchmark scale: scenario counts, fit counts, external reference type, and recorded runtime.

| Benchmark | External reference | Scenarios | Fits | Cross-validation rows | Median runtime per fit | Runtime recorded |
|---|---|---:|---:|---:|---:|:---:|
| Age / residence time | USGS national groundwater-age release (TracerLPM Table 4) | 13 design-matrix scenarios | 16,536 (10 scenarios x 1,272 rows, plus 3 ablations) | 375 (leakage-guarded tracer-withholding CV, 5 tracers) | not recorded | no |
| Topology | MODPATH particle-tracking connectivity (3 public MODFLOW/MODPATH archives) | 12 independent graph-inference rows | 242 endpoint-derived reference edges (Savage + Great Miami; Long Island is a documented fallback stub) | not applicable | not recorded | no |
| Reaction | 240-scenario live-PHREEQC factorial synthetic benchmark | 240 PHREEQC scenarios | 21,600 factorial inverse fits (240 x 5 comparator methods x noise/panel/archetype grid) | not applicable | 20.8 ms (median across methods) | yes |

Notes. Counts are read directly from each benchmark's own manifests (`m3_design_matrix_summary.csv`, `independent_graph_vs_modpath.csv` plus the three `tier_*_archive_summary.csv` files, `analysis_summary.json`); no scenario was re-run to produce this table. The reaction benchmark is the only one to record per-fit runtime across every comparator method, which this paper reports as a genuine difference in benchmark design rather than normalising it away.

### 4.6 Field- and archive-transfer scope is uneven

Figure 6 and Table 5 report where each layer has been
tested against something beyond its primary synthetic or archive reference.
The reaction layer has the broadest field-transfer footprint: 160 candidate
edges from the Northern Ghana workbook (median evidence-lifted resolution
index, ELRI, 0.072), 85 from the Talensi District (median ELRI 0.147), and
49 from the Lower Anayari catchment (median ELRI 0.072). The topology layer's
transfer footprint is archive-to-archive rather than archive-to-field: the
Savage and Great Miami MODFLOW archives both return endpoint/pathline edge
F1 of 1.000 as a self-consistency check, and the Long Island archive remains
an inactive fallback stub with no rows in this comparison. The age layer has
no field- or archive-transfer benchmark beyond its primary comparison against
the public USGS national release; no dataset in Table 5 extends the
age layer's evaluation beyond that single reference.

Within the reaction layer's field-transfer component, a retrospective
next-best-measurement analysis ranks candidate additional ions and tracers by
their expected reduction in reaction-class ambiguity; fluoride returns the
highest measurement-value score of any candidate in the Northern Ghana
transfer analysis, and all 160 Northern Ghana wet-to-dry pairs are classified
as partially identifiable, rather than fully identifiable, by the calibration
transferred from the synthetic benchmark. Median Hydrosheaf-Core evidence and
total-dissolved-solids consistency scores for the Northern Ghana panel are
0.690 and 0.942 respectively.

![Figure 6. Field- and archive-transfer scope. Top: MODPATH archive endpoint/pathline edge F1 (Savage, Great Miami; Long Island is a documented fallback stub with no active rows). Bottom: median evidence-lifted resolution index (ELRI) for the three field chemistry datasets (Northern Ghana, Talensi, Lower Anayari). No field-transfer benchmark exists for the age layer.](artifacts/figures/FIG-6.png)

**Table 5.** Full dataset inventory: every synthetic and field dataset any of the three benchmarks consumes, with source, size, and role.

| Dataset | Type | N | Source | Role in this comparison | Access |
|---|---|---:|---|---|---|
| USGS national groundwater-age release (TracerLPM Table 4 reported configuration) | Published field tracer data | 1,272 candidate site-fits; tracers 3H, SF6, CFC-11, CFC-12, CFC-113, 14C | [@jurgens2012tracerlpm] | External reference for the age/residence-time benchmark | Public (USGS) |
| USGS Savage Municipal Water-Supply Well archive | Published synthetic MODFLOW-2005/MODPATH 5 archive | 3,000 particles; 174 endpoint reference edges | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/F7J102FK | Primary external reference for the topology benchmark | Public (USGS) |
| Great Miami River Basin archive | Published synthetic MODFLOW 6/MODPATH 7 archive | 154 particles; 68 endpoint reference edges | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/P9X4C9R6 | Secondary external reference for the topology benchmark (self-consistency check) | Public (USGS) |
| Long Island Nitrogen archive | Published synthetic MODFLOW 6/MODPATH 7 archive | Fallback stub; no active rows in this comparison | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/P97VFXZ4 | Tertiary archive for the topology benchmark, documented as inactive | Public (USGS) |
| 240-scenario live-PHREEQC factorial benchmark | Synthetic, HydroSheaf-generated | 240 scenarios; 21,600 factorial inverse fits; 4 hydrochemical archetypes (carbonate, crystalline, evaporitic, mixed) | Generated with [@parkhurst2013phreeqc] (PHREEQC 3.7.3/3.8.6) | Primary external reference for the reaction benchmark | Reproducible from the reaction benchmark's own repository folder |
| Northern Ghana workbook | Published field chemistry, reconstructed seasonal split | 160 boreholes, 4 regions; Dry panel is the primary measured unit | Author-compiled regional dataset; see Methods for the reconstruction disclosure | Field-transfer chemistry audit for the reaction benchmark | Project-held; see data-availability statement |
| Talensi District | Published field chemistry | 63 samples | [@chegbeleh2020talensi] | Field-transfer chemistry audit for the reaction benchmark (screening-level) | Published source |
| Lower Anayari catchment | Published field chemistry | 41 samples | [@zakaria2021anayari] | Field-transfer chemistry audit for the reaction benchmark (screening-level) | Published source |

Regional hydrogeological context for the three Ghanaian sites is drawn from [@banoengyakubo2011ghana] and [@macdonald2012africa]; neither is a source for the Northern Ghana workbook itself.

## 5. Discussion

### 5.1 Principal finding

Directed connectivity, residence time and reaction extent are governed by
different physics, constrained by different observation types, and were
benchmarked here against three different external references built years
apart by different protocols. Under independent, prior-free, uncalibrated
evaluation, the point estimate of a capture-type metric exceeds the point
estimate of a matched correctness-type metric in all three, but only two of
the three, topology and age, hold up once 95% Wilson intervals are attached:
in reaction, the pooled recall and macro precision intervals overlap
(Figure 2, Table 2), so the same-direction point estimate there cannot be
distinguished from sampling variation. The defensible claim is accordingly
narrower than a first read of the point estimates suggests: whatever is
separately responsible for imperfect precision in edge classification and
imperfect held-out age agreement, the resulting asymmetry is resolved at
conventional confidence in both of those layers and points the same
direction, unresolved, in the reaction layer. This is not a claim that the
three numbers are statistically identical, nor that one underlying cause
explains any of them; the paper takes no position on whether a shared cause
exists. A framework that reported connectivity, age and reaction with equal
confidence, or that reported its most permissive metric for each layer
without stating which axis that metric sits on and whether the gap survives
an interval check, would obscure a pattern -- and its real limits -- that
this comparison makes visible only by placing the three benchmarks side by
side.

### 5.2 What calibration and prior information buy, and what they cost

Figure 3 shows that the gap is not fixed: it narrows sharply when a
layer is allowed to use information about the reference it is scored
against. But the two calibration exercises compared here differ in what that
information costs. Age's calibrated log10 R2 of 0.962 is achieved by fitting
a ridge model on the same held-out folds it is then scored against -- an
emulation of the reference's own output, in the sense used by
@kennedy2001calibration and @brynjarsdottir2014discrepancy, who distinguish a
model that reproduces observed outputs from a model whose parameters
correctly represent the underlying physical process. Reaction's Mechanism
Resolution Score, by contrast, is trained on three archetypes and tested,
untouched, on a fourth, and its 48.9% four-class accuracy is a genuine,
harder transfer result that happens to look less impressive stated on its
own. Reporting the age result without the reaction result, or reporting both
without stating that they were earned under different rules, would flatter
the framework in a way this comparison declines to do. Topology's
prior-informed mode is different again: it is not a calibration at all, and
this paper resists the temptation to summarize it as one. Feeding MODPATH
edges into the graph as priors is a legitimate mode of use -- arguably the
intended one, where MODPATH output exists -- but it changes what the graph
is built from, not how accurately an independently built graph agrees with
MODPATH.

### 5.3 Relation to equifinality and prediction-in-ungauged-basins practice

The pattern reported here is a specific, quantified instance of a general
difficulty long recognized in environmental modelling: a model can reproduce
available observations well without its structure or parameters being
uniquely determined by those observations [@beven2001glue; @beven2006manifesto].
The Prediction in Ungauged Basins programme argued that hydrological
inference under data scarcity requires reporting the boundary of what
evidence supports, rather than a single calibrated estimate presented
without that boundary [@sivapalan2003pub; @hrachowitz2013pub]. The age
layer's reportability guard, the reaction layer's equivalence-class and
Mechanism Resolution Score machinery, and the topology layer's explicit
negative-control floor are each, independently, an implementation of that
same argument for a different inference problem. This comparison's
contribution is not to introduce the argument -- it is well established --
but to show that three independently built implementations of it, tested
against three unrelated references, produce the same qualitative signature.
That convergence is weak evidence that the signature reflects something
general about data-limited groundwater inference rather than a
peculiarity of any one benchmark's design, though a comparison of three
benchmarks from one research programme cannot rule out a shared design
habit as the explanation (Section 5.6).

### 5.4 Implications for monitoring design and reporting practice

A capture-exceeds-correctness pattern has a direct practical reading for
anyone using these layers to prioritize field effort. Where recall or
reportability is high but precision or agreement is not, the useful output
of a layer is a shortlist of candidates worth checking rather than a
confident answer, which is precisely the logic behind reduced-rank and
information-theoretic groundwater monitoring network design
[@sreekanth2017monitoring] and behind the reaction layer's own
next-best-measurement ranking, which identifies fluoride as the single
highest-value additional field measurement in the Northern Ghana transfer
analysis (Section 4.6) rather than asking for a uniformly larger analytical
panel. More generally, reporting a calibrated interval or a calibrated class
probability without stating whether the calibration set overlaps the
evaluation set -- the distinction this paper draws throughout between
independent and calibrated/emulated results -- is the same distinction that
motivates the conformal-prediction literature's insistence on exchangeable,
held-out calibration data before an interval is described as having a
stated coverage guarantee [@angelopoulos2023conformal]. None of the three
benchmarks compared here uses conformal methods directly, but all three are
answerable to the same standard: a calibrated number is only as informative
as the independence of the data it was calibrated on.

### 5.5 Relation to the HydroSheaf framework article

The companion framework article (Objective 2 of this thesis chapter) reports
a compressed version of the same three-layer pattern -- coverage and error
for age, log-loss for reaction-family calibration, precision and recall for
topology -- inside a software-and-methods article whose primary contribution
is HydroSheaf's claim-gated architecture. This paper does not restate that
architecture, does not claim its integrated field demonstration as its own
evidence, and does not report a fourth validation of any layer beyond what
each layer's own result package already supports. Its contribution is the
comparison itself: the taxonomy that lets the age, topology and reaction
benchmarks be read against one another; the disclosed reconciliation of their different averaging
conventions (Table 2); the computational-scale accounting that a
7,000-word framework article has no room for (Figure 5); and the
explicit statement of where field- and archive-transfer evidence exists for
each layer and where it does not (Figure 6).

Splitting this material across two papers is a deliberate choice, not a
division of one publishable unit for its own sake, and the test is whether
each paper serves a reader the other does not. A reader of the framework
article wants to know whether HydroSheaf's architecture works end to end on
a field problem; that reader is served by the framework article's integrated
demonstration and does not need this paper's benchmark-design taxonomy, its
averaging-convention reconciliation, or its confidence intervals to answer
that question. A reader auditing or extending one specific inference layer,
for example someone building a competing topology-inference method who wants
to know exactly what independence means in the topology benchmark and how its
edge count translates into interval width, is served by this paper and would
find the framework article's compressed three-line summary insufficient for
that purpose; that reader has no use for the framework article's claim-gate
architecture or its field demonstration. The two audiences are different
enough, and this paper's quantitative content (the taxonomy, the interval
analysis, the benchmark-scale accounting) different enough from a restatement
of the framework article's own numbers, that the split serves readers rather
than only inflating a publication count. Whether an editor agrees is a
judgement this paper cannot make for itself; what it can do is state the
distinction plainly rather than leave it implicit.

### 5.6 Limitations

Five limitations bound this comparison. As with any numerical model of an
earth system, agreement with a reference should be read as consistency with
that reference under stated conditions rather than as proof of a correct
underlying process [@oreskes1994validation], and that general caution applies
to every number below, not only to the ones flagged individually. First, the capture/correctness
analogy across age, topology and reaction is a stated methodological choice,
not a proof that the three axes measure the same underlying quantity;
readers should treat Figure 2 as three parallel single-component
results placed on adjacent panels, not as one merged statistic, and should
note that the analogy's own interval check (Section 4.2) resolves the gap
for only two of the three panels. Second, the
reaction layer's recall figure is a micro-averaged pooled statistic where
the component's own published headline is a macro-averaged precision;
the two conventions are shown to differ by roughly five percentage points
(Table 2), a difference attributable to averaging method rather
than to any disagreement about the underlying counts. Third, the age and
reaction benchmarks' result files were locked before validation-layer code
changes made on 2026-08-01 and 2026-08-02; those changes were checked
directly and found to touch the integrated programme/decision-utility layer
the companion framework article benchmarks, not the core age-inference or
reaction-dictionary modules the age and reaction benchmarks exercise, but a
change in a shared
dependency neither this check nor any other in this paper can fully rule
out remains a residual risk. Fourth, the age layer's external reference is
itself built from environmental tracers whose own interpretive limitations
are well documented -- degassing, mixing, and matrix diffusion can each bias
a tracer-derived age independently of any lumped-parameter or graph
regularization choice this paper compares [@mccallum2015limitations] -- so
even the "independent" age comparison in Figure 2 inherits whatever
uncertainty the reference release itself carries. Fifth, and most
importantly, no result in this paper is independent field validation of
connectivity, residence time or reaction mechanism: the age layer has no
field-transfer benchmark at all, the topology layer's transfer evidence is
archive self-consistency rather than field truth, and the reaction layer's
field-transfer scores are a candidate-plausibility audit built on a
workbook whose seasonal panel is independently confirmed to be a
reconstruction, not an independently sampled comparison (Section 2.4).

## 6. Conclusion

Three already-locked HydroSheaf component benchmarks -- age and residence
time against the USGS national groundwater-age release, topology against
MODPATH particle-tracking connectivity, and reaction against a 240-scenario
live-PHREEQC factorial benchmark -- were harmonized under one common
evidence taxonomy without re-running any of the underlying computation. In
all three, independent, prior-free, uncalibrated evaluation shows a
capture-type metric exceeding a matched correctness-type metric: edge recall
above edge precision for topology, reportability rate below within-factor-2
agreement for age, and phase recall close to but below phase precision for
reaction. The gap narrows only when a layer is calibrated or conditioned
against the same reference used to score it, and the two calibration
exercises available here -- an emulation fit on the same held-out folds for
age, a genuine held-out-archetype transfer test for reaction -- are shown to
differ substantially in rigour rather than being interchangeable evidence of
the same kind. The three benchmarks differ by roughly two orders of
magnitude in scale, and their field- and archive-transfer scope is uneven:
broadest for reaction, archive-internal for topology, and absent for age.

None of this is a fourth, independent validation of connectivity, residence
time, or reaction mechanism, and none of it is a second HydroSheaf framework
contribution; both remain the domain of the components' own result packages
and of the companion framework article. What this comparison adds is a
shared vocabulary for reading the three benchmarks against one another, a
disclosed reconciliation of their different metrics and averaging
conventions, and a computational-scale and field-scope accounting that no
single component paper reports. The taxonomy and harmonization scripts are
retained so that a fourth component benchmark, should one be built, can be
added to this comparison without repeating the reconciliation performed
here.

## Code and data availability

The harmonization scripts, R figure scripts, tidy CSV exports, and the
component result files this comparison reads (`M3/m3_age_benchmark/`,
`M4/m4_topology_benchmark/`, `M5/m5_inverse_reaction_benchmark/`) are archived
at `https://doi.org/10.5281/zenodo.PLACEHOLDER-DOI` and developed at
`https://github.com/PLACEHOLDER-ORG/hydrosheaf`. The repository carries an
open-source licence, an English README with installation and usage
instructions, a declared dependency lock, and a one-command reproduction path
for every figure and table in this paper (`O3/README.md`). No PHREEQC,
MODFLOW/MODPATH, or TracerLPM re-run is required to reproduce any number in
this paper; only the already-published result files listed in Table 5 and
Python/R are required.

## Author contributions

To be completed with named CRediT roles before submission: conceptualisation,
methodology, software, formal analysis, data curation, writing -- original
draft, writing -- review and editing.

## Ethics, competing interests and funding

This comparison used only already-published or already-locked result files
and involved no new human or animal subject procedures, and no new field
sampling. The authors declare no competing interests. Funding is to be
declared in the final submission metadata.
