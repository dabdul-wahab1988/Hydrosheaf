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
