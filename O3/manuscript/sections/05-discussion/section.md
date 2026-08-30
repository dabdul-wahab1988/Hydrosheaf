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
