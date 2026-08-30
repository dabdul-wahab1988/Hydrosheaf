# Conditional Value of Graph Priors in Multi-Tracer Groundwater-Age Inference: A Hydrosheaf Benchmark of Public Aquifer Datasets

## Abstract

How old is the water coming out of a well? The answer usually comes from a handful of chemical tracers, and different tracer combinations can support very different answers for the same sample, a long-standing problem in groundwater dating. One proposed fix is to let nearby wells inform one another's estimates, on the assumption that age should vary smoothly along a flow path. The objective here was to test that assumption using Hydrosheaf, an open framework, applied to 1,272 wells from a public United States Geological Survey (USGS) release; published USGS ages, being model outputs rather than directly observed truths, were treated only as a comparison target. After screening for trustworthy fits, 329 of 1,272 wells qualified, and Hydrosheaf's ages tracked the published values closely (87.5% within a factor of two). The main finding is that letting neighbouring wells influence one another did not help: the best-performing connection pattern nudged one error measure in the right direction negligibly while making a practical accuracy measure worse, and two deliberately wrong connection patterns, reversed and randomised, made things clearly worse, as intended; an independent method reached the same conclusion. A further, unresolved finding emerged along the way: for roughly one well in four, no single age distribution could reconcile all of its tracers within measurement uncertainty, and none of seven candidate explanations tested could account for it. The main contribution and practical implication is that Hydrosheaf currently supports pressure-testing and rejecting groundwater-age assumptions, not a replacement for existing tracer-dating methods, and not field-validated proof of any specific flow path.

## 1. Introduction

Knowing how long water has been underground tells you a great deal about it: whether a well is likely to pick up a pesticide applied last decade, whether a contamination problem will clear up quickly once its source is removed, and how far a chemical signal has travelled through the rock. This residence time, or groundwater age, links the physical structure of an aquifer to the water-quality problems that show up at the tap [@cook2000determining; @jurgens2022over].

The trouble is that age is never measured directly. It is inferred from environmental tracers: trace gases and isotopes that entered the water at recharge and have been decaying, accumulating, or otherwise changing ever since. Tritium (^3^H) and the helium it decays into mark water recharged since the 1950s nuclear-testing era. Sulfur hexafluoride (SF~6~) and chlorofluorocarbons (CFCs) are industrial gases with well-documented atmospheric histories, useful for the same few decades. Carbon-14 extends the record to thousands of years, and helium-4 accumulates over tens of thousands more. Each tracer requires its own correction assumptions [@darling2012practicalities; @han2016review], and a real sample rarely carries a full, unambiguous set of them. As a result, the same measured concentrations can often be explained by more than one story about how the water arrived: a young pulse, an old and dispersed mixture, or some blend of both [@massoudieh2014bayesian; @mccallum2014bias]. This is age equifinality, and it is the central difficulty this study addresses.

The standard response, formalised in tools such as the U.S. Geological Survey's TracerLPM workbook [@jurgens2012tracerlpm], is to fit a small family of idealised age distributions to each well's tracers independently. This works well and remains the field's reference method, but it treats every well as its own island. It does not ask whether the fitted age at one well is consistent with the fitted age at its neighbour, even when nearby wells plausibly sit on the same flow path.

Hydrosheaf is a broader computational framework for graph-based hydrogeochemical inference, of which this paper tests one component [@abdulwahab2026hydrosheaf]. The idea being tested is simple to state: represent wells as nodes in a network, connect nodes that plausibly lie along the same flow path, and let that network gently penalise implausible age reversals, such as a downstream well appearing much younger than the well upstream of it. If the network reflects real hydrogeology, this should sharpen noisy age estimates. If it does not, it can just as easily impose false patterns on wells that have nothing to do with one another. Graph-based smoothing of this kind is well established in machine learning [@belkin2006manifold], but its value in any particular application depends entirely on whether the graph encodes real structure, and that must be demonstrated, not assumed.

This paper treats that as an open, testable question rather than a premise. A benchmark was built that compares single-well age estimates against several network variants, including two deliberately wrong ones, a reversed-direction network and a randomised network, that share the same structure as the real candidates but none of the underlying logic. A genuinely useful network should beat the single-well baseline and comfortably beat both deliberately wrong versions; if it does not, the network is not earning its place in the analysis. A further test addresses something the reference-age comparison alone cannot: whether the network can predict a tracer that was deliberately hidden from it, which rules out the network simply re-deriving an answer it was quietly given [@oreskes1994verification].

Since M3 was first written, an independent and methodologically distinct check of the same question became available: rather than fitting one best-guess age distribution per well, it asks which age distributions, if any, remain mathematically consistent with the tracers once measurement uncertainty is taken into account, and whether a candidate network is compatible with those possibilities. This diagnostic remains at a development, exploratory stage, but it allows the same graph-benefit question to be tested a second way, and, in building it, a separate and unresolved finding emerged about how often a well's own tracers fail to agree with each other at all, independent of any network. Both are reported, clearly labelled as exploratory, alongside the main benchmark.

The objectives are, first, to establish which fitted ages can be trusted at all, given that tracer scarcity and multiple candidate age models make some fits unreliable by the framework's own criteria; second, to test, with negative controls, whether a flow-path network improves on treating each well independently; third, to test the same question again through target-tracer withholding, which does not depend on the published ages being correct; and fourth, to report, honestly and without overclaiming, what the newer exploratory check adds. Throughout, the published USGS ages are treated as a comparison target derived from a well-established model, not as an independently verified truth, and no candidate network is claimed to represent a confirmed groundwater flow path.

## 2. Materials and methods

### 2.1 Design philosophy

The benchmark was designed to test, rather than assume, whether letting nearby wells inform one another's age estimates helps. A candidate flow-path network was therefore always evaluated against a single-well baseline and against two negative controls built to share the same size and shape as the real network but none of its hydrogeological logic: a network with every connection pointing the wrong way, and a network connecting wells at random. A network result was accepted as a genuine improvement only when it beat the single-well baseline on both a statistical-agreement measure and a practical accuracy measure at once, and only when the negative controls did not show a comparable effect. Published USGS ages were used throughout as a comparison target for evaluating agreement, not as directly observed true ages, because they are themselves the product of a modelling workflow [@jurgens2012tracerlpm].

### 2.2 Data

The primary dataset was the USGS national public-supply aquifer release for the continental United States, 2004-2017 [@jurgens2022data], covering 1,279 wells across 21 principal aquifers, with tracer measurements, dissolved-gas model outputs, published age estimates, and site metadata. After harmonising sample identifiers, dates, coordinates, and tracer fields across the release's component tables, 1,272 rows remained (Figure 1; Supplementary Table S1). Two further public releases were examined but not pooled into the main comparison: a Western Principal Aquifers release [@faulkner2019data] that lacks the coordinates needed for the network benchmark, and a Mississippi River Valley alluvial-aquifer release [@gratzer2025groundwater] that lacks the published age table needed for a like-for-like comparison.

The tracer suite spanned modern to very old water: tritium (^3^H) and tritiogenic helium-3 for water recharged since roughly 1950, sulfur hexafluoride and chlorofluorocarbons for the same young-water window, carbon-14 for water up to several tens of thousands of years old, and helium-4 as a supporting old-water indicator (Table 1). The tritium input history used to interpret ^3^H was built from precipitation isotope records (Figure 2; Supplementary Methods S1).

### 2.3 Age model

Groundwater age was represented by a transit-time distribution, a mathematical description of the mixture of travel times contributing to a sample, drawn from a small family of shapes already standard in tracer hydrology (piston-flow, exponential, dispersion, gamma, and related mixtures) [@jurgens2012tracerlpm]. For a given transit-time distribution and a tracer's known atmospheric or decay history, the concentration expected at a well followed

<!-- EQ:EQ-1 -->
$$
\widehat{C}_j(t_s, \theta) = \int_0^{\tau_{\max}} I_j(t_s - \tau)\, \exp(-\lambda_j \tau)\, g(\tau \mid \theta)\, \mathrm{d}\tau,
$$

where $\widehat{C}_j$ is the predicted value of tracer $j$, $t_s$ is the sample year, $I_j$ is that tracer's input history, $\lambda_j$ is its radioactive decay constant (zero for stable tracers), $g(\tau \mid \theta)$ is the transit-time distribution with parameters $\theta$, and $\tau_{\max}$ is the oldest age considered. All tracer observations available for a well were fitted jointly, so that a single set of parameters had to satisfy every measured tracer at once, weighted by each tracer's measurement uncertainty. Full derivations, the age-fraction and graph-prior penalty terms, and the discretisation used to evaluate the integral are given in Supplementary Methods S1-S4.

Because several plausible transit-time distributions can fit the same tracers almost equally well, every fit was passed through a reference-free reportability guard before being used in any accuracy comparison: a fit had to converge, use a distribution family actually supported by the reported USGS workflow where one was requested, be constrained by more tracer observations than free parameters, and leave no comparably good alternative fit more than roughly a factor of three away in age. Fits failing this guard were retained in the dataset but excluded from every reference-agreement result; published ages played no part in the guard itself (Supplementary Methods S2).

### 2.4 Flow-path network

Wells were represented as nodes and hypothesised flow connections as directed edges, built from metadata available before any comparison: geographic coordinates, aquifer and study-area membership, well depth, and a hydraulic-gradient proxy, at increasingly restrictive levels of physical plausibility. A directed edge encoded an expected upstream-to-downstream or shallow-to-deep ordering, not a confirmed flow path. Where an edge's expected ordering was violated, the network gently pulled the downstream age toward the upstream age plus an expected travel increment, at one of four pre-declared strengths (none, weak, medium, strong; Supplementary Methods S4). The wrong-direction control reversed every real edge; the randomised control connected the same number of same-sized groups of wells by chance (fixed random seed, for exact reproducibility). Edge distance and vertical separation varied systematically by construction rule, from a median of roughly 14 km for the hydraulic-proxy network to roughly 284 km for the randomised control (Supplementary Figure S2), which is expected given how each network family was built and was not itself treated as evidence for or against any candidate.

### 2.5 Scenario matrix and metrics

The full comparison matrix (Table 2) included the reported-configuration and reported-age-fraction emulation scenarios, an automatic model-selection scenario, several tracer and correction ablations, and the network variants above, applied to all 1,272 wells with the reportability guard applied throughout. Direct comparisons between scenarios used only the rows both scenarios could evaluate in common; results on each scenario's full available support are reported separately and are not directly rankable against each other because their supports differ (Table 3).

Accuracy against the published ages was assessed on a base-ten logarithmic scale, because groundwater ages span several orders of magnitude, using the median absolute log-error, the root-mean-square log-error, and the share of estimates within a factor of two or ten of the published value (Supplementary Methods S5). A network result required a lower root-mean-square error and unchanged-or-better factor-of-two agreement simultaneously to be accepted as an improvement.

A second, more demanding test avoided the published ages altogether. For every eligible well, one tracer at a time was removed before model fitting and before network construction, so that neither the single-well fit nor any connected neighbour could see the withheld measurement, and the withheld tracer's own concentration was then predicted and compared with its true measured value. This target-withholding test cannot be satisfied by a network simply reproducing information it was quietly given.

### 2.6 An independent, exploratory cross-check

A separate, more recently developed diagnostic approached the same network question from a different angle: rather than fitting one best age distribution per well, it asked which age distributions, if any, are mathematically consistent with a well's tracers once measurement uncertainty is taken into account, and tested whether a candidate flow-path edge was compatible with the possibilities at both of its wells. This diagnostic is explicitly at a development and exploratory stage; it was not used to revise any scalar result above, and it is reported separately with that status stated throughout (Supplementary Methods S6-S7).

### 2.7 Reproducibility

The full pipeline was rerun end to end from the current codebase for this revision, using a fixed 90-step age grid, a fixed random seed for the randomised control, and machine-readable manifests recording software versions and file checksums for every input and output. Every regenerated result was compared line by line against the previously locked benchmark; the outcome of that comparison, and the small number of cases showing ordinary floating-point-level variation rather than a substantive change, are reported in Supplementary Methods S8.

## 3. Results

### 3.1 How many wells could be trusted at all

Of the 1,272 harmonised wells, 329 passed the reportability guard under the scenario that most closely follows the USGS reported workflow (Table 2). A further 289 passed under a scenario that additionally used the USGS-reported age-fraction breakdown, and 309 passed when Hydrosheaf was left to select its own best-fitting model rather than following the reported configuration. These are three different, overlapping subsets, not three views of the same 329 wells, and because the age-fraction scenario reuses numbers from the same USGS release as the comparison ages, its closer agreement reflects shared bookkeeping rather than an independently stronger result. Roughly a quarter of all wells were unreachable by any scenario tested, most often for lack of enough tracers to outnumber the free parameters in even the simplest age model. Carbon-14 and helium-4, the two old-water indicators, agreed with each other at 176 wells, conflicted outright at 3, and showed only partial tension at 8; the remaining wells carried only one of the two or neither, so no old-water constraint could be cross-checked at all (Supplementary Table S3). These categories describe how often the available old-water evidence agrees with itself, not an independent check on whether the resulting ages are correct.

### 3.2 Agreement with the published ages

Across scenarios, reference-workflow agreement and reportability varied together (Figure 3; Table 2). On its 329 reportable wells, the reported-configuration scenario matched the published USGS ages with a median error of a factor of 1.07 (0.0279 on the log scale) and 87.5% of estimates within a factor of two (Figure 4). Agreement was closest for old and intermediate-age water and weakest for modern water under 50 years old, where several tracers can point in different directions at once (Supplementary Table S2). Letting Hydrosheaf select its own model, rather than following the reported configuration, performed worse on its own reportable subset (67.3% within a factor of two), which is expected: automatic selection was not built to outperform a manually curated reference workflow, and the two subsets are not the same wells, so this is not a fair one-to-one contest. On the 40 wells common to both the reported-configuration and age-fraction scenarios, the age-fraction scenario's closer agreement was statistically distinguishable from the reported-configuration scenario by both a paired Wilcoxon signed-rank test and a paired t-test (p < 0.05 for each), and a bootstrap estimate of the difference in median absolute error excluded zero (Table 4). This confirms the two scenarios are measurably different from one another; it does not upgrade either one from reference-workflow agreement to independent validation, since both draw on the same USGS release.

### 3.3 Did the flow-path network help?

No. Across four pre-declared connection strengths and seven network variants, the joint improvement test, requiring both a lower error and unchanged-or-better practical accuracy, was met by exactly none of the physically motivated candidates (Figure 5; Table 5). The closest any candidate came was the hydraulic-gradient-informed network at its weakest setting: it nudged the error measure very slightly in the right direction, by about one part in a thousand on the logarithmic scale, but the share of wells within a factor of two of the true value fell at the same time, from 87.5% to 86.9%, so it still failed the joint test. Every other candidate network made the error measure worse outright, by margins fifty to five hundred times larger than that one near-miss.

The two negative controls behaved as they were designed to: the wrong-direction network raised the error measure by roughly 0.10 on the log scale, and the randomised network by roughly 0.65, both far more damaging than any real candidate. That the model responded strongly and consistently to deliberately wrong topology is itself informative: it confirms the network mechanism does something rather than nothing, while also showing that "the network changed the answer" is not, by itself, evidence that a candidate network is correct.

### 3.4 Predicting a hidden tracer

Reference-age agreement can be inflated by scenarios that share bookkeeping with the comparison target, so a second test withheld one tracer at a time and asked whether a network could predict it. Reportable, target-withheld fits were obtained for 121 of 794 eligible tritium measurements, 75 of 262 sulfur-hexafluoride measurements, and 169 of 1,103 carbon-14 measurements (Figure 6). For all three tracers, every candidate network's prediction error was statistically indistinguishable from the single-well baseline, differing by at most a few tenths of one percent; the randomised control, by contrast, was measurably worse for tritium and carbon-14. Chlorofluorocarbon-11 and -12 could not be evaluated this way at all: only four and six wells respectively passed the reportability guard for those tracers, and none of those wells was connected by an eligible network edge, so the comparison is undefined rather than a demonstration that networks and baselines perform equally (Supplementary Figure S3).

### 3.5 An independent check, and an unresolved puzzle

A separate, exploratory diagnostic re-asked the network question using a different method: instead of fitting one age distribution per well, it characterised the full range of transit-time distributions still consistent with a well's tracers under measurement uncertainty, and tested whether a candidate network edge was compatible with that range at both endpoints. Consistent with the main result above, this check likewise found no support for a candidate network sharpening the possibilities beyond what each well already implied on its own.

Building this check surfaced a separate finding that has nothing directly to do with networks: for 975 of 3,501 evaluable well-tracer combinations (27.85%), no non-negative transit-time distribution existed that could reconcile the well's tracers within their stated uncertainty at all (Figure 7). This infeasibility rate rose sharply with the number of tracers being reconciled at once, from 1.5% with a single tracer to 75.4% with five, and concentrated overwhelmingly on specific pairs of tracers rather than being spread evenly: the three chlorofluorocarbons conflicted with one another at 19-33%, and tritium conflicted with its own decay-derived helium partner at 17.9%, while the two most independent tracers, carbon-14 and tritium, almost never conflicted (0.8%). Seven candidate explanations were tested and each was rejected in turn: the effect was not explained by the number of independent constraints alone, by uniformly under-stated measurement uncertainty, by individual measurements lying outside their tracers' physically achievable range, by multiple distinct flow paths mixing in one sample, by a refittable local correction (the relevant correction is supplied upstream by the source release, not fitted locally), by a single shared correction error moving all three chlorofluorocarbons together, or by any one chlorofluorocarbon acting alone. This diagnostic is exploratory and development-stage; it does not revise any scalar result reported above, and no cause for the conflict is asserted here (Supplementary Information S6-S7).

### 3.6 Data-limited experiments

Several planned comparisons could not be evaluated at all with the available data. A test comparing corrected against raw chlorofluorocarbon and sulfur-hexafluoride measurements found zero wells with both values usable side by side, so no gas-correction effect could be estimated. An experiment using only tracers for water older than 1,000 years produced no reportable fits whatsoever, while the corresponding young-tracer-only experiment remained usable but noticeably less accurate (53.8% within a factor of two on 65 wells) than the full reported configuration. Removing helium-4 from an otherwise complete tracer set changed the error measure by an amount indistinguishable from zero on the 65 wells where the comparison was possible, and a paired comparison of a raw-carbon-14 ablation against the reported correction likewise showed no detectable effect on its 61 common wells (Table 6). These results describe the limits of what this particular dataset can support; they do not establish that any tracer is generally unimportant.

A controlled synthetic check, using simulated wells with known true ages rather than the USGS data, confirmed separately that the network mechanism can, in principle, resolve a specific kind of ambiguity, a tritium bomb-pulse measurement consistent with two different recharge dates, when flow ordering is known exactly (Supplementary Figure S4). This demonstrates the mechanism functions as intended; it is not evidence that the same benefit carries over to real, imperfectly known networks, which is precisely the question the main benchmark above addresses.

### 3.7 Summary

Six results, taken together, define what this benchmark supports and does not. First, only about a quarter of wells produced a fit trustworthy enough to use in any accuracy comparison. Second, on those wells, Hydrosheaf's ages tracked the published USGS values closely, understood as agreement with a modelled reference rather than as proof of true age. Third, no flow-path network met the pre-declared joint improvement test, while both deliberately wrong networks were clearly more harmful, confirming the network mechanism is doing something without validating any one candidate. Fourth, predicting a tracer deliberately hidden from the network told the same story. Fifth, an independently built, exploratory check reached the same network conclusion by a different method. Sixth, that same check surfaced an unresolved case where a quarter of tracer combinations cannot be reconciled by any transit-time distribution at all, for a reason seven tested explanations failed to identify.

## 4. Discussion

### 4.1 What a flow-path network adds, and what it does not

Well-by-well age-dating tools such as TracerLPM remain valuable precisely because they are interpretable and economical with data [@jurgens2012tracerlpm]. Their limitation is not that they are wrong, but that they never ask whether a well's fitted age makes sense next to its neighbour's. This study built that missing check and, having built it, found that it does not currently pay off on this dataset: no candidate network met the pre-declared bar for a genuine improvement, on either the reference-agreement test or the harder hidden-tracer test, and an independently built exploratory method reached the same verdict by a different route. The most useful reading of that result is not "networks do not work," but "this particular way of proposing connections between these particular wells did not contain information the tracers had not already supplied." Coordinate proximity, aquifer membership, and depth are reasonable starting hypotheses about which wells share a flow path, but reasonable is not the same as correct, and the results here suggest they were not correct often enough, on this dataset, to earn a place in the analysis.

The negative controls make that reading more, not less, informative. A reversed-direction network and a randomised network are structurally identical to the real candidates and share none of their hydrogeological reasoning, and both made age estimates measurably worse, the randomised control especially so. A method that could not tell a plausible network from a nonsensical one would be useless as a diagnostic; this one clearly can. What it cannot yet do, on this dataset, is turn a plausible network into a demonstrably better one. That distinction matters for how the framework should be read: as a tool that can reject bad topology with some confidence, not yet as one that can confirm good topology.

### 4.2 Why the near-miss matters

The hydraulic-gradient-informed network deserves a second look precisely because it came closest to succeeding. Its effect on the error measure was real but tiny, and its failure came entirely from the accompanying drop in practical accuracy. A less careful analysis, using only the error measure and without a pre-declared joint rule, could have reported this as a positive result. Requiring both measures to move favourably at once, and checking that finding against negative controls before interpreting it, is what prevented that. This is offered less as a claim about this specific network than as an illustration of why a benchmark needs a rule fixed in advance, rather than a rule discovered after seeing which measure looks best.

### 4.3 An unresolved tracer conflict, honestly reported

The clearest new result in this revision is not about networks at all. Roughly one well-tracer combination in four cannot be reconciled by any physically sensible transit-time distribution, and the pattern is concentrated: the three chlorofluorocarbons disagree with one another far more than any of them disagrees with an unrelated tracer such as carbon-14. Seven explanations were tested for this and every one failed, including the two that looked most promising going in, a shared correction error moving all three chlorofluorocarbons together, and one chlorofluorocarbon behaving badly on its own. A natural next step, not available with the data at hand, would be to check whether the conflict tracks reducing chemical conditions, since chlorofluorocarbon-11 is known to degrade under low-oxygen conditions; that check needs dissolved-oxygen or related measurements at the same wells and same times as the tracer measurements, which this release does not provide. Reporting a well-documented, thoroughly narrowed-down puzzle without a resolution is a deliberate choice here: incompletely explained results are common in tracer hydrology, and stating plainly what was ruled out is more useful to the next researcher than quietly setting the finding aside.

### 4.4 Implications for aquifers with sparse data

For an aquifer where wells are few and tracers are incomplete, a defensible reading of these results is that Hydrosheaf currently earns its keep as a screening and rejection tool rather than as a stand-alone dating method: it flags fits too poorly constrained to trust, and it can reject a proposed flow-path network with reasonable confidence when the network performs no better than an obviously wrong one. Whether it can go further and confirm a genuinely useful network is still open, and the present results argue for testing that question locally with the same discipline used here, rather than assuming a network of convenience will help simply because it is hydrogeologically plausible on paper. This is especially relevant for aquifers in semi-arid and lower-resource settings, including parts of West Africa, where monitoring networks tend to be sparse and residence-time estimates can inform long-term water-supply decisions; no data from such a setting were evaluated here, however, and extending this benchmark to one is a natural next step rather than a demonstrated result.

### 4.5 Where this sits within the wider Hydrosheaf programme

Hydrosheaf, more broadly, is a modular framework for inferring groundwater flow and chemistry from sparse tracer data, of which the residence-time and network component evaluated here is one part; companion work addresses reaction chemistry, transport, and uncertainty separately. Nothing here should be read as validating those other components, and nothing in those other components was drawn on to support the results reported here.

### 4.6 Limitations

The most important limitation is that "agreement with the published USGS ages" is not the same as "agreement with the true age of the water," because the published ages are themselves model outputs. The target-withholding test and the exploratory set-valued check were included specifically to reduce dependence on that comparison, and both told the same story as the reference-agreement test, which increases confidence in the overall conclusion without removing the underlying limitation. A second limitation is coverage: every well examined here comes from one national USGS release, and nothing about semi-arid, fractured-rock, or data-sparse African aquifers was tested directly. A third is that the unresolved tracer conflict remains unresolved; the explanations ruled out narrow the problem considerably but do not solve it. A fourth is that some ablation experiments, particularly around helium-4 and old-water tracers, had too few usable wells to support strong conclusions either way.

### 4.7 Recommendations

Future sampling for a network-based approach like this one should prioritise tracer combinations that remain informative even when one tracer is later held back for validation: at least one reliable young-water tracer, one reliable old-water tracer, and the metadata needed to judge whether a proposed network connection is physically sensible in the first place. Future validation of any specific network should include the same negative-control discipline used here, applied to that network's own setting, before the network is trusted. Finally, given the unresolved chlorofluorocarbon conflict identified here, any study leaning heavily on chlorofluorocarbon-derived ages in this kind of release would benefit from an independent redox indicator collected at the same time and place as the tracer sample.

## 5. Conclusions

This study set out to test, rather than assume, whether letting nearby wells inform one another's groundwater-age estimates improves on treating each well independently. On 329 reportable wells from a public USGS release, it did not: no candidate flow-path network met a pre-declared, jointly applied bar for improvement over a single-well baseline, while two deliberately wrong networks were clearly more harmful, and predicting a tracer deliberately withheld from the network told the same story. An independently built, exploratory check, using a different method, reached the same conclusion, and in the process surfaced an unresolved finding of its own: for about a quarter of evaluable tracer combinations, no transit-time distribution can reconcile the measurements at all, for a reason that seven tested explanations failed to identify.

The contribution of this work is therefore a disciplined procedure for testing and, where warranted, rejecting groundwater-age assumptions, not a demonstration that networks of this kind are broadly useful, nor a claim to have identified true groundwater ages, confirmed a flow path, or resolved the tracer conflict it surfaced. Hydrosheaf, evaluated on this narrow but consequential slice of its capability, currently supports a screening and falsification role in multi-tracer groundwater dating. Whether a genuinely informative network exists for a given aquifer remains an empirical question to be asked locally, with the same negative controls and pre-declared decision rule used here, rather than assumed from hydrogeological plausibility alone.

## Author contributions

Dickson Abdul-Wahab: conceptualisation, methodology, software, formal analysis, visualisation, writing - original draft. Dickson Adomako: supervision, validation, hydrogeological interpretation, and writing - review and editing. Gibrilla Abass: water-resources interpretation, validation, and writing - review and editing. Ebenezer Aquisman Asare: nuclear-isotope methodology, supervision, validation, and writing - review and editing. Samuel Ganyaglo: water-resources supervision, interpretation, validation, and writing - review and editing. All authors reviewed and approved the final paper.

## Declaration of competing interest

The authors declare no competing interests related to this paper.

## Data and code availability

The benchmark uses public USGS data releases cited in the reference list. Tritium precipitation isotope input data were obtained from the IAEA/WMO Global Network of Isotopes in Precipitation through WISER [@iaeawmo2026gnip]. The Hydrosheaf implementation used for this study is organised under the public repository referenced in [@abdulwahab2026hydrosheaf]. The groundwater-age benchmark workflow is located under `M3.1/m3_age_benchmark` and includes scripts for design-matrix analysis, gas-correction auditing, tracer-withholding cross-validation, the real-USGS graph benchmark, and publication table and figure generation; figures and the study-area map were produced in R under `M3.1/analysis/r`. Randomised network controls used a fixed seed of 42.

## Tables

**Table 1.** Nuclear and dissolved-gas tracer properties used to benchmark groundwater residence-time inference.

| Tracer | Decay or process scale | Target age range (yr) | Unit | Benchmark role |
| --- | --- | --- | --- | --- |
| Tritium (3H) | 12.32 | 0-60 | TU | bomb-pulse young-water tracer |
| 3H/3He | 12.32 parent | 0-50 | TU-equivalent | closed-system ingrowth apparent age |
| Carbon-14 (14C) | 5730 | 1000-30000 | pmC | old-groundwater radiocarbon constraint |
| Helium-4 (4He) | accumulation | 1000-100000 | ccSTP/g | radiogenic accumulation screening constraint |
| SF6 | stable | 0-50 | pptv | atmospheric-equivalent gas tracer |
| CFC-11/12/113 | stable | 0-60 | pptv | atmospheric-equivalent gas tracers |

**Table 2.** Active design-matrix performance on the national release, regenerated from current HEAD.

| Scenario | Total wells | Reportable N | Median \|log10 error\| | log10 RMSE | log10 R2 | Within factor 2 (%) | Within factor 10 (%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ablation: no 4He | 1272 | 65 | 0.131 | 0.3981 | -2.997 | 69.2 | 100.0 |
| Ablation: raw 14C | 1272 | 61 | 0.1469 | 0.419 | -3.196 | 67.2 | 100.0 |
| Ablation: raw gases | 1272 | 47 |  |  |  |  |  |
| Hydrosheaf model selection | 1272 | 309 | 0.1304 | 0.6098 | 0.764 | 67.3 | 91.6 |
| Old-water 14C ensemble | 1272 | 59 | 0.1469 | 0.4138 | -4.244 | 67.8 | 100.0 |
| Old-water ensemble 4He uncertainty | 1272 | 59 | 0.1469 | 0.4138 | -4.244 | 67.8 | 100.0 |
| Old-water 4He uncertainty | 1272 | 66 | 0.1361 | 0.3955 | -2.618 | 69.7 | 100.0 |
| Reported-model parity | 1272 | 66 | 0.1361 | 0.3955 | -2.618 | 69.7 | 100.0 |
| Screened young-gas correction | 1272 | 66 | 0.1361 | 0.3955 | -2.618 | 69.7 | 100.0 |
| Old tracers only | 1272 | 0 |  |  |  |  |  |
| Young tracers only | 1272 | 67 | 0.1922 | 0.6872 | -0.1 | 53.8 | 92.3 |
| Reported-output-constrained fraction sensitivity | 1272 | 289 | 0.0216 | 0.1964 | 0.97 | 91.7 | 99.7 |
| Strict reported-configuration parity | 1272 | 329 | 0.0279 | 0.2769 | 0.937 | 87.5 | 98.8 |

**Table 3.** Modelling-mode comparison on unequal full-available supports and the N = 40 common support.

| Support | Mode | N | Median \|log10 error\| | log10 RMSE | log10 R2 | Within factor 2 (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Full Available | USGS reported-configuration parity emulation | 329 | 0.0279 | 0.2769 | 0.937 | 87.5 |
| Full Available | Reported-output-constrained sensitivity (not independent validation) | 289 | 0.0216 | 0.1964 | 0.97 | 91.7 |
| Full Available | Hydrosheaf model selection | 309 | 0.1304 | 0.6098 | 0.764 | 67.3 |
| Full Available | Screened young-gas correction | 66 | 0.1361 | 0.3955 | -2.618 | 69.7 |
| Full Available | Reported-model parity | 66 | 0.1361 | 0.3955 | -2.618 | 69.7 |
| High-N Common Support (N=40) | USGS reported-configuration parity emulation | 40 | 0.0752 | 0.304 | -1.016 | 75.0 |
| High-N Common Support (N=40) | Reported-output-constrained sensitivity (not independent validation) | 40 | 0.0653 | 0.2854 | -0.777 | 75.0 |
| High-N Common Support (N=40) | Reported-model parity | 40 | 0.0801 | 0.3079 | -1.069 | 77.5 |
| Design Common Support (N=40) | USGS reported-configuration parity emulation | 40 | 0.0752 | 0.304 | -1.016 | 75.0 |
| Design Common Support (N=40) | Reported-output-constrained sensitivity (not independent validation) | 40 | 0.0653 | 0.2854 | -0.777 | 75.0 |
| Design Common Support (N=40) | Hydrosheaf model selection | 40 | 0.0938 | 0.3607 | -1.838 | 77.5 |
| Design Common Support (N=40) | Screened young-gas correction | 40 | 0.0801 | 0.3079 | -1.069 | 77.5 |
| Design Common Support (N=40) | Reported-model parity | 40 | 0.0801 | 0.3079 | -1.069 | 77.5 |

**Table 4.** Paired comparison of reported-configuration emulation and reported-output-constrained fraction sensitivity.

| Test or metric | Comparison | Statistic | p-value | 95% CI lower | 95% CI upper | Interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| Paired Wilcoxon signed-rank test | Reported-Configuration Emulation vs Reported-Output-Constrained Sensitivity | 5881.0 | < 0.001 |  |  | Significant difference in absolute error distributions (p < 0.05). |
| Paired t-test | Reported-Configuration Emulation vs Reported-Output-Constrained Sensitivity | 2.5389 | 0.012 |  |  | Significant difference in mean absolute errors (p < 0.05). |
| Bootstrap Difference in MAE (strict - agefractions) | Reported-Configuration Emulation vs Reported-Output-Constrained Sensitivity | 0.0047 |  | 0.0012 | 0.0081 | Positive difference indicates age-fraction constraints reduce median absolute error (MAE). |
| Bootstrap Difference in RMSE (strict - agefractions) | Reported-Configuration Emulation vs Reported-Output-Constrained Sensitivity | 0.0336 |  | 0.007 | 0.0777 | Positive difference indicates age-fraction constraints reduce RMSE. |

**Table 5.** Graph-age benchmark on the 329 strict reportable nodes.

| Graph family | Prior strength | Nodes | Edges | log10 RMSE (single) | log10 RMSE (graph) | Delta log10 RMSE | Within factor 2, single (%) | Within factor 2, graph (%) | Violations before | Violations after | Meets robust-improvement rule |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Coordinate (global) | None | 329 | 329 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 127 | 127 | No |
| Coordinate (global) | Weak | 329 | 329 | 0.2769 | 0.8532 | 0.57632 | 87.5 | 69.6 | 127 | 110 | No |
| Coordinate (global) | Medium | 329 | 329 | 0.2769 | 0.9237 | 0.64677 | 87.5 | 66.9 | 127 | 109 | No |
| Coordinate (global) | Strong | 329 | 329 | 0.2769 | 0.9395 | 0.66266 | 87.5 | 66.9 | 127 | 109 | No |
| Coordinate (study unit) | None | 329 | 319 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 123 | 123 | No |
| Coordinate (study unit) | Weak | 329 | 319 | 0.2769 | 0.7516 | 0.47472 | 87.5 | 71.4 | 123 | 103 | No |
| Coordinate (study unit) | Medium | 329 | 319 | 0.2769 | 0.8079 | 0.53102 | 87.5 | 68.1 | 123 | 102 | No |
| Coordinate (study unit) | Strong | 329 | 319 | 0.2769 | 0.8178 | 0.54094 | 87.5 | 68.1 | 123 | 102 | No |
| Depth constrained | None | 329 | 218 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 77 | 77 | No |
| Depth constrained | Weak | 329 | 218 | 0.2769 | 0.5531 | 0.27618 | 87.5 | 77.2 | 77 | 66 | No |
| Depth constrained | Medium | 329 | 218 | 0.2769 | 0.5836 | 0.3067 | 87.5 | 75.4 | 77 | 65 | No |
| Depth constrained | Strong | 329 | 218 | 0.2769 | 0.587 | 0.31008 | 87.5 | 75.1 | 77 | 65 | No |
| Hydraulic proxy | None | 329 | 30 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 12 | 12 | No |
| Hydraulic proxy | Weak | 329 | 30 | 0.2769 | 0.2755 | -0.00135 | 87.5 | 86.9 | 12 | 11 | Yes |
| Hydraulic proxy | Medium | 329 | 30 | 0.2769 | 0.2766 | -0.0003 | 87.5 | 86.9 | 12 | 11 | Yes |
| Hydraulic proxy | Strong | 329 | 30 | 0.2769 | 0.2766 | -0.00028 | 87.5 | 86.9 | 12 | 11 | Yes |
| Parameter smoothing | None | 329 | 30 | 0.2769 | 0.2783 | 0.00142 | 87.5 | 87.2 | 12 | 12 | No |
| Parameter smoothing | Weak | 329 | 30 | 0.2769 | 0.2815 | 0.00462 | 87.5 | 86.9 | 12 | 8 | No |
| Parameter smoothing | Medium | 329 | 30 | 0.2769 | 0.2823 | 0.00545 | 87.5 | 86.9 | 12 | 9 | No |
| Parameter smoothing | Strong | 329 | 30 | 0.2769 | 0.2851 | 0.00825 | 87.5 | 86.6 | 12 | 8 | No |
| Wrong-direction control | None | 329 | 30 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 18 | 18 | No |
| Wrong-direction control | Weak | 329 | 30 | 0.2769 | 0.3788 | 0.10193 | 87.5 | 84.2 | 18 | 18 | No |
| Wrong-direction control | Medium | 329 | 30 | 0.2769 | 0.3858 | 0.10887 | 87.5 | 84.2 | 18 | 18 | No |
| Wrong-direction control | Strong | 329 | 30 | 0.2769 | 0.3859 | 0.10903 | 87.5 | 84.2 | 18 | 18 | No |
| Randomised control | None | 329 | 319 | 0.2769 | 0.2769 | 0.0 | 87.5 | 87.5 | 157 | 157 | No |
| Randomised control | Weak | 329 | 319 | 0.2769 | 0.9934 | 0.71652 | 87.5 | 55.9 | 157 | 160 | No |
| Randomised control | Medium | 329 | 319 | 0.2769 | 1.0948 | 0.81791 | 87.5 | 52.0 | 157 | 158 | No |
| Randomised control | Strong | 329 | 319 | 0.2769 | 1.1142 | 0.83729 | 87.5 | 51.1 | 157 | 159 | No |

**Table 6.** Paired ablation and sensitivity effects on reportable common rows.

| Scenario | Paired N | Median delta \|log10 error\| | Mean delta \|log10 error\| | Improved fraction (%) | Gained factor-2 | Lost factor-2 |
| --- | --- | --- | --- | --- | --- | --- |
| Ablation: no 4He | 65 | 0.0 | 0.0 | 0.0 | 0 | 0 |
| Ablation: raw 14C | 61 | 0.0 | 0.0054 | 6.6 | 0 | 0 |
| Hydrosheaf model selection | 64 | -0.0017 | -0.0446 | 53.1 | 11 | 5 |
| Old-water 14C ensemble | 59 | 0.0 | -0.0003 | 16.9 | 0 | 0 |
| Old-water ensemble 4He uncertainty | 59 | 0.0 | -0.0003 | 16.9 | 0 | 0 |
| Old-water 4He uncertainty | 66 | 0.0 | -0.0 | 1.5 | 0 | 0 |
| Screened young-gas correction | 66 | 0.0 | 0.0 | 0.0 | 0 | 0 |
| Young tracers only | 38 | 0.0 | -0.0054 | 44.7 | 0 | 0 |
| Reported-output-constrained fraction sensitivity | 40 | 0.0 | -0.012 | 35.0 | 0 | 1 |
| Strict reported-configuration parity | 44 | 0.0 | 0.0023 | 18.2 | 0 | 1 |

## Figures

![Figure 1. Geographic distribution of the 1,272 USGS national public-supply sites used in the M3.1 benchmark, coloured by aquifer lithology group, with the N = 329 strict reportable subset marked.](manuscript/artifacts/figures/FIG-1_site_map.png)

**Figure 1.** Geographic distribution of the 1,272 USGS national public-supply sites used in the M3.1 benchmark, coloured by aquifer lithology group, with the N = 329 strict reportable subset marked.

![Figure 2. Atmospheric input histories of tritium, sulfur hexafluoride and chlorofluorocarbon-12 used for young-groundwater residence-time interpretation.](manuscript/artifacts/figures/FIG-2_atmospheric_histories.png)

**Figure 2.** Atmospheric input histories of tritium, sulfur hexafluoride and chlorofluorocarbon-12 used for young-groundwater residence-time interpretation.

![Figure 3. Active design-matrix reference agreement and reportability across scenarios, regenerated from current HEAD.](manuscript/artifacts/figures/FIG-3_design_matrix_performance.png)

**Figure 3.** Active design-matrix reference agreement and reportability across scenarios, regenerated from current HEAD.

![Figure 4. Strict reported-configuration emulation against model-derived USGS reference ages, N = 329 reportable subset.](manuscript/artifacts/figures/FIG-4_strict_parity_scatter.png)

**Figure 4.** Strict reported-configuration emulation against model-derived USGS reference ages, N = 329 reportable subset.

![Figure 5. Graph-age effects on the N = 329 strict reportable nodes: change in RMSE against change in within-factor-2 agreement, one point per candidate family or negative control at each non-zero prior strength.](manuscript/artifacts/figures/FIG-5_graph_benchmark.png)

**Figure 5.** Graph-age effects on the N = 329 strict reportable nodes: change in RMSE against change in within-factor-2 agreement, one point per candidate family or negative control at each non-zero prior strength.

![Figure 6. Leakage-guarded target-tracer withholding for 3H, SF6 and 14C: RMSE for the single-node baseline against three candidate/negative-control graphs.](manuscript/artifacts/figures/FIG-6_tracer_withholding.png)

**Figure 6.** Leakage-guarded target-tracer withholding for 3H, SF6 and 14C: RMSE for the single-node baseline against three candidate/negative-control graphs.

![Figure 7. Exploratory set-valued compatibility audit (development-stage): pairwise tracer infeasibility rates and the rise of local infeasibility with conditioning-set size.](manuscript/artifacts/figures/FIG-7_infeasibility_audit.png)

**Figure 7.** Exploratory set-valued compatibility audit (development-stage): pairwise tracer infeasibility rates and the rise of local infeasibility with conditioning-set size.
