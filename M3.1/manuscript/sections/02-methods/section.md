## 2. Materials and methods

### 2.1 Design philosophy

The benchmark was designed to test, rather than assume, whether letting nearby wells inform one another's age estimates helps. A candidate flow-path network was therefore always evaluated against a single-well baseline and against two negative controls built to share the same size and shape as the real network but none of its hydrogeological logic: a network with every connection pointing the wrong way, and a network connecting wells at random. A network result was accepted as a genuine improvement only when it beat the single-well baseline on both a statistical-agreement measure and a practical accuracy measure at once, and only when the negative controls did not show a comparable effect. Published USGS ages were used throughout as a comparison target for evaluating agreement, not as directly observed true ages, because they are themselves the product of a modelling workflow [@jurgens2012tracerlpm].

### 2.2 Data

The primary dataset was the USGS national public-supply aquifer release for the continental United States, 2004-2017 [@jurgens2022data], covering 1,279 wells across 21 principal aquifers, with tracer measurements, dissolved-gas model outputs, published age estimates, and site metadata. After harmonising sample identifiers, dates, coordinates, and tracer fields across the release's component tables, 1,272 rows remained ([[FIGREF:FIG-1]]; Supplementary Table S1). Two further public releases were examined but not pooled into the main comparison: a Western Principal Aquifers release [@faulkner2019data] that lacks the coordinates needed for the network benchmark, and a Mississippi River Valley alluvial-aquifer release [@gratzer2025groundwater] that lacks the published age table needed for a like-for-like comparison.

The tracer suite spanned modern to very old water: tritium (^3^H) and tritiogenic helium-3 for water recharged since roughly 1950, sulfur hexafluoride and chlorofluorocarbons for the same young-water window, carbon-14 for water up to several tens of thousands of years old, and helium-4 as a supporting old-water indicator ([[TABREF:TAB-1]]). The tritium input history used to interpret ^3^H was built from precipitation isotope records ([[FIGREF:FIG-2]]; Supplementary Methods S1).

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

The full comparison matrix ([[TABREF:TAB-2]]) included the reported-configuration and reported-age-fraction emulation scenarios, an automatic model-selection scenario, several tracer and correction ablations, and the network variants above, applied to all 1,272 wells with the reportability guard applied throughout. Direct comparisons between scenarios used only the rows both scenarios could evaluate in common; results on each scenario's full available support are reported separately and are not directly rankable against each other because their supports differ ([[TABREF:TAB-5]]).

Accuracy against the published ages was assessed on a base-ten logarithmic scale, because groundwater ages span several orders of magnitude, using the median absolute log-error, the root-mean-square log-error, and the share of estimates within a factor of two or ten of the published value (Supplementary Methods S5). A network result required a lower root-mean-square error and unchanged-or-better factor-of-two agreement simultaneously to be accepted as an improvement.

A second, more demanding test avoided the published ages altogether. For every eligible well, one tracer at a time was removed before model fitting and before network construction, so that neither the single-well fit nor any connected neighbour could see the withheld measurement, and the withheld tracer's own concentration was then predicted and compared with its true measured value. This target-withholding test cannot be satisfied by a network simply reproducing information it was quietly given.

### 2.6 An independent, exploratory cross-check

A separate, more recently developed diagnostic approached the same network question from a different angle: rather than fitting one best age distribution per well, it asked which age distributions, if any, are mathematically consistent with a well's tracers once measurement uncertainty is taken into account, and tested whether a candidate flow-path edge was compatible with the possibilities at both of its wells. This diagnostic is explicitly at a development and exploratory stage; it was not used to revise any scalar result above, and it is reported separately with that status stated throughout (Supplementary Methods S6-S7).

### 2.7 Reproducibility

The full pipeline was rerun end to end from the current codebase for this revision, using a fixed 90-step age grid, a fixed random seed for the randomised control, and machine-readable manifests recording software versions and file checksums for every input and output. Every regenerated result was compared line by line against the previously locked benchmark; the outcome of that comparison, and the small number of cases showing ordinary floating-point-level variation rather than a substantive change, are reported in Supplementary Methods S8.
