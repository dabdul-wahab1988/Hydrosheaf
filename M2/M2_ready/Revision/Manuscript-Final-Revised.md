# Hydrosheaf: A Reproducible Computational Framework for Inferring Directed Groundwater Hydrochemical Evolution Networks from Sparse Multitracer Data

**Manuscript CAGEO-D-26-00847 (revised), Computers and Geosciences (Elsevier)**

Dickson Abdul-Wahab^a\*^, Ebenezer Aquisman Asare^b^, Dickson Adomako^a^, Gibrilla Abass^c^, Samuel Ganyaglo^a,d^

^a^ Department of Nuclear Science and Applications, School of Nuclear and Allied Sciences, University of Ghana, Atomic-Kwabenya, Accra, Ghana

^b^ Nuclear Chemistry and Environmental Research Centre, National Nuclear Research Institute (NNRI), Ghana Atomic Energy Commission (GAEC), Box LG 80, Legon-Accra, Ghana

^c^ Department of Water Resources and Aquaculture Management, School of Sustainable Development, University of Environment and Sustainable Development, Somanya, Ghana

^d^ Water Resources Research Centre, National Nuclear Research Institute, GAEC, Box LG 80, Kwabenya-Accra, Ghana

\* Corresponding author: Dickson Abdul-Wahab (dabdul_wahab11@yahoo.com; ORCID: https://orcid.org/0000-0001-7446-5909)


## Abstract

Sparse groundwater datasets combine incomplete major-ion chemistry, irregular spatial coverage, optional isotope measurements, and limited hydraulic-head information, hindering reproducible reconstruction of groundwater evolution pathways. Hydrosheaf is an open Python framework for inferring directed hydrochemical evolution networks from such data. Sampling locations are graph nodes; candidate pathways are directed edges constrained by hydraulic, spatial, conservative-tracer, isotope, residence-time, and thermodynamic evidence. Each retained edge is decomposed into non-reactive transport and residual chemical change before sparse inverse fitting from a configurable reaction dictionary, saturation-index gating, and process-stability ranking (PSI). Validation is a multi-tier benchmark: synthetic reaction and transport recovery, MODPATH-referenced topology recovery, age benchmarks, PHREEQC-compatible forward checks, and a data-limited field demonstration in Northern Ghana. Independent no-prior topology inference achieves F1 = 0.62 (precision 0.49, recall 0.84) against 174 MODPATH reference edges but overconnects, so edges are screening-level hypotheses; the prior-assisted F1 = 1.00 is an ingestion check. Reaction-extent recovery is moderate and limited by reaction-dictionary degeneracy (active-reaction R2 = 0.23, MAE = 0.37 mmol/L, 54% false activations). Synthetic ages recover with log10 corr-R2 = 0.98 and median absolute error 17.9 y, while public USGS parity is identifiability-gated: 356 of 1272 fits identifiable, median |log10 error| = 0.02, 92% within a factor of two. PHREEQC-compatible forward checks (median RMSE = 0.02 mmol/L, NSE = 1.00, 88% feasible) are a proxy, not live kinetic simulation. On 258 retained field edges (from 572 candidates), median chemistry R2 = 0.70 (Lower Anayari 0.53, Talensi 0.82) with median PSI = 0.97 and no independent process-truth labels; the field cohomology diagnostic reports homogeneous nullities H0 = 0 and H0 = 10, with positive affine obstruction energies of 0.295 and 0.128 at Lower Anayari and Talensi, respectively. Thus neither field network has an exact affine global section, and the nullities should not be read as counts of exact section assignments. Hydrosheaf thus provides a reusable, auditable workflow yielding screening-level, uncertainty-ranked hypotheses for groundwater process discovery.

## Keywords

hydrogeochemistry; directed graph inference; sheaf-inspired consistency scoring; inverse modelling; process stability index; residence time; PHREEQC; uncertainty quantification; reproducible geoscience software

## Highlights

- Hydrosheaf infers directed groundwater hydrochemical evolution networks
- Sheaf-style scoring combines hydraulic, chemical, isotope, and age evidence
- Sparse inverse fitting with thermodynamic gates screens reaction hypotheses
- Process-stability index ranks process hypotheses under input uncertainty
- Multi-tier validation: synthetic, MODPATH, age, PHREEQC proxy, and field



## 1. Introduction

Groundwater datasets are characteristically sparse, multivariate, spatially irregular, and temporally incomplete. A typical field investigation yields hydrochemical concentrations at tens to hundreds of wells, isotope and tracer measurements at a subset of those wells, hydraulic-head records at further reduced coverage, and spatial metadata of variable completeness. The scientific challenge is not merely to interpret these data in isolation but to reconstruct computationally plausible, directed groundwater evolution pathways from incomplete, heterogeneous evidence. This is a geoscientific inference problem that requires the simultaneous integration of spatial topology, conservative transport, reactive geochemistry, age information, and thermodynamic feasibility, a combination that exceeds the scope of any single existing analytical tool. Groundwater quality is not static but evolves continuously as water moves through the subsurface, accumulating dissolved constituents through mineral weathering, ion exchange, redox transformation, and mixing with waters of different ages and origins (Ma et al., 2019; Yang et al., 2018). This evolution is governed by three interacting dimensions: the residence time of water in the aquifer, which determines the duration of water–rock contact; the connectivity structure of the flow system, which determines which minerals and rock surfaces water encounters; and the geochemical reactions that transform initial recharge chemistry into the observed groundwater composition (Casillas-Trasvina et al., 2022; Villegas et al., 2018). The problem is particularly acute in data-limited aquifer settings, where semi-arid climates, sparse monitoring infrastructure, and complex subsurface heterogeneity combine to produce maximum interpretive ambiguity from minimum field data (Borzí, 2025). In such settings, the analyst is caught between over-interpreting sparse chemistry and abandoning quantitative pathway reconstruction altogether; neither outcome serves groundwater management or monitoring design. Reproducible, open computational workflows are therefore not merely convenient but scientifically necessary, because site-specific manual interpretation is neither transferable nor independently verifiable.

Conventional groundwater interpretation relies on complementary but operationally disconnected tools. Piper diagrams, Gibbs plots, ion-ratio analyses, and hydrochemical facies classifications identify water types and qualitatively characterise geochemical evolution, but they operate on individual samples and do not compute directed pathways, reaction extents, or network-level connectivity (Roy et al., 2020). Conventional interpretation treats each sampled location as an independent observation whose chemistry is attributed to processes at that point or along a flow path selected by the analyst on the basis of geological intuition, spatial proximity, or hydrochemical similarity (Teeple et al., 2021; Wang et al., 2025). This point-based paradigm has two well-recognised limitations. First, a single groundwater sample typically represents a mixture of waters that have followed different flow paths and carry different ages, such that neither the apparent age nor the chemistry can be unambiguously attributed to a single evolutionary trajectory (Casillas-Trasvina et al., 2022; Newman et al., 2021); multi-component mixing, between modern shallow recharge and older confined water, for example, can produce intermediate chemistry, yielding an apparent lumped-parameter-model (LPM) age that reflects neither end-member correctly, as documented in the Belgian Neogene aquifer (Casillas-Trasvina et al., 2022), fault-connected Algerian systems (Ghiglieri et al., 2021), and aquifer complexes with inter-unit leakage (Murgulet et al., 2016). Second, point-based hydrochemical classifications such as Piper diagrams, hydrochemical facies, and bivariate ion ratios identify water types but do not identify the spatial or temporal sequence that generated them (Díaz-Puga et al., 2016; Song et al., 2024); two samples of the same facies can arise from entirely different evolutionary pathways in contrasting lithological and climatic settings (Wang et al., 2025). Inverse hydrogeochemical codes, principally PHREEQC and NETPATH, provide rigorous mass-balance modelling but depend entirely on the analyst's prior selection of an initial and a final water pair assumed to represent a connected up-gradient/down-gradient pathway (Parkhurst & Appelo, 2013). When flow paths are uncertain, as they are in virtually every data-limited aquifer, different assumed pairs yield different but equally valid reaction solutions, so non-uniqueness is structural rather than merely statistical (Manu et al., 2023; Slimani et al., 2017). The same Ca-HCO₃ water type, for instance, may be produced by carbonate dissolution, silicate weathering with CO₂ uptake, mixing with Ca-rich end-members, ion exchange, or combinations of these processes (Roy et al., 2020; Li et al., 2019; Zhi et al., 2021). Dissolution, precipitation, and exchange reactions can offset each other, for example, gypsum dissolution adding Ca²⁺ and SO₄²⁻ while calcite precipitation removes Ca²⁺, generating net budgets indistinguishable from alternative reaction combinations (Slimani et al., 2017; Roy et al., 2020). In the Qaidam Basin, Guo et al. (2025) demonstrated that halite and gypsum dissolution, silicate weathering, carbonate precipitation, and cation exchange were all necessary to explain observed hydrochemical evolution, but that numerous phase-extent combinations reproduced the same final water chemistry within analytical precision. In the Hancheng Mining Area, Gibbs and Piper plots indicated rock-weathering control, yet the interpretation still required separate attribution to dolomite, calcite, gypsum, halite, and pyrite oxidation, highlighting the inadequacy of graphical methods alone (Kou et al., 2024). Phase selection strongly controls the solution: Manu et al. (2023) explicitly used combinatorial inverse modelling in the Ghanaian Pra Basin because different mineral assemblages produced chemically valid but geologically contrasting pathway interpretations, and multiple tracers in the Pahute Mesa groundwater study were required to reveal mixing zones and long-term regional flow that no single method could identify (Kwicklis et al., 2021). Combinatorial inverse modelling (Manu et al., 2023) systematises this procedure by enumerating many candidate sample pairs, yet it still starts from analyst-chosen pairs, does not score the resulting edges jointly against hydraulic, tracer, isotope, and age evidence, does not separate transport effects from reactions, and provides no uncertainty ranking of the competing solutions. The practical consequence is that combinatorial search accelerates the arithmetic but not the interpretation: the analyst still chooses the search space, and every enumerated pair is treated as equally admissible regardless of hydraulic, tracer, or age support. In contrast, Hydrosheaf runs inverse models only on topology-supported edges, so pathway selection is tested against independent evidence before any reaction fit is attempted. The problem is compounded by uncertainty in the flow-path assumption itself. If the analyst incorrectly pairs wells that are not hydraulically connected, the inferred reaction extents are artefacts rather than records of genuine geochemical evolution (Manu et al., 2023; Teeple et al., 2021). Seltzer et al. (2021) further showed that anthropogenic carbonate inputs can distort ¹⁴C-based age estimates by thousands of years, illustrating how geochemical-model uncertainty propagates directly into age estimates and thence into assessments of groundwater renewability. Separately, lumped-parameter residence-time models use tracers such as tritium and radiocarbon to estimate groundwater age at individual nodes, but these single-node age estimates do not constrain reaction magnitudes, do not test topological consistency, and cannot independently validate or reject proposed flow connections (Thiros et al., 2023; Rädle et al., 2022). Single-tracer inversion conflates mean residence time with the shape of the transit-time distribution and tracer-input-function uncertainty (Ma et al., 2019; Rädle et al., 2022), with tritium-based ages relying on the ³H half-life of 12.32 yr (Lucas & Unterweger, 2000), and when multiple tracers disagree, as is common in mixed-age systems, the analyst must choose between competing age models without a formal arbitration framework (Casillas-Trasvina et al., 2022). In heterogeneous hard-rock aquifers, nearby wells can sample entirely different hydraulic domains whose age and chemistry patterns are only interpretable when aquifer architecture is considered (Comte et al., 2018; Cao et al., 2020). Point samples can also miss vertical connectivity and leakage: isotopic differences between shallow and deep alluvial groundwater in Beijing were invisible to point-chemistry classification but clearly interpretable in a depth-connected network context (He et al., 2021; Joshi et al., 2018). River–aquifer exchange further illustrates that isolated well sampling can miss discharge-integrated flow-path signals apparent only when surface-water and groundwater nodes are treated as parts of a connected network (Newman et al., 2021; Smerdon & Gardner, 2022; Binet et al., 2017; Sun et al., 2025), and dynamic groundwater flow tracking has quantified time lags of nitrate, chloride, and tritium in lowland stream networks (Kaandorp et al., 2021). Eastern Dahomey Basin analyses further demonstrated that aquifer heterogeneity and anthropogenic influences obscure the simple spatial trends that point-based interpretation assumes (Aladejana et al., 2020). Applied sequentially rather than jointly, these three strands, hydrochemical classification, inverse reaction modelling, and residence-time estimation, produce locally plausible but globally unverified interpretations that cannot be cross-examined against one another (Casillas-Trasvina et al., 2022). Multi-tracer field studies across contrasting settings such as arid-basin alluvial aquifers in northwestern China (Xiao et al., 2018; Yang et al., 2018), confined coastal aquifers in Colombia (Villegas et al., 2018), and carbonate karst systems in Morocco (Roubil et al., 2022) and the Xishan karst aquifer near Beijing (Qin et al., 2017) consistently demonstrate that groundwater evolves from fresh, recently recharged waters toward more mineralised, mixed, or palaeorecharge compositions through carbonate and silicate dissolution, evaporite dissolution, cation exchange, redox reactions, evaporation, and anthropogenic inputs (Caschetto et al., 2016). In arid and semi-arid aquifers, groundwater age increases with distance from recharge zones and depth while chemistry shifts toward higher salinity, as documented in loess-basin systems (Ling et al., 2022) and basin-scale confined aquifers (Purtschert et al., 2022). In mining-influenced and urban environments, age tracers combined with redox chemistry have revealed that modern water fractions transport arsenic, nitrate, and sulphate into shallow and deeper aquifer zones via coupled flow and reaction mechanisms (Richards et al., 2022; Cheng et al., 2022; Lapworth et al., 2017). Seasonal pumping-induced recharge further complicates the coupling between surface-water chemistry and groundwater composition (Guo et al., 2018), while semi-arid zoning models confirm that hydrochemical processes change along flow systems from leaching through cation exchange, evaporation, and concentration (Li et al., 2019). Stable-isotope interpretation of evaporative enrichment rests on the equilibrium fractionation of ¹⁸O and ²H between water and vapour (Majoube, 1971), and reviews of radioactive-isotope and other residence-time tracers have mapped the possibilities, challenges, and limitations of recharge estimation (Cartwright et al., 2017).

The absence this paper addresses is not another groundwater case study but a reusable computational framework capable of jointly handling directed graph construction, transport adjustment, sparse inverse reaction fitting, residence-time constraints, thermodynamic gating, and uncertainty diagnostics. Existing approaches cover individual tasks: PHREEQC solves reactions given assumed paths; MODFLOW/MODPATH simulates physical flow given calibrated parameters; graph-based network methods infer connectivity without reaction attribution; and machine-learning approaches predict hydrochemical patterns without attributing chemical change to specific processes. None of these strands integrates into a common inference structure in which the evidence streams can veto one another. Representing groundwater sampling locations as nodes in a directed graph and inferred flow connections as directed edges provides a natural and flexible structure for this integration (Schiavo et al., 2022; Taccari et al., 2024; Borzí, 2025). Empirical evidence for the network character of groundwater systems is strong: in the Lez karst system, δ¹⁸O, δ²H, and ⁸⁷Sr/⁸⁶Sr indicated interaction between deep storage and the main Jurassic aquifer through fracture-controlled connections (Bicalho et al., 2017); in the Galilee Basin, hydraulic-head patterns and isotopic evidence jointly supported connected spring discharge from shallow and deep regional flow paths (Keegan-Treloar et al., 2024); in fractured gneiss in Malawi, isotopes and major ions identified the Ntcheu Fault as a conduit linking recharge and discharge zones along a directed flow network (Kambuku et al., 2018); and in the Yuquan Mountain karst system, multiple evidence streams identified distinct recharge zones and fault-controlled pathways undetectable from point chemistry (Sun et al., 2025). Within a directed graph, each edge is a testable hypothesis scored, tested, and either retained or discarded on the basis of physical and geochemical evidence (Borzí, 2025). Graph constraints improve reaction interpretation in two complementary ways: they restrict the set of node pairs submitted to the inverse solver to those for which connectivity is physically admissible, eliminating reaction solutions that contradict elevation, head, or tracer evidence; and they enable network-level consistency checks (that is, age must increase monotonically along directed edges, chemistry must evolve in the direction of the inferred reaction) impossible to apply when samples are interpreted individually (Roubil et al., 2022; Schiavo et al., 2022). Sheaf-theoretic consistency motivates this: the requirement that local observations on edges must agree globally across the network parallels the requirement in cellular sheaf theory that local sections must be consistent across the covering structure (Robinson, 2020; Hansen & Ghrist, 2019), and while formal sheaf implementations have not previously been demonstrated for groundwater networks (Borzí, 2025), the underlying principle is directly implementable as a scoring mechanism for directed-edge plausibility. We use that principle to motivate a diagnostic affine coboundary layer while keeping edge selection as an explicit weighted multi-criteria score. The implementation does not form or spectrally decompose a graph or sheaf Laplacian; it reports homogeneous nullity, affine obstruction energy, and leave-one-edge-out leverage on the retained network (Section 2.4; Supplementary Method S1). Graph-based methods have also been applied directly to groundwater systems: physical conduit-network graphs represent hydraulic connections as directed edges between conduit intersections and spring outlets in karst systems (Fandel et al., 2022; Malenica et al., 2018); probabilistic preferential-flow networks use aquifer-bottom topographic gradients and Monte Carlo sampling to generate directed edge-confidence maps for porous aquifers (Schiavo et al., 2022); and graph-based flow approximations in discrete fracture networks preserve breakthrough-curve structure at reduced computational cost (Karra et al., 2018; Viswanathan et al., 2018). Graph neural networks and data-driven connectivity models learn statistical associations in groundwater level time series but do not produce geochemically interpretable mass transfers or age-constrained flow paths (Bai & Tahmasebi, 2023; Taccari et al., 2024), and deep-learning surrogates have been proposed for groundwater flow equations (Zhang et al., 2022). Stochastic particle-tracking methods have been used to characterise flow-path and capture-zone uncertainty in urban groundwater systems (Alberti et al., 2018; Colombo et al., 2021), and ensemble particle-tracking studies have shown that connections below minimum probability thresholds produce internally inconsistent age-order fields (Gonzalez et al., 2020; Juckem & Starn, 2021). Where calibrated flow models exist, published MODPATH particle-tracking archives provide physically based flow-path references (Harte, 2021), and particle-tracking travel times are sensitive to the underlying flow model and porosity assumptions (Meyer et al., 2018; Baker et al., 2025). Public USGS tracer-age compilations provide independent residence-time benchmarks for such comparisons (Jurgens et al., 2022). Python scripting could, in principle, combine PHREEQC wrappers such as phreeqpy or iPhreeqc with graph-analysis libraries such as NetworkX to assemble parts of this workflow. To our knowledge, however, no published, versioned workflow integrates transport decomposition, sparse inverse fitting, thermodynamic gating, uncertainty diagnostics, and a multi-tier validation suite under one configurable pipeline; ad hoc combinations lack the consistency checks, provenance logging, and benchmark infrastructure needed to make the results auditable and transferable. Assembling such a pipeline from scratch also leaves the decisive choices implicit: which pairs to test, how to score them against independent evidence, and how to propagate input uncertainty into the final hypothesis ranking. Reproducibility is an additional requirement: methods applied in isolation are difficult to compare, audit, or transfer because they depend on undocumented analyst choices for flow-path pairs, mineral-phase lists, and age-correction assumptions. A reproducible, open workflow transforms these subjective choices into configurable, logged, and testable algorithmic decisions.

Hydrosheaf is a modular Python framework designed to infer directed groundwater hydrochemical evolution networks from sparse multivariate field data. Sampling sites are represented as nodes in a directed graph, and inferred evolution pathways as directed edges, each treated as a testable hypothesis that water can evolve from an up-gradient node to a down-gradient node under physically plausible flow conditions. Along each retained edge, a two-stage inverse solver first removes non-reactive transport effects attributable to evaporation or dilution mixing and then fits the minimum set of geochemical reactions capable of explaining the residual hydrochemical change. The pipeline is organised modularly, so each component can be tested in isolation and exchanged without disrupting the remainder of the workflow. Inputs comprise the eight principal major ions, spatial coordinates, and optionally hydraulic head, stable isotopes, and nuclear tracers. Outputs include edge-level reaction extents, transport parameters, residence-time estimates, thermodynamic feasibility flags, process-stability indices (PSI), and structured provenance records. The framework accepts the minimum data available in data-limited settings and activates higher-level modules whenever optional tracers or head data are present, making it transferable beyond its Northern Ghana demonstration sites. The Lower Anayari demonstration catchment, whose hydrogeochemistry and recharge processes have been characterised previously (Abdul-Wahab et al., 2021), exemplifies this transferability.

Hydrosheaf makes five algorithmic contributions. First, groundwater hydrochemical evolution is formulated as a directed graph inference problem rather than a manual flow-path selection task. Second, sheaf-inspired consistency scoring evaluates whether proposed edges are simultaneously supported by hydraulic gradients, conservative-tracer ordering, stable-isotope consistency, and age-order coherence, providing a joint-evidence edge filter that also localises subregions of inconsistent inference. Third, sparse inverse reaction fitting with L1 regularisation selects the most parsimonious reaction set from a pre-compiled stoichiometric dictionary, and PHREEQC-derived saturation-index gates enforce thermodynamic feasibility on the inferred reaction directions. Fourth, the process-stability index (PSI) quantifies the probability that each inferred reaction remains active under Monte Carlo input perturbation, converting best-fit outputs into robustness-ranked process hypotheses. Fifth, the transport-reaction decomposition and its uncertainty diagnostics, Monte Carlo propagation, bootstrap resampling, and variance decomposition, make the separation of physical from chemical change explicit and auditable. In addition, the framework is validated through a multi-tier benchmark strategy covering synthetic reaction and transport recovery, MODPATH-referenced topology recovery, public and synthetic age benchmarks, PHREEQC-compatible forward checks, and a data-limited field demonstration; we present this benchmark suite as a validation strategy rather than as an algorithmic contribution in itself.

The aim of this paper is to present, validate, and demonstrate Hydrosheaf as a reusable computational geoscience tool. Section 2 describes the computational framework and algorithm design; Section 3 describes the multi-tier validation strategy; Section 4 reports benchmark and field results; and Section 5 discusses the contribution, its limitations, and its practical scope. The Northern Ghana datasets serve as field demonstrations of transferability, not as the primary scientific target of the paper.





## 2. Computational Framework and Algorithm Design

### 2.1 Framework Overview and Design Philosophy

Hydrosheaf is designed as a modular, end-to-end computational workflow that converts sparse, multivariate groundwater observations into directed hydrochemical evolution networks with quantified uncertainty. The framework addresses the fundamental inference problem that groundwater sampling sites, hydrochemical observations, isotope measurements, and spatial metadata rarely arrive in a form that directly supports flow-path reconstruction or reaction attribution; the computational challenge is to impose physically defensible structure on this heterogeneous evidence simultaneously rather than sequentially. The workflow is organised into eight processing layers and approximately 41 functional modules (Figure 1). Data enter through command-line, Python API, or automated workflow entry points and pass through a preprocessing stage before reaching the central inference engine. The inference engine integrates seven topology and prior-inference modules with five geochemical constraint modules, whose outputs are consumed by a six-stage inverse solver. The final output layer generates edge result tables, network summaries, diagnostic plots, and structured interpretation reports. A module-disabling mechanism ensures that any module dependent on absent optional data is automatically bypassed, preserving pipeline continuity under data-limited conditions. The modular architecture allows each component to be tested in isolation and exchanged without disrupting the remainder of the workflow, which is a prerequisite for reproducible scientific software (Table 1).

### 2.2 Data Model and Preprocessing

The framework operates on a tabular sample schema in which each row represents one groundwater sample and columns correspond to field identifiers, spatial metadata, major-ion concentrations, optional field parameters, optional stable isotopes, and optional nuclear tracers (Table S1). Required fields comprise spatial coordinates (latitude/longitude or projected x/y), sample and site identifiers, and concentrations of the eight principal major ions — Ca²⁺, Mg²⁺, Na⁺, K⁺, HCO₃⁻, Cl⁻, SO₄²⁻ and NO₃⁻ — expressed in mmol/L or mg/L. Hydraulic head or ground-surface elevation, pH, electrical conductivity (EC) and temperature are recommended. Stable isotopes (δ¹⁸O, δ²H) and nuclear tracers (³H, ¹⁴C) are optional but activate higher-order modules when present. Concentrations are converted to a common pipeline unit basis of mmol/L using stored conversion factors; all downstream modules operate on this basis.

Charge balance is evaluated for each sample using the charge balance error (CBE). Because charge balance requires comparison of equivalent charges, the CBE is computed on equivalent concentrations in meq/L, obtained from the mmol/L basis as meq/L = |z_i| × mmol/L:


$$
\mathrm{CBE} = \frac{\sum_{i} z_i^{+}\, C_i^{+} \;-\; \sum_{j} |z_j^{-}|\, C_j^{-}}{\sum_{i} z_i^{+}\, C_i^{+} \;+\; \sum_{j} |z_j^{-}|\, C_j^{-}}
\tag{1}
$$

Equation 1. Charge balance error, where *C*~i~ is the equivalent concentration of ion *i* in meq/L (meq/L = |z~i~| × mmol/L), z~i~ is the ionic valence, and the superscripts + and − denote cations and anions respectively. Samples with |CBE| > 0.10 are flagged but not discarded. Detection-limit substitution is applied to censored values, and negative concentrations are replaced with a small positive floor and flagged accordingly. Missing optional fields are recorded in a data availability matrix that governs downstream module activation. All preprocessing operations are logged to a provenance record retained in the output archive.

### 2.3 Directed Graph Construction

Each groundwater sampling site is represented as a node v~i~ in the node set *V* of a directed graph *G* = (*V*, *E*). Candidate directed edges e~ij~ = (v~i~ → v~j~) are generated for each node in four steps, and all four are exposed as configuration fields (Table S3). First, admissible partners are restricted to those within a search radius edge_radius_km (default 5 km), measured as haversine distance; the radius is an a priori, reproducibility-fixed benchmark default rather than a site-calibrated hydraulic length scale. A graph-stage sensitivity run at 3, 5, 7.5 and 10 km (Supplementary Table S14) changed the retained topology from 117–123 edges at Lower Anayari (Manu) and 85–162 at Talensi; Jaccard overlap with the 5-km graph ranged from 0.44–0.67 and 0.18–0.29, respectively. Field topology and downstream metrics are therefore conditional on the selected radius and should be re-run with a site-specific scale when such evidence is available. Second, each surviving pair is assigned a directional confidence

$$
p_{ij} = \Phi\!\left( \frac{h_i - h_j}{\sigma_{\Delta h}} \right)
$$

where Φ is the standard normal CDF, h~i~ and h~j~ are the measured or proxied hydraulic heads, and σ~Δh~ is the standard deviation of their difference, taken from the head covariance matrix where one is available and otherwise as √(σ~i~² + σ~j~²) from the per-node head uncertainties. Third, candidates are ranked by descending p~ij~ with haversine distance as the tie-break, and at most `edge_max_neighbors` primary edges are retained per node (default 3). Fourth, an edge is kept only if p~ij~ ≥ `edge_p_min` (default 0.75) and the hydraulic gradient |h~i~ − h~j~| / d~ij~ is at least `edge_gradient_min` (default 10⁻⁴); pairs that fall below the gradient floor but lie inside the search radius are retained separately as lateral, dispersive mixing candidates with p~ij~ shrunk toward 0.5 in proportion to the gradient deficit, and pairs whose screen depths differ by more than `edge_depth_mismatch` (default 20 m) are flagged. The maximum-neighbour bound is a deliberate compromise: too many neighbours inflate the false-positive rate through long-range edges that are rarely supported by hydrochemical evidence, whereas too few neighbours miss valid connections in sparse networks, where the true downstream pathway may skip over intermediate sites that were not sampled.

Direction assignment is governed by a ranked evidence hierarchy: (1) the measured or interpolated hydraulic head gradient, constituting the primary control; (2) the elevation difference, used as a head proxy where direct head measurements are unavailable; and (3) spatial proximity, applied only as a fallback when neither head nor elevation data are present. The tier actually used for each node is recorded on the edge, so head-derived and elevation-proxy edges remain distinguishable in the output. Candidate edges supported by user-specified external priors, such as MODPATH particle-tracking endpoint pairs, receive an additive prior weight drawn from a configurable physics-prior vector. An optional three-dimensional probabilistic builder is also provided for layered systems; it replaces the hard radius cut-off with a Gaussian distance decay, P(d) = exp(−d²/2r²), whose length scale r is the same `edge_radius_km` (default 5 km), so that P falls to 0.61 at one radius and 0.14 at two. All results reported in this paper use the two-dimensional builder described above, in which distance enters as the radius cut-off and the ranking tie-break only, and not as a multiplicative decay term. The computational output of this stage is a weighted directed multigraph whose edge set represents the full set of physically plausible pathway hypotheses available for subsequent refinement and fitting.

### 2.4 Sheaf-Inspired Edge Refinement and Consistency Scoring

Following initial construction, the candidate graph is refined using a sheaf-inspired consistency framework (Hansen and Ghrist, 2019; Robinson, 2020). The construction is defined precisely as follows, The term sheaf refers to the node stalks, edge restriction maps and coboundary diagnostic; edge selection itself is a weighted multi-criteria operation, and no graph-Laplacian eigendecomposition is used. Each node v~i~ carries a stalk s~i~ ∈ ℝ^d^, a local observation vector whose entries are the measured hydrochemical concentrations, isotope values, and hydraulic head (or elevation) at the sampling site. Each directed edge e~ij~ carries a restriction map from the upstream stalk to the downstream stalk, which encodes the expected transformation of the local observations under the proposed flow connection; the restriction map is the affine predicted-downstream map α~ij~s~i~ + δ~ij~ used in the edge-level consistency residual:


$$
\rho(e_{ij}) = \sqrt{\sum_{d} w_{d}\left(\alpha_{ij}\, s_{i,d} + \delta_{ij,d} - s_{j,d}\right)^{2}}
\tag{2}
$$

where w~d~ is the weight assigned to the *d*-th observation dimension, α~ij~ is an edge-scale parameter, and δ~ij,d~ is the expected offset along dimension *d* (for example, the sign-consistent change in hydraulic head or age). A directed edge is considered consistent if the observed downstream vector s~j~ is close to the value predicted by applying the restriction map to the upstream observation s~i~; the residual ρ(e~ij~) measures the mismatch. The global consistency metric for the network is obtained by aggregating the edge residuals over the retained edges.

Four evidence streams are evaluated within this construction: (i) hydraulic direction consistency, verified against head or elevation; (ii) conservative-tracer ordering, checked using Cl⁻ or EC as a dilution indicator; (iii) stable-isotope consistency, verified against the local meteoric water line (LMWL) and the expected direction of evaporative enrichment; and (iv) age-ordering coherence, assessed by the criterion that the downstream apparent age must not be significantly younger than the upstream apparent age along a proposed directed edge (§2.5). Each candidate edge receives a weighted multi-criteria score from these streams, and the retained set is chosen per node by a soft selection (inverse-temperature `sheaf_soft_beta`, default 2.0) capped at `edge_max_neighbors`; edges that lose this selection are removed. This selection step is a weighted multi-criteria score, and we describe it as such rather than as a sheaf-theoretic operation.

The sheaf structure enters as a separate diagnostic layer computed on the retained graph. Stacking the edge relations gives the affine coboundary operator *D*, with one row per edge and observation dimension,

$$
D_{(e,d),(u,d)} = \sqrt{w_e}\,\alpha_e, \qquad D_{(e,d),(v,d)} = -\sqrt{w_e},
\qquad b_{(e,d)} = -\sqrt{w_e}\,\delta_{e,d}
$$

so that D has shape (|E|·d) × (|V|·d) and the affine edge relations are represented by Dx = b. This is a cellular-sheaf coboundary matrix with d-dimensional node stalks and an upstream restriction map alpha_e I. The reported quantities are the homogeneous nullity H0 = dim ker D, the first cohomology dimension H1 = |E|·d − rank D, the affine obstruction energy min_x ||Dx − b||², the per-edge leverage (the increase in obstruction energy attributable to each edge), and the maximum cycle obstruction. H0 is therefore the dimension of the homogeneous nullspace, not by itself the dimension of an affine global-section set. An exact affine global section exists only when the obstruction energy is zero within numerical tolerance; if the system is consistent, its affine solution set has dimension H0, whereas a positive obstruction energy means that the affine solution set is empty. In particular, H0 = 0 means that the homogeneous system has only the zero solution, not that no affine node assignment exists. We do not compute a spectral decomposition of any Laplacian: localisation is performed by the leave-one-edge-out leverage statistic, not by eigenvectors.

Expressed plainly for a hydrogeological readership, the separation between the two layers matters. Edge selection is a weighted multi-criteria score, and we make no claim that it outperforms a well-designed score of the same four evidence streams. The coboundary layer adds a network-level closure diagnostic that no per-edge score can provide: it tests whether the affine right-hand side b is compatible with the retained relations, measures the residual obstruction when it is not, and localises that obstruction with leave-one-edge-out leverage. H0 describes the dimension of the homogeneous nullspace; it is not a yes/no test for affine solvability. A network can consist entirely of individually plausible edges and still have a positive affine obstruction because the inconsistency lives in cycles rather than in any single edge.

### 2.5 Residence Time and Age Constraint Module

When tracer data are available, each node is assigned a posterior residence-time estimate using a lumped-parameter model (LPM) inversion. For tritium (³H; half-life 12.32 yr, Lucas & Unterweger, 2000), apparent ages are estimated from measured ³H activities and the reconstructed atmospheric input function using an exponential or piston-flow model. For radiocarbon (¹⁴C), a dispersion-corrected inversion with dead-carbon adjustment is applied to yield mean residence times (MRT) in the range 10² to 10⁴ yr. Where ³H/³He data are available, tritiogenic helium ingrowth is used to compute apparent ages without reliance on historical input reconstruction (Supplementary Method S3). Posterior age estimates are expressed as probability intervals obtained from Monte Carlo propagation of analytical uncertainty or from a Bayesian network update (Supplementary Method S4).

The network Bayesian update propagates age-ordering constraints from directed edges into individual node posteriors. The posterior age at a downstream node v~j~ is conditioned on the ordering criterion

$$\tau_{j} \ge \tau_{i} + \varepsilon_{\min},$$

where epsilon_min is the minimum downstream increase (default min_downstream_increase = 0.0 yr); with the default, the downstream posterior must be no younger than the upstream posterior. The criterion is evaluated on the posterior intervals rather than on point estimates, which is essential for the wide, overlapping age ranges that older waters typically produce (Purtschert et al., 2022). An edge whose intervals overlap is retained and flagged as "not resolved at stated uncertainty"; such overlap cases are not counted as severe violations. Only non-overlapping reversals that exceed the severe-violation threshold (0.3 in log10 age units) receive the severe age-coherence flag. The current audit returns these overlap and severe flags; it does not apply a continuous overlap-proportional weight or silently convert the flags into an edge confidence. Any exclusion or penalty based on these flags is a downstream policy choice, not part of the field results reported here. Analysis of the synthetic benchmark age network (Section 4.4) shows that 84% of the violations of the ordering criterion are interval-overlap cases of this kind, so the audit does not treat imprecise age evidence as a hard rejection. When tracer data are entirely absent, the age module is disabled and topology refinement proceeds without age-ordering constraints, relying exclusively on head, chemistry and isotope evidence.

### 2.6 Transport and Inverse Reaction Fitting Along Edges

For each retained directed edge e~ij~, Hydrosheaf solves a two-stage decomposition that separates non-reactive transport effects from geochemical reactions before fitting reaction extents.

**Stage 1: transport decomposition.** The upstream concentration vector x~u~ ∈ ℝ^n^ (where *n* is the number of ions in the configured ion order) is transformed by one of two transport operators. Under the evaporation model, concentrations are scaled by an evaporation factor γ:


$$
x_{T} = \gamma\, x_{u}, \qquad \gamma \ge 1
\tag{3}
$$

Under the mixing model, the upstream vector is blended with an endmember vector x~e~:


$$
x_{T} = (1 - f)\, x_{u} + f\, x_{e}, \qquad 0 \le f \le 1
\tag{4}
$$

Both γ and *f* are optimised by a preliminary unconstrained least-squares fit using EC or TDS and the conservative ion signals (Cl⁻, Br⁻) as anchor constraints. The transport-corrected residual, representing the chemistry attributable to reactive processes alone, is:


$$
r = x_{v} - x_{T}
\tag{5}
$$

where x~v~ is the observed downstream concentration vector.

The choice between the evaporation and mixing operators is made algorithmically per edge: both models are fitted, and the better-fitting candidate is selected as the one minimising the combined weighted objective — the transport and residual-chemistry mismatch plus the L1, EC/TDS, isotope, and kinetic penalties. The transport fit itself is weighted towards the conservative anchor signals (EC/TDS, Cl⁻, Br⁻), and the transport-model probabilities reported with each edge are derived from the information scores of the candidate fits. This selection criterion is part of the configuration (Table S3) and requires no analyst intervention on individual edges.

Chloride and EC are treated as conservative anchors only where no halite dissolution, agricultural chloride input, or anthropogenic salinity is indicated; the framework reports transport-model fit diagnostics for every edge so that departures from conservative behaviour are visible rather than assumed away. The limits of this treatment are quantified by a multi-endmember stress test (Supplementary Method S2; Section 3.2): when the generating process couples two-endmember mixing with simultaneous reactions, the single-endmember transport model absorbs most of the reaction signal (median recovered halite extent 0.04 mmol/L against a true 0.40 mmol/L, and calcite 0.00 mmol/L against a true 0.20 mmol/L) while the chemistry fit remains high (R² = 0.999), demonstrating the mixing–reaction coupling limit. The distinguishability of halite dissolution, evapoconcentration and mixing under such conditions is discussed in Section 5.

**Stage 2: sparse reaction fitting.** The residual *r* is expressed as a linear combination of reaction stoichiometry columns drawn from a pre-compiled reaction dictionary *S* ∈ ℝ^(n×m)^, where *m* is the number of candidate reactions (Table S2). The dictionary encompasses carbonate dissolution (calcite, dolomite), evaporite dissolution (gypsum, halite), silicate weathering (albite), cation exchange (Ca–Na, Mg–Na), redox reactions (pyrite oxidation, denitrification), and nitrate-source inputs, with each column encoding the stoichiometric change in mmol/L per unit reaction extent. The reaction extent vector ξ ∈ ℝ^m^ is estimated by minimising a penalised weighted least-squares objective:


$$
\hat{\xi} = \arg\min_{\xi} \left[ \left\| W^{1/2}(S\xi - r) \right\|_2^2 + \lambda_1 \left\| P\xi \right\|_1 + \lambda_2 \left\| \xi \right\|_2^2 \right]
\tag{6}
$$

where *W* = diag(w~1~, …, w~n~) is a diagonal ion-weight matrix assigning higher weight to conservative tracers, *P* = diag(p~1~, …, p~m~) is a reaction-specific penalty-scale matrix incorporating geological bias, λ₁ is the L1 sparsity hyperparameter (default 0.002), and λ₂ is a Tikhonov regularisation coefficient. λ₂ defaults to 0.0; a numerical ridge floor at the 10⁻¹⁰ scale is added to the diagonal of the normal equations purely for the stability of the linear solve, and λ₂ is not optimised jointly with λ₁. The L1 penalty induces sparsity, selecting the smallest reaction set capable of explaining the transport-corrected residual and thereby reducing the non-uniqueness inherent in inverse hydrogeochemical modelling. The optimal λ₁ is identified over a grid sweep via the Akaike Information Criterion corrected for small samples (AICc). The matrix rank and condition number κ(*S*) are reported in the edge diagnostics, with systems having κ ≥ 100 flagged as ill-conditioned.

### 2.7 PHREEQC Thermodynamic Gates and Geological Bias Constraints

Chemical mass balance alone is insufficient to guarantee that inferred reaction extents are geochemically feasible: a statistically optimal ξ̂ may invoke mineral dissolution at supersaturation or precipitation under undersaturation. To prevent such solutions, PHREEQC-derived saturation indices (SI) are used as thermodynamic gates. For each candidate mineral *k*, forward speciation is run on the upstream chemistry to obtain SI~k~^(u)^ and on the downstream chemistry to obtain SI~k~^(v)^. The gate logic enforces: (i) if ξ̂~k~ > 0 (net dissolution inferred), then SI~k~^(u)^ < τ~SI~ (upstream undersaturated); and (ii) if ξ̂~k~ < 0 (net precipitation inferred), then SI~k~^(u)^ > τ~SI~ (upstream supersaturated). The diagnostic band half-width is set to τ~SI~ = 0.2 (Table S3). Reactions whose inferred extents violate these conditions are flagged with a thermodynamic bound-violation code, and PHREEQC-derived bounds are imposed as inequality constraints in a constrained re-fit. Concentration logic gates additionally suppress reactions that imply measurable changes in ions falling below analytical detection limits at the sampled site. Geological bias penalties p~k~ in the *P* matrix are assigned from a configurable mineralogy-prior dictionary that encodes the expected prevalence of each phase in the lithological setting: reactions inconsistent with the known mineralogy receive elevated p~k~ values that effectively require stronger chemical evidence before selection by the sparse solver. Forward feasibility of the surviving reaction sets is then evaluated with a PHREEQC-compatible mass-balance proxy (Section 4.5): the inferred extents are applied to reproduce the downstream composition using locked PHREEQC-consistent saturation fields, yielding edge-level RMSE and Nash–Sutcliffe efficiency (NSE) values. Live PHREEQC kinetic execution is not part of the current release, so the forward step is a screening-level feasibility check rather than a kinetic simulation.

### 2.8 Uncertainty Diagnostics and Process Stability Index

Uncertainty is characterised through three complementary modes. Monte Carlo propagation introduces aleatory uncertainty by perturbing input concentrations with Gaussian noise scaled to the configured analytical relative sigma — defaults of 4% for major ions and 0.5‰ for stable isotopes — and re-running the full edge-fitting pipeline across *N* realisations (default 100). These sigma values are configurable defaults representing routine analytical precision for major-ion and stable-isotope measurement, not universal analytical constants; they can be set to laboratory-specific values for individual datasets (Table S3). Bootstrap resampling introduces epistemic uncertainty by drawing with replacement from the available ion-concentration observations and re-fitting whilst holding mean inputs fixed. Variance decomposition assigns the total parameter variance of each output (γ, *f*, or ξ̂~k~) to its aleatory and epistemic components. A fully Bayesian alternative is also supported (Supplementary Method S5).

The process stability index (PSI) for reaction *k* on edge e~ij~ is defined as:


$$
\Psi_k(e_{ij}) = \frac{1}{N} \sum_{n=1}^{N} \mathbb{1}\left( |\hat{\xi}_{k,n}| \ge \varepsilon \right)
\tag{7}
$$

where 1[·] is the indicator function and ε is the activation threshold. A PSI value near unity indicates that reaction *k* is robustly activated across the full input-uncertainty envelope; values below 0.5 indicate that the inferred process is sensitive to analytical noise and should be treated as a provisional hypothesis requiring further confirmation. Sensitivity indices computed by finite-difference perturbation identify which measured species most strongly determine each inferred reaction extent, thereby guiding targeted re-sampling strategies. Missing-ion sensitivity analysis re-fits the reaction model after systematically zeroing the weights of each ion or ion group, and flags solutions that change qualitatively under ion loss as missing-ion sensitive.

### 2.9 Reproducible Implementation

Hydrosheaf is implemented in Python and organised as a structured package with command-line entry points, a Python API, and an automated workflow driver. Configuration is managed through a hierarchical settings object that exposes all hyperparameters — including ion order, transport-model selection, the λ₁ grid, the SI threshold, and the PSI trial count — as user-settable fields with documented defaults (Table S3). All outputs are written to a structured archive containing edge result tables, network summary files, diagnostic figures, and a JSON provenance manifest; the manifest schema and a worked example from the canonical runs are provided in the Supplementary Information. Validation scripts for the synthetic benchmark, the MODPATH topology test, and the field demonstration are supplied with the code repository, together with an automated test suite and a continuous-integration workflow, enabling independent replication of all reported results.

The revision is versioned: the source code is pinned to commit `463e1ce` (full SHA `463e1ceed6f7cb5fbe90cf96877d88cea181e6be`), with an exported environment file and installation instructions provided in the repository, and an archived, DOI-minted release is planned (minting and publication are pending the authors' explicit request). Candidate-edge construction scales approximately linearly in the number of nodes — about 3*n* candidate edges under maximum-neighbour pruning rather than the O(n²) of all-pairs enumeration — with a measured construction time of 0.06 s for a 320-node network. The source code, example workflows, validation scripts, and reproducibility materials are publicly available at https://github.com/dabdul-wahab1988/Hydrosheaf.

## 3. Benchmark and Validation Design

### 3.1 Validation Rationale

Hydrosheaf combines topology inference, age constraints, inverse reactions, thermodynamic screening, and uncertainty diagnostics; these components must be tested separately to avoid treating a successful field example as full method validation. The validation suite is therefore organised as a multi-tier benchmark that evaluates topology, reaction recovery, residence-time inference, PHREEQC feasibility, and field transferability as distinct evidence streams, with each tier reporting its own reference and metrics (Table 2).

### 3.2 Synthetic Reaction and Transport Recovery

Synthetic benchmarks use known transport perturbations and reaction extents to test whether Hydrosheaf can recover active processes under noise and missing-data scenarios. The benchmark comprises 100 locked Monte Carlo realisations in which known evaporation fractions, mixing fractions, and reaction extents are imposed on simulated chemistry, and the pipeline is evaluated on its ability to recover the generating parameters. These controlled tests assess identifiability and sparse model selection where field ground truth is unavailable. Because the synthetic tests deliberately share the reaction dictionary between data generation and inversion, the tier also reports identifiability diagnostics rather than relying on recovery scores alone: the dictionary (14 candidate reactions over 11 ions) has rank 8 of 11 with an effectively infinite condition number and 10 column pairs with |cosine| > 0.7, and leave-one-out dictionary-sensitivity runs quantify how single-reaction removals change recovery (Supplementary Method S2). A multi-endmember stress scenario couples two-endmember mixing with simultaneous halite and calcite dissolution to test the mixing–reaction coupling limit discussed in §2.6 (Supplementary Method S2).

### 3.3 MODPATH Topology Validation

Directed-edge recovery is evaluated against a published MODPATH particle-tracking archive — the USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 model archive (Harte, 2021; USGS data release, DOI 10.5066/F7J102FK). The tier serves two distinct roles. First, it is a physics-prior ingestion-fidelity check: when the MODPATH endpoint pairs are supplied as priors, the pipeline recovers exactly the 174 reference directed edges (precision = recall = F1 = 1.00), demonstrating that the prior machinery ingests external physics without loss. This perfect result is a pipeline-integrity check, not an inference test. Second, and primarily, the tier is an independent no-prior inference test: with MODPATH information withheld, topology is constructed from head gradients (with elevation as the head proxy) and hydrochemical evidence alone, yielding F1 = 0.62 (precision 0.49, recall 0.84; 147 true positives, 155 false positives, 27 false negatives) against the same 174-edge reference. The no-prior case is the primary inference test of the framework.

Agreement with MODPATH in either mode demonstrates recovery of well-to-well endpoint connectivity, not the pathline geometry or travel time that MODPATH produces: an inferred edge here is a directed connection supported by head or elevation gradients and chemical ordering, not a particle trajectory traced through a calibrated flow field, and travel-time-related outputs are known to be sensitive to the underlying flow model and porosity assumptions (Meyer et al., 2018; Baker et al., 2025). Three baseline graph-construction rules are compared against the same reference in the Supplementary Information: all-pairs elevation-drop construction (F1 = 0.47), proximity k-nearest-neighbour construction (F1 = 0.00), and conservative-tracer ordering, which could not be evaluated on the Savage archive because it contains no hydrochemistry. This separates physically admissible flow connectivity from chemical similarity, reducing the risk that a good reaction fit is mistaken for an independently supported groundwater pathway.

### 3.4 Residence-Time and PHREEQC Checks

Synthetic and public tracer-age benchmarks evaluate temporal consistency. The public tier uses the identifiability-gated subset of a published USGS tracer-age compilation (Jurgens et al., 2022; USGS data release, DOI 10.5066/P9W7T0DN), and the design explicitly separates identifiable from non-identifiable rows so that parity statistics are not inflated by unconstrained estimates. PHREEQC forward checks test whether inferred reactions remain saturation-index feasible, and accepted edges are required to be chemically, temporally, and thermodynamically plausible where supporting data exist.

### 3.5 Ghana Field Demonstration

Lower Anayari and Talensi are used as data-limited field demonstrations, not as the sole target of the study. They test whether the open workflow remains interpretable when optional tracers, head data, and geological priors are incomplete. Both sites are processed through the same graph construction as every other tier — the probabilistic builder of §2.3 followed by the refinement of §2.4 — with elevation substituting for hydraulic head, so the field topology is subject to the same directional confidence, gradient and neighbour limits as the benchmarked topology, and the same sheaf-cohomology diagnostics are reported for it. The role of this tier is to demonstrate workflow transferability under field data sparsity; it is not an independent validation of the hydrogeological process interpretation, and the site-specific context together with the independent evidence that is and is not available are stated in Section 4.6 and the Supplementary Information.


## 4. Results

All metrics reported in this section and in Tables 1–6 were regenerated from a
single locked pipeline (the canonical benchmark re-runs, the identifiability-gated
public-age benchmark, and the field/PSI chain) and audited for exact agreement
across text, tables and figures; the audit record is maintained in
metric_audit.md. Where a value differs from an earlier report, the regenerated
value is authoritative.

## 4.1 Workflow Execution and Reproducible Outputs

The Hydrosheaf pipeline executed end to end across every benchmark and field
scenario configured in the validation suite (Table 1), producing the full set of
outputs described in the framework design. For each scenario the
data-preprocessing module ingested sample chemistry, applied charge-balance
screening after unit conversion to a common millimolar basis, built the
data-availability matrix, and logged all preprocessing decisions to a JSON
provenance manifest. The graph-construction module generated candidate directed
edge sets, assigned confidence weights using the ranked evidence hierarchy, and
ingested MODPATH-derived physics priors where an archive was supplied. The
inverse solver exported per-edge result tables containing reaction extents,
transport parameters (γ for evaporation, f for mixing), residual norms,
thermodynamic feasibility flags, AICc values and chemistry R² statistics;
residence-time posteriors, PSI probability vectors and forward-validation
metrics were written as structured comma-separated files alongside
publication-ready figures.

Missing optional data were handled gracefully: modules dependent on isotope or
tracer fields that were absent were automatically disabled without interrupting
the run, confirming the graceful-degradation behaviour specified in the
framework design. Candidate-edge construction scaled near-linearly in the number
of samples once k-nearest-neighbour pruning was applied, generating
approximately three candidate edges per node and completing a 320-node graph in
0.06 s, consistent with the expected reduction from quadratic to linear scaling
(Supplementary Information).

## 4.2 Synthetic Benchmark Performance

The synthetic benchmark was evaluated over 100 locked Monte Carlo realisations
containing known evaporation fractions, mixing fractions and reaction extents.

**Transport recovery.** The transport stage recovered the generating model
reliably. The correct transport model (evaporation or mixing) was detected on
91.7% of edges overall — 75.5% for evaporation-dominated edges and 99.7% for
mixing-dominated edges — with a median absolute relative bias of 15.8% and a
median absolute parameter error of 0.058 on the recovered transport parameter.
Transport recovery is therefore the strongest claim the synthetic tier supports.

**Reaction-extent recovery.** Reaction-extent recovery was substantially more
limited, and is reported here with the diagnostics the identifiability analysis
provides. Across 2,100 active reaction rows (|true extent| > 0.01 mmol/L; 4,900
inactive rows), the correlation R² between recovered and true extents was 0.23,
with MAE = 0.37 mmol/L and RMSE = 0.62 mmol/L (Figure 3, panel B). The
false-activation rate, defined as the proportion of genuinely inactive reactions
inferred with extents exceeding 0.05 mmol/L, was 54.1%, indicating that the L1
penalty does not fully suppress spurious activations at this threshold. Median
chemistry R² across fitted edges remained 1.00 (Table 2): the inversion
reproduces observed compositions almost exactly, but it does so through extent
combinations that are not uniquely identified.

The gap between perfect chemical closure and weak extent-level recovery is
explained by reaction-dictionary degeneracy. The 14-reaction × 11-ion dictionary
has rank 8 of 11 (deficiency 3) and an effectively infinite condition number, so
the stoichiometric system is underdetermined; ten of the dictionary column pairs
exhibit |cosine| > 0.7, with calcite–anorthite and CaNa_exch–NaCa_exch collinear
at |cosine| = 1.00. Leave-one-out dictionary sensitivity over five realisations
confirms that no single reaction removal resolves the degeneracy: removing
albite improves active-reaction R² from 0.09 to 0.17 and lowers MAE from 0.40 to
0.34 mmol/L, but the remaining dictionary stays rank-deficient and the extent
recovery improves only marginally.

**Process-family recovery.** Aggregating to process-family level — the level at
which stability is ultimately reported through PSI — recovery remained limited
but quantifiable: family-level R² = 0.16, a sign-match rate of 36.5% and a
dominant-family hit rate of 48%. The supported claim class for reaction
inference is therefore process-family stability as quantified by PSI, not unique
extent-level attribution; extent-level claims are restricted to the limits
documented here and in the identifiability diagnostics of the Supplementary
Information.

**Isotope recovery.** Isotope-shift recovery reproduced the expected
point-versus-network contrast. Pointwise recovery was effectively null
(R² = −0.03, MAE = 0.58‰), but edge-mean aggregation recovered the signal
reliably (R² = 0.99) against an isotope-difference noise amplitude of
σ_Δ = 0.71‰ — that is, √2 times the 0.5‰ analytical sigma applied to each of the
two nodes forming the edge (Figure 3, panel D) —
indicating that isotope signals are recovered at the network scale even where
per-sample noise dominates.

**PHREEQC-compatible proxy.** Forward feasibility was evaluated with a
PHREEQC-compatible proxy that applies locked PHREEQC-consistent saturation
fields; live PHREEQC kinetic simulation is not executed in this version, and the
proxy results are labelled as such. The proxy achieved median RMSE = 0.02 mmol/L
and NSE = 1.00 across the validation suite, with a feasibility fraction of 0.88,
i.e. 88% of proposed reactions passed the saturation-index gates.

**Multi-endmember mixing–reaction stress test.** To probe the limits of the
single-endmember transport stage under coupled mixing and reaction, a dedicated
scenario generated two-endmember mixing (f1 = 0.20, f2 = 0.15; true combined
fraction 0.35) with simultaneous halite (0.40 mmol/L) and calcite (0.20 mmol/L)
extents. The single-endmember transport stage selected the mixing model on 73%
of edges, with a median recovered fraction of f = 0.24 against the true total of
0.35, and the mixing step absorbed most of the reaction signal: median recovered
halite extent fell to 0.04 mmol/L and calcite to 0.00 mmol/L, while median
chemistry R² remained at 0.999. The scenario demonstrates that mixing and
reaction are coupled in the framework's target settings and quantifies the
limits of distinguishing halite dissolution, evapoconcentration and mixing when
chloride cannot be treated as strictly conservative (per-edge results in the
Supplementary Information).

## 4.3 Directed Topology Recovery

Directed-edge recovery was evaluated against the USGS Savage Municipal
Water-Supply Well MODFLOW-2005/MODPATH5 archive, whose 174 particle-tracking
endpoint pairs define the reference connectivity (Harte, 2021; USGS data
release, DOI 10.5066/F7J102FK). The independent, no-prior mode — head-gradient
construction using elevation as the head proxy, downhill-directed edges and
k = 2 nearest-neighbour pruning — generated 302 candidate edges, of which 147
were true positives, 155 false positives and 27 false negatives (precision =
0.49, recall = 0.84, F1 = 0.62; Figure 2; Table 5). This is screening-level
performance with systematic overconnection (false positives exceed true
positives), which is consistent with the benchmark's information limit: the
Savage flow field is outlet-convergent, so true and false candidate edges
converge on the same outlet cells, and the FP:TP ratio largely reflects how many
up-gradient neighbours each outlet attracts under k = 2 pruning rather than a
systematic directional error.

When MODPATH endpoint pairs were supplied as physics priors, the graph
construction module reproduced all 174 reference edges exactly (174 true
positives, 0 false positives, 0 false negatives; precision = recall = F1 = 1.00;
Table 5). This prior-assisted result is a pipeline-integrity and ingestion check
— it confirms that supplied priors are ingested and propagated without loss —
and is not evidence of independent inference; the no-prior result carries the
inference claim. The same MODPATH archive therefore plays two distinct roles,
benchmark and prior, and the circularity this creates is bounded by reporting
the independent mode as the primary result.

Three baseline construction rules were compared against the same 174-edge
reference (Supplementary Information): all-pairs elevation-drop construction
achieved F1 = 0.47 (precision = 0.33, recall = 0.84 — the same recall as k = 2
pruning with roughly twice the false positives); proximity kNN with k = 2
produced F1 = 0.00 (306 inferred edges, zero true positives); and
conservative-tracer ordering could not be evaluated because the Savage archive
contains no hydrochemistry. Finally, the agreement demonstrated here is
well-to-well endpoint connectivity: an inferred edge records that two wells are
plausibly connected in the direction indicated by the hydraulic and chemical
evidence, not the pathline geometry, travel time or porosity-dependent transport
that MODPATH itself produces.

## 4.4 Residence-Time and Age-Order Performance

Pooled across the synthetic age classes, recovery of log10(age) achieved
correlation R² = 0.98 with a median absolute error of 17.86 y over 1,400 pooled
comparisons (Figure 5a; Table 2). Per-class accuracy follows the expected
information gradient (Table 4): the young/recharge class was recovered tightly
(R² = 0.97, MAE = 0.81 y), the intermediate class at R² = 0.85 with MAE = 15.75
y, and the old/deep class at R² = 0.82 with MAE = 188.49 y; the mixed and fossil
classes have no meaningful R² (NA) with MAEs of 103.23 y and 1,069.81 y
respectively, reflecting the wide posterior ranges that very old and mixed
waters impose without 39Ar or 81Kr constraints (Purtschert et al., 2022). The
age-order consistency index — the fraction of directed edges satisfying the
downstream ageing criterion τj ≥ τi − ε — was 0.85.

Wide posterior intervals feed back directly into the age-ordering constraint,
and the framework handles them explicitly. Across the canonical age-network
output, 15.2% of directed edges violated the ordering criterion; of these, 84.3%
were interval-overlap cases in which the posterior ranges of the two nodes
overlap but the order is not established, and these edges are retained with an
explicit "unresolved at stated uncertainty" flag rather than rejected. Only 2.4%
of edges were severe, non-overlapping reversals, which receive an
age-coherence-failure flag. The current audit does not apply a continuous
overlap-proportional weight or an automatic refinement penalty.

Against the public USGS tracer-age benchmark (the M3 tracerlpm parity
age-fractions dataset; Jurgens et al., 2022; USGS data release, DOI
10.5066/P9W7T0DN), fits were screened for identifiability before evaluation:
356 of the 1,272 fits were identifiable and the remaining 916 were flagged as
non-identifiable and excluded from the metric. On the identifiable subset the
median |log10 error| was 0.022, log10 RMSE was 0.235, 91.9% of estimates fell
within a factor of two and 98.6% within a factor of ten, with log10 R² = 0.960
(Figure 5b; Supplementary Figure S1). Because the USGS ages and age fractions
are reported model outputs rather than measured true ages, this is
screening-level parity with published model outputs, not true-age validation.

## 4.5 Regularisation, Thermodynamic Gating and Forward Validation

AICc-guided regularisation selected an L1 penalty of λ = 0.0483 at the AICc
minimum for both Lower Anayari and the Talensi Mining Area (Figure 6),
indicating that a single regularisation strength balances residual error and
parsimony across the two datasets. Thermodynamic gating applied PHREEQC
saturation-index checks with a diagnostic band half-width τSI = 0.2 (Table S3),
excluding dissolution extents proposed for supersaturated phases and
precipitation extents for undersaturated ones. Forward feasibility of the
surviving reaction sets was then evaluated with the PHREEQC-compatible proxy
described in Section 4.2: median RMSE = 0.02 mmol/L and NSE = 1.00 across the
validation suite, with a feasibility fraction of 0.88. These values confirm
that the inverse-derived reaction sets remain thermodynamically plausible at
screening level, subject to the proxy caveat stated explicitly in Section 4.2.

## 4.6 Field Demonstration and Process-Stability Results

The field demonstration on the Lower Anayari and Talensi datasets produced
interpretable directed hydrochemical evolution networks from sparse major-ion
and stable-isotope data (Figure 4; Table 6). The graph stage generated 572
candidate edges (422 Lower Anayari, 150 Talensi) and the refinement retained 258
— 121 for Lower Anayari and 137 for Talensi — rejecting 314. Because no
hydraulic-head records exist at either site, the directional confidence p~ij~
was evaluated on the elevation proxy with the elevation-tier sigma of 1.0 m
(Table S3); the median retained confidence was 0.93 at Talensi and 0.50 at
Lower Anayari, the latter because station elevations there are recorded as
constants, so most Lower Anayari edges fall below the gradient floor and are
retained explicitly as lateral, dispersive mixing candidates rather than as
gradient-supported directed edges (§2.3). No retained edge runs against the
elevation proxy at either site. The transport stage assigned mixing on 172 edges
and evaporation on 86. Median chemistry R² across all retained edges was 0.70
overall (0.53 for Lower Anayari, 0.82 for Talensi; Table 2; Supplementary Figure
S3), and median PSI was 0.97. Fourteen edges returned negative chemistry R²,
indicating that the fitted reaction set reproduces the downstream composition
worse than that composition's own mean; these are reported rather than removed.
The PSI family distribution was Evaporites 107, Conservative 84, Redox 55,
Carbonates 7 and Plagioclase 5, consistent with the evaporite-bearing,
mining-influenced crystalline-basement setting of the Talensi area (Song et al.,
2024) and the exchange-dominated basement aquifers of Lower Anayari
(Abdul-Wahab et al., 2021).

The sheaf-cohomology diagnostics separate the two networks in a way no per-edge statistic does. Lower Anayari returns homogeneous nullity H0 = 0 with H1 = 800 and an affine obstruction energy of 0.295 over 56 independent cycles. Because the obstruction energy is positive, the affine edge system has no exact global section at the reported numerical tolerance; H0 = 0 itself means that the homogeneous coboundary has only the zero solution. Talensi returns H0 = 10 with H1 = 760 and a lower but still positive obstruction energy of 0.128 over 74 cycles. Its H0 = 10 is a ten-dimensional homogeneous nullspace, not a ten-dimensional family of exact affine assignments, so the Talensi affine system also has no exact global section at the reported tolerance. Per-edge obstruction leverage localises the failure without any spectral computation: the most influential single edge accounts for 0.065 of the Lower Anayari obstruction against a median of 0.0018, and 0.017 against a median of 0.0007 at Talensi. The higher obstruction at Lower Anayari is consistent with its weaker directional evidence, since edges retained on chemistry alone are not constrained to close around cycles.

The seven zero-closure edges reported at the previous revision — all terminating on the sink wells NAB14 and THDW1, with chemistry R2 = 0.0 and all reaction extents approximately 0 — belonged to the superseded unfiltered as-run candidate graph. They are not part of the canonical 572-to-258 output, which contains zero unresolved null edges. Table 6 lists the top-ranked explained edges (all at the AICc-selected lambda = 0.0483; Figure 6), covering redox processes and nitrate input on Talensi edges. On one of those six edges the largest fitted extent and the most stable PSI family name the same process family; on the other five they disagree, and Table 6 reports the disagreement rather than resolving it, because it is the edge-level signature of the dictionary degeneracy quantified in Section 4.2.

PSI behaviour on the degenerate mineral pairs determines whether the stability
metric separates them in practice. At the site level, CaNa_exch PSI exceeds
NaCa_exch PSI at both sites — 0.73 versus 0.46 at Lower Anayari and 0.83 versus
0.42 at Talensi (Figure 7; Supplementary Table S10) — so the perturbation
envelope consistently prefers one exchange direction per edge. For calcite and
dolomite, by contrast, PSI values were both near zero (0.05 and 0.002 at Lower
Anayari; 0.06 and 0.04 at Talensi), so the metric does not falsely promote
either member of that pair when carbonates are not stably active (Carbonates
family: 7 of 258 edges). We emphasise that PSI measures stability to input
perturbation, not the correctness of process attribution.

Consistent with the validation design (Section 3.5), the field demonstration
tests workflow transferability and internal hydrochemical closure under
data-limited conditions, not independent process truth: no independent
mineralogical, hydraulic-head or tracer-truth dataset exists for either site
against which the inferred processes could be verified, and the results are
therefore reported as screening-level process hypotheses.

## Tables

**Table 1.** Hydrosheaf module architecture and validation links.

| Module                    | Scientific role                    | Main inputs                        | Main outputs          | Core method                 | Validation link |
| ------------------------- | ---------------------------------- | ---------------------------------- | --------------------- | --------------------------- | --------------- |
| Data preprocessing        | Harmonise field and chemistry data | ions, isotopes, coordinates        | cleaned dataset       | unit conversion, QC         | Table S1        |
| Graph construction        | Infer candidate flow topology      | coordinates, elevation, head proxy | directed edges        | spatial/head/chemical rules | Fig. 2          |
| Residence-time inversion  | Estimate age/MRT class             | tracers, isotope data              | posterior age classes | LPM/Bayesian update         | Fig. 5          |
| Inverse hydrogeochemistry | Fit reactions along edges          | major ions, minerals               | reaction extents      | logic gates, sparse fit     | Fig. 3, Fig. 6  |
| PHREEQC validation        | Check thermodynamic feasibility    | chemistry, mineral phases          | SI constraints        | forward diagnostics         | Fig. S2         |
| Robustness module         | Quantify process stability         | bootstrap/MC outputs               | PSI probabilities     | uncertainty propagation     | Fig. 7          |

**Table 2.** Multi-tier validation suite used to evaluate Hydrosheaf.

| Validation tier             | Dataset/source                                                                             | What is tested                  | Reference/target          | Main metric                                                                                              | Related figure  |
| --------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------- | ------------------------- | -------------------------------------------------------------------------------------------------------- | --------------- |
| Synthetic benchmark         | simulated benchmark suite                                                                  | transport and reaction recovery | known truth               | median chemistry R2=1.00 (reaction recovery in Fig. 3B)                                                  | Fig. 3          |
| MODPATH topology            | particle-tracking reference                                                                | directed-edge recovery          | MODPATH edges             | no-prior F1=0.62 (P=0.49, R=0.84); prior-assisted F1=1.00 (ingestion check)                              | Fig. 2          |
| Residence-time benchmarking | synthetic + M3 identifiability-gated public USGS benchmark (tracerlpm_parity_agefractions) | age agreement                   | known MRT/public age      | synthetic R2=0.98, median AE=17.86 y; public M3 (n=356 identifiable): median \|log10\|=0.02, 92% within 2x | Fig. 5, Fig. S1 |
| PHREEQC validation          | geochemical forward check                                                                  | reaction feasibility            | SI/forward model          | RMSE=0.02, NSE=1.00                                                                                      | Fig. S2         |
| Ghana field demonstration   | Lower Anayari/Talensi                                                                      | field process discovery         | hydrochemical consistency | median R2=0.70, median PSI=0.97                                                                          | Fig. 4, Fig. 7  |

**Table 3.** Performance and diagnostic metrics used in Hydrosheaf validation.

| Metric           | Formula/meaning                            | Used for                    | Interpretation                             |
| ---------------- | ------------------------------------------ | --------------------------- | ------------------------------------------ |
| Information gain | percent reduction in posterior uncertainty | temporal inversion          | higher means stronger network constraint   |
| PSI (edge)       | edge-level process inclusion probability   | spatial discovery           | probability path is robust to input noise  |
| PSI (region)     | regional phase stability index             | geologic province           | probability process exists in aquifer type |
| R2               | explained variance                         | age/chemistry recovery      | closer to 1 is stronger agreement          |
| RMSE             | root mean squared residual                 | chemistry/age               | lower values indicate better fit           |
| NSE              | Nash-Sutcliffe efficiency                  | forward validation          | >0.5 indicates useful predictive skill     |
| AICc             | small-sample model-selection criterion     | regularisation/model choice | lower favours parsimonious model           |
| F1-score         | harmonic mean of precision and recall      | topology inference          | balanced connectivity performance          |

**Table 4.** Residence-time benchmark accuracy by validation group.

| Validation group      | Reference age range (y) | Hydrosheaf inferred range (y) | R2   | MAE (y)             | Age-order consistency                                                          | Interpretation                                         |
| --------------------- | ----------------------- | ----------------------------- | ---- | ------------------- | ------------------------------------------------------------------------------ | ------------------------------------------------------ |
| Fossil                | 4500.0 - 4500.0         | 2255.0 - 9949.8               | NA   | 1069.81             | 0.85                                                                           | synthetic age-class recovery                           |
| Intermediate          | 60.0 - 260.0            | 38.5 - 433.3                  | 0.85 | 15.75               | 0.85                                                                           | synthetic age-class recovery                           |
| Mixed                 | 300.0 - 300.0           | 77.4 - 1047.9                 | NA   | 103.23              | 0.85                                                                           | synthetic age-class recovery                           |
| Old/Deep              | 520.0 - 2500.0          | 305.6 - 7370.6                | 0.82 | 188.49              | 0.85                                                                           | synthetic age-class recovery                           |
| Young/Recharge        | 5.0 - 25.0              | 3.7 - 34.4                    | 0.97 | 0.81                | 0.85                                                                           | synthetic age-class recovery                           |
| Public USGS screening | 5.5 - 49900.0           | 1.4 - 50000.0                 | 0.96 | median \|log10\| 0.02 | M3 identifiability-gated public USGS benchmark (tracerlpm_parity_agefractions) | n=356 identifiable fits; 92% within 2x, 99% within 10x |

**Table 5.** MODPATH topology agreement statistics.

| Site/model                                  | Candidate edges | MODPATH reference edges | True positives | False positives | False negatives | Precision | Recall | F1-score | Mode                                            |
| ------------------------------------------- | --------------- | ----------------------- | -------------- | --------------- | --------------- | --------- | ------ | -------- | ----------------------------------------------- |
| USGS MODPATH benchmark (no-prior inference) | 302.00          | 174.00                  | 147.00         | 155.00          | 27.00           | 0.49      | 0.84   | 0.62     | independent heuristic inference (head-gradient) |
| USGS MODPATH benchmark (prior-assisted)     | 174.00          | 174.00                  | 174.00         | 0.00            | 0.00            | 1.00      | 1.00   | 1.00     | physics-prior ingestion fidelity check          |

**Table 6.** Field discovery and process-stability results.

| Site    | Flow path/edge       | Dominant process | Reaction extent (mmol/L) | Selected lambda | RMSE/NSE                 | PSI probability | Interpretation                                                                                  |
| ------- | -------------------- | ---------------- | ------------------------ | --------------- | ------------------------ | --------------- | ----------------------------------------------------------------------------------------------- |
| Talensi | Talensi_TGW49->TGW36 | redox process    | 2.45                     | 0.0483          | objective 0.42 / R2 1.00 | 1.00            | provisional redox process; most stable PSI family: Evaporites (extent and stability disagree)   |
| Talensi | Talensi_TGW34->TGW19 | nitrate input    | 2.09                     | 0.0483          | objective 0.56 / R2 1.00 | 1.00            | provisional nitrate input; most stable PSI family: Conservative (extent and stability disagree) |
| Talensi | Talensi_TGW42->TGW45 | nitrate input    | 1.22                     | 0.0483          | objective 0.45 / R2 1.00 | 1.00            | provisional nitrate input; most stable PSI family: Redox (extent and stability disagree)        |
| Talensi | Talensi_TGW3->TGW8   | redox process    | 1.94                     | 0.0483          | objective 0.10 / R2 1.00 | 1.00            | redox process signal (PSI-supported)                                                            |
| Talensi | Talensi_TGW32->TGW12 | nitrate input    | 2.14                     | 0.0483          | objective 0.21 / R2 1.00 | 1.00            | provisional nitrate input; most stable PSI family: Conservative (extent and stability disagree) |
| Talensi | Talensi_TGW44->TGW45 | nitrate input    | 1.45                     | 0.0483          | objective 0.51 / R2 1.00 | 1.00            | provisional nitrate input; most stable PSI family: Redox (extent and stability disagree)        |

Note to Table 6. The two process columns answer different questions and are not expected to coincide. "Dominant process" is the reaction carrying the largest fitted extent on that edge; "most stable PSI family" is the process family with the highest Monte Carlo inclusion probability on the same edge. They agree on one of the six edges and disagree on the other five. That disagreement is a reported diagnostic, not a labelling error: it is the edge-level expression of the dictionary degeneracy quantified in Section 4.2 and Tables S6-S8b, and it is precisely why extent-level attribution is not claimed. The six edges listed are the highest-ranked explained edges by rank score (chemistry R2 x PSI), and are identical to those in Supplementary Table S5. The superseded seven null edges are absent from the canonical retained graph, rather than being hidden by this table.

## Figure captions

**Figure 1.** Hydrosheaf computational architecture for graph-constrained inverse
hydrogeochemical modelling. The workflow links data entry, hydrochemical and
physical input layers, data conditioning, topology/prior inference,
geochemical constraints, sparse inverse reaction fitting, optional uncertainty
validation, and final network-based outputs.

**Figure 2.** Topology validation of Hydrosheaf against MODPATH-derived
reference flow paths (USGS Savage Municipal Water-Supply Well archive, 174
reference edges). In the independent no-prior mode (head-gradient construction,
elevation-as-head, downhill, k = 2), Hydrosheaf recovered 147 true positives,
155 false positives and 27 false negatives (precision = 0.49, recall = 0.84,
F1 = 0.62). In prior-assisted mode all 174 reference edges were reproduced
exactly (174/0/0; F1 = 1.00); this mode is an ingestion-fidelity check of
MODPATH-derived priors, not independent inference. Agreement demonstrates
well-to-well endpoint connectivity, not pathline geometry or travel time.

**Figure 3.** Synthetic benchmark validation of Hydrosheaf inversion performance
across 100 locked realisations. Panel B shows active-reaction extent recovery
(correlation R² = 0.23, MAE = 0.37 mmol/L, RMSE = 0.62 mmol/L); the transport
stage detected the generating model on 91.7% of edges (evaporation 75.5%,
mixing 99.7%). Panel C shows robustness under degraded input conditions. Panel D
shows isotope-shift recovery: pointwise R² = −0.03 (MAE = 0.58‰) and edge-mean
R² = 0.99 against an isotope-difference noise amplitude σ_Δ = 0.71‰ (= √2 × the
0.5‰ per-node analytical sigma).

**Figure 4.** Spatiotemporal hydrogeochemical process-discovery networks for
Lower Anayari and Talensi. The figure shows all 258 retained directed
flow/process edges (121 Lower Anayari, 137 Talensi) with process-identification
probability networks, illustrating spatially variable groundwater evolution
pathways across elevation gradients and geochemical process families.

**Figure 5.** Residence-time benchmarking of Hydrosheaf age inference against
synthetic truth and the public USGS tracer-age reference. Panel a: synthetic
recovery, pooled log10 correlation R² = 0.98 with median absolute error 17.86 y.
Panel b: identifiability-gated public USGS parity (n = 356 identifiable fits of
1,272): median |log10 error| = 0.022, log10 RMSE = 0.235, 91.9% within a factor
of two, 98.6% within a factor of ten, log10 R² = 0.960 — screening-level parity
with reported USGS model outputs, not true-age validation.

**Figure 6.** Optimal Hydrosheaf model selection using AICc-minimum
regularisation. The AICc minimum selects λ = 0.0483 for both Lower Anayari and
the Talensi Mining Area, balancing model residual error and parsimony.

Figure 7. Process Stability Index (PSI) matrix for hydrogeochemical process identification under perturbation testing. The heatmap compares site-level stability summaries for inferred geochemical processes between Lower Anayari basement aquifers and the Talensi mining/agriculture aquifer system; its cells are not an edge-level median. Across retained field edges, the median PSI is 0.97. The PSI family distribution is Evaporites 107, Conservative 84, Redox 55, Carbonates 7 and Plagioclase 5.

**Supplementary Figure S1.** Public USGS residence-time parity check using the
M3 tracerlpm age-fraction-constrained benchmark, identifiability-gated. Of the
1,272 fitted rows, 356 are identifiable and 916 are flagged non-identifiable and
excluded; on the identifiable subset the median |log10 error| = 0.022, log10
RMSE = 0.235, 91.9% of estimates are within a factor of two and 98.6% within a
factor of ten, with log10 R² = 0.960. The USGS ages and age fractions are
reported model outputs, so this is screening-level parity with published model
outputs rather than true-age validation.

**Supplementary Figure S2.** PHREEQC-compatible forward-feasibility proxy. Panel
(a) shows the RMSE distribution (median = 0.02 mmol/L); panel (b) relates RMSE
to Nash–Sutcliffe efficiency (median NSE = 1.00); and panel (c) reports the
feasibility fraction (0.88) of proposed reactions passing the saturation-index
gates. The simulator uses locked PHREEQC-consistent saturation fields; this is a
PHREEQC-compatible proxy, and live PHREEQC kinetic execution is not part of the
current release.

**Supplementary Figure S3.** Field-scale chemical discovery fit for the Lower
Anayari and Talensi groundwater networks. The histograms show edge-level
chemistry R² for the 258 retained edges: median R² = 0.53 for Lower Anayari
(121 edges) and 0.82 for Talensi (137 edges), with an overall median of 0.70.
These fits describe hydrochemical closure of a field demonstration and do not
provide independent process-truth validation.





## 5. Discussion

### 5.1 Hydrosheaf as a Computational Geoscience Contribution

The question this study set out to answer was whether a single reproducible computational workflow
could convert sparse groundwater observations into directed, testable, uncertainty-ranked
hydrochemical evolution networks. The validation suite answers that question at two distinct
levels: directed-connectivity inference and the stability ranking of process-family hypotheses are supported across the
synthetic, benchmark and field tiers, whereas extent-level attribution of individual reactions is
not uniquely recoverable with the current reaction dictionary (Section 5.5). Within that claim
boundary, Hydrosheaf contributes a reusable framework whose combined inference structure is the
methodological novelty: directed graph construction defines candidate groundwater evolution
pathways; sheaf-style consistency scoring evaluates whether hydraulic, chemical, isotope and
residence-time evidence support each pathway; sparse inverse reaction fitting estimates parsimonious
reaction sets along retained edges; PHREEQC-derived saturation-index gates screen thermodynamic
feasibility; and the process-stability index (PSI) quantifies process robustness under input
uncertainty. Rather than replacing established tools such as PHREEQC or MODPATH, Hydrosheaf links
their respective strengths within a reproducible workflow designed for sparse multitracer datasets.
The integration converts manual decisions about flow-path selection, reaction dictionaries,
regularisation settings and thermodynamic screening into documented, configurable and repeatable
computational steps. The framework is therefore positioned not as a site-specific groundwater case
study but as a transferable method for graph-constrained inverse hydrogeochemical interpretation in
data-limited aquifers.

For a hydrogeological readership, the practical value of the sheaf construction (Section 2.4; Supplementary Method S1) can be stated plainly. Each node carries a stalk holding the locally observed chemical, isotope and head evidence at that sampling site, and the restriction map of a proposed directed edge predicts the downstream observation from the upstream one under that flow connection; the consistency residual is the weighted mismatch between prediction and observation, aggregated over the retained network. Relative to a simpler weighted score of the same four evidence streams, the construction adds a network-level closure diagnostic and a way to localise its obstruction through leave-one-edge-out leverage. It does not create an edge-level classifier that is claimed to outperform a carefully designed weighted score. On the field networks, Lower Anayari returns homogeneous nullity H0 = 0 and affine obstruction energy 0.295, while Talensi returns H0 = 10 and obstruction energy 0.128 over 137 edges (Section 4.6). Both positive obstruction energies mean that neither retained field network has an exact affine global section at the reported numerical tolerance; the H0 values describe the homogeneous nullspaces, not exact assignment counts. Edge selection itself remains a weighted multi-criteria score, and the cohomology is a diagnostic on the result, not a better classifier. Neither layer resolves reaction-level non-uniqueness, which is governed by the stoichiometric dictionary (Section 5.5).

### 5.2 Distinction from Existing Tools and Approaches

PHREEQC and NETPATH solve geochemical mass-balance problems but require the analyst to supply a
pre-selected initial–final water pair representing an assumed flow path (Parkhurst & Appelo, 2013);
Hydrosheaf replaces this assumption with principled directed-graph inference, running inverse models
only on topology-supported edge pairs and thereby removing subjective pathway selection as a source
of non-uniqueness (Manu et al., 2023; Slimani et al., 2017). Combinatorial inverse modelling, as
implemented by Manu et al. (2023) for the Pra Basin, partially automates this selection by
enumerating candidate water pairs and testing them exhaustively. The practical distinction is that
combinatorial enumeration still scores pairs against the chemistry alone; it does not score candidate
connections against hydraulic, tracer, isotope and age-order evidence jointly, nor does it separate
non-reactive transport from residual reactions with uncertainty ranking as a routine output.
Hydrosheaf embeds the pairwise-selection step inside the graph stage, where candidate edges already
carry directional, spatial and age evidence, and couples that selection with the transport–reaction
decomposition and the PSI diagnostics. The automation is therefore not simply the enumeration step
moved into code; it is the joint scoring of connectivity and chemistry that the graph structure
provides.

MODFLOW/MODPATH provides physically based pathline inference through calibrated flow fields, but it
demands a calibrated hydraulic-conductivity field, boundary conditions and effective-porosity
assumptions that are usually unavailable in the data-limited settings Hydrosheaf targets; that demand
is precisely what motivates Hydrosheaf's design, which operates from basic chemistry, coordinates and,
where available, head or elevation data. Where a MODPATH archive already exists, it serves a dual role
in this study. First, it is a source of physics priors: in the prior-assisted mode, the pipeline
reproduced all 174 reference directed edges exactly (precision = recall = F1 = 1.00). This result is
an ingestion-fidelity check — it confirms that externally supplied connectivity information is
ingested without loss — and it is not an independent validation of the inference. Second, the same
archive provides an independent connectivity benchmark: in the no-prior mode, in which edges are
inferred from head/elevation gradients and chemical ordering alone, Hydrosheaf achieves F1 = 0.62
(precision 0.49, recall 0.84) against the 174-edge reference (Section 4.3; Figure 2; Table 5). The
circularity risk is stated explicitly: the prior-assisted result must not be cited as independent
validation, because the reference information entered the pipeline as input; the independent no-prior
result carries the inference claim. The claim itself should also be qualified in scope. A directed
edge inferred from a head or elevation gradient and chemical ordering is conceptually weaker than an
advective pathline traced through a calibrated flow field; agreement with MODPATH endpoints
demonstrates recovery of well-to-well connectivity, not of pathline geometry, travel time or
porosity-dependent transport, all of which are sensitive to model and porosity assumptions (Meyer et
al., 2018; Baker et al., 2025). The reference archive is the published USGS Savage Municipal
Water-Supply Well MODFLOW-2005/MODPATH5 model (Harte, 2021; USGS data release,
DOI 10.5066/F7J102FK).

Hydrochemical diagram methods such as Piper, Gibbs and facies classification are interpretively
useful but do not compute reaction extents, propagate uncertainty or produce testable network
outputs. Graph-based and data-driven connectivity models (Bai & Tahmasebi, 2023; Taccari et al.,
2024) learn statistical associations in groundwater level time series but do not produce
geochemically interpretable mass balances or age-constrained flow paths. Python wrapper combinations —
for example, PHREEQC bindings such as phreeqpy or iPhreeqc used together with a graph library such as
NetworkX — could in principle partially replicate components of the framework, and the literature
reviewed here documents no published, versioned workflow that integrates transport decomposition,
sparse inverse fitting, thermodynamic gating, uncertainty diagnostics and a multi-tier validation
suite under one configurable pipeline. The integration gap is therefore not the existence of the
individual components but their joint, versioned assembly with documented defaults and reproducible
outputs, which is what Hydrosheaf provides. The combination of graph topology inference,
transport–reaction decomposition, thermodynamic gating, age constraints and PSI diagnostics is not
replicated by any single existing tool or workflow documented in the literature reviewed here
(Viswanathan et al., 2018; Schiavo et al., 2022; Rädle et al., 2022).

### 5.3 Interpretation of the Validation Results

The multi-tier validation is organised as a validation strategy rather than a single case study: each
tier isolates one claim class, and the results must be read at the claim class they support. At the
extent level, the synthetic benchmark shows that recovery of individual reaction extents is limited
under the locked pipeline: active-reaction recovery gives a correlation R² = 0.23 with mean absolute
error 0.37 mmol/L and RMSE 0.62 mmol/L, and 54.1% of genuinely inactive reactions are assigned
extents above the 0.05 mmol/L activation threshold (Figure 3, Panel B; Table 2). These values
supersede those reported at initial submission, which were drawn from unarchived analysis snapshots;
they are regenerated end to end from the pinned repository revision (Code and data availability). The
difference reflects a correction to the transport-selection logic in the canonical pipeline, and it
is reported because the revised framework no longer claims near-perfect extent recovery. Because the
benchmark generates its chemistry from the same reaction dictionary used for inversion, it is best
interpreted as an identifiability stress test under idealised conditions rather than as independent
evidence of field-level recovery. At the transport level, the same benchmark recovers the generating
model in 91.7% of edges (evaporation 75.5%, mixing 99.7%) with median absolute relative bias 15.8%,
and median edge-level chemistry R² is 1.00, confirming that mass-balance closure is achievable along
retained edges (Figure 3, Panel A). Isotope shifts are not recoverable at the point scale
(R² = −0.03, MAE 0.58‰) but are recovered at the network scale (edge-mean R² = 0.99, noise amplitude
σ_Δ = 0.71‰) (Figure 3, Panel D).

The topology tier separates the inference claim from the ingestion check (Section 5.2): the no-prior
mode achieves F1 = 0.62 against the MODPATH reference (302 candidates, 147 true positives, 155 false
positives, 27 false negatives), while the prior-assisted mode achieves F1 = 1.00 as an
ingestion-fidelity result (Figure 2; Table 5). Baseline graph-construction rules confirm the value of
the head-gradient evidence hierarchy: all-pairs elevation-drop construction yields F1 = 0.47 (the
same recall as the k2 rule but twice the false positives), proximity-based k-nearest-neighbour
construction yields F1 = 0.00, and conservative-tracer ordering cannot be evaluated on the Savage
archive, which contains no hydrochemistry (Supplementary Table S9). Residence-time
recovery is strong at the network level: pooled log10 correlation R² = 0.98 with median absolute
error 17.86 years across age classes (Figure 5; Table 4), an age-order consistency index of 0.85, and
— after identifiability gating — a public USGS parity result of median |log10 error| 0.022, log10
RMSE 0.235, 91.9% of estimates within a factor of two and 98.6% within a factor of ten
(log10 R² = 0.960; 356 identifiable rows of 1,272, with 916 non-identifiable rows flagged and
excluded; Figure 5, Panel b; Supplementary Figure S1). The PHREEQC forward checks confirm the
thermodynamic feasibility of the inferred reaction sets (median RMSE 0.02 mmol/L, median NSE 1.00,
feasibility fraction 0.88) (Table 2).

The Ghana field demonstration shows that the workflow remains interpretable under real data-limited
conditions: 258 of 572 candidate edges were retained across both sites with median chemistry R² = 0.70
(Lower Anayari 0.53, Talensi 0.82), median PSI = 0.97, and process-family counts of Evaporites 107,
Conservative 84, Redox 55, Carbonates 7 and Plagioclase 5 (Figure 4; Table 6). The seven zero-closure edges
from the previous revision belonged to the superseded unfiltered candidate graph, not to the canonical
572-to-258 output, which contains zero unresolved null edges. Table 6 therefore retains the disagreement
between the dominant extent label and the most stable PSI family on five of six edges; only one edge
agrees, and no extent-level attribution is claimed. Two claim classes are therefore supported by the validation suite:
the framework ranks the stability of process-family hypotheses under the specified perturbation model — the present benchmarks do not establish reliable identification of the true reaction family — while extent-level attribution is not uniquely
recoverable with the current dictionary (Section 5.5). What remains unvalidated is the absolute
accuracy of inferred reaction extents against independent field mineralogical measurements, and the
generalisability of the topology inference to settings with strongly transient or seasonally
reversing flow directions.

### 5.4 Practical Relevance and Transferability

Hydrosheaf is designed for contexts where dense monitoring networks and calibrated numerical flow
models are unavailable but basic chemistry, coordinates, elevation and optionally isotope or tracer
data exist — conditions common across the semi-arid tropics, data-sparse basement aquifer systems and
early-stage investigation settings worldwide (Borzí, 2025; Newman et al., 2021). The minimum viable
input schema requires only major-ion concentrations and spatial coordinates, so the framework can
produce directed process networks even in the most data-limited scenarios.

The Northern Ghana applications illustrate the intended use, and their evidentiary scope should be read precisely. Both sites lie in the crystalline basement terrain of the Upper East Region: Lower Anayari is a crystalline/alluvial setting with agricultural land use (Abdul-Wahab et al., 2021), while Talensi is a mining-influenced crystalline basement catchment (Song et al., 2024). The available evidence comprises major-ion chemistry and stable water isotopes (delta18O and delta2H); no nuclear tracers (3H, 14C) and no hydraulic-head data were available at either site, so age-ordering constraints were disabled and the directional confidence pij was evaluated on the elevation proxy with the elevation-tier sigma (Table S3). That substitution is weaker at Lower Anayari than at Talensi: station elevations at Lower Anayari are recorded as constants, so most retained edges there fall below the gradient floor and are kept as lateral mixing candidates at pij = 0.50 rather than as gradient-supported directed edges, and the affine cohomology diagnostic returns homogeneous nullity H0 = 0 with obstruction energy 0.295 (Section 4.6); this is a trivial homogeneous nullspace with positive affine inconsistency, not an independent claim that no node assignment exists. Independent process-truth labels, mineralogical confirmation of the inferred phases, and hydraulic-head records — the evidence that would permit independent validation of the inferred flow paths and reactions — are absent. The field demonstration is therefore exactly that: a demonstration that the workflow transfers to real sparse data and produces internally consistent, PSI-ranked process networks, not an independent validation of the hydrogeological interpretation. High chemical reconstruction performance does not by itself confirm the inferred flow paths or reactions.

### 5.5 Limitations and Risk Control

Several limitations constrain the scope of Hydrosheaf's outputs, and the principal ones are now
quantified. The reaction dictionary (14 reactions × 11 ions) is rank-deficient: it has rank 8 of a
possible 11, an infinite condition number (rank deficiency 3), and 10 column pairs with |cosine| >
0.7, including calcite–anorthite and CaNa-exchange–NaCa-exchange at |cosine| = 1.00. This structure
directly limits extent-level recovery: family-level recovery achieves R² = 0.16 with 36.5% sign
match and a 48% dominant-family hit rate, and leave-one-out dictionary sensitivity tests show that no
single reaction removal resolves the degeneracy (removing albite improves active-recovery R² from
0.09 to 0.17 and MAE from 0.40 to 0.34 mmol/L, but the residual collinearity persists). The PSI
metric behaves informatively on the degenerate pairs: it separates CaNa-exchange from NaCa-exchange
on the field edges (CaNa_exch PSI 0.73 vs NaCa_exch 0.46 at Lower Anayari and 0.83 vs 0.42 at Talensi; Figure 7), whereas both
calcite and dolomite PSI are near zero on the sampled field edges, so carbonate extents are not
robustly attributable on these data. We emphasise that PSI measures stability to input perturbation, not the correctness of process attribution: a high PSI indicates that the selected hypothesis is robust under the imposed noise model, not that it is the true generating process. This behaviour echoes the non-uniqueness documented for
Ca-HCO3⁻ water types, where carbonate dissolution, silicate weathering and mixing can generate
indistinguishable signatures (Zhi et al., 2021). The implication for interpretation is explicit:
Hydrosheaf outputs are screening-level, uncertainty-ranked process hypotheses; PSI ranks how stable
each process-family hypothesis is under the specified perturbation model, but the present benchmarks
do not establish reliable identification of the true reaction family, and the framework cannot
uniquely attribute extents for degenerate pairs. Neither the family ranking nor the extents should be
used beyond screening purposes without independent confirmation.

A second, distinct limitation concerns the coupling of mixing and reaction. In the multi-endmember
stress test, two-endmember mixing was combined with simultaneous halite dissolution (0.40 mmol/L)
and calcite dissolution (0.20 mmol/L); because the transport stage assumes a single endmember, the
misspecified transport model absorbs most of the reaction signal — recovered halite extents were
0.04 mmol/L and calcite extents 0.00 mmol/L — while the chemistry fit remained effectively perfect
(R² = 0.999). Chloride and EC therefore cannot always be treated as conservative anchors in
evaporite-bearing, mining-impacted or agricultural systems, where halite dissolution,
evapoconcentration, agricultural chloride input and mixing can confound the transport correction.
The framework flags this condition through transport-fit diagnostics and through low PSI on reactions
whose signal is absorbed by the transport stage; where such flags appear, independent tracer
sampling — for example 87Sr/86Sr, boron isotopes or nitrate isotopes — is recommended to confirm the
processes before interpretation.

Other limitations follow from the data and design assumptions. Sparse major-ion chemistry, particularly in the absence of trace elements such as fluoride, iron and phosphate, limits the number of independently resolvable reactions, and the condition-number diagnostics flag ill-conditioned systems in which the stoichiometric matrix cannot support unique reaction attribution. Missing tracer data prevent age-ordering constraints from filtering topologically inconsistent edges, which may raise the false-positive edge rate in datasets lacking 3H, 14C or stable isotopes. The graph-inference module performs poorly under strongly transient or reversing flow conditions, because the evidence hierarchy assumes steady or time-averaged directional consistency. Where age posteriors are wide — as for fossil and mixed waters, whose reference ages are poorly constrained without argon-39 or krypton-81 tracers — the age-ordering criterion is evaluated on posterior intervals: edges violating the ordering beyond tolerance receive an age-coherence failure flag, while overlapping intervals are retained with an explicit "not resolved at stated uncertainty" flag. The current audit does not apply a continuous overlap-proportional weight; any exclusion or penalty based on these flags is a downstream policy choice rather than part of the reported field results (Section 2.5). Residual non-uniqueness therefore persists even after L1 regularisation and PHREEQC gating, and low-PSI inferences should be treated as provisional process hypotheses rather than confirmed geochemical interpretations; they identify candidate processes worthy of confirmation by repeat sampling, mineralogical analysis or independent tracer tests. Hydrosheaf reduces non-uniqueness substantially relative to unconstrained manual inverse modelling, but it does not eliminate it, and its outputs are best described as transparent, reproducible, uncertainty-ranked process hypotheses suited to screening-level diagnosis rather than to replacing calibrated reactive transport simulation.

## Code and data availability

The Hydrosheaf source code, documentation and general benchmark infrastructure are available in the public GitHub repository: https://github.com/dabdul-wahab1988/Hydrosheaf. The exact analysis snapshot used in this manuscript is identified by commit 463e1ce on branch codex/m3-correctness, but that revision and the revised M2 field package were not present on the public origin when this package was frozen. The revision materials therefore include a pre-release reproducibility manifest with SHA-256 hashes, the Python 3.14.6 environment export, the locked configurations, the benchmark scripts and the document-generation inputs. The package also identifies the automated test suite (tests/) and the continuous-integration workflow (.github/workflows/ci.yml). Benchmark scripts reproducing the figures and tables are provided under M2/m2_benchmark/scripts/ (synthetic benchmark, MODPATH topology validation, PHREEQC-compatible forward checks, field demonstration and radius sensitivity) and M3/m3_age_benchmark/scripts/ (public residence-time benchmark). Test datasets are included in the repository; third-party datasets are cited by their persistent identifiers — the USGS public-supply groundwater-age data release (Jurgens et al., 2022; DOI 10.5066/P9W7T0DN) and the USGS Savage Municipal Water-Supply Well MODPATH data release (Harte, 2021; DOI 10.5066/F7J102FK). Runtime and scalability are reported in the Supplementary Information: candidate-edge construction is pruned by the maximum-neighbour rule and scales approximately linearly with node count (approximately 3n candidates; 0.06 s at 320 nodes, sub-second for typical field networks), and the full 100-realisation benchmark runtime is reported in Table S13. Every pipeline run writes a JSON provenance manifest recording inputs, configuration, module decisions and outputs. The authors must publish the exact analysis snapshot as a versioned release and add a persistent DOI before resubmission; until then, the revised field package is a pre-release artifact and the public repository alone does not provide complete reproduction of the reported field results.

## References

Abdul-Wahab, D., Adomako, D., Abass, G., Adotey, D. K., Anornu, G., & Ganyaglo, S. (2021). Hydrogeochemical and isotopic assessment for characterizing groundwater quality and recharge processes in the Lower Anayari catchment of the Upper East Region, Ghana. Environment, Development and Sustainability, 23(4), 5297–5315. https://doi.org/10.1007/s10668-020-00815-w

Aladejana, J., Kalin, R., Hassan, I., Sentenac, P., & Tijani, M. (2020). Origin and residence time of groundwater in the shallow coastal aquifer of Eastern Dahomey Basin, Southwestern Nigeria, using δ18O and δD isotopes. Applied Sciences, 10, 7980. https://doi.org/10.3390/app10227980

Alberti, L., Colombo, L., & Formentin, G. (2018). Null-space Monte Carlo particle tracking to assess groundwater PCE (tetrachloroethene) diffuse pollution in north-eastern Milan functional urban area. Science of the Total Environment, 621, 326–339. https://doi.org/10.1016/j.scitotenv.2017.11.253

Bai, T., & Tahmasebi, P. (2023). Graph neural network for groundwater level forecasting. Journal of Hydrology, 616, 128792. https://doi.org/10.1016/j.jhydrol.2022.128792

Baker, E. A., Juckem, P., Feinstein, D., & Hart, D. (2025). A regional model comparison between MODPATH and MT3D of groundwater travel time distributions. Groundwater, 63(6), 861–873. https://doi.org/10.1111/gwat.70024

Bicalho, C. C., Batiot-Guilhe, C., Seidel, J. L., Tavignot, C., & Jourde, H. (2017). Geochemical investigation of groundwater dynamics and flow paths in a karst aquifer (Lez spring, southern France). Journal of Hydrology, 554, 789–805. https://doi.org/10.1016/j.jhydrol.2017.09.025

Binet, S., Joigneaux, E., Pauwels, H., Albéric, P., Fléhoc, C., & Bruand, A. (2017). Water exchange, mixing and transient storage between a saturated karstic conduit and the surrounding aquifer: Groundwater flow modeling and inputs from stable water isotopes. Journal of Hydrology, 544, 278–289. https://doi.org/10.1016/j.jhydrol.2016.11.042

Borzí, I. (2025). Modeling groundwater resources in data-scarce regions for sustainable management: Methodologies and limits. Hydrology, 12, 11. https://doi.org/10.3390/hydrology12010011

Cao, F., Jaunat, J., Vergnaud-Ayraud, V., Devau, N., Labasque, T., Guillou, A., Guillaneuf, A., Hubert, J., Aquilina, L., & Ollivier, P. (2020). Heterogeneous behaviour of unconfined chalk aquifers inferred from combination of groundwater residence time, hydrochemistry and hydrodynamic tools. Journal of Hydrology, 581, 124433. https://doi.org/10.1016/j.jhydrol.2019.124433

Cartwright, I., Cendón, D., Currell, M., & Meredith, K. (2017). A review of radioactive isotopes and other residence time tracers in understanding groundwater recharge: Possibilities, challenges, and limitations. Journal of Hydrology, 555, 797–811. https://doi.org/10.1016/j.jhydrol.2017.10.053

Caschetto, M., Colombani, N., Mastrocicco, M., Petitta, M., & Aravena, R. (2016). Estimating groundwater residence time and recharge patterns in a saline coastal aquifer. Hydrological Processes, 30, 4202–4213. https://doi.org/10.1002/hyp.10942

Casillas-Trasvina, A., Rogiers, B., Beerten, K., Pärn, J., Wouters, L., & Walraevens, K. (2022). Using helium-4, tritium, carbon-14 and other hydrogeochemical evidence to evaluate the groundwater age distribution: The case of the Neogene aquifer, Belgium. Journal of Hydrology X, 16, 100132. https://doi.org/10.1016/j.hydroa.2022.100132

Cheng, L., Jiang, C., Li, C., & Zheng, L. (2022). Tracing sulfate source and transformation in the groundwater of the Linhuan Coal Mining Area, Huaibei Coalfield, China. International Journal of Environmental Research and Public Health, 19, 14434. https://doi.org/10.3390/ijerph192114434

Colombo, L., Gzyl, G., Mazzon, P., Łabaj, P., Frączek, R., & Alberti, L. (2021). Stochastic particle tracking application in different urban areas in Central Europe: The Milano (IT) and Jaworzno (PL) case study to secure the drinking water resources. Sustainability, 13(18), 10291. https://doi.org/10.3390/su131810291

Comte, J.-C., Cassidy, R., Nitsche, J., Ofterdinger, U., Flynn, R., & Allen, A. R. (2018). Catchment-scale residence-time structure in weathered/fractured hard-rock aquifers. Water Resources Research, 54, 8966–8987. https://doi.org/10.1029/2018WR022974

Díaz-Puga, M. A., Vallejos, Á., Sola, F., Daniele, L., Molina, L., & Pulido-Bosch, A. (2016). Groundwater flow and residence time in a karst aquifer using ion and isotope characterization. International Journal of Environmental Science and Technology, 13, 2579–2596. https://doi.org/10.1007/s13762-016-1094-0

Fandel, C., Miville, F., Ferré, T., Goldscheider, N., & Renard, P. (2022). The stochastic simulation of karst conduit network structure using anisotropic fast marching, and its application to a geologically complex alpine karst system. Hydrogeology Journal, 30, 927–946. https://doi.org/10.1007/s10040-022-02464-x

Ghiglieri, G., Buttau, C., Arras, C., Funedda, A., Soler, A., Barbieri, M., Carrey, R., Domènech, C., Torrentó, C., Otero, N., & Carletti, A. (2021). Using a multi-disciplinary approach to characterize groundwater systems in arid and semi-arid environments: The case of Biskra and Batna regions (NE Algeria). Science of the Total Environment, 757, 143797. https://doi.org/10.1016/j.scitotenv.2020.143797

Gonzalez, D., Janardhanan, S., Pagendam, D. E., & Gladish, D. W. (2020). Probabilistic groundwater flow, particle tracking and uncertainty analysis for environmental receptor vulnerability assessment of a coal seam gas project. Water, 12(11), 3177. https://doi.org/10.3390/w12113177

Guo, L., Ding, Y., Fang, H., An, C., Jiang, W., & Yang, N. (2025). Integrating inverse modeling to investigate hydrochemical evolution in arid endorheic watersheds: A case study from the Qaidam Basin, Northwestern China. Water, 17, 2074. https://doi.org/10.3390/w17142074

Guo, X.-R., Zuo, R., Wang, J., Meng, L., Teng, Y., Shi, R., Gao, X., & Ding, F. (2018). Hydrogeochemical evolution of interaction between surface water and groundwater affected by exploitation. Groundwater, 57. https://doi.org/10.1111/gwat.12805

Hansen, J., & Ghrist, R. (2019). Toward a spectral theory of cellular sheaves. Journal of Applied and Computational Topology, 3, 315–358. https://doi.org/10.1007/s41468-019-00038-7

Harte, P. T. (2021). MODFLOW-2005, MODPATH, and MOC3D used for groundwater flow simulation, pathlines analysis, and solute transport in the crystalline-rock aquifer in the vicinity of the Savage Municipal Water-Supply Well Superfund Site, Milford, New Hampshire. U.S. Geological Survey data release. https://doi.org/10.5066/F7J102FK

He, Z., Han, D., Song, X., Yang, L., Zhang, Y., Ma, Y., Bu, H., Li, B., & Yang, S. (2021). Variations of groundwater dynamics in alluvial aquifers with reclaimed water restoring the overlying river, Beijing, China. Water, 13, 806. https://doi.org/10.3390/w13060806

Joshi, S. K., Rai, S. P., Sinha, R., Gupta, S., Densmore, A. L., Rawat, Y. S., & Shekhar, S. (2018). Tracing groundwater recharge sources and processes in the northwestern Indian alluvial aquifer using environmental isotopes (δ18O, δD, and 3H). Journal of Hydrology, 559, 835–847. https://doi.org/10.1016/j.jhydrol.2018.02.056

Juckem, P. F., & Starn, J. J. (2021). Re-purposing groundwater flow models for age assessments: Important characteristics. Groundwater, 59(5), 710–727. https://doi.org/10.1111/gwat.13088

Jurgens, B. C., et al. (2022). Data for distribution of groundwater age in aquifers used for public supply in the continental United States, 2004–2017 (Version 1.1). U.S. Geological Survey data release. https://doi.org/10.5066/P9W7T0DN

Kaandorp, V. P., Broers, H. P., van der Velde, Y., Rozemeijer, J., & de Louw, P. G. B. (2021). Time lags of nitrate, chloride, and tritium in streams assessed by dynamic groundwater flow tracking in a lowland landscape. Hydrology and Earth System Sciences, 25, 3691–3711. https://doi.org/10.5194/hess-25-3691-2021

Kambuku, D., Tsujimura, M., & Kagawa, S. (2018). Groundwater recharge and flow processes as revealed by stable isotopes and geochemistry in fractured Hornblende-biotite-gneiss, Rivirivi Catchment, Malawi. African Journal of Environmental Science and Technology, 12(1), 1–14. https://doi.org/10.5897/AJEST2017.2406

Karra, S., O'Malley, D., Hyman, J. D., Viswanathan, H. S., & Srinivasan, G. (2018). Modeling flow and transport in fracture networks using graphs. Physical Review E, 97, 033304. https://doi.org/10.1103/PhysRevE.97.033304

Keegan-Treloar, R., Banks, E. W., Cartwright, I., Irvine, D. J., Webb, J. A., Werner, A. D., & Currell, M. J. (2024). Using major ions and stable isotopes to improve conceptualisation of a spring-aquifer system in the Galilee Basin, Australia. Hydrogeology Journal, 32, 1211–1228. https://doi.org/10.1007/s10040-024-02777-z

Kou, X., Zhao, Z., Duan, L., & Sun, Y. (2024). Hydrogeochemical behavior of shallow groundwater around Hancheng Mining Area, Guanzhong Basin, China. Water, 16, 660. https://doi.org/10.3390/w16050660

Kwicklis, E., Farnham, I., Hershey, R. L., Visser, A., & Hoaglund, J. (2021). Understanding long-term groundwater flow at Pahute Mesa and vicinity, Nevada National Security Site, USA, from naturally occurring geochemical and isotopic tracers. Hydrogeology Journal, 29, 2725–2749. https://doi.org/10.1007/s10040-021-02397-x

Lapworth, D. J., Krishan, G., MacDonald, A. M., & Rao, M. S. (2017). Groundwater quality in the alluvial aquifer system of northwest India: New evidence of the extent of anthropogenic and geogenic contamination. Science of the Total Environment, 599–600, 1433–1444. https://doi.org/10.1016/j.scitotenv.2017.04.223

Li, M., Liang, X., Xiao, C., Cao, Y., & Hu, S. (2019). Hydrochemical evolution of groundwater in a typical semi-arid groundwater storage basin using a zoning model. Water, 11, 1334. https://doi.org/10.3390/w11071334

Ling, X., Ma, J., Chen, P., Liu, C., & Horita, J. (2022). Isotope implications of groundwater recharge, residence time and hydrogeochemical evolution of the Longdong Loess Basin, Northwest China. Journal of Arid Land, 14, 34–55. https://doi.org/10.1007/s40333-022-0051-7

Lucas, L. L., & Unterweger, M. P. (2000). Comprehensive review and critical evaluation of the half-life of tritium. Journal of Research of the National Institute of Standards and Technology, 105(4), 541–549. https://doi.org/10.6028/jres.105.043

Ma, B., Jin, M., Liang, X., & Li, J. (2019). Application of environmental tracers for investigation of groundwater mean residence time and aquifer recharge in fault-influenced hydraulic drop alluvium aquifers. Hydrology and Earth System Sciences, 23, 427–446. https://doi.org/10.5194/hess-23-427-2019

Majoube, M. (1971). Oxygen-18 and deuterium fractionation between water and steam. Journal de Chimie Physique et de Physico-Chimie Biologique, 68(10), 1423–1436.

Malenica, L., Gotovac, H., Kamber, G., Simunovic, S., Allu, S., & Divic, V. (2018). Groundwater flow modeling in karst aquifers: Coupling 3D matrix and 1D conduit flow via control volume isogeometric analysis—Experimental verification with a 3D physical model. Water, 10(12), 1787. https://doi.org/10.3390/w10121787

Manu, E., De Lucia, M., & Kühn, M. (2023). Water–rock interactions driving groundwater composition in the Pra Basin (Ghana) identified by combinatorial inverse geochemical modelling. Minerals, 13, 899. https://doi.org/10.3390/min13070899

Meyer, R., Engesgaard, P., Hinsby, K., Piotrowski, J. A., & Sonnenborg, T. O. (2018). Estimation of effective porosity in large-scale groundwater models by combining particle tracking, auto-calibration and 14C dating. Hydrology and Earth System Sciences, 22, 4843–4865. https://doi.org/10.5194/hess-22-4843-2018

Murgulet, D., Cook, M., & Murgulet, V. (2016). Groundwater mixing between different aquifer types in a complex structural setting discerned by elemental and stable isotope geochemistry. Hydrological Processes, 30, 410–423. https://doi.org/10.1002/hyp.10589

Newman, C., Paschke, S., & Keith, G. (2021). Natural and anthropogenic geochemical tracers to investigate residence times and groundwater–surface-water interactions in an urban alluvial aquifer. Water, 13, 871. https://doi.org/10.3390/w13060871

Parkhurst, D. L., & Appelo, C. A. J. (2013). Description of input and examples for PHREEQC version 3: A computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations. U.S. Geological Survey Techniques and Methods, 6-A43. https://doi.org/10.3133/tm6A43

Purtschert, R., Love, A., Jiang, W., Lu, Z.-T., Yang, G.-M., Fulton, S., Wohling, D. L., Shand, P., Aeschbach, W., Broder, L., Müller, P., & Tosaki, Y. (2022). Residence times of groundwater along a flow path in the Great Artesian Basin determined by 81Kr, 36Cl and 4He: Implications for palaeo hydrogeology. Science of the Total Environment, 159886. https://doi.org/10.2139/ssrn.4204475

Qin, D., Zhao, Z., Han, D., Qian, H., & Liang, X. (2017). Isotope and hydrogeochemical constraints on groundwater recharge and mixing in the Xishan karst aquifer, Beijing. Hydrogeology Journal, 25, 1735–1750.

Richards, L., Kumari, R., Parashar, N., Kumar, A., Lu, C., Wilson, G., Lapworth, D., Niasar, V. J., Ghosh, A., Chakravorty, B., Krause, S., Polya, D., & Gooddy, D. (2022). Environmental tracers and groundwater residence time indicators reveal controls of arsenic accumulation rates beneath a rapidly developing urban area in Patna, India. Journal of Contaminant Hydrology, 249, 104043. https://doi.org/10.1016/j.jconhyd.2022.104043

Robinson, M. (2020). Assignments to sheaves of pseudometric spaces. Compositionality, 2. https://doi.org/10.32408/compositionality-2-2

Roubil, A., El Ouali, A., Bülbül, A., Lahrach, A., Mudry, J., Mamouch, Y., Essahlaoui, A., El Hmaidi, A., & El Ouali, A. (2022). Groundwater hydrochemical and isotopic evolution from High Atlas Jurassic limestones to Errachidia Cretaceous Basin (Southeastern Morocco). Water, 14, 1747. https://doi.org/10.3390/w14111747

Roy, A., Keesari, T., Mohokar, H., Pant, D., Sinha, U. K., & Mendhekar, G. N. (2020). Geochemical evolution of groundwater in hard-rock aquifers of South India using statistical and modelling techniques. Hydrological Sciences Journal, 65, 951–968. https://doi.org/10.1080/02626667.2019.1708914

Rädle, V., Kersting, A., Schmidt, M., Ringena, L., Robertz, J., Aeschbach, W., & Oberthaler, M. (2022). Multi-tracer groundwater dating in southern Oman using Bayesian modeling. Water Resources Research, 58. https://doi.org/10.1029/2021WR031776

Schiavo, M., Riva, M., Guadagnini, L., Zehe, E., & Guadagnini, A. (2022). Probabilistic identification of preferential groundwater networks. Journal of Hydrology, 610, 127906. https://doi.org/10.1016/j.jhydrol.2022.127906

Seltzer, A. M., Bekaert, D. V., Barry, P. H., Durkin, K., Mace, E. K., Aalseth, C. E., Zappala, J. C., Mueller, P., Jurgens, B. C., & Kulongoski, J. T. (2021). Groundwater residence time estimates obscured by anthropogenic carbonate. Science Advances, 7. https://doi.org/10.1126/sciadv.abf3503

Slimani, R., Guendouz, A., Trolard, F., Moulla, A. S., Hamdi-Aïssa, B., & Bourrié, G. (2017). Identification of dominant hydrogeochemical processes for groundwaters in the Algerian Sahara supported by inverse modeling of chemical and isotopic data. Hydrology and Earth System Sciences, 21, 1669–1691. https://doi.org/10.5194/hess-21-1669-2017

Smerdon, B. D., & Gardner, W. P. (2022). Characterizing groundwater flow paths in an undeveloped region through synoptic river sampling for environmental tracers. Hydrological Processes, 36(1), e14464. https://doi.org/10.1002/hyp.14464

Song, Y., Guo, J., Li, F., Wang, J., Ma, F., Wu, G., & Li, G. (2024). Investigation into factors controlling groundwater evolution in mining areas with an integrated approach. Heliyon, 10. https://doi.org/10.1016/j.heliyon.2024.e38860

Sun, Y., Wang, L., Zhang, Q., & Dong, Y. (2025). Recharge sources and flow pathways of karst groundwater in the Yuquan Mountain Spring catchment area, Beijing: A synthesis based on isotope, tracers, and geophysical evidence. Water, 17, 2292. https://doi.org/10.3390/w17152292

Taccari, M., Wang, H., Nuttall, J., Chen, X., & Jimack, P. (2024). Spatial-temporal graph neural networks for groundwater data. Scientific Reports, 14. https://doi.org/10.1038/s41598-024-75385-2

Teeple, A. P., Ging, P. B., Thomas, J. V., Wallace, D. S., & Payne, J. D. (2021). Hydrogeologic framework, geochemistry, groundwater-flow system, and aquifer hydraulic properties used in the development of a conceptual model of the Ogallala, Edwards-Trinity (High Plains), and Dockum aquifers in and near Gaines, Terry, and Yoakum Counties, Texas. U.S. Geological Survey Scientific Investigations Report. https://doi.org/10.3133/sir20215009

Thiros, N. E., Siirila-Woodburn, E. R., Dennedy-Frank, P. J., Williams, K. H., & Gardner, W. P. (2023). Constraining bedrock groundwater residence times in a mountain system with environmental tracer observations and Bayesian uncertainty quantification. Water Resources Research, 59, e2022WR033282. https://doi.org/10.1029/2022WR033282

Villegas, P., Paredes, V., Betancur, T., Taupin, J., & Toro, L. E. (2018). Groundwater evolution and mean water age inferred from hydrochemical and isotopic tracers in a tropical confined aquifer. Hydrological Processes, 32, 2158–2175. https://doi.org/10.1002/hyp.13160

Viswanathan, H. S., Hyman, J. D., Karra, S., O'Malley, D., Srinivasan, S., Hagberg, A., & Srinivasan, G. (2018). Advancing graph-based algorithms for predicting flow and transport in fractured rock. Water Resources Research, 54(9), 6085–6099. https://doi.org/10.1029/2017WR022368

Wang, M., Fang, J., Feng, F., Yao, T., Shan, Y., & Su, W. (2025). Geological and anthropogenic factors jointly influence hydrochemical interactions between groundwater and surface water in the middle and lower Yellow River. ACS Omega, 10, 48034–48050. https://doi.org/10.1021/acsomega.5c03965

Xiao, Y., Shao, J., Frape, S., Cui, Y., Dang, X., Wang, S., & Ji, Y. (2018). Groundwater origin, flow regime and geochemical evolution in arid endorheic watersheds: A case study from the Qaidam Basin, northwestern China. Hydrology and Earth System Sciences, 22, 4381–4400. https://doi.org/10.5194/hess-22-4381-2018

Yang, N., Wang, G., Shi, Z., Zhao, D., Jiang, W., Guo, L., Liao, F., & Zhou, P. (2018). Application of multiple approaches to investigate the hydrochemistry evolution of groundwater in an arid region: Nomhon, Northwestern China. Water, 10, 1667. https://doi.org/10.3390/w10111667

Zhang, X., Zhu, Y., Wang, J., Ju, L., Qian, Y., Ye, M., & Yang, J. (2022). GW-PINN: A deep learning algorithm for solving groundwater flow equations. Advances in Water Resources, 165, 104243. https://doi.org/10.1016/j.advwatres.2022.104243

Zhi, W., Li, L., Dong, W., Brown, W., Kaye, J., Steefel, C., & Williams, K. H. (2021). Non-uniqueness in identifying processes generating Ca-HCO3 water types: Carbonate dissolution, silicate weathering, and mixing. Geochimica et Cosmochimica Acta, 300, 120–138.

## 6. Conclusions

Hydrosheaf introduces a reproducible computational framework for inferring directed groundwater hydrochemical evolution networks from sparse multitracer observations. Its contribution is integrative: directed graph construction defines candidate pathways; sheaf-style consistency scoring tests those pathways jointly against hydraulic, chemical, isotope, and age evidence; a two-stage solver separates transport from residual chemical change; sparse inverse fitting with thermodynamic gates estimates parsimonious reaction sets; and Monte Carlo, bootstrap, and variance diagnostics quantify uncertainty, culminating in the process-stability index (PSI), which ranks reactions by their robustness to input perturbation. Because inverse models are run only on topology-supported edges, the framework reduces reliance on manual pathway selection and makes the remaining topology and reaction non-uniqueness explicit rather than claiming to remove it.

The results define the claim scope precisely. Independent no-prior topology inference recovered well-to-well connectivity with F1 = 0.62 (precision 0.49, recall 0.84) against 174 MODPATH reference edges; because false positives exceed true positives, inferred edges are screening-level connectivity hypotheses, and the prior-assisted F1 = 1.00 is an ingestion check rather than an inference result. Synthetic active-reaction recovery was moderate (R2 = 0.23, MAE = 0.37 mmol/L, false-activation rate 54%), limited by degeneracy of the reaction dictionary, with the rank, condition number, and collinearity diagnostics quantified in Sections 4.2 and 5.5. Synthetic network-age recovery was strong (log10 corr-R2 = 0.98, median absolute error 17.9 y), while public USGS parity was identifiability-gated, with 356 of 1272 fits identifiable and 92% of those within a factor of two. Field chemistry fits reached median R2 = 0.70 across 258 retained edges (Lower Anayari 0.53, Talensi 0.82) with median PSI = 0.97, but no independent process-truth labels exist for these sites.

The supported claim class is therefore the stability ranking of process-family hypotheses under the specified perturbation model, not reliable identification of the true reaction family. PSI distinguishes hypotheses that survive input uncertainty from noise-sensitive ones, and the outputs are best described as screening-level, uncertainty-ranked hypotheses for groundwater process discovery, not confirmed process truth. Hydrosheaf's transferable contribution is the reusable, auditable workflow itself, which converts sparse groundwater observations into testable network hypotheses and provides a transparent alternative to manual flow-path selection in data-limited aquifer investigations.

## Acknowledgement

The lead author acknowledges the GetFund Scholarship secretariat for the PhD tuition fee sponsorship that made this research possible. We would also like to thank the anonymous peer reviewers whose critical feedback and insightful suggestions significantly improved the clarity and scientific rigour of this article. Finally, we are grateful to the administrative and technical staff at the University of Ghana for their logistical assistance and support throughout the project.

## Funding

The lead author received no specific funding for this research project.

## Declaration of competing interest

The authors declare no known competing financial interests or personal relationships that could have influenced this work.

## CRediT authorship contribution statement

DA-W: Conceptualization; Methodology; Software; Writing – Original Draft (primary synthesis). EAA: Writing – Original Draft; Software (graph topology). DA: Writing – Review & Editing (academic validation and structural review). GA: Writing – Review & Editing. SG: Supervision; Writing – Review & Editing.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During preparation of this work, the authors used language-editing tools to improve clarity. The authors reviewed and edited the content and take full responsibility for the submitted manuscript.
