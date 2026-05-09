<!-- Markdown companion generated from M2outline.docx for version-controlled manuscript planning. -->

# A) Hydrosheaf: An Integrated Validation Framework for Coupled Groundwater Residence-Time, Topology and Inverse Hydrogeochemical Modelling

# C) Outline with word counts

## 1. Introduction — [1,500] words

### 1.1 Groundwater evolution as a coupled age-transport-reaction problem — [300] words

### 1.2 Limitations of point-based hydrochemical and residence-time interpretation — [300] words

### 1.3 Flow-path uncertainty and non-uniqueness in inverse geochemical modelling — [300] words

### 1.4 Need for graph-constrained coupling of residence time and geochemical reactions — [300] words

### 1.5 Objectives, novelty and structure of the article — [300] words

## 2. Materials and Methods — [4,000] words

### 2.1 Hydrosheaf framework overview — [450] words

#### 2.1.1 Software architecture and main workflow — [200] words

#### 2.1.2 Conceptual model: nodes, edges, tracers and reactions — [250] words

### 2.2 Data model and preprocessing — [450] words

#### 2.2.1 Required sample, chemistry, isotope and coordinate fields — [150] words

#### 2.2.2 Unit harmonisation and ion-order configuration — [150] words

#### 2.2.3 Missing-value, detection-limit and quality-control handling — [150] words

### 2.3 Graph construction and edge inference — [650] words

#### 2.3.1 Spatial, elevation and head-proxy edge generation — [200] words

#### 2.3.2 Primary, lateral and 3D edge structures — [150] words

#### 2.3.3 Physics priors from MODPATH and external models — [150] words

#### 2.3.4 Sheaf-style topology refinement using chemical and age consistency — [150] words

### 2.4 Residence-time and nuclear-tracer modelling — [700] words

#### 2.4.1 Nuclide definitions and tracer input functions — [150] words

#### 2.4.2 Lumped parameter models: piston, exponential and dispersion assumptions — [200] words

#### 2.4.3 Single-node tracer-age inversion — [150] words

#### 2.4.4 Network-enhanced Bayesian age inference — [200] words

### 2.5 Inverse hydrogeochemical modelling — [850] words

#### 2.5.1 Transport fitting for evaporation and mixing — [150] words

#### 2.5.2 Conservative-tracer weighting — [100] words

#### 2.5.3 Sparse reaction fitting and L1 regularisation — [200] words

#### 2.5.4 Mineral/reaction dictionary and process terms — [150] words

#### 2.5.5 PHREEQC saturation-index constraints — [150] words

#### 2.5.6 Objective functions, residuals and model selection — [100] words

### 2.6 Uncertainty, validation and outputs — [600] words

#### 2.6.1 Bayesian, bootstrap and Monte Carlo uncertainty modes — [150] words

#### 2.6.2 PHREEQC kinetic forward validation — [150] words

#### 2.6.3 MODFLOW/MODPATH topology-validation design and prior-mode separation — [150] words

#### 2.6.4 Tables, diagnostic plots and reproducible exports — [150] words

### 2.7 Implementation, availability and reproducibility — [300] words

#### 2.7.1 Package structure, tests and command-line interface — [150] words

#### 2.7.2 Reproducible analysis scripts and data-access assumptions — [150] words

## 3. Results — [3,000] words

### 3.1 Framework implementation results — [600] words

#### 3.1.1 Verified software modules and workflow execution — [200] words

#### 3.1.2 Graph-edge and physics-prior outputs — [150] words

#### 3.1.3 Residence-time output fields and uncertainty diagnostics — [150] words

#### 3.1.4 Hydrogeochemical reaction output fields — [100] words

### 3.2 Synthetic benchmark recovery — [700] words

#### 3.2.1 Recovery of known transport processes — [150] words

#### 3.2.2 Recovery of known reaction extents — [200] words

#### 3.2.3 Sensitivity to missing data and ion weighting — [150] words

#### 3.2.4 Robustness to edge-density and topology uncertainty — [200] words

### 3.3 Nuclear-tracer and residence-time validation — [600] words

#### 3.3.1 Agreement with public LPM/TracerLPM-style age estimates — [200] words

#### 3.3.2 Behaviour under mixed-age and ambiguous tracer conditions — [200] words

#### 3.3.3 Network-age constraints and downstream ageing consistency — [200] words

### 3.4 MODFLOW/MODPATH topology validation — [450] words

#### 3.4.1 Conversion of MODPATH endpoints to graph priors — [150] words

#### 3.4.2 Agreement between reduced-order graph edges and MODPATH particle-tracking connectivity — [200] words

#### 3.4.3 Failure modes and topology uncertainty — [100] words

### 3.5 PHREEQC-constrained hydrogeochemical validation — [450] words

#### 3.5.1 Saturation-index constraints and reaction feasibility — [150] words

#### 3.5.2 Effect of thermodynamic constraints on sparse reaction solutions — [150] words

#### 3.5.3 Forward kinetic validation diagnostics — [150] words

### 3.6 Demonstration on a data-limited aquifer dataset — [200] words

#### 3.6.1 End-to-end workflow demonstration — [100] words

#### 3.6.2 Diagnostic value for process-pathway interpretation — [100] words

## 4. Discussion — [2,500] words

### 4.1 Scientific contribution of Hydrosheaf — [500] words

#### 4.1.1 Coupling residence time, topology and inverse reactions — [200] words

#### 4.1.2 Moving from point interpretation to process networks — [150] words

#### 4.1.3 Contribution to Nuclear Earth Science and isotope hydrology — [150] words

### 4.2 Comparison with existing approaches — [600] words

#### 4.2.1 PHREEQC/NETPATH and manual inverse modelling — [150] words

#### 4.2.2 MODFLOW/MODPATH and physical particle tracking — [150] words

#### 4.2.3 VISHMOD, GraphFlow and virtual-sampling approaches — [150] words

#### 4.2.4 PINNs, GNNs and data-driven surrogates — [150] words

### 4.3 Interpretation of validation outcomes — [450] words

#### 4.3.1 What synthetic benchmarks prove — [150] words

#### 4.3.2 What public residence-time datasets prove — [150] words

#### 4.3.3 What MODPATH and PHREEQC validations prove — [150] words

### 4.4 Practical relevance for data-limited aquifers — [400] words

#### 4.4.1 Screening-level diagnosis where calibrated flow models are absent — [150] words

#### 4.4.2 Prioritising isotope sampling and monitoring nodes — [150] words

#### 4.4.3 Water-quality management implications — [100] words

### 4.5 Limitations — [350] words

#### 4.5.1 Dependence on input chemistry, tracer availability and graph assumptions — [120] words

#### 4.5.2 Non-uniqueness and identifiability under sparse data — [120] words

#### 4.5.3 Limits of benchmark transferability — [110] words

### 4.6 Future development — [200] words

#### 4.6.1 Full regional MODFLOW/MODPATH coupling — [80] words

#### 4.6.2 Expanded nuclear-tracer and isotope datasets — [70] words

#### 4.6.3 Community benchmarking and open reproducibility — [50] words

## 5. Conclusions — [700] words

### 5.1 Main methodological conclusions — [200] words

### 5.2 Validation conclusions — [200] words

### 5.3 Implications for groundwater residence-time and hydrogeochemical diagnosis — [200] words

### 5.4 Final contribution statement — [100] words

## 6. Data and Code Availability — [200] words

### 6.1 Software repository, version and licence — [80] words

### 6.2 Public datasets and benchmark resources — [80] words

### 6.3 Reproducibility package — [40] words

## 7. Author Contributions, Competing Interests and Acknowledgements — [100] words

- Total estimated length — [12,000] words
- D) Proposed Tables and Figures
- • Table 1: Hydrosheaf module architecture — Purpose: summarise the computational components of the framework — Key variables: module, function, input, output, validation role.
- • Table 2: Required and optional data fields — Purpose: clarify minimum and enhanced data requirements — Key variables: site ID, coordinates, ions, stable isotopes, nuclear tracers, hydraulic heads, MODPATH priors, PHREEQC inputs.
- • Table 3: Residence-time modelling options in Hydrosheaf — Purpose: compare nuclear and non-nuclear residence-time methods — Key variables: tracer, model type, assumptions, age range, uncertainty output, limitation.
- • Table 4: Inverse hydrogeochemical reaction dictionary — Purpose: document process terms available for sparse fitting — Key variables: reaction label, stoichiometry, process interpretation, PHREEQC constraint, relevance.
- • Table 5: Validation design and benchmark datasets — Purpose: show how the framework is evaluated — Key variables: benchmark type, dataset source, target variable, validation metric, expected evidence.
- • Table 6: Comparison of Hydrosheaf with existing approaches — Purpose: position novelty against established methods — Key variables: method, topology handling, residence-time support, thermodynamic rigour, data requirement, limitation.
- • Figure 1: Hydrosheaf conceptual architecture — Plot type: workflow diagram — What it demonstrates: integration of samples, graph topology, residence-time inference, inverse reaction fitting and validation outputs.
- • Figure 2: Graph-constrained groundwater process network — Plot type: schematic directed graph — What it demonstrates: nodes, primary edges, lateral edges, physics priors, age constraints and reaction extents.
- • Figure 3: Residence-time inference workflow — Plot type: flowchart — What it demonstrates: tracer input histories, LPM selection, single-node inversion, network Bayesian inference and uncertainty output.
- • Figure 4: Sparse inverse hydrogeochemical fitting workflow — Plot type: algorithm diagram — What it demonstrates: transport fitting, conservative-tracer weighting, reaction residual fitting, PHREEQC constraints and objective scoring.
- • Figure 5: Benchmark validation design — Plot type: multi-panel validation schematic — What it demonstrates: synthetic recovery, public tracer-age validation, MODPATH topology comparison and PHREEQC forward validation.
- • Figure 6: 3D Process-network output — Plot type: 3D network map or edge-coloured graph — What it demonstrates: Spatially connected reaction pathways in 3D, dominant processes and uncertainty indicators mapped to depth.
- • Figure 7: Sensitivity and uncertainty summary — Plot type: tornado plot or uncertainty cascade — What it demonstrates: effects of graph parameters, tracer assumptions, reaction dictionary and input uncertainty on model outputs.
- Implementation alignment update for the revised M2 manuscript
- The outline should now be written around the committed curated benchmark package in M2/m2_benchmark. That package contains synthetic recovery tests, public tracer-age validation, MODPATH-to-graph topology validation, PHREEQC forward diagnostics, DGMETA dissolved-gas validation and a Northern Ghana data-limited demonstration.
- Raw public input archives and large local source files are intentionally excluded from Git. Reproducibility should be described through the committed scripts, generated result tables, figures, source manifests and README files.
- Required wording guardrail: present M2 as an integrated validation package with tier-specific evidence. Do not collapse all tiers into a single claim of full age, topology and reactive-transport equivalence.
