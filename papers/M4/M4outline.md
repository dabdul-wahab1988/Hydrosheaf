<!-- Markdown companion generated from M4outline.docx for version-controlled manuscript planning. -->

# A) Validating Reduced-Order Hydrosheaf Groundwater Topology Against MODPATH Advective Connectivity Under Controlled Benchmarks

# C) Outline with word counts

## 1. Introduction — [1,500] words

### 1.1 Groundwater connectivity as a controlling uncertainty in aquifer diagnosis — [250] words

### 1.2 MODFLOW/MODPATH particle tracking as a numerical reference for advective connectivity — [250] words

### 1.3 Graph-constrained inference as a reduced-order alternative for data-limited aquifers — [250] words

### 1.4 Problem: unvalidated graph topology can impose false flow pathways — [250] words

### 1.5 Research gap: limited benchmarking of inferred groundwater graphs against particle-tracking outputs — [250] words

### 1.6 Aim, objectives, hypotheses and article contribution — [250] words

## 2. Materials and Methods — [4,200] words

### 2.1 Study design and benchmarking logic — [300] words

### 2.2 Controlled MODPATH benchmark design and source-archive scope — [650] words

#### 2.2.1 Primary curated MODPATH reference archive and inclusion criteria — [125] words

#### 2.2.2 Great Miami/Wright-Patterson archive as future external-validation candidate — [125] words

#### 2.2.3 Long Island MODFLOW 6/MODPATH 7 archive as future external-validation candidate — [125] words

#### 2.2.4 Coastal Connecticut/Long Island Sound archive as future external-validation candidate — [125] words

#### 2.2.5 Optional generic MODPATH models for controlled sensitivity testing — [150] words

### 2.3 Data acquisition, model archive screening and inclusion criteria — [350] words

### 2.4 MODFLOW/MODPATH reference processing — [550] words

#### 2.4.1 Extraction of model grids, heads, boundaries and hydraulic stresses — [125] words

#### 2.4.2 Extraction of MODPATH pathlines, endpoints and travel times — [125] words

#### 2.4.3 Conversion of particle trajectories into directed connectivity graphs — [150] words

#### 2.4.4 Treatment of transient versus steady-state simulations — [75] words

#### 2.4.5 Quality control of disconnected, terminated or weakly constrained particles — [75] words

### 2.5 Hydrosheaf graph-constrained connectivity inference — [650] words

#### 2.5.1 Node definition from wells, cells, discharge zones and particle release locations — [125] words

#### 2.5.2 Candidate-edge construction from distance, head, elevation, depth and aquifer membership — [150] words

#### 2.5.3 Directional weighting and probabilistic edge confidence — [125] words

#### 2.5.4 MODPATH-informed graph priors used inside Hydrosheaf — [125] words

#### 2.5.5 Strict separation of independent validation from prior-informed Hydrosheaf runs — [125] words

### 2.6 Benchmark scenarios — [550] words

#### 2.6.1 Spatial-only graph inference — [100] words

#### 2.6.2 Head-gradient and depth-constrained graph inference — [100] words

#### 2.6.3 Hydrostratigraphic graph inference — [100] words

#### 2.6.4 Sparse-data emulation with reduced input variables — [100] words

#### 2.6.5 Physics-prior Hydrosheaf run using MODPATH-derived priors — [100] words

#### 2.6.6 Negative-control scenarios using perturbed or randomised topology — [50] words

### 2.7 Graph-to-particle comparison metrics — [650] words

#### 2.7.1 Directed edge precision, recall and F1 score — [125] words

#### 2.7.2 Pathline-to-graph path similarity — [125] words

#### 2.7.3 Endpoint capture-zone agreement — [125] words

#### 2.7.4 Travel-time rank correlation and edge-weight consistency — [100] words

#### 2.7.5 Topological distance, betweenness and source-sink preservation — [100] words

#### 2.7.6 False-positive and false-negative connectivity diagnostics — [75] words

### 2.8 Uncertainty, sensitivity and robustness analysis — [350] words

### 2.9 Software implementation and reproducibility protocol — [400] words

### 2.10 Statistical analysis and interpretation thresholds — [300] words

## 3. Results — [3,000] words

### 3.1 Characteristics of selected MODFLOW/MODPATH benchmark models — [350] words

### 3.2 Reference particle-tracking connectivity structures — [400] words

### 3.3 Hydrosheaf graph structures under alternative inference scenarios — [450] words

### 3.4 Independent graph inference compared with MODPATH reference connectivity — [600] words

#### 3.4.1 Directed edge-level agreement — [150] words

#### 3.4.2 Source-to-discharge path agreement — [150] words

#### 3.4.3 Capture-zone and endpoint agreement — [150] words

#### 3.4.4 Travel-time and pathway-order agreement — [150] words

### 3.5 MODPATH-informed prior-mode results and why they are not independent validation — [400] words

### 3.6 Failure modes and hydrogeologic causes of mismatch — [400] words

### 3.7 Computational performance relative to full particle tracking — [250] words

### 3.8 Summary of benchmark outcomes — [150] words

## 4. Discussion — [2,300] words

### 4.1 Main finding: when graph-constrained connectivity reproduces particle-tracking structure — [300] words

### 4.2 What MODPATH benchmarking demonstrates, and what it cannot prove — [300] words

### 4.3 Contribution to Hydrosheaf as a reduced-order hydrogeologic inference framework — [300] words

### 4.4 Implications for data-limited semi-arid aquifers where full numerical models are unavailable — [350] words

### 4.5 Role of MODFLOW/MODPATH validation in a Nuclear Earth Science PhD thesis — [250] words

### 4.6 Methodological risks — [400] words

#### 4.6.1 False connectivity from sparse observations — [100] words

#### 4.6.2 Overfitting to numerical-model assumptions — [100] words

#### 4.6.3 Scale mismatch between wells, grid cells and hydrostratigraphic units — [100] words

#### 4.6.4 Limits of transferring USGS benchmarks to Ghanaian aquifers — [100] words

### 4.7 Recommendations for future Hydrosheaf-MODFLOW/MODPATH integration — [250] words

### 4.8 Broader contribution to groundwater model validation and decision support — [150] words

## 5. Conclusions — [500] words

### 5.1 Principal findings — [150] words

### 5.2 Methodological contribution — [125] words

### 5.3 Practical relevance for data-limited aquifer studies — [125] words

### 5.4 Final thesis-relevant conclusion — [100] words

## 6. Data and Code Availability — [150] words

- Public USGS model archives, Hydrosheaf scripts, processed graph outputs and reproducibility files — [150] words

## 7. Declarations — [350] words

### 7.1 Author contributions — [75] words

### 7.2 Funding and acknowledgements — [75] words

### 7.3 Competing interests — [50] words

### 7.4 Software and version statement — [75] words

### 7.5 Limitations of benchmark reuse — [75] words

- Total estimated length: [12,000] words, excluding references and supplementary material.
- D) Proposed Tables and Figures
- • Table 1: Controlled MODPATH benchmark source and validation scope — Purpose — Summarises the curated reference archive, hydrogeologic setting, MODFLOW/MODPATH versions, particle-tracking outputs, inclusion criteria and validation limits.
- • Table 2: Benchmark inclusion and exclusion criteria — Purpose — Defines which public models are suitable for Hydrosheaf validation; key variables: model availability, pathline availability, endpoint availability, spatial metadata, heads, aquifer units, travel times.
- • Table 3: Hydrosheaf graph-inference scenarios — Purpose — Documents spatial-only, head-constrained, hydrostratigraphic, data-limited and physics-prior graph configurations; key variables: node type, edge rule, direction rule, prior strength, threshold.
- • Table 4: Graph-to-MODPATH validation metrics — Purpose — Defines quantitative comparison metrics; key variables: precision, recall, F1 score, path similarity, capture-zone overlap, travel-time rank correlation, false-positive rate.
- • Table 5: Benchmark performance by model archive — Purpose — Reports Hydrosheaf agreement with MODPATH for each public benchmark; key variables: dataset, number of nodes, number of MODPATH paths, edge agreement, path agreement, endpoint agreement.
- • Table 6: Failure-mode classification — Purpose — Identifies why inferred graphs disagree with particle tracking; key variables: mismatch type, hydrogeologic cause, data cause, model-scale cause, recommended correction.
- • Table 7: Reproducibility and software configuration — Purpose — Records Hydrosheaf, FloPy, MODFLOW, MODPATH and Python versions; key variables: package version, executable version, random seed, input archive, output folder.
- • Figure 1: Benchmarking workflow — Flow diagram — Shows conversion from MODFLOW/MODPATH archives to reference graphs, Hydrosheaf-inferred graphs and validation metrics.
- • Figure 2: Controlled benchmark domain and future external-validation archives — Multi-panel map — Shows the controlled benchmark domain and identifies additional public archives only as future external-validation targets.
- • Figure 3: MODPATH pathlines converted to directed graph edges — Network/pathline overlay — Demonstrates how particle trajectories become reference connectivity structures.
- • Figure 4: Hydrosheaf versus MODPATH adjacency agreement — Heatmap/confusion matrix — Shows true-positive, false-positive and false-negative directed connections.
- • Figure 5: Capture-zone agreement — Map or polygon-overlap plot — Demonstrates overlap between Hydrosheaf-inferred source areas and MODPATH-derived contributing areas.
- • Figure 6: Path-similarity performance across inference scenarios — Box plot — Compares spatial-only, head-constrained, hydrostratigraphic and physics-prior graph inference.
- • Figure 7: Travel-time rank agreement — Scatter plot — Shows relation between MODPATH travel time and Hydrosheaf edge/path weighting.
- • Figure 8: Sensitivity to graph-prior strength — Line plot — Demonstrates how connectivity accuracy changes as the graph prior is weakened or strengthened.
- • Figure 9: Failure-mode examples — Multi-panel diagnostic figure — Shows representative mismatches caused by pumping, boundary effects, anisotropy, vertical gradients and sparse nodes.
- Source anchors used for this outline include USGS documentation for [MODFLOW 6](https://www.usgs.gov/software/modflow-6-usgs-modular-hydrologic-model), [MODPATH](https://www.usgs.gov/software/modpath-a-particle-tracking-model-modflow), the curated benchmark archive used in Hydrosheaf, and additional public model archives treated as future external-validation candidates until fully processed.
- Implementation alignment update for the revised M4 manuscript
- The outline should now cite hydrosheaf.validation.topology and M4/m4_topology_benchmark. Implemented capabilities include directed-edge confusion diagnostics, false-positive and false-negative reporting, false-positive rate, false-negative rate, scale-mismatch diagnostics, independent validation against MODPATH reference edges and MODPATH-informed prior construction.
- The manuscript must separate two modes: independent graph inference compared against MODPATH, and MODPATH-informed graph priors used inside Hydrosheaf. These are scientifically different and should not be merged in results or claims.
- Required wording guardrail: M4 validates reduced-order topology reproduction under controlled benchmark conditions. It should not claim broad external-archive benchmarking unless those archives have been fully processed and reported.
