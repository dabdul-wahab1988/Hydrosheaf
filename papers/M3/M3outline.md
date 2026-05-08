<!-- Markdown companion generated from M3outline.docx for version-controlled manuscript planning. -->

# A) Testing Graph-Regularized Multi-Tracer Groundwater Age Inference with Tritium, Helium, SF6, CFCs and Carbon-14 in Hydrosheaf

# C) Outline with word counts

## Abstract — [300] words

- Study rationale, datasets, benchmarking design, principal findings and contribution — [300] words

## Keywords — [50] words

- Groundwater age; tritium; helium; sulfur hexafluoride; carbon-14; graph topology; Hydrosheaf; lumped-parameter models — [50] words

## 1. Introduction — [1,500] words

### 1.1 Background: groundwater residence time as a control on aquifer vulnerability and geochemical evolution — [300] words

### 1.2 Nuclear and atmospheric tracers in groundwater-age assessment — [300] words

### 1.3 Current reliance on well-wise lumped-parameter age models — [250] words

### 1.4 Limitation: weak use of spatial connectivity and aquifer network structure — [250] words

### 1.5 Rationale for testing graph regularization under controlled conditions — [200] words

### 1.6 Aim, objectives and research questions — [200] words

## 2. Materials and Methods — [4,100] words

### 2.1 Study design and benchmarking philosophy — [300] words

### 2.2 Public USGS groundwater-age datasets — [450] words

#### 2.2.1 National public-supply aquifer age dataset, 2004-2017 — [150] words

#### 2.2.2 Western Principal Aquifers groundwater-age dataset, 2004-2018 — [150] words

#### 2.2.3 Mississippi River Valley alluvial aquifer tracer-age dataset, 2018-2020 — [150] words

### 2.3 Tracer variables and hydrogeologic attributes — [350] words

#### 2.3.1 Tritium and tritiogenic helium-3 — [100] words

#### 2.3.2 SF6 and dissolved-gas corrections — [100] words

#### 2.3.3 Carbon-14 and radiogenic helium-4 indicators of older groundwater — [100] words

#### 2.3.4 Ancillary well, aquifer, depth and chemistry variables — [50] words

### 2.4 Data harmonisation and quality-control protocol — [450] words

#### 2.4.1 Unit standardisation and censoring treatment — [100] words

#### 2.4.2 Missing-data rules and tracer-availability classes — [100] words

#### 2.4.3 Screening for SF6 contamination and excess-air uncertainty — [100] words

#### 2.4.4 Screening for helium terrigenic contribution — [75] words

#### 2.4.5 Carbon-14 interpretive constraints and dead-carbon uncertainty — [75] words

### 2.5 Baseline groundwater-age models — [500] words

#### 2.5.1 Piston-flow model — [100] words

#### 2.5.2 Exponential mixing model — [100] words

#### 2.5.3 Dispersion model — [100] words

#### 2.5.4 Binary young-old mixture model — [100] words

#### 2.5.5 USGS-reported age estimates as benchmark reference values — [100] words

### 2.6 Hydrosheaf multi-tracer age-inference and graph-regularization framework — [750] words

#### 2.6.1 Graph construction from spatial, depth and aquifer constraints — [150] words

#### 2.6.2 Directional connectivity and probable recharge-to-discharge ordering — [150] words

#### 2.6.3 Sheaf-style consistency constraints across connected wells — [150] words

#### 2.6.4 Graph-prior settings, randomized negative controls and regularization diagnostics — [150] words

#### 2.6.5 Separation of tracer evidence from graph prior influence — [150] words

### 2.7 Benchmark scenarios — [450] words

#### 2.7.1 Single-well tracer-only inference — [100] words

#### 2.7.2 Multi-tracer well-wise inference — [100] words

#### 2.7.3 Single-node versus graph-regularized age inference — [100] words

#### 2.7.4 Graph-constrained inference with hydrochemical consistency checks — [100] words

#### 2.7.5 Data-limited emulation by tracer removal experiments — [50] words

### 2.8 Performance metrics and statistical evaluation — [450] words

#### 2.8.1 Bias, RMSE, MAE and rank correlation against reference ages — [100] words

#### 2.8.2 Classification skill for modern, mixed and premodern groundwater — [100] words

#### 2.8.3 Posterior uncertainty width and calibration — [100] words

#### 2.8.4 Spatial smoothness versus over-regularisation diagnostics — [75] words

#### 2.8.5 Cross-validation by aquifer, region and tracer availability — [75] words

### 2.9 Reproducibility, software implementation and data availability — [450] words

#### 2.9.1 Hydrosheaf configuration and computational workflow — [150] words

#### 2.9.2 Public data retrieval and processing scripts — [100] words

#### 2.9.3 Version control, parameter logging and random seeds — [100] words

#### 2.9.4 Limitations of public-data benchmarking — [100] words

## 3. Results — [3,000] words

### 3.1 Dataset coverage and tracer availability — [400] words

#### 3.1.1 Spatial distribution of wells and aquifer systems — [100] words

#### 3.1.2 Distribution of tracer combinations — [100] words

#### 3.1.3 Age-class distribution across aquifers — [100] words

#### 3.1.4 Missingness and censoring structure — [100] words

### 3.2 Baseline age-model performance — [450] words

#### 3.2.1 Tritium-only classification performance — [100] words

#### 3.2.2 3H/3He apparent-age performance — [100] words

#### 3.2.3 SF6-derived modern-age performance — [100] words

#### 3.2.4 Carbon-14 older-groundwater performance — [100] words

#### 3.2.5 Failure modes under tracer disagreement — [50] words

### 3.3 Graph topology and inferred hydrogeologic connectivity — [450] words

#### 3.3.1 Edge-density and connected-component structure — [100] words

#### 3.3.2 Depth-ordered and aquifer-constrained connectivity patterns — [100] words

#### 3.3.3 Recharge-to-discharge age gradients along graph paths — [100] words

#### 3.3.4 Identification of topological inconsistencies — [100] words

#### 3.3.5 Sensitivity to graph-prior strength — [50] words

### 3.4 When graph regularization improves or degrades age inference — [550] words

#### 3.4.1 Improvement cases relative to tracer-only models — [125] words

#### 3.4.2 Degradation cases under wrong topology, tracer disagreement or over-regularization — [125] words

#### 3.4.3 Performance by aquifer type and hydrogeologic setting — [125] words

#### 3.4.4 Posterior uncertainty reduction and calibration — [125] words

#### 3.4.5 Cases where graph constraints degrade performance — [50] words

### 3.5 Tracer-specific diagnostic results — [450] words

#### 3.5.1 Tritium and helium sensitivity to modern recharge — [100] words

#### 3.5.2 SF6 vulnerability to contamination and excess-air assumptions — [100] words

#### 3.5.3 Carbon-14 correction uncertainty in older groundwater — [100] words

#### 3.5.4 Radiogenic helium-4 as qualitative old-water support — [75] words

#### 3.5.5 Multi-tracer conflict resolution by network constraints — [75] words

### 3.6 Data-limited emulation experiments — [450] words

#### 3.6.1 Loss of performance after removing helium — [100] words

#### 3.6.2 Loss of performance after removing SF6 — [100] words

#### 3.6.3 Loss of performance after removing carbon-14 — [100] words

#### 3.6.4 Minimum tracer set required for defensible age classification — [100] words

#### 3.6.5 Implications for semi-arid and low-resource aquifers — [50] words

### 3.7 Summary of key results — [250] words

## 4. Discussion — [2,400] words

### 4.1 What network enhancement adds beyond conventional LPMs — [350] words

### 4.2 Hydrogeologic meaning of graph-constrained residence-time structure — [300] words

### 4.3 Contribution to nuclear-isotope groundwater dating — [300] words

### 4.4 Implications for data-limited aquifers — [350] words

### 4.5 Relevance to Hydrosheaf development and PhD thesis integration — [250] words

### 4.6 Methodological risks and interpretive boundaries — [350] words

#### 4.6.1 Risk of imposing false connectivity — [100] words

#### 4.6.2 Risk of masking tracer contamination — [100] words

#### 4.6.3 Risk of overconfident posterior ages — [75] words

#### 4.6.4 Limits of transferring USGS benchmarks to Ghanaian aquifers — [75] words

### 4.7 Comparison with existing tracer-age modelling approaches — [250] words

### 4.8 Recommendations for future sampling and model validation — [250] words

## 5. Conclusions — [500] words

### 5.1 Main findings — [150] words

### 5.2 Methodological contribution — [125] words

### 5.3 Practical contribution for data-limited aquifers — [125] words

### 5.4 Final thesis-relevant conclusion — [100] words

## 6. Data, Code and Reproducibility Statement — [150] words

- Public USGS data, Hydrosheaf scripts, parameter files and archived outputs — [150] words
- Total estimated length: [12,000] words, excluding references and supplementary material.
- D) Proposed Tables and Figures
- • Table 1: Public USGS groundwater-age datasets used in the benchmark — Purpose — Documents dataset source, period, aquifer coverage, sample size, tracer coverage and DOI/URL; key variables: dataset name, years, wells, aquifers, 3H, 3He, SF6, 14C, 4He, reference age fields.
- • Table 2: Tracer interpretation framework — Purpose — Summarises tracer age ranges, assumptions and failure modes; key variables: tracer, age range, required corrections, contamination risk, interpretive class.
- • Table 3: Benchmark model scenarios — Purpose — Defines baseline and Hydrosheaf model variants; key variables: model ID, tracer inputs, graph prior, chemistry constraint, validation metric.
- • Table 4: Model-performance summary — Purpose — Compares single-node, graph-regularized and randomized-graph controls; key variables: RMSE, MAE, bias, Spearman correlation, classification accuracy, posterior interval width and uncertainty coverage.
- • Table 5: Cross-validation results by aquifer system — Purpose — Tests transferability across hydrogeologic settings; key variables: aquifer, number of wells, tracer completeness, baseline error, Hydrosheaf error, improvement percentage.
- • Table 6: Data-limited emulation results — Purpose — Quantifies performance loss when selected tracers are removed; key variables: removed tracer, retained tracer set, age-class accuracy, uncertainty increase, failure mode.
- • Figure 1: Spatial distribution of USGS benchmark wells — Map — Shows geographic coverage and tracer availability across principal aquifers.
- • Figure 2: Hydrosheaf network-age inference workflow — Conceptual workflow diagram — Demonstrates data ingestion, QC, graph construction, tracer-age modelling, posterior inference and benchmarking.
- • Figure 3: Tracer-age distributions by aquifer and tracer combination — Violin/box plots — Demonstrates how age structure varies across aquifer systems and tracer availability classes.
- • Figure 4: Single-node versus graph-regularized age estimates — Scatter plot with 1:1 line — Shows agreement, bias, improvement cases and degradation cases relative to USGS/reference age estimates.
- • Figure 5: Graph topology and residence-time gradients — Network plot/map — Demonstrates whether inferred graph paths preserve plausible young-to-old groundwater ordering.
- • Figure 6: Model-error decomposition by tracer type — Bar or heatmap — Shows which tracers and tracer combinations produce the largest uncertainty or disagreement.
- • Figure 7: Data-limited tracer-removal experiment — Line plot or grouped bars — Demonstrates how age-inference performance changes when 3He, SF6 or 14C is unavailable.
- • Figure 8: Posterior uncertainty calibration — Reliability plot — Demonstrates whether Hydrosheaf uncertainty intervals are statistically credible across age classes.
- Useful public data anchors: USGS national public-supply groundwater-age dataset [2004-2017](https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004), USGS Western Principal Aquifers groundwater-age dataset [2004-2018](https://www.usgs.gov/data/data-groundwater-age-western-principal-aquifers-2004-2018), and USGS Mississippi River Valley alluvial aquifer tracer-age dataset [2018-2020](https://data.usgs.gov/datacatalog/data/USGS%3A60900920d34e93746a710344).
- Implementation alignment update for the revised M3 manuscript
- The outline should now cite the committed Hydrosheaf nuclear modules and M3/m3_age_benchmark package. Implemented capabilities include dissolved-gas correction, tracer-input handling, multi-tracer age estimates, joint LPM fitting, tracer-disagreement diagnostics, LPM identifiability diagnostics, tracer-removal sensitivity and graph age-coherence auditing.
- The core M3 test is no longer whether graph information is generally beneficial. It is: when does graph regularization improve age inference, when does it degrade inference, and which tracer combinations remain robust?
- Required wording guardrail: claim graph-regularization benefit only where a benchmark row shows lower error than single-node inference and randomized negative-control graphs do not show comparable gains.
