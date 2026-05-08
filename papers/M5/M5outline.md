<!-- Markdown companion generated from M5outline.docx for version-controlled manuscript planning. -->

# A) Sparse Linear Inverse Reaction Modelling Screened by PHREEQC Thermodynamic and Forward-Validation Diagnostics in Hydrosheaf

# C) Outline with word counts

## 1. Introduction — [1,500] words

### 1.1 Groundwater evolution as an inverse hydrogeochemical problem — [250] words

### 1.2 Non-uniqueness in conventional inverse mass-balance modelling — [250] words

### 1.3 Sparse reaction fitting as a parsimonious alternative — [250] words

### 1.4 Need for thermodynamic screening of inferred reactions — [250] words

### 1.5 Research gap: limited validation of sparse inverse reaction models against PHREEQC constraints — [250] words

### 1.6 Aim, objectives, hypotheses and article contribution — [250] words

## 2. Materials and Methods — [4,200] words

### 2.1 Study design and validation philosophy — [300] words

### 2.2 Hydrosheaf sparse inverse reaction framework — [650] words

#### 2.2.1 Edge-wise concentration-change formulation — [125] words

#### 2.2.2 Separation of transport effects from reaction residuals — [125] words

#### 2.2.3 Sparse reaction fitting using L1-regularised coordinate descent — [150] words

#### 2.2.4 Reaction extents, residual vectors and goodness-of-fit metrics — [125] words

#### 2.2.5 Handling missing ions and reduced ion-order configurations — [125] words

### 2.3 Reaction dictionary and hydrogeochemical process representation — [500] words

#### 2.3.1 Carbonate, evaporite, silicate and fluoride-bearing mineral reactions — [125] words

#### 2.3.2 Ion-exchange reactions — [100] words

#### 2.3.3 Nitrate input and denitrification terms — [100] words

#### 2.3.4 Signed versus non-negative reaction extents — [75] words

#### 2.3.5 Stoichiometric limitations of linear reaction dictionaries — [100] words

### 2.4 PHREEQC thermodynamic screening and saturation-index bound system — [650] words

#### 2.4.1 PHREEQC speciation and saturation-index calculation — [125] words

#### 2.4.2 Translation of saturation indices into reaction bounds — [150] words

#### 2.4.3 Thermodynamic admissibility of dissolution and precipitation — [125] words

#### 2.4.4 Treatment of unavailable PHREEQC results — [100] words

#### 2.4.5 Redox and nitrate-related override constraints — [75] words

#### 2.4.6 Difference between equilibrium screening and kinetic validation — [75] words

### 2.5 Benchmark datasets and test cases — [650] words

#### 2.5.1 PHREEQC-generated synthetic reaction-path benchmarks — [150] words

#### 2.5.2 Controlled mineral-mixture experiments in silico — [125] words

#### 2.5.3 Public groundwater-chemistry datasets with major ions — [125] words

#### 2.5.4 Northern Ghana hydrochemical dataset as data-limited semi-arid case application — [125] words

#### 2.5.5 Data-limited experiments using ion removal and noise addition — [125] words

### 2.6 Validation scenarios — [600] words

#### 2.6.1 Unconstrained sparse fitting — [100] words

#### 2.6.2 PHREEQC saturation-index constrained fitting — [100] words

#### 2.6.3 Sparse fitting with ion-exchange enabled — [100] words

#### 2.6.4 Sparse fitting with nitrate and denitrification terms — [100] words

#### 2.6.5 PHREEQC kinetic forward checks of inferred sparse extents — [100] words

#### 2.6.6 Negative-control scenarios with chemically implausible pathways — [100] words

### 2.7 Evaluation metrics — [500] words

#### 2.7.1 Reaction-selection precision, recall and F1 score — [100] words

#### 2.7.2 Concentration reconstruction error — [100] words

#### 2.7.3 Extent-estimation bias and uncertainty — [100] words

#### 2.7.4 Thermodynamic consistency score — [100] words

#### 2.7.5 Kinetic forward-validation error — [100] words

### 2.8 Sensitivity and uncertainty analysis — [250] words

### 2.9 Software implementation and reproducibility — [300] words

### 2.10 Statistical comparison and interpretation thresholds — [200] words

## 3. Results — [3,000] words

### 3.1 Benchmark chemistry and reaction-path characteristics — [350] words

### 3.2 Performance of unconstrained sparse reaction fitting — [400] words

### 3.3 Effect of PHREEQC thermodynamic constraints — [500] words

#### 3.3.1 Reduction of thermodynamically inadmissible reactions — [125] words

#### 3.3.2 Improvement in reaction-selection specificity — [125] words

#### 3.3.3 Effect on concentration reconstruction error — [125] words

#### 3.3.4 Cases where PHREEQC constraints reduce fit flexibility — [125] words

### 3.4 Recovery of known reaction extents in synthetic benchmarks — [450] words

### 3.5 Kinetic PHREEQC forward-check diagnostics — [400] words

### 3.6 Sensitivity to L1 penalty, ion weights and missing ions — [400] words

### 3.7 Application to semi-arid groundwater chemistry — [350] words

#### 3.7.1 Dominant fitted reactions and process-network structure — [100] words

#### 3.7.2 Carbonate, silicate, evaporite and fluoride-related signatures — [100] words

#### 3.7.3 Nitrate and denitrification indications — [75] words

#### 3.7.4 Thermodynamic plausibility of inferred field processes — [75] words

### 3.8 Summary of key validation outcomes — [150] words

## 4. Discussion — [2,300] words

### 4.1 Main finding: thermodynamic constraints improve interpretability more than raw numerical fit — [300] words

### 4.2 Contribution to inverse hydrogeochemical modelling — [300] words

### 4.3 Why sparse reaction fitting is useful in data-limited aquifers — [300] words

### 4.4 Role of PHREEQC in screening and stress-testing Hydrosheaf reaction pathways — [300] words

### 4.5 Implications for semi-arid aquifers and groundwater-quality diagnosis — [250] words

### 4.6 Methodological limitations — [500] words

#### 4.6.1 Linear stoichiometric approximation of nonlinear geochemical systems — [125] words

#### 4.6.2 Dependence on reaction dictionary completeness — [125] words

#### 4.6.3 Ambiguity between mixing, evaporation and reaction residuals — [125] words

#### 4.6.4 Limits of field validation without mineralogical and isotopic confirmation — [125] words

### 4.7 Relevance to the Hydrosheaf PhD thesis — [200] words

### 4.8 Future development: stronger PHREEQC coupling and field validation — [150] words

## 5. Conclusions — [500] words

### 5.1 Principal findings — [150] words

### 5.2 Methodological contribution — [125] words

### 5.3 Practical contribution for data-limited groundwater studies — [125] words

### 5.4 Final thesis-relevant conclusion — [100] words

## 6. Data and Code Availability — [150] words

- Hydrosheaf code, PHREEQC inputs, benchmark scripts, processed outputs and reproducibility archive — [150] words

## 7. Declarations — [350] words

### 7.1 Author contributions — [75] words

### 7.2 Funding and acknowledgements — [75] words

### 7.3 Competing interests — [50] words

### 7.4 Software versions and PHREEQC database statement — [75] words

### 7.5 Reuse limitations of public and synthetic datasets — [75] words

- Total estimated length: [12,000] words, excluding references and supplementary material.
- D) Proposed Tables and Figures
- • Table 1: Hydrosheaf reaction dictionary — Purpose — Documents reactions fitted by the sparse inverse model; key variables: reaction label, process class, stoichiometric vector, expected hydrochemical effect, signed/non-signed status.
- • Table 2: PHREEQC constraint rules — Purpose — Shows how saturation indices and redox indicators are translated into reaction bounds; key variables: mineral, SI threshold, dissolution rule, precipitation rule, bound type.
- • Table 3: Benchmark scenarios — Purpose — Compares unconstrained, PHREEQC-constrained, exchange-enabled, nitrate-enabled and kinetic-validation configurations; key variables: scenario ID, constraints, active reactions, L1 penalty, validation target.
- • Table 4: Synthetic benchmark recovery performance — Purpose — Quantifies whether Hydrosheaf recovers known reactions; key variables: true reaction, inferred reaction, extent error, precision, recall, F1 score.
- • Table 5: Thermodynamic consistency results — Purpose — Summarises improvement after PHREEQC screening; key variables: inadmissible reactions, SI conflicts, consistency score, residual error.
- • Table 6: Kinetic PHREEQC forward-validation results — Purpose — Tests whether inferred extents reproduce downstream chemistry under forward simulation; key variables: reaction extents, predicted concentrations, observed concentrations, RMSE, ion-specific errors.
- • Table 7: Semi-arid field-data application summary — Purpose — Reports process diagnosis from the Northern Ghana hydrochemical dataset; key variables: pathway, dominant reactions, nitrate term, fluoride-related term, fit residual, thermodynamic flag.
- • Figure 1: Hydrosheaf sparse reaction fitting workflow — Flow diagram — Demonstrates transport correction, reaction residual calculation, sparse fitting, PHREEQC screening and validation.
- • Figure 2: Reaction dictionary stoichiometric matrix — Heatmap — Shows which ions are affected by each mineral, exchange, nitrate and denitrification reaction.
- • Figure 3: Synthetic benchmark recovery — Bar/scatter plot — Compares true and inferred reaction extents across controlled PHREEQC-generated cases.
- • Figure 4: Unconstrained versus PHREEQC-constrained reaction selection — Paired bar plot — Shows reduction of chemically implausible reactions after thermodynamic screening.
- • Figure 5: Saturation-index consistency plot — SI-versus-extent plot — Demonstrates whether inferred dissolution/precipitation directions are thermodynamically plausible.
- • Figure 6: Concentration reconstruction accuracy — Observed-versus-predicted scatter plot — Shows ion-wise fit quality for major ions and nitrate.
- • Figure 7: Sensitivity to sparse penalty — Line plot — Shows how L1 regularisation controls reaction parsimony, reconstruction error and thermodynamic consistency.
- • Figure 8: Field application process network — Directed graph — Shows inferred groundwater evolution pathways and dominant reaction classes in the semi-arid case study.
- • Figure 9: Failure-mode diagnostics — Multi-panel figure — Shows examples of ambiguous mixing, missing-ion bias, reaction-dictionary incompleteness and PHREEQC constraint conflicts.
- Implementation alignment update for the revised M5 manuscript
- The outline should now cite hydrosheaf.validation.reaction and M5/m5_inverse_reaction_benchmark. Implemented capabilities include sparse reaction fitting checks, L1 penalty sensitivity, missing-ion sensitivity, transport/reaction residual separation, PHREEQC saturation-index bound violation checks, and linearized reaction-dictionary rank/condition diagnostics.
- The strongest contribution is sparse linear inverse reaction fitting constrained, screened and stress-tested using thermodynamic and forward-validation diagnostics.
- Required wording guardrail: describe Hydrosheaf as a sparse linear inverse reaction model screened and stress-tested using PHREEQC thermodynamic and forward-validation diagnostics, not as a nonlinear coupled PHREEQC inversion engine.
