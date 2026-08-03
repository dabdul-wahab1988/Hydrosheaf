# M2.3 outline and word contract

## Target

Computers & Geosciences, research article. The journal requires a substantive
contribution in both the computational and the geoscientific dimension, and
mandates a public, documented repository at submission.

## Positioning

M2.3 is a software and methods article. Its contribution is an implemented,
auditable evidence-to-claim architecture for groundwater inference from
incomplete hydraulic, tracer and hydrochemical evidence, together with an honest
account of which inference targets that architecture can and cannot reach.

The paper is deliberately not a universal benchmark of MODFLOW, MODPATH,
TracerLPM or PHREEQC, and it is not a field-validation study.

## Central argument

The paper's spine is a single contrast that runs through every results section:

> Where an inference target is structurally identifiable from the available
> evidence, the framework recovers it; where it is not, the framework's value is
> that it says so rather than returning a confident number.

The evidence for that argument, in order:

1. Transport parameters and network-constrained residence time are recoverable
   (median absolute error 0.058; log10 R2 0.983 against 0.922 for single-node).
2. Point reaction extents are not (active-term R2 -0.023; 54.1% of truth-zero
   terms leak above the 0.05 mmol/L threshold).
3. Recast as calibrated reaction-family probabilities with abstention, the same
   evidence supports a gated claim (log loss 0.896 against 2.852; false
   commitment 0.038).
4. The same logic governs kinetics (prediction passes, parameter identification
   abstains under k-times-A confounding) and topology (recall 0.845, precision
   0.487).
5. Against external references, agreement is screening-level (held-out
   uncalibrated log10 R2 0.482), and calibrating to the reference raises R2 to
   0.962 while measuring emulation rather than agreement.
6. On real data, the framework's useful output is the evidence ceiling itself.

## Word budget: 7,000 main-text prose words

Prose only; excludes title, headings, references, table content, figure captions
and equations.

| Section | Target | Purpose |
|---|---:|---|
| Abstract | 250 | Problem, contribution, what is and is not established |
| 1. Introduction | 1,000 | Inverse problem, equifinality, gap, contribution, objectives |
| 2. Data | 1,100 | Full reader-facing description of the three datasets and their ceiling |
| 3. Software and methods | 1,500 | Architecture, components, validation design, claim gates |
| 4. Results | 1,900 | Recovery, gates, external comparison, field application |
| 5. Discussion | 1,000 | What the pattern means, comparison with existing tools, limitations |
| 6. Conclusions | 250 | Exact claim and next stage |
| **Total** | **7,000** | |

Supplementary Methods: 3,000 words.

## Display items

| ID | Type | Content |
|---|---|---|
| FIG-1 | Figure | Evidence-to-claim architecture |
| FIG-2 | Figure | Field setting: locations, charge balance, facies, isotopes |
| FIG-3 | Figure | Variable availability and evidence ceiling |
| FIG-4 | Figure | Synthetic recovery: transport, residence time, reaction extents |
| FIG-5 | Figure | Locked programme gates and prospective decision utility |
| FIG-6 | Figure | External comparison: published ages and particle-tracking topology |
| FIG-7 | Figure | Field application: in-sample closure and retained alternatives |
| TAB-1 | Table | Dataset inventory and supportable use |
| TAB-2 | Table | Package capability inventory by module |
| TAB-3 | Table | Locked gate results against predeclared thresholds |
| TAB-4 | Table | Recomputed values that disagree with earlier internal reports |

## Claims that must not appear

- Field validation of groundwater age, directed connectivity or reaction mechanism.
- Measured seasonal change in the Northern Ghana panel.
- Universal superiority over any specialist forward model.
- Calibrated emulation of a reference described as independent agreement.
- Candidate edges described as measured flow paths.
- A formal sheaf theorem.
