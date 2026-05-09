# E4c Northern Ghana Field-Hydrochemistry Validation Report

Run timestamp: 2026-05-08T22:21:51Z

## Source

- Workbook: `data/NorthenGhana/NorthernGhana.xlsx`.
- Sheets used: `Dry` and `Wet`.
- Supplementary information: `data/NorthenGhana/SI.pdf` if retained for provenance.
- Public repository/DOI status: not embedded in the corrected workbook; cite the final public source in manuscript text if required.

## Dataset Profile

- Wells: 160.
- Hydrochemical records: 320.
- Quantitative CBE records used for fitting: 294.
- Source graph edges supplied: 0.
- Generated graph edges used: 1019.
- Candidate edge fits run: 753.
- Top-20 minimum chemistry R2: 0.999.
- Median chemistry R2: 0.979.

## Top 20 Inferred Connectivity

| season   | edge_id              | transport_model   |   objective_score |   chemistry_r2 |
|:---------|:---------------------|:------------------|------------------:|---------------:|
| Wet      | Wet:NGW-125->NGW-155 | mix               |       1.48382e-05 |       0.999996 |
| Wet      | Wet:NGW-083->NGW-093 | mix               |       5.21127e-05 |       0.999991 |
| Wet      | Wet:NGW-146->NGW-149 | mix               |       7.5385e-05  |       0.999986 |
| Wet      | Wet:NGW-119->NGW-081 | mix               |       0.00018675  |       0.999972 |
| Wet      | Wet:NGW-058->NGW-051 | mix               |       0.00027123  |       0.999933 |
| Wet      | Wet:NGW-156->NGW-149 | mix               |       0.000467232 |       0.999914 |
| Wet      | Wet:NGW-058->NGW-049 | mix               |       0.000857003 |       0.999831 |
| Wet      | Wet:NGW-043->NGW-051 | mix               |       0.00100001  |       0.999754 |
| Dry      | Dry:NGW-080->NGW-065 | mix               |       0.00132063  |       0.999875 |
| Wet      | Wet:NGW-139->NGW-149 | mix               |       0.00157756  |       0.999711 |
| Dry      | Dry:NGW-104->NGW-111 | mix               |       0.00173774  |       0.999845 |
| Dry      | Dry:NGW-042->NGW-065 | mix               |       0.00182438  |       0.999827 |
| Dry      | Dry:NGW-053->NGW-065 | mix               |       0.0018472   |       0.999825 |
| Wet      | Wet:NGW-104->NGW-111 | mix               |       0.00192982  |       0.999684 |
| Dry      | Dry:NGW-119->NGW-107 | mix               |       0.00207837  |       0.999882 |
| Wet      | Wet:NGW-125->NGW-127 | mix               |       0.00222588  |       0.999487 |
| Dry      | Dry:NGW-100->NGW-107 | mix               |       0.00229366  |       0.99987  |
| Wet      | Wet:NGW-137->NGW-130 | mix               |       0.0023694   |       0.999516 |
| Wet      | Wet:NGW-150->NGW-149 | mix               |       0.00241698  |       0.999557 |
| Wet      | Wet:NGW-140->NGW-142 | mix               |       0.00243905  |       0.999603 |

## Reviewer Interpretation

This Northern Ghana run replaces the local `manu.xlsx` pilot for M2. It validates Hydrosheaf's ability to ingest the corrected field workbook, compute charge-balance screening, harmonize major ions to mmol/L, infer graph priors with Hydrosheaf's probabilistic graph mechanism, and fit sparse hydrochemical process edges under wet/dry seasonal conditions.

Guardrail: this is field-hydrochemistry and data-limited workflow evidence. The corrected workbook does not contain an independent process-truth graph, tracer-age output, MODPATH pathline truth, or external PHREEQC benchmark solution, so it should not be used to replace E1, E2, or E3.
