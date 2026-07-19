# M7.2 Strong Integration — Final result

## Outcome

M7.2 now provides a reproducible and substantially stronger integration
benchmark, but it does **not** establish positive incremental topology value
from age in the present generator. That negative result must remain visible.

The initial run and the fresh-seed confirmation are both retained. The
confirmatory protocol and code were committed as `3d886b9` before any
4101-series case was generated.

## Fresh-seed confirmation

The confirmation used development seeds 2101–2106 and untouched test seeds
4101–4112. Candidate recall was 0.9537.

### Incremental age

The no-age hydraulics/chemistry model achieved:

- PR-AUC 0.4640;
- ROC-AUC 0.8942;
- Brier score 0.1618; and
- F1 0.5068.

The development-frozen, uncertainty-aware age compatibility gate achieved:

- PR-AUC 0.4623;
- ROC-AUC 0.8893;
- Brier score 0.1602; and
- F1 0.5068.

The paired case-block contrasts were:

- PR-AUC difference -0.00171, 95% CI [-0.00557, +0.00007];
- F1 difference 0.000, 95% CI [0.000, 0.000].

The gate retained 102 of 103 candidate-contained true edges (0.9903) and
rejected 29.2% of false candidates. This shows valid age-based incompatibility
screening, but the rejected false edges were not influential enough to improve
the primary ranking metric. Comparisons with permuted age also failed to
exclude zero.

Therefore the confirmatory incremental-age endpoint failed. No third
fresh-seed analysis was created.

### Bayesian ages

Age estimation itself performed well:

- mean MAE 2.76 years;
- mean bias -0.98 years;
- mean 95% interval coverage 0.917;
- every case met R-hat and ESS rules; and
- zero divergences.

This distinction matters: HydroSheaf dates the synthetic water accurately, but
those ages are largely redundant with hydraulics/chemistry for this particular
edge-ranking task.

### Topology convergence

The topology repair succeeded. All 12 fresh cases met the strict every-edge
rule with:

- maximum edge R-hat 1.00645;
- minimum edge ESS 939;
- maximum graph-size R-hat 1.00309; and
- minimum graph-size bulk ESS 2,303.

The posterior MAP graph remained moderately overconnected: case F1 ranged
roughly 0.39–0.67. Convergence makes this uncertainty/result trustworthy; it
does not make topology accuracy high.

### PHREEQC and reaction coverage

PHREEQC succeeded for all samples and direction constraints were active on
every fitted candidate edge. They materially changed an average of 56.75 edge
objectives per case.

Constrained dominant-family accuracy was 0.602, versus 0.563 unconstrained:

- denitrification 1.00;
- sulfate reduction 1.00;
- silicate weathering 0.917;
- iron reduction 0.75;
- carbonate weathering 0.00;
- carbonate precipitation 0.00.

Thus the constraints are demonstrably active and improve overall family
recovery, but the carbonate dictionary remains non-identifiable/misclassified.

### Field prequential evidence

The field result is unchanged because the data and protocol are unchanged:
graph ridge clearly beats persistence but is not distinguishable from the
expanding mean-delta baseline. This supports within-campaign chemistry
hold-forward only, not field topology, age, or reaction truth.

## What the result means for the paper

The implementation is now much more defensible than M7.1:

- official external flow/pathline simulators;
- a generator that imports no HydroSheaf code;
- separated truth and observables;
- accurate, converged multi-tracer Bayesian ages;
- fully converged constrained topology chains;
- active PHREEQC bounds;
- six-process reaction coverage;
- adverse age controls;
- case/well-block uncertainty; and
- a strict field claim guardrail.

Scientifically, however, the current model is not yet a clean positive
integration story. The strongest honest claims are:

1. Bayesian ages are accurate and computationally converged.
2. Age safely rejects about 29% of false candidate edges at 99% true-edge
   retention, but does not improve overall edge ranking.
3. PHREEQC constraints are operationally active and improve aggregate reaction
   family recovery.
4. Topology posterior convergence is solved, while MAP accuracy remains
   moderate.
5. Carbonate process separation and field external validation remain open.

A Q1 submission can use this as a rigorous benchmark/negative-result section,
but should not claim that age improves topology accuracy. A stronger positive
integration paper would require a scientifically justified setting in which
age contributes information not already encoded by heads and chemistry, plus
multi-year or independent-basin field validation.
