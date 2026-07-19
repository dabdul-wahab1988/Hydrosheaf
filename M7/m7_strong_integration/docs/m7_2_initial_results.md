# M7.2 Initial locked run — audit result

This document records the first full run on development seeds 2101–2106 and
locked seeds 3101–3112. It is retained as a failed/mixed primary analysis and
is not overwritten by the confirmatory amendment.

## Main integration result

Candidate recall was 0.9815. Hydraulics plus constrained chemistry achieved
PR-AUC 0.5127 and F1 0.5568. Adding the original age score produced PR-AUC
0.5214 (difference +0.0087) and F1 0.5512 (difference -0.0055).

The case-block bootstrap did not establish an improvement over the no-age
model:

- PR-AUC difference +0.0070, 95% CI [-0.0108, +0.0258];
- F1 difference -0.0057, 95% CI [-0.0252, +0.0161].

The full model did outperform the within-case permuted-age adverse control:

- PR-AUC difference +0.1470, 95% CI [+0.0562, +0.2264];
- F1 difference +0.0380, 95% CI [+0.0006, +0.0749].

However, the fitted age coefficient had the physically wrong sign. The cause
was the universal-velocity travel residual: slow true MODPATH segments could
receive larger mismatch costs than shorter false candidates. Therefore this
run does not demonstrate defensible positive incremental age value, despite
the positive PR-AUC point estimate.

## Bayesian ages

All 12 cases met the sampler rules. Mean age MAE was 2.99 years, mean bias
-0.88 years, and mean 95% interval coverage 0.875. Individual maximum R-hat
values were 1.0011–1.0040, with minimum bulk ESS above 1,668 and zero
divergences. Numerical convergence and tracer identifiability were both
successful.

## Topology posterior

Graph-size diagnostics cleared the rules in every case. Eight of 12 cases
cleared the stricter every-edge rule. Cases 3107, 3108, 3110, and 3111 had
maximum edge R-hat values 1.0108–1.0129 despite minimum edge ESS values above
607. They are correctly marked non-converged under the preregistered 1.01
limit.

## PHREEQC and reactions

PHREEQC succeeded for all 144 samples. Direction constraints were active on
every fitted candidate edge and materially changed an average of 55.9 edge
objectives per case. Constrained dominant-family accuracy was 0.613, compared
with 0.509 without constraints.

Recovery was strong for sulfate reduction (1.00), silicate weathering (0.958),
denitrification (0.955), and iron reduction (0.75). Carbonate weathering and
precipitation were not separated correctly (both 0.00). This remains an
important model limitation rather than a hidden success.

## Field prequential result

Across 138 wells and 20 August issue batches, graph ridge had MAE 0.1744 in
log1p concentration units, versus 0.3425 for persistence and 0.1802 for an
expanding mean-delta baseline. The paired well-block differences were:

- graph ridge minus persistence -0.1681, 95% CI [-0.1847, -0.1522];
- graph ridge minus expanding mean-delta -0.0058, 95% CI
  [-0.0130, +0.0012].

Thus graph ridge clearly beats persistence but is not distinguishable from the
strong non-graph baseline. The claim remains within-campaign chemistry
hold-forward only.

## Decision

The initial run passes Bayesian aging and active-constraint checks, but fails
the stated goals of universal topology convergence and defensible incremental
age value. The follow-up protocol is fixed in
`m7_2_confirmatory_protocol.md`; it uses new seeds and preserves this result.

