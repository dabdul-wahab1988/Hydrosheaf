# Frontier active-learning development decision

Decision date: 2026-07-28  
Development run: `DEV-M8-FRONTIER-AL-20260728-01`  
Planned untouched run: `RUN-M8-FRONTIER-AL-20260728-01`

## Purpose

This record separates inspected development evidence from the untouched locked
test. Seeds 7401--7408 were used to fit the frozen measurement and prior
calibration models. Seeds 7451--7458 were used for method development. Locked
seeds 7601--7624 had not been executed or inspected when this decision was
written.

## Development findings

The original topology-uncertainty heuristic was not an adequate scientific
active-learning method: it ranked edges but did not represent measurement
modalities, likelihoods, resource costs, joint hypotheses, decision utility,
model discrepancy, batch redundancy, abstention, or posterior validation.

The implemented method closes those software and methodological gaps with
scenario-robust Bayesian value of information. Development exposed two
failures that were corrected before locking:

1. direct quadrature over 256 particles was too slow; exact predictive
   equivalence-class collapse retained the same scalar EIG while making the
   benchmark practical;
2. an uncalibrated hydraulic prior and unshrunk class-conditional likelihoods
   produced poor Brier scores. A shared development-only logistic prior
   calibration and coefficient shrinkage of 0.25 corrected this without giving
   any strategy privileged information.

Development tuning compared decision weights 0.75, 0.85, 0.95 and 1.00 and
likelihood-shrinkage values 0, 0.25, 0.50, 0.75 and 1.00. The frozen choices are
0.95 decision weight and 0.25 likelihood shrinkage.

On the eight development-tuning cases, the proposed robust policy had median
Brier score 0.04319 and median entropy reduction per cost 0.06067, versus
0.06866 and 0.01155 for random acquisition. Against the strong legacy policy,
the median paired differences were +0.00041 in Brier score and +0.00089 nats
per cost in entropy reduction; their intervals did not establish superiority.

## Frozen claim logic

Because development did not establish superiority to the strong legacy policy,
the confirmatory claim is deliberately narrower. The untouched test must show:

- Brier-score and information-efficiency superiority to random acquisition;
- Brier-score noninferiority to legacy within +0.01;
- entropy-reduction-per-cost noninferiority to legacy within -0.01 nats per
  relative-cost unit;
- no material PR-AUC harm, candidate recall at least 0.80, actionability at
  least 0.90, and byte-identical regeneration.

The two absolute noninferiority margins are deliberately small relative to the
development median legacy values (0.04331 Brier and 0.05616 nats per cost).
They prevent the additional multimodal, cost-aware and robustness capabilities
from being accepted at the price of material predictive or information harm.

If any locked gate fails, the method will receive a mixed or negative claim.
No setting, seed, comparator, metric, interval or margin may change after the
locked-test outcomes are inspected.
