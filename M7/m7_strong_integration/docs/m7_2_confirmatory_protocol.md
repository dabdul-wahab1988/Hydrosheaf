# M7.2 Fresh-seed confirmatory amendment

## Why an amendment is required

The initial locked analysis exposed two method-level defects without using
hidden truth inside inference:

1. a single-velocity travel residual caused a physically reversed fitted age
   coefficient; and
2. 16 transitions per retained topology draw left four cases just above the
   strict every-edge R-hat rule.

The initial outputs under `results/m7_2_strong` and their negative/mixed
interpretation are frozen. Seeds 3101–3112 will not be reused to claim the
amendment works.

## Frozen amendment

Development remains seeds 2101–2106. The new, previously unexecuted
confirmatory seeds are 4101–4112.

Age evidence is restricted to the one-sided posterior probability that
downstream water is older than upstream water, with dating uncertainty
propagated. The universal-velocity travel mismatch remains available as a
diagnostic flag but its fusion cost weight is zero.

The development-fixed compatibility gate is:

```text
age_cost <= 0.05
```

This threshold retained all development true edges and rejected 29.1% of
development false candidates. It is not refit on confirmatory cases. The
no-age baseline probability is retained for compatible candidates and capped
at `1e-6` for incompatible candidates. This enforces the physical sign of age
evidence and prevents an unconstrained fusion coefficient from reversing it.

The adverse control permutes age costs within each aquifer before applying the
same frozen gate.

Topology sampling uses the same target posterior and hard constraints, with 32
transition updates between retained draws instead of 16. This changes
effective exploration only; it does not use truth or alter the posterior.

All other generator, chemistry, PHREEQC, field, feature, threshold-selection,
and bootstrap procedures remain as specified in `m7_2_protocol.md`.

## Confirmatory endpoints

The age amendment is considered successful only if:

- the full model has a positive paired case-block PR-AUC difference from the
  no-age baseline;
- its 95% interval excludes zero;
- the paired F1 difference is non-inferior, with lower 95% limit above -0.02;
- at least 95% of candidate-contained true edges pass the age gate; and
- the full model outperforms permuted age.

Topology convergence succeeds only if all 12 cases satisfy the unchanged
R-hat <= 1.01 and ESS >= 400 rules for graph size and every edge indicator.

If any criterion fails, the failure is retained and no third analysis is
created from these confirmatory seeds.

## Replay

```powershell
.venv\Scripts\python.exe M7\m7_strong_integration\scripts\run_m7_2_strong_integration.py --confirmatory
```

