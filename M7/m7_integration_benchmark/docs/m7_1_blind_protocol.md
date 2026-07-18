# M7.1 Blind Replicated Integration Protocol

## Purpose and claim boundary

M7.1 tests whether HydroSheaf can recover topology and reactions from observables
when the joint truth is known only to an external evaluator. It is a synthetic
capability and integration test. It is not field validation, a calibrated field
digital twin, or evidence that every HydroSheaf subsystem improves every aquifer.

The original M7 remains locked. M7.1 writes only under
`results/m7_1_blind/`.

## Truth–inference firewall

- `m7_1_truth_model.py` generates observations and hidden truth and imports no
  HydroSheaf code.
- `m7_1_inference.py` imports HydroSheaf but never imports the truth generator.
  Its public contract accepts observations, a random seed, frozen fusion
  coefficients, and a frozen decision threshold.
- The runner is the only evaluation layer that sees both objects.
- A truth-poisoning test reverses edges, corrupts ages, and replaces process
  labels while keeping observations fixed. Serialized inference outputs must
  remain identical.
- Randomness uses NumPy generators with explicit integer seeds; Python's
  process-randomized `hash()` is not used.

## Anti-inverse-crime generator

The generator uses response vectors that are not copied from HydroSheaf's
reaction dictionary. It adds nonlinear concentration response, heteroscedastic
measurement error, variable dilution/recharge mixing, multi-parent mixing, and
occasional unmodelled nitrate pulses. Three topology archetypes are replicated:
branching, convergent, and leaky. Every fourth seed is shifted into an old-water
regime (>70 years), where tritium is intentionally uninformative; the remaining
seeds occupy the tracer-informative regime.

Mg, sulfate, and Fe are excluded from inverse fitting. They are used only after
fitting to score reaction-delta generalization. A process absent from the fitted
dictionary is reported as unsupported rather than relabelled as a successful
recovery.

## Predeclared tiers

### Tier A: replicated held-out topology/chemistry test

- Development aquifers: seeds 10000–10019.
- Test aquifers: seeds 20000–20099.
- Candidate generation is deliberately broad. Candidate-universe recall is
  reported separately from selector performance.
- Actual calls are made to `infer_edges`, `refine_edges_with_sheaf`, and
  `fit_network`.
- Five direct scores are retained: hydraulic-spatial, age-only,
  sheaf-multievidence, chemistry, and equal-weight geometric fusion. The
  sheaf-multievidence score intentionally overlaps hydraulic, isotope, chloride,
  age, and global chemistry evidence and is not treated as an independent
  ablation. The age-only configuration zeros head, isotope, chloride, and global
  chemistry weights.
- A sixth score is a regularized logistic fusion fitted only on development
  candidates. Its coefficients and all method thresholds are frozen before any
  test metric is computed.
- All ablations use the identical candidate universe and independent copies of
  mutable edge objects.

Primary topology metrics are aquifer-level F1, MCC, and PR-AUC. Secondary metrics
are precision, recall, Brier score, 10-bin expected calibration error, structural
Hamming distance, and reachability Jaccard. Means and 95% aquifer-bootstrap
intervals use 2,000 resamples. Paired method-minus-hydraulic intervals are
computed on the same 100 aquifers.

### Tier B: heavy module and diagnostic audit

- Seeds 30000–30007.
- Actual calls are made to `infer_topology_map_edges`, `run_phreeqc`, and
  `infer_network_ages_bayesian`.
- PHREEQC outputs are passed into a second `fit_network` call so thermodynamic
  bounds affect reaction fits. The posterior MAP graph drives Bayesian network
  aging, and posterior mean ages are injected into a second age-only sheaf pass.
  This is a sequential coupling audit, not a converged iterative digital twin.
- The topology posterior uses three chains, 5,000 retained samples, 1,500 burn-in
  samples, a DAG constraint, weak-connectivity constraint, minimum edge count,
  and maximum out-degree. Constrained chains start from the same feasible
  observable-derived graph and follow independent burn-in trajectories.
- Bayesian aging uses four chains and 1,000 retained draws per chain. Reported
  diagnostics include age R-hat, bulk/tail ESS, MCSE, and divergences.
- Reported topology diagnostics include acceptance, edge-count R-hat/ESS, worst edge
  R-hat, minimum edge ESS, PHREEQC success, Bayesian age MAE, and 95% age-interval
  coverage.

Tier B has only eight aquifers and is an execution/convergence audit. It is not
used for superiority claims. Posterior and Bayesian-age MAE/coverage are
withheld from confirmatory interpretation unless R-hat ≤ 1.01, required ESS ≥
400, all expected nodes are present, and divergences are zero. Module execution,
PHREEQC success, posterior convergence, age convergence, and sequential coupling
are separate status fields; a completed call is not labelled converged.

## Reproducibility and failure policy

The runner writes `replicate_metrics.csv`, `heavy_module_audit.csv`,
`summary.json`, and `manifest.json`. The manifest records seeds, frozen
coefficients/thresholds, package versions, configuration, runtime, module call
counts, and every heavy-tier failure. It captures the Git commit and source-tree
status before execution; hashes the complete HydroSheaf runtime tree, all M7
scripts, dependency metadata/lockfile, and databases; and verifies the commit
and every runtime-input hash again at the end. No failed seed is replaced. A
truth-poison failure, heavy-tier exception, commit change, or runtime-input
change makes the command fail.

Reproduce from the repository root after installing `.[dev]`:

```bash
python M7/m7_integration_benchmark/scripts/run_m7_1_blind_benchmark.py
```

Permitted claim: development-trained logistic integration improved selected synthetic
topology metrics on held-out aquifer realizations under this generator.

Prohibited claims: field validation; universal superiority; calibrated
operational digital twin; unique reaction identification; or valid posterior
probabilities when convergence diagnostics fail.
