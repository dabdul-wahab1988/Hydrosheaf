## Data and code availability

### Reproducibility

All results in this manuscript are derived from immutable locked runs whose
manifests record environment, seed families, source hashes and artifact
checksums:

- `RUN-M8-CONFIRM-20260728-01` (matched-model transport and kinetics),
  protocol locked 2026-07-27T17:22:33Z before execution; seed families
  2026072801 (data) and 2026072802 (bootstrap);
- `RUN-M8-INDEPENDENT-20260728-01` (independent-model robustness and
  active-learning portability); development seed family 2026072810,
  locked-test family 2026072820;
- `RUN-M8-FRONTIER-AL-20260728-01` (topology active learning), protocol locked
  2026-07-28T12:45:27Z before locked-test seeds 7601-7624 were executed;
  development seeds 7401-7408 and tuning seeds 7451-7458.

Lock files, run manifests, raw outputs and the byte-identical reproduction of
the frontier run are stored under `provenance/` in the benchmark package.
Registered artifacts (tables, the confirmatory figure, claim decisions and
the scientific-workflow qualification record) are stored under
`manuscript/artifacts/` and `provenance/runs/`.

### Software

The HydroSheaf package (v0.5.1) and the benchmark scripts are available in the
project repository; the benchmark package is `M8/m8_calibration_benchmark`
with generators `scripts/run_m8_confirmatory.py`,
`scripts/run_m8_independent_active.py`,
`scripts/run_m8_frontier_active_learning.py` and a bundled
`generation_script.py` entry point. The frontier benchmark requires MODFLOW
6.7.0 and MODPATH 7.2.001 binaries, hash-locked in
`m8_frontier_active_learning_protocol.lock.json`. The independent numerical
truth solver is implemented inside the benchmark script and shares no code
with the calibration forward model.

### Data

No field data are used in this study; all results are controlled synthetic
outputs generated with fixed seeds as described in Section 2. Ground-truth
synthetic realisations and per-replicate calibration records are preserved in
the locked run directories.

### Author contributions and competing interests

[AUTHOR CONTRIBUTIONS: to be completed by the authors.]

The authors declare no competing interests.
