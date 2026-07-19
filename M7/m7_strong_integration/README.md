# M7.2 — Strong integration validation

M7.2 is the stronger, truth-blind integration experiment. It leaves the frozen
M7.1 outputs unchanged and asks whether age, hydraulics, constrained reaction
inversion, and topology uncertainty add defensible value when the truth is
generated outside HydroSheaf.

The benchmark has three deliberately separate evidence sources:

1. Official MODFLOW 6 and MODPATH 7 executables generate heterogeneous flow and
   particle travel times. A standalone script, which imports no HydroSheaf code,
   generates nonlinear chemistry and tracer observations.
2. HydroSheaf receives observables only. Six development cases fit and freeze
   the fusion model and decision thresholds; twelve different seeds are the
   locked test set.
3. The Northern Ghana workbook supplies strict within-campaign prequential
   hydrochemistry prediction. It is not represented as field topology, age, or
   reaction-truth validation.

See [docs/m7_2_protocol.md](docs/m7_2_protocol.md) for the estimands,
guardrails, and convergence criteria. The generated result interpretation is
written to [docs/m7_2_results.md](docs/m7_2_results.md) after the locked run.

## Replay

The official executables are intentionally not committed. Place `mf6.exe` and
`mp7.exe` in `.codex_work/modflow-bin`, then run:

```powershell
.venv\Scripts\python.exe M7\m7_strong_integration\scripts\run_m7_2_strong_integration.py
```

For a non-claim-bearing smoke run:

```powershell
.venv\Scripts\python.exe M7\m7_strong_integration\scripts\run_m7_2_strong_integration.py --quick
```

The full run writes `results/m7_2_strong`. Each case retains blind
observations, separate MODPATH/truth tables, simulator provenance, and
diagnostics. Executable hashes and versions are recorded in every case.

