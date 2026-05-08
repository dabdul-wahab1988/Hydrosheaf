# M5 Inverse Reaction Benchmark

This folder contains the M5 hydrogeochemical inverse-modelling validation
benchmark. The benchmark supports the defensible Hydrosheaf framing:

> a sparse linear inverse reaction model screened and stress-tested using
> PHREEQC thermodynamic and forward-validation diagnostics.

It does not frame Hydrosheaf as a fully coupled nonlinear PHREEQC inverse
solver.

Run from the repository root:

```powershell
python M5/m5_inverse_reaction_benchmark/scripts/run_m5_inverse_reaction_benchmark.py
```

Outputs are written to `results/`, `tables/`, and `docs/`.
