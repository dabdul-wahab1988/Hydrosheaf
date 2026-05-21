# M4 Topology Benchmark

This folder contains a focused topology-validation benchmark for Hydrosheaf.
The benchmark asks:

> Can Hydrosheaf's reduced-order graph construction reproduce reference
> advective connectivity from MODPATH particle tracking under controlled
> benchmark conditions?

The benchmark keeps two modes separate:

- independent graph inference compared against MODPATH reference connectivity;
- MODPATH-informed graph priors used inside Hydrosheaf.

Run from the repository root:

```powershell
python M4/m4_topology_benchmark/scripts/run_m4_topology_benchmark.py
```

Outputs are written to `results/`, `tables/`, and `docs/`.

For the manuscript-ready analysis package, run:

```powershell
python M4/m4_topology_benchmark/scripts/run_m4_manuscript_analysis.py
```

This command also ingests the M2 public USGS Savage MODFLOW/MODPATH archive
validation into M4 when the M2 external outputs are present.

To refresh the endpoint-derived MODPATH validation from the selected M2
MODPATH output files as well, run:

```powershell
python M4/m4_topology_benchmark/scripts/run_m4_manuscript_analysis.py --run-modpath-validation
```

Manuscript-ready outputs are written to:

- `figures/Manuscript_Ready/`
- `tables/Manuscript_Ready/`
- `docs/02_figures.md`
- `docs/03_tables.md`

The figure and table registers define the allowed claim and required guardrail
for each artefact.
