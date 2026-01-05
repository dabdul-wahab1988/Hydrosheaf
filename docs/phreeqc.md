# PHREEQC Integration

## Overview
Hydrosheaf can call PHREEQC to compute saturation indices (SI), ionic strength, and charge balance error. These outputs drive sign constraints on mineral reaction extents.

## Execution Modes
- `phreeqc_mode="phreeqpython"`: uses the `phreeqpython` bindings.
- `phreeqc_mode="subprocess"`: calls an external PHREEQC executable.

If PHREEQC is unavailable, results are marked with `phreeqc_ok=false` and constraints are disabled for affected edges.

## Input Template
Hydrosheaf generates one `SOLUTION` block per sample and a shared `SELECTED_OUTPUT` block:

```phreeqc
SOLUTION 1
temp      25.0
pH        7.0
units     mg/L
Ca        40.0
Mg        24.0
Na        23.0
Cl        35.0
S(6)      96.0   as SO4
N(5)      62.0   as NO3
F         19.0
Alkalinity 61.0  as HCO3
END
```

Selected output (see `hydrosheaf/phreeqc/templates.py`) requests ionic strength, charge balance, and SI for calcite, dolomite, gypsum, anhydrite, halite, and fluorite.

## Database Selection
If `phreeqc_database` is set, Hydrosheaf prepends a `DATABASE <path>` line to the input file.
The default database ships in `hydrosheaf/databases/phreeqc.dat`.
When using `phreeqpython`, the VIPhreeqc engine can only load databases from a directory
and may reject some PHREEQC `.dat` files. If your custom database fails to load, point
`phreeqc_database` to the `phreeqpython`-bundled `phreeqc.dat` or switch to
`phreeqc_mode="subprocess"` with an external executable.

## SI → Bounds Logic
With threshold `tau` (default 0.2):
- If `SI_u < -tau` and `SI_v < -tau`: dissolution-only (`z >= 0`).
- If `SI_v > tau`: precipitation-only (`z <= 0`).
- Otherwise: free (`z` unbounded).

Non-mineral reactions default to dissolution-only unless explicitly marked signed.

## Assumptions
- HCO3 is treated as alkalinity (as HCO3) when constructing the solution.
- Temperature defaults to 25 °C if `temp_c` is not provided.
