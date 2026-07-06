# M5 R Figure Layer

These scripts are the preferred R plotting layer for Q1-journal polishing.
They read from `results/m5_results.duckdb` when the R `duckdb` package is
available and never recompute the benchmark analysis. If the R `duckdb` package
is unavailable, the script can still render from the generated CSV/JSON mirrors
written by the same analysis run.

Required R packages:

```r
install.packages(c(
  "DBI", "duckdb", "dplyr", "ggplot2", "patchwork",
  "scales", "viridis", "ragg", "ggrepel"
))
```

Run from `M5/m5_inverse_reaction_benchmark`:

```powershell
Rscript r_figures/plot_m5_publication_figures.R
```

Outputs are written to `figures/r_publication/` and
`figures/r_publication/supplementary/` as PNG, TIFF, and PDF. The Python
database-backed figures in `figures/publication/` remain a reproducible fallback
and quality-control display layer, but the R outputs are the intended
Nature-style manuscript figures.
