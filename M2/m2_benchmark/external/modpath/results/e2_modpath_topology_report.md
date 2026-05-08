# E2 MODPATH Topology Validation Report

Run timestamp UTC: 2026-05-05T07:47:07.975276+00:00

Source: USGS Savage MODFLOW/MODPATH archive, DOI `10.5066/F7J102FK`.

## Outputs

- Graph priors: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\modpath\results\modpath_graph_priors.csv`
- Topology agreement: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\modpath\results\modpath_topology_agreement.csv`
- Particle pathline summary: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\modpath\results\modpath_pathline_particles.csv`
- Figure S2: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\figures\figure_s2_modpath_to_graph_prior_real_archive.png`

## Summary

|   n_endpoint_edges |   n_pathline_edges |   true_positive_edges |   false_positive_edges |   false_negative_edges |   edge_precision |   edge_recall |   edge_f1 |   direction_agreement_rate |   mean_source_receptor_overlap |   pathline_elapsed_time_median |   pathline_elapsed_time_p90 |   edge_error_rate |
|-------------------:|-------------------:|----------------------:|-----------------------:|-----------------------:|-----------------:|--------------:|----------:|---------------------------:|-------------------------------:|-------------------------------:|----------------------------:|------------------:|
|                174 |                174 |                   174 |                      0 |                      0 |                1 |             1 |         1 |                          1 |                              1 |                        293.331 |                     550.655 |                 0 |

## Interpretation

Hydrosheaf converts compact MODPATH5 endpoints into directed cell-to-cell graph priors. Compact MODPATH5 pathlines are used as an independent check that endpoint-derived source-receptor edges preserve particle-tracking direction and overlap. This validates the physics-prior ingestion layer; it does not validate geochemical edge inference because the MODPATH archive does not contain paired hydrochemical samples for the same graph nodes.
