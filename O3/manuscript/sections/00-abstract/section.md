## Abstract

Groundwater connectivity, residence time and reaction pathways are usually
benchmarked as separate problems, each against its own reference, leaving it
unclear whether a result reflects the inference problem or the benchmark's
design. The objective of this paper is to harmonize three already-locked
HydroSheaf benchmarks -- age against the USGS national groundwater-age
release (1,272 candidate fits), topology against MODPATH connectivity from
three public MODFLOW archives, and reaction inversion against a
240-scenario live-PHREEQC benchmark (21,600 fits) -- under one evidence
taxonomy, without re-running any of the three. The method is a disclosed,
auditable arithmetic pass over already-published result files: grouping,
pooling confusion counts, and computing ratios and confidence intervals,
with every aggregation discrepancy recorded rather than resolved silently.
The result is a capture-type metric (edge recall 0.845; age-fit
reportability 0.259; reaction-phase recall 0.599) exceeding its matched
correctness-type metric (edge precision 0.487; held-out age agreement within
a factor of 2, 0.564; reaction-phase precision 0.639) under independent,
prior-free, uncalibrated evaluation in all three layers. Attaching 95%
Wilson intervals to each proportion shows this gap is resolved,
non-overlapping, for topology and age, but not for reaction, whose interval
overlap means the point-estimate ordering there is not distinguishable from
sampling variation. Calibrating each resolved layer against its own scoring
reference narrows the gap but differs substantially in rigour between
layers, and the three benchmarks differ by two orders of magnitude in
scale. The contribution is this interval-aware cross-component comparison,
unavailable from any one underlying benchmark alone, offered as a companion
to the HydroSheaf framework article rather than a fourth validation of any
one component.
