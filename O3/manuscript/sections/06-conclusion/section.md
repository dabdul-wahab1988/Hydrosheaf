## 6. Conclusion

Three already-locked HydroSheaf component benchmarks -- age and residence
time against the USGS national groundwater-age release, topology against
MODPATH particle-tracking connectivity, and reaction against a 240-scenario
live-PHREEQC factorial benchmark -- were harmonized under one common
evidence taxonomy without re-running any of the underlying computation. In
all three, independent, prior-free, uncalibrated evaluation shows a
capture-type metric exceeding a matched correctness-type metric: edge recall
above edge precision for topology, reportability rate below within-factor-2
agreement for age, and phase recall close to but below phase precision for
reaction. The gap narrows only when a layer is calibrated or conditioned
against the same reference used to score it, and the two calibration
exercises available here -- an emulation fit on the same held-out folds for
age, a genuine held-out-archetype transfer test for reaction -- are shown to
differ substantially in rigour rather than being interchangeable evidence of
the same kind. The three benchmarks differ by roughly two orders of
magnitude in scale, and their field- and archive-transfer scope is uneven:
broadest for reaction, archive-internal for topology, and absent for age.

None of this is a fourth, independent validation of connectivity, residence
time, or reaction mechanism, and none of it is a second HydroSheaf framework
contribution; both remain the domain of the components' own result packages
and of the companion framework article. What this comparison adds is a
shared vocabulary for reading the three benchmarks against one another, a
disclosed reconciliation of their different metrics and averaging
conventions, and a computational-scale and field-scope accounting that no
single component paper reports. The taxonomy and harmonization scripts are
retained so that a fourth component benchmark, should one be built, can be
added to this comparison without repeating the reconciliation performed
here.
