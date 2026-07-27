## Open Research

**Author contributions.** [To be completed by the submitting authors using
a CRediT-style statement, for example: conceptualisation, methodology,
software, formal analysis, investigation, data curation, writing (original
draft), writing (review and editing), and visualisation, attributed to
each named author individually rather than as "all authors contributed
equally."]

**Competing interests.** [To be completed by the submitting authors. State
any financial or non-financial competing interests, or state that none are
declared.]

**Data availability.** The locked benchmark result tables (evidence panel
and case-block contrast summaries, topology-conditioned age sensitivity
tables, reaction-family support tables, the Northern Ghana data-scope
audit record, and the multiplicity-correction robustness check introduced
in this revision) and the per-case simulation provenance records
(executable identity, checksums, and version strings for every generated
case) are held in the authors' version-controlled project repository at
`https://github.com/dabdul-wahab1988/Hydrosheaf`. Prior to submission, these result
tables and provenance records should be deposited in a persistent,
publicly accessible repository (for example Zenodo or a comparable
AGU-recommended repository) and the resulting dataset DOI substituted for
this placeholder. The Northern Ghana workbook audited in Experiment 4
(`data/FieldData/NorthenGhana/NorthernGhana.xlsx`) is this project's own
field dataset and is held in the same repository; it should be deposited
alongside the result tables above and cited by the same dataset DOI. An
earlier revision instead audited a different, antecedent study's own
derived reprocessing of the same boreholes
(`Aquifers_Dataset_Mendeley.xlsx`, documented in `data/FieldData/
NorthenGhana/SI.pdf` as the supplementary information for "Graph-inverted
Ghanaian aquifers under aridification"); that workbook is not this
project's data, is not redistributed by this manuscript, and has been
removed from the analysis (DECISIONS.md).

**Code availability.** The independent synthetic-truth generator, the
evidence-fusion and topology-conditioning analysis scripts, the reaction
bootstrap procedure, and the multiplicity-correction script used for the
robustness check reported in this revision are held in the same
version-controlled repository referenced above (protocol-freeze commit
d336e87). This repository should be archived at a persistent, versioned
release (for example via a Zenodo-linked GitHub release) and the resulting
software DOI substituted for this placeholder before submission. Software
versions used to generate the reported results were MODFLOW 6 (build
6.7.0), MODPATH 7 (version 7.2.001), FloPy 3.10.0, PHREEQC via the
IPhreeqc interface, scikit-learn 1.9.0, numpy 2.4.6, scipy 1.17.1, and
pandas 2.3.3, run under Python 3.14.6, as detailed in Supplementary
Methods.
