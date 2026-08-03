# Resolution of O3 review issues

Source review: `O3_manuscript_review.md` (score 60/100, return for major
revisions, no fabrication or critical-flaw override triggered). The issues
below were resolved with data already in hand in the same session; the
remainder are left as genuine outstanding work, listed at the end.

## Resolved

**Methodology Issue 3 (Major) -- no uncertainty quantification.** 95% Wilson
score confidence intervals were computed for every capture-type and
correctness-type proportion in Table 2 (`derive_headline_metrics.py`, using
the underlying counts already present in each component's own result files:
147/174 and 147/302 for topology, 329/1272 and 381/675 for age, 408/681 for
reaction recall; the reaction precision interval was inverted from `M5`'s
own published false-discovery-rate 95% CI rather than recomputed). This
produced a genuine, previously unknown finding: the reaction layer's capture
and correctness intervals overlap (0.562-0.635 vs 0.598-0.681), so its
point-estimate gap is not statistically resolved, while topology's and
age's are. Every section touching the central claim (Abstract, Introduction,
Results 4.2, Discussion 5.1 and 5.6, Table 2, Figure 2) was rewritten to
report this qualification rather than the unqualified point-estimate
pattern. This is a strengthening correction, not a cosmetic fix: it narrows
the paper's central claim to what the evidence actually supports.

**Title and Abstract Issue 2 (Major) -- contribution framed as bookkeeping.**
Abstract rewritten to lead with the empirical, interval-aware finding as the
contribution rather than "vocabulary and accounting."

**Introduction Issue 2 (Major) -- incommensurability disclosed too late.**
The Introduction's statement of the central finding (Section 1, closing
paragraph) now states the interval-resolution qualification and the
non-identical-statistics caveat at first mention, not only in Methods
Section 3.4.

**Results and Discussion Issue 1 (Major) -- "largest separation" ranking
buried the incommensurability caveat.** Results Section 4.2 restructured so
the interval-resolution finding precedes the magnitude comparison, and the
magnitude comparison is now stated only for the two layers (topology, age)
where the gap is resolved.

**Cross-Section Consistency Audit finding -- Methods/Results mismatch on
graph-regularization strengths.** Methods Section 3.2 no longer promises a
weak/medium/strong three-way breakdown that Results does not deliver;
reworded to match what Section 4.4 actually reports.

**Tables and Figures Issues 1-2 -- non-comparable panel scales, undocumented
diagnostic-variant duplication.** Figure 4's caption now states explicitly
that its three panels are not on a comparable scale. Table 1 gained a note
identifying the Hodge-pruned and projected-gradient topology rows as
diagnostic variants of the Head-gradient row, not independent replications.

**References Issue 1 (Major) -- companion manuscripts uncitable.** Added
`hydrosheaf_m2_3_inprep`, `hydrosheaf_m3_inprep`, `hydrosheaf_m4_inprep`, and
`hydrosheaf_m5_inprep` as `@unpublished` bibliography entries and cited each
at its first mention in the Introduction, so a reader has a citable handle
on the source material.

**Novelty Issue 2 (Major) -- salami-slicing risk under-addressed.**
Discussion Section 5.5 expanded with a concrete, falsifiable test for why
this paper and the companion framework article serve different readers,
rather than asserting non-duplication in one sentence.

## Verified, not changed

The citation spot-check (Beven and Freer 2001, Pollock 2016, Oreskes,
Shrader-Frechette and Belitz 1994) found no misrepresentation; no bibliography
change was needed beyond the additions above. The Research Integrity Red
Flag Scan found no fabrication indicators; no change was needed.

## Outstanding, not addressed in this session

1. **Methodology Issue 1 (Critical-adjacent, not fixed by CI addition alone).**
   The three capture/correctness pairs remain three different statistics
   computed three different ways; adding confidence intervals quantifies the
   uncertainty within each pair but does not make the three pairs formally
   commensurable with each other. This is disclosed at every point the claim
   is made (Introduction, Methods 3.4, Results 4.2, Discussion 5.1 and 5.6)
   but not structurally resolved. A future revision could restrict the
   headline cross-component figure to topology and reaction only (the two
   layers with a literal precision/recall pair) and report age separately.
2. **Methodology Issue 2.** The pooled reaction-recall statistic remains a
   new number not published by `M5` itself; disclosed in Methods, Table 2,
   and this ledger, not removed.
3. **Methodology Issue 4.** No empirical re-run confirming `M3`/`M4` results
   are unchanged under current code was performed; the module-dependency
   argument in Section 3.1/5.6 stands as written.
4. **Data Availability.** The Zenodo DOI and GitHub URL remain placeholders
   by the author's explicit instruction; this is the single largest blocker
   to actual submission and is disclosed as such throughout.
5. **Supplementary Methods.** Not yet drafted as a separate delivered file;
   `Outline.md` specifies its target length and required subsections.
6. **CRediT roles and funding statement.** Remain placeholders pending
   actual submission metadata.
7. **Results and Discussion Issue 3 (MRS class-level breakdown), Issue 2
   (external benchmark comparison for the "something general" claim), and
   Issue 4 (connect the age-fraction leakage result back to the central
   thesis).** Not addressed; left for a future revision pass.
