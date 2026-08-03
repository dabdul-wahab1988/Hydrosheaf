#!/usr/bin/env python3
"""Write an evidence-checked second-round review of the M7 revision."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TASKS = json.loads((ROOT / "manuscript/revision/tasks.json").read_text(encoding="utf-8"))
CHANGES = {
    item["comment_id"]: item
    for item in json.loads((ROOT / "manuscript/revision/change_log.json").read_text(encoding="utf-8"))
}


FOUND = {
    "R01": "The title now reads 'Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks.' The aim names 'two independent controlled-synthetic generator systems,' and the Methods has a dedicated non-transfer subsection.",
    "R02": "The abstract uses five labelled paragraphs. The primary paragraph states that M7.4 'failed the prespecified complete gate,' and the follow-up paragraph uses the same wording for M7.5.",
    "R03": "The Introduction now says 'In a targeted, non-systematic search' and names the searched source types. The versioned search JSON records the 30 July 2026 targeted queries and expressly disclaims an exhaustive absence claim.",
    "R04": "Main Table 1 contains seven audit rows and the requested columns for generator, development/test counts, comparator, metrics, lock or multiplicity family, and permitted claim.",
    "R05": "The phrase 'worth collecting' is absent. The revised text says the experiments 'do not establish a value-of-information threshold or a measurement-cost recommendation.'",
    "R06": "The active runner and its archived copy both hash to a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191. This equals the confirmatory-lock value. The recovery manifest identifies the archived task-session source and records that no locked test was rerun.",
    "R07": "The availability statement now correctly says that commit 53beb460... does not contain the later files and is not their release identifier. The complete files are present locally with run IDs and hashes, but the text also states that the new commit, release, and persistent DOI remain required.",
    "R08": "The Methods subsection 'Generator construct validity and non-transfer boundary' states what M7.3 can test, what the scalar M7.4/M7.5 generator can test, and that neither transfers to the other or to field validity.",
    "R09": "The revision openly states that no practical margins or prospective power analysis were prespecified. It then labels a 20,000-replicate simulation as post-review planning and reports probabilities of excluding zero and clearing the planning margins in the main Results and Table S12.",
    "R10": "The Methods defines separate complete families of 120 M7.4 and 560 M7.5 contrasts, with 10,000 shared case-block resamples and two-sided max-z simultaneous 95% intervals. Tables S10 and S11 provide the adjusted results.",
    "R11": "The paper now defines the reaction estimand as discrimination among the six planted archetypes under the specified stoichiometry, mineral assemblage, noise model, and two panels. It explicitly says that out-of-dictionary reactions were not tested.",
    "R12": "The Methods reports M7.3 prevalence as 54/414 in development and 108/827 in locked test, and one-third for the scalar generator. It describes unweighted fitting as a prespecified choice for the generated per-candidate probability estimand, not as the only honest method.",
    "R13": "The Results reports the restored runner hash, states that it matches the confirmatory lock, and says the stored locked test was not rerun. This is supported by the content-addressed archive manifest.",
    "R14": "The Results now reports relative effects and margins. For example, the M7.5 PR-AUC mean difference is given as 3.09% of the edge-local mean and is said not to exceed the post-review 0.02 planning margin. The topology-age effects are treated the same way.",
    "R15": "After the full-family correction, the paper retains the incompatible-cycle conflict-localisation result but withdraws support for noisy/missing overall gains and M7.5 scenario ranking gains. The abstract and conclusion follow those corrected decisions.",
    "R16": "The Results and Discussion heading is now exactly 'Carbonate reactions were not recovered under either tested indicator panel.' The supporting paragraph confines the result to six planted archetypes, 3% noise, and the two tested panels.",
    "R17": "Figure 5 was regenerated as 'Northern Ghana evidence and claim boundary using M7 evidence only.' The M6 panel is absent; the four panels now cover evidence availability, defensible claims, workbook coverage, and the truth-free hold-forward audit.",
    "R18": "The Results says all 198 generated candidates were retained in every arm, with 53 TP, 145 FP, and 0 conditional FN. It reports conditional F1 0.4223, end-to-end F1 0.4206, one pre-scoring missed truth edge, and recall 0.9815. Table S13 gives the arm-level thresholds and counts.",
    "R19": "The Conclusion contains the requested sentence stating that the experiments support the sheaf layer as a prespecified model of non-identity relations and a global-compatibility diagnostic under the tested scalar scenarios, with global fallback for missing endpoints.",
    "R20": "The revised Figure 5 is M7-only and the figure-source manifest identifies it as an M7 field-supportability boundary. There is no M6 evidence inside the result figure.",
    "R21": "The revised main manuscript contains seven figures and four tables. The design map, compact M7.3 decision table, M7.4 means, and M7.5 means remain in the main text; detailed results are in thirteen supplementary tables.",
    "R22": "The revised figure-source manifest records 7.08 inches and a minimum 8 point label size for Figures 6 and 7. Both figures are legible in the Word and LibreOffice renderings.",
    "R23": "The Introduction and Discussion add recent groundwater tracer, hydrochemistry, age-model, and machine-learning sources from 2022 to 2025. They also discuss non-sheaf alternatives, including graph regularisation, Gaussian-process smoothing, flow-network inversion, and residual diagnostics.",
    "R24": "The Methods and Supplement now say that precision-recall curves 'can give a more informative view under class imbalance' and explicitly retain ROC-AUC rather than claiming unconditional scalar superiority.",
    "R25": "The declarations are correctly ordered and the MIT licence is stated, but author roles, funding, competing interests, dataset DOI, software DOI, and final release remain explicit pre-submission requirements. These details are not yet supplied.",
    "R26": "The Methods now says: 'The reaction solver was the inference method under evaluation and was not part of the synthetic generator.'",
    "R27": "The section is titled 'Declarations and open research' and appears in the requested order: author contributions, funding, competing interests, data availability, and code availability.",
    "R28": "Clean DOCX files were rebuilt. Word produced complete 33-page and 22-page PDFs, while LibreOffice produced complete 31-page and 21-page PDFs without the former libpng failure. All 55 Word-rendered pages were inspected, including the final references and Table S13.",
    "R29": "The abstract states that the experiments do not establish general predictive superiority or field validity. The Conclusion says: 'They do not validate HydroSheaf as a whole or establish general predictive superiority over weighted graphs.'",
    "R30": "The Results, Discussion, Conclusion, Figure 7 caption, and Table 4 caption call the selected method local-first/global-fallback and report selected local weight 1.0. The paper says this is not evidence that a general local/global blend is superior.",
    "R31": "The Limitations identify a different generator family or field data with independently measured connectivity as the next claim-bearing step. They explicitly exclude temporal, three-dimensional, vadose-zone, vector-stalk, and active-learning performance from the M7 evidence.",
    "R32": "The abstract, Methods, design table, and Limitations consistently distinguish the 6/12-case M7.3 system from the 32/64-case M7.4 and 64/128-case M7.5 scalar designs.",
    "R33": "The Introduction enumerates seven linked audits, Table 1 contains seven audit rows, and the Conclusion begins 'Across seven linked audits.' No conflicting six-audit statement remains.",
    "R34": "The exact runner source is restored, the current SHA-256 matches the lock, and the Results states this verification. Unit tests exercise the recovered runner without executing the locked test set.",
    "R35": "The manuscript no longer presents d336e87 or 53beb460... as a valid release identifier for all experiments. It gives the M7.3 freeze separately and accurately describes the local M7.4/M7.5 run records, but no new versioned release identifier exists yet.",
    "R36": "Figure 7 says the selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Table 4 gives the same qualification. Internal machine arm names remain only in provenance-oriented supplementary material.",
}


REASONS = {
    "R01": "The revision fixes both the title-level error and the deeper scope problem. A reader can now tell which generator supports each claim.",
    "R02": "The hierarchy and gate failures are now unmistakable in the abstract.",
    "R03": "The authors selected the reviewer-approved claim-restriction route and no longer imply a systematic or exhaustive review.",
    "R04": "The new table contains every requested design field and resolves the confirmatory-versus-diagnostic ambiguity.",
    "R05": "The practical claim is now aligned with the absence of measurement costs and field decision loss.",
    "R06": "The exact locked source is recoverable and content-addressed, while the locked outputs remain untouched.",
    "R08": "The non-transfer boundary is explicit and repeated where the results are interpreted.",
    "R09": "The requested planning analysis is quantitative and is correctly labelled as post-review rather than retrospective preregistration.",
    "R10": "The multiplicity problem is resolved for the entire published contrast families, not only a selected subset.",
    "R11": "The narrower claim is scientifically defensible and directly follows the alternative offered by the reviewer.",
    "R12": "The estimand, prevalence, and design rationale are now precise without making an absolute methodological claim.",
    "R13": "The provenance sentence is now independently verifiable from the restored bytes and lock file.",
    "R14": "Statistical intervals are now paired with relative magnitude and explicit, honestly post-review planning margins.",
    "R15": "The manuscript's supported scenario claims now match the simultaneous intervals.",
    "R16": "The universal non-identifiability claim has been fully withdrawn.",
    "R17": "The M7 field figure no longer imports evidence from another module.",
    "R18": "The identical F1 values and candidate-generation conditioning are now fully explained.",
    "R19": "The evaluative language was replaced with a test-bounded conclusion.",
    "R20": "The visual and its provenance now describe one M7 audit only.",
    "R21": "The main-text table burden is reduced exactly as requested without removing complete supplementary results.",
    "R22": "Final-size legibility is now measurable and documented.",
    "R23": "The positioning is materially broader and more current while remaining honest about the targeted search scope.",
    "R24": "The citation now supports the precise claim made.",
    "R26": "The grammatical error and the inference-generator independence statement are both corrected.",
    "R27": "The requested title and ordering are present. Missing declaration content is assessed separately under R25.",
    "R28": "The former interoperability defect is absent in two independent rendering applications.",
    "R29": "The contribution is now a conditional representation and diagnostic result, not a framework-level validation claim.",
    "R30": "The nomenclature now matches what the selected weight actually did.",
    "R31": "The limitation and next validation step are concrete, and unsupported capability extrapolation is excluded.",
    "R32": "The previously conflicting generator and test-size descriptions are consistent.",
    "R33": "The audit count is now internally consistent.",
    "R34": "The critical source-lock inconsistency has been resolved without a new confirmatory run.",
    "R36": "The explanatory labels prevent the reader from mistaking weight 1.0 for a fitted two-source blend.",
}


PARTIAL = {
    "R07": (
        "The factual availability statement and local technical package are corrected, but the comment also required a committed, tagged, persistently deposited release. That external release has not yet been created, so a reader cannot cite or retrieve an immutable public M7.4/M7.5 package.",
        "Create one versioned repository commit containing the complete M7.4/M7.5 protocols, runners, tests, manifests, and immutable outputs. Tag the release, deposit it in a persistent repository, obtain the software and data identifiers, and replace the pending language with those exact identifiers before submission.",
    ),
    "R25": (
        "The revision appropriately refuses to invent declarations, but the original submission blockers remain unresolved. Explicit notes that metadata are missing are accurate interim text, not final declarations.",
        "The named authors must supply a final CRediT statement, funding declaration, and competing-interest declaration. After the repository deposits are made, insert the dataset DOI, software DOI, licence, and versioned release, then remove the interim instructions to authors.",
    ),
    "R35": (
        "The false commit claim is gone, which is an important correction, but the requested accurate experiment-specific release identifier does not yet exist for M7.4/M7.5. A local run ID and a file hash do not replace a retrievable versioned release.",
        "After publishing the complete M7.4/M7.5 package, cite its commit, tag, and persistent identifier separately from the M7.3 d336e87 freeze. Verify that a clean checkout of that release can regenerate or validate every stated artifact.",
    ),
}


def sid(task: dict) -> str:
    return task["id"].rsplit("-", 1)[-1]


def main() -> None:
    lines = [
        "# Re-Review Report: Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks",
        "",
        "## Part A: Comment-by-comment assessment",
        "",
    ]
    adequate = 0
    partial = 0
    for task in TASKS:
        short = sid(task)
        change = CHANGES[task["id"]]
        lines += [
            f"### {task['id']}",
            "",
            "**Original comment.**",
            "",
            task["comment"],
            "",
            "**Authors' claimed response.**",
            "",
            change["after_summary"],
            "",
            "**What was found in the revised manuscript.**",
            "",
            FOUND[short],
            "",
        ]
        if short in PARTIAL:
            partial += 1
            reason, guidance = PARTIAL[short]
            lines += [
                "**Assessment: PARTIALLY ADDRESSED.**",
                "",
                reason,
                "",
                "**Guidance for further revision.**",
                "",
                guidance,
                "",
            ]
        else:
            adequate += 1
            lines += [
                "**Assessment: ADEQUATELY ADDRESSED.**",
                "",
                REASONS[short],
                "",
            ]

    ratio = 100.0 * adequate / len(TASKS)
    lines += [
        "## Part B: Current assessment summary",
        "",
        "The authors engaged seriously with the technical review. The revised paper is not a cosmetic rewrite. It recovers and archives the exact M7.5 runner, adds full-family simultaneous inference for M7.4 and M7.5, supplies an honestly labelled post-review precision analysis, reduces the main table burden, rebuilds the Ghana figure without M6 evidence, and narrows every major scientific claim to what the controlled generators can establish. Most importantly, the revised paper now answers the ordinary-weighted-graph question directly. The answer is conditional: the sheaf layer encodes non-identity maps and provides global compatibility and conflict localisation, but it does not outperform the strong edge-local comparator overall under the complete gate.",
        "",
        "The statistical interpretation is much stronger. The full 120- and 560-contrast families prevent selective scenario reporting, and the adverse heterogeneous-affine and null overall estimator results remain visible. The planning margins are not presented as if they had been prespecified. The distinction between the MODFLOW/MODPATH generator and the scalar affine generator is also clear throughout, which prevents capability tests from being mistaken for field validation or cross-generator replication.",
        "",
        "The remaining issues are operational but publication-critical. M7.4 and M7.5 still lack a retrievable versioned repository release and persistent identifier. The declarations also still require author-supplied roles, funding, competing interests, and final data and software identifiers. These omissions do not invalidate the numerical findings, but they prevent the package from meeting the reproducibility and submission requirements set by the original critical comments.",
        "",
        "The trajectory is strongly positive. Once the release and declaration work is completed and checked from a clean checkout, no further scientific reanalysis appears necessary for this review round.",
        "",
        "## Part C: New issues",
        "",
        "The revisions did not introduce a new scientific or statistical contradiction. The added precision and practical-magnitude analysis is clearly labelled as post-review planning, and the manuscript does not use it to recast the locked tests as adequately powered. The different Word and LibreOffice page counts reflect layout pagination; both render complete text, figures, tables, and references.",
        "",
        "## Part D: Structured recommendation",
        "",
        f"**Adequacy summary.** {adequate} of {len(TASKS)} comments ({ratio:.1f}%) were adequately addressed. {partial} were partially addressed and none were inadequately addressed. Two overlapping critical reproducibility comments remain partially open because the versioned M7.4/M7.5 release does not yet exist; one major declarations comment also remains partially open.",
        "",
        "**1. Recommendation: Return for major revisions.** The scientific and statistical revisions are sufficient, but the unresolved critical release requirements trigger one more major round under the review decision framework. The required work is bounded: publish and verify the complete versioned M7.4/M7.5 package, mint the persistent identifiers, and replace the interim declaration text with final author-supplied statements.",
        "",
        "**2. Study design and evidential support: Yes.** The design is appropriate for the bounded controlled-synthetic questions, the adverse controls are suitable, and the conclusions now remain within the evidence supported by each generator.",
        "",
        "**3. Methods reproducibility: No, but it can be addressed with revision.** The methods and local artifacts are sufficiently detailed, and the exact locked runner is recovered. Repetition by an independent reader is still impeded by the absence of a complete versioned public release for M7.4/M7.5.",
        "",
        "**4. Statistics and uncertainty: Yes.** The revision uses case-block resampling, full-family simultaneous intervals, transparent post-review planning, and practical-magnitude reporting. The retained claims match those analyses.",
        "",
        "**5. Guidance on overstated claims: This was not needed further.** The earlier overstatements have been rewritten. The current conclusion gives a suitably bounded conditional representation claim.",
        "",
        "**6. Presentation clarity: Yes.** The paper is dense but coherent, the main table count is manageable, the figures are readable at the recorded dimensions, and both DOCX files render completely in Word and LibreOffice.",
        "",
    ]
    out = ROOT / "manuscript/review/M7_SECOND_ROUND_REVIEW_2026-07-30.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out} with {adequate}/{len(TASKS)} adequate")


if __name__ == "__main__":
    main()
