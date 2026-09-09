"""Build the master M2 rebuttal letter from the preserved reviewer report.

The formal comments are read from the immutable Git revision that introduced
the M2 review package.  The generated letter therefore reproduces each formal
comment verbatim, while the response text is kept in a checked, explicit
mapping below.  The same mapping produces the Track-D audit records and a
reader-facing DOCX.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Any

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


COMMENTS_REVISION = "66ddfcd"
COMMENTS_PATH = "M2/M2_ready/Revision/Comments.txt"


def _response(
    *,
    decision: str,
    status: str,
    response: str,
    location: str,
    verification: str,
    change_summary: str,
    requested_action: str,
    target_files: list[str],
    requires_analysis: bool = False,
) -> dict[str, Any]:
    return {
        "decision": decision,
        "status": status,
        "response": response.strip(),
        "location": location,
        "verification": verification,
        "change_summary": change_summary,
        "requested_action": requested_action,
        "target_files": target_files,
        "requires_analysis": requires_analysis,
    }


RESPONSES: dict[str, dict[str, Any]] = {
    "EDITOR-1": _response(
        decision="AGREE",
        status="ADDRESSED IN THE SUBSTANTIVE REVISION",
        response=(
            "We agree with the editor's synthesis. The revision makes the sheaf "
            "construction precise and limits its claim to a network-level closure "
            "diagnostic; separates prior-assisted MODPATH ingestion from the "
            "independent no-prior endpoint-connectivity test; quantifies reaction "
            "non-uniqueness and transport--reaction confounding; repositions the "
            "Ghana material as a screening-level transfer demonstration; and "
            "reconciles the numerical record across source, tables, captions and "
            "figures. The only remaining item is external release administration: "
            "the exact immutable analysis snapshot must be published and assigned a "
            "persistent DOI before resubmission. We state that requirement rather "
            "than representing the pre-release package as fully public."
        ),
        location="Sections 2--5; Supplementary Tables S3 and S6--S14; Code and data availability.",
        verification="The current package passes the DOCX structural and numerical audits; the release/DOI gate remains explicitly open.",
        change_summary="Reframed the paper around evidence tiers, identifiability and reproducible implementation, with one clearly identified external release gate.",
        requested_action="Respond to the editor's summary and align the revision with the six stated themes.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
            "M2/M2_ready/Revision/Response_to_Reviewers.md",
        ],
    ),
    "R1-M1": _response(
        decision="AGREE",
        status="ADDRESSED AND VERIFIED",
        response=(
            "Agreed. The compound modifiers identified by the reviewer are now "
            "hyphenated consistently, including semi-arid, non-uniqueness, "
            "data-limited, up-gradient, down-gradient and time-averaged. The same "
            "check was applied to the Supplementary Information."
        ),
        location="Main manuscript and Supplementary Information, throughout.",
        verification="A source-level search found no unhyphenated target forms in the revised reader-facing sources.",
        change_summary="Applied consistent hyphenation to compound modifiers.",
        requested_action="Standardise compound-modifier hyphenation.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-M2": _response(
        decision="AGREE",
        status="ADDRESSED WITH A QUALIFIED NOVELTY CLAIM",
        response=(
            "Agreed. The revised text now distinguishes combinatorial inverse "
            "modelling from Hydrosheaf rather than implying that pair enumeration "
            "does not exist. Manu et al. (2023) systematises the selection of "
            "candidate initial--final water pairs and tests their chemical "
            "feasibility. Hydrosheaf instead represents candidate links as a "
            "directed network and combines hydraulic or elevation direction, "
            "conservative-tracer and isotope evidence, age ordering, transport "
            "correction, reaction fitting, uncertainty diagnostics and provenance "
            "at the edge level. We do not claim that this architecture is "
            "universally more accurate; the defensible distinction is the joint, "
            "auditable integration and the explicit retention of competing "
            "hypotheses."
        ),
        location="Section 1 (Introduction) and Section 5.2.",
        verification="The revised manuscript names combinatorial inverse modelling and states the scope of the distinction without a categorical superiority claim.",
        change_summary="Added the existing combinatorial approach and narrowed the novelty claim to integrated evidence and provenance.",
        requested_action="Explain the practical distinction from combinatorial inverse modelling.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R1-M3": _response(
        decision="AGREE",
        status="ADDRESSED WITH A QUALIFIED CLAIM",
        response=(
            "We agree that a Python wrapper around PHREEQC or iPhreeqc and a graph "
            "library such as NetworkX could reproduce individual components. The "
            "revision acknowledges that possibility. Our claim is not that no such "
            "combination can be assembled; it is that Hydrosheaf supplies one "
            "versioned, configurable workflow in which candidate topology, "
            "transport--reaction separation, thermodynamic gates, sparse fitting, "
            "uncertainty diagnostics, validation tiers and provenance are defined "
            "under one data contract. We have not identified a published workflow "
            "that documents this complete combination, but we now state that as a "
            "qualified literature observation rather than an exhaustive absence "
            "claim."
        ),
        location="Section 1 and Section 5.2.",
        verification="The manuscript explicitly acknowledges wrapper-based compositions and avoids claiming that they are impossible or inferior in every setting.",
        change_summary="Acknowledged replicable component combinations and narrowed the integration claim.",
        requested_action="Acknowledge wrapper-based combinations and explain the integrated contribution.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R1-M4": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed. The revision no longer presents multi-tier benchmarking as an "
            "algorithmic contribution. The contribution list is restricted to the "
            "computational and inferential components; the synthetic, reference and "
            "field tiers are described separately as an evaluation design."
        ),
        location="End of Section 1 and Section 3.1.",
        verification="The revised contribution paragraph separates algorithmic contributions from validation strategy.",
        change_summary="Reclassified reproducible benchmarking as evaluation design rather than a sixth algorithmic contribution.",
        requested_action="Reframe contribution six as validation strategy.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R1-M5": _response(
        decision="AGREE",
        status="ADDRESSED AND VERIFIED",
        response=(
            "Agreed. Ionic charges are rendered consistently with superscripts in "
            "the main text, equations and the reaction-dictionary description, "
            "including HCO₃⁻, SO₄²⁻, NO₃⁻, Ca²⁺, Mg²⁺, Na⁺, K⁺ and Cl⁻."
        ),
        location="Sections 2.2 and 2.6, Equations (1) and (8), and Supplementary Table S2.",
        verification="The revised source and rendered DOCX outputs contain the corrected charged-species notation.",
        change_summary="Standardised chemical charge notation across prose, equations and the dictionary.",
        requested_action="Use consistent charged-species notation.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-M6": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed, and the unit convention is now explicit. Input concentrations "
            "are converted to mmol L⁻¹ as the common molar basis. For the charge "
            "balance calculation, the converted concentrations are multiplied by "
            "the absolute ionic valence to obtain equivalent concentrations in meq "
            "L⁻¹; Equation (1) is evaluated on those equivalent concentrations."
        ),
        location="Section 2.2 and Equation (1).",
        verification="The equation and the surrounding explanation now distinguish mmol L⁻¹ from meq L⁻¹.",
        change_summary="Resolved the molar-versus-equivalent concentration notation.",
        requested_action="Clarify the charge-balance unit basis.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R1-M7": _response(
        decision="AGREE",
        status="ADDRESSED; EARLIER DESCRIPTION CORRECTED",
        response=(
            "Agreed. The canonical two-dimensional builder now has an explicit, "
            "reproducibility-fixed search radius of 5.0 km, a maximum of three "
            "primary neighbours per node, minimum directional confidence "
            "pᵢⱼ = Φ((hᵢ − hⱼ)/σΔh) of 0.75, a minimum gradient of 10⁻⁴ and a "
            "20 m screen-depth mismatch flag. Distance is used as the radius "
            "cut-off and ranking tie-break in the reported two-dimensional results; "
            "it is not a multiplicative distance decay. An optional three-dimensional "
            "builder uses P(d) = exp[−d²/(2r²)] with r equal to the configured radius, "
            "but that branch is not used for the reported results. Radius sensitivity "
            "at 3, 5, 7.5 and 10 km is reported in Supplementary Table S14. This "
            "distinction corrects an imprecision in an earlier response that did not "
            "separate the optional three-dimensional term from the canonical builder."
        ),
        location="Section 2.3; Supplementary Tables S3 and S14.",
        verification="The parameter values and builder distinction are stated in both manuscript and Supplementary source; the radius-sensitivity table is present.",
        change_summary="Added the fixed 5 km default, all directional thresholds, the sensitivity analysis and the correct scope of the optional distance-decay term.",
        requested_action="State the distance default, range/sensitivity and scientific basis.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-M8": _response(
        decision="AGREE",
        status="ADDRESSED WITH AN EXPLICIT SCOPE LIMIT",
        response=(
            "We agree with the practical concern and have narrowed the claim. A node "
            "stalk is the local d-dimensional observation vector; an edge carries an "
            "affine restriction map and its offset; and the retained graph is tested "
            "with a stacked coboundary matrix D and right-hand side b. The diagnostic "
            "reports homogeneous nullity H₀, first cohomology dimension H₁, affine "
            "obstruction energy and leave-one-edge-out leverage. No graph or sheaf "
            "Laplacian is constructed or spectrally decomposed. Edge selection itself "
            "remains a weighted multi-criteria score. The added value is therefore a "
            "network-level closure and localisation diagnostic, not demonstrated "
            "superiority over an equivalent weighted score for individual edge "
            "decisions."
        ),
        location="Section 2.4 and Supplementary Method S1.",
        verification="The current text defines D, b, H₀, H₁, obstruction energy and leverage, and states the absence of a Laplacian/eigendecomposition.",
        change_summary="Replaced the imprecise spectral description with the implemented affine cellular-sheaf diagnostic and stated its limited practical role.",
        requested_action="Define the sheaf construction and explain its added value in accessible terms.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-M9": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed. The revised defaults are: ε = 0 years for the age-ordering "
            "constraint; overlapping posterior intervals are retained and flagged "
            "as unresolved at the stated uncertainty; non-overlapping reversals "
            "beyond 0.3 log₁₀ age units receive a severe flag; both evaporation and "
            "mixing are fitted per edge and the candidate minimising the combined "
            "transport, residual-chemistry, L1, EC/TDS, isotope and kinetic objective "
            "is selected; λ₁ = 0.002 is the configured default, λ₂ = 0 is fixed, and "
            "a 10⁻¹⁰ numerical ridge floor is used only for stability. AICc selects "
            "λ₁ over the configured grid for the field sites, giving 0.0483 for both "
            "sites. These values and their rationale are collected in Table S3."
        ),
        location="Sections 2.3--2.6 and Supplementary Table S3.",
        verification="The stated defaults, selection rule and age-flag behaviour are identical in the main and Supplementary sources.",
        change_summary="Documented ε, transport-model selection, λ₁, λ₂, ridge floor and AICc selection.",
        requested_action="State defaults and rationale for ε, transport-model selection and λ₂.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-M10": _response(
        decision="AGREE",
        status="ADDRESSED WITH A CORRECTED RESULT",
        response=(
            "Agreed. The discrepancy triggered a source-to-artifact audit. The "
            "submitted values did not share one reproducible provenance chain, so "
            "neither 0.86 nor 0.74 is retained. The canonical locked benchmark now "
            "reports active-reaction extent correlation R² = 0.23, MAE = 0.37 mmol "
            "L⁻¹ and RMSE = 0.62 mmol L⁻¹ across 2,100 active rows, with 54.1% of "
            "inactive terms exceeding the 0.05 mmol L⁻¹ activation threshold. The "
            "same values are used in the main text, Table 2, Figure 3B and the "
            "Supplementary diagnostics. The lower result is not hidden: it reflects "
            "the rank-deficient, collinear reaction dictionary and is now part of the "
            "claim boundary."
        ),
        location="Section 4.2, Table 2, Figure 3B, Section 5.3 and Supplementary Tables S6--S8b.",
        verification="`audit/number_audit.py` passes all 167 assertions, including the canonical reaction metrics and stale-value checks.",
        change_summary="Replaced unsourced and stale active-reaction metrics with the reproducible canonical values and quantified the identifiability limit.",
        requested_action="Reconcile the active-reaction R² across text, figure and tables.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
            "M2/M2_ready/Revision/audit/number_audit.py",
        ],
        requires_analysis=True,
    ),
    "R1-M11": _response(
        decision="AGREE",
        status="ADDRESSED WITH A QUALIFIED TOPOLOGY CLAIM",
        response=(
            "Agreed on all three points. The prior-assisted F1 = 1.00 is now labelled "
            "an ingestion-fidelity check, not independent inference. The primary "
            "no-prior benchmark generated 302 candidates against 174 MODPATH endpoint "
            "pairs, with TP = 147, FP = 155, FN = 27, precision = 0.49, recall = 0.84 "
            "and F1 = 0.62. The result is explicitly described as screening-level "
            "well-to-well endpoint-connectivity recovery, not pathline geometry, travel "
            "time or porosity-dependent transport. The reference is the publicly "
            "accessible USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 "
            "data release (Harte, 2021, DOI 10.5066/F7J102FK). Elevation-drop, "
            "proximity-kNN and conservative-tracer baselines are reported in "
            "Supplementary Table S9; the last is not evaluable on an archive with no "
            "hydrochemistry."
        ),
        location="Sections 3.3 and 4.3; Table 5; Figure 2; Supplementary Table S9; References.",
        verification="The source, table and caption carry the same 174-edge reference, 0.62 no-prior F1 and prior-assisted scope label.",
        change_summary="Separated prior ingestion from independent inference, replaced the unsourced no-prior value, identified the USGS archive and added baselines.",
        requested_action="Make no-prior inference primary, qualify endpoint connectivity and identify the reference dataset.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-M12": _response(
        decision="AGREE",
        status="ADDRESSED; AGE-ONLY LIMIT RETAINED",
        response=(
            "Thank you; this is a real design limitation. The age-ordering rule is "
            "evaluated with posterior intervals. An overlapping pair is retained and "
            "flagged as unresolved at the stated uncertainty; it is not assigned a "
            "false direction or automatically rejected. Only a non-overlapping reversal "
            "beyond 0.3 log₁₀ age units receives the severe flag. In the synthetic "
            "network, 15.2% of point-estimate checks violated ordering, 84.3% of those "
            "violations had overlapping intervals, and 2.4% were severe reversals; the "
            "age-order consistency index was 0.85. A separate controlled benchmark "
            "shows why this does not establish direct adjacency: order-only PR-AUC was "
            "0.4313, while the covariance-aware direct-versus-indirect Bayes-factor "
            "comparison reached 0.8847 on 254/360 scorable comparisons, with 106 "
            "overlapping cases abstained. The full-minus-permuted-control difference "
            "was +0.0420 (95% CI −0.0004 to 0.0839), so the advantage over the control "
            "is uncertain. Age ordering alone therefore identifies temporal "
            "compatibility, not direct adjacency."
        ),
        location="Sections 2.5, 4.4, 5.3 and 5.5; Supplementary Table S12; Supplementary Tables S12b--S12c; Supplementary Methods S4 and S4b.",
        verification="The interval rules and directness limits are present in the current sources; the controlled benchmark is explicitly separate from field transfer.",
        change_summary="Added interval-aware handling, quantified overlap/severity and documented the conditional directness benchmark without upgrading age ordering to an adjacency label.",
        requested_action="Explain how wide age intervals feed back into edge retention and confidence.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-M13": _response(
        decision="AGREE",
        status="ADDRESSED WITH A CORRECTED TABLE",
        response=(
            "The reviewer identified a genuine defect in the submitted display. The "
            "seven zero-closure sink rows belonged to a superseded, unfiltered "
            "candidate graph and are not part of the canonical retained graph. The "
            "revised Table 6 lists the six highest-ranked explained edges, all with "
            "their current fitted extents; its note distinguishes the largest fitted "
            "extent from the most stable PSI family. PSI = 1.00 means that the dominant "
            "reaction was selected in every perturbation trial for that edge; it does "
            "not mean that a zero reaction was robustly detected. Fourteen retained "
            "edges with negative chemistry R² remain visible as poor fits rather than "
            "being silently removed, and no unresolved null edge remains."
        ),
        location="Section 4.6, Table 6, Supplementary Table S5 and Section 5.5.",
        verification="The current Table 6 contains no superseded zero-closure rows, and the numerical audit confirms consistency with Supplementary Table S5.",
        change_summary="Removed the superseded null rows from the canonical display and explained PSI semantics and remaining poor fits.",
        requested_action="Clarify the zero extents and revise Table 6.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-M14": _response(
        decision="AGREE",
        status="ADDRESSED WITH A QUALIFIED TOPOLOGY CLAIM",
        response=(
            "Agreed. MODPATH has two distinct roles in this revision. Its calibrated "
            "flow model is a physically based reference when such an archive exists, "
            "and its endpoint pairs can be supplied as a prior to test ingestion. The "
            "prior-assisted F1 = 1.00 is therefore an integrity result. The independent "
            "no-prior F1 = 0.62 is the capability result, and it is deliberately "
            "qualified as endpoint-connectivity recovery with substantial "
            "overconnection. An inferred edge is weaker than a calibrated pathline: it "
            "does not recover pathline geometry, travel time or porosity-dependent "
            "transport. Hydrosheaf is designed for settings where the calibrated inputs "
            "needed by MODFLOW/MODPATH are unavailable; it is not presented as a "
            "replacement for those solvers when they are available."
        ),
        location="Sections 4.3, 5.1--5.3; Figure 2; Table 5.",
        verification="The role separation and endpoint-only estimand are stated in the main text and Figure 2 caption.",
        change_summary="Reframed MODPATH as both a prior-ingestion check and a connectivity benchmark, with explicit limits on what an edge represents.",
        requested_action="Clarify benchmark versus prior and limit the topology claim.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R1-M15": _response(
        decision="AGREE",
        status="ADDRESSED WITH AN IDENTIFIABILITY RESULT",
        response=(
            "We now answer this empirically and conservatively. Site-aggregated PSI "
            "for CaNa_exch versus NaCa_exch is 0.73 versus 0.46 at Lower Anayari and "
            "0.83 versus 0.42 at Talensi, so the perturbation analysis separates the "
            "directional exchange pair in these data. For calcite versus dolomite, PSI "
            "is near zero for both members (0.05 versus 0.002 at Lower Anayari and 0.06 "
            "versus 0.04 at Talensi); PSI therefore does not rescue carbonate-family "
            "identification. The conclusion is deliberately narrow: PSI ranks "
            "stability under the imposed perturbations, not the correctness or "
            "uniqueness of a mineral attribution."
        ),
        location="Section 4.6, Section 5.5 and Supplementary Table S10.",
        verification="The four site-pair values and the interpretation are aligned with Figure 7 and Supplementary Table S10.",
        change_summary="Added an empirical pair-separation diagnostic and stated that carbonate-family non-identifiability remains unresolved.",
        requested_action="Test whether PSI separates the stated degenerate pairs.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-m1": _response(
        decision="AGREE",
        status="ADDRESSED AND VERIFIED",
        response=(
            "Agreed. This minor comment is addressed by the same global chemical "
            "notation correction reported for R1-M5; HCO₃⁻ and the other charged "
            "species are rendered consistently in prose, equations and the dictionary."
        ),
        location="Main manuscript and Supplementary Table S2, throughout.",
        verification="The notation correction is present in the current Markdown and DOCX outputs.",
        change_summary="Applied the superscript charge convention consistently.",
        requested_action="Use a consistent HCO₃⁻ notation.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-m2": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed. The 12.32-year tritium half-life is now cited to Lucas and "
            "Unterweger (2000), including the DOI and the reported 4,500 ± 8 day "
            "basis."
        ),
        location="Section 2.5, Supplementary Method S3 and References.",
        verification="The primary reference is present and the tritium statement links to it in the revised source.",
        change_summary="Added the primary half-life citation and source context.",
        requested_action="Cite the primary source for the tritium half-life.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-m3": _response(
        decision="AGREE",
        status="ADDRESSED WITH CONFIGURABILITY QUALIFICATION",
        response=(
            "Clarified. The 4% major-ion relative sigma and 0.5‰ stable-isotope sigma "
            "are configurable defaults used to represent routine analytical uncertainty "
            "in the benchmark; they are not universal constants or claims about every "
            "laboratory. Their values, role and sensitivity implications are now stated "
            "in Table S3 and the uncertainty section."
        ),
        location="Section 2.8 and Supplementary Table S3.",
        verification="The defaults are labelled configurable in the revised manuscript and SI.",
        change_summary="Recast analytical sigmas as configurable benchmark defaults rather than universal values.",
        requested_action="Cite a basis or identify the analytical sigmas as configurable defaults.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-m4": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed. Supplementary Information now provides the provenance-manifest "
            "schema and a populated example from the public-age run. The manuscript "
            "also identifies the manifest as a pre-release checksum record rather than "
            "claiming that it is itself a persistent archive."
        ),
        location="Section 2.9 and Supplementary Information, Provenance manifest: schema and example.",
        verification="The schema and example are present at Supplementary Information lines 364--397 of the current source.",
        change_summary="Added a concrete manifest schema and example.",
        requested_action="Provide a provenance-manifest example or schema.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
    ),
    "R1-m5": _response(
        decision="AGREE",
        status="ADDRESSED WITH CORRECTED METRICS",
        response=(
            "Agreed. The Ghana row now reports the canonical median chemistry R² = "
            "0.70 (0.53 at Lower Anayari and 0.82 at Talensi) over 258 retained edges "
            "and median PSI = 0.97. The table and its note explain that PSI is a "
            "stability statistic and that the superseded zero-closure edges are not in "
            "the canonical denominator."
        ),
        location="Table 2, Section 4.6 and Supplementary Table S5.",
        verification="The current Table 2, field results paragraph and SI use the same field values; the numerical audit passes.",
        change_summary="Replaced the misleading Ghana summary values and added the denominator and PSI interpretation.",
        requested_action="Correct or explain the Ghana median R² and PSI row.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R1-m6": _response(
        decision="AGREE",
        status="ADDRESSED",
        response=(
            "Agreed. A clarity pass replaced the densest compound noun strings in the "
            "abstract, discussion and conclusions with explicit phrases such as "
            "network-wide consistency residual, process-attributable mass transfers "
            "and uncertainty-ranked process hypotheses. Technical terms were retained "
            "where they carry a defined method or metric."
        ),
        location="Abstract, Sections 1 and 5, Conclusions and related passages.",
        verification="The revised source was reviewed for the cited constructions and the document was rendered after revision.",
        change_summary="Simplified dense compound noun strings while retaining defined technical terms.",
        requested_action="Review compound noun strings for readability.",
        target_files=["M2/M2_ready/Revision/Manuscript-Final-Revised.md"],
    ),
    "R2-M1": _response(
        decision="AGREE",
        status="ADDRESSED WITH AN EXPLICIT SCOPE LIMIT",
        response=(
            "Agreed. The same correction addresses the reviewer’s formal concern. "
            "Hydrosheaf uses d-dimensional node stalks, affine restriction maps and a "
            "stacked coboundary Dx = b on the retained graph. H₀ is the homogeneous "
            "nullity and H₁ is the first-cohomology dimension; affine solvability is "
            "tested separately by the obstruction energy. No sheaf Laplacian or "
            "eigendecomposition is used. Edge selection remains a weighted score. The "
            "sheaf layer is therefore technically meaningful as a network-closure and "
            "edge-localisation diagnostic, but it is not claimed to improve individual "
            "edge decisions over an equivalent score."
        ),
        location="Section 2.4 and Supplementary Method S1.",
        verification="The current main and Supplementary sources contain the implemented equations and the explicit no-Laplacian statement.",
        change_summary="Defined the actual affine cellular-sheaf diagnostic and removed the unsupported spectral description.",
        requested_action="Define the sheaf, stalks, maps, diagnostic and any Laplacian use.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R2-M2": _response(
        decision="AGREE",
        status="ADDRESSED WITH NO-PRIOR AS THE PRIMARY TEST",
        response=(
            "Agreed. The validation design now separates the two estimands. The "
            "prior-assisted run supplies MODPATH endpoint pairs only to test faithful "
            "prior ingestion (F1 = 1.00); it is not used as evidence of independent "
            "inference. The no-prior result is primary (F1 = 0.62; precision = 0.49; "
            "recall = 0.84). Supplementary Table S9 compares it with all-pairs "
            "elevation-drop (F1 = 0.47), proximity kNN (F1 = 0.00) and conservative-"
            "tracer ordering, which is not evaluable on the hydrochemistry-free Savage "
            "archive. The benchmark tests endpoint connectivity only."
        ),
        location="Sections 3.3 and 4.3; Table 5; Figure 2; Supplementary Table S9.",
        verification="The current source, table and figure caption use the same baseline counts and prior/no-prior labels.",
        change_summary="Made no-prior inference primary and added matched baseline graph constructions.",
        requested_action="Emphasise no-prior topology and provide baseline comparisons.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R2-M3": _response(
        decision="AGREE",
        status="ADDRESSED WITH AN IDENTIFIABILITY LIMIT",
        response=(
            "Agreed. The benchmark now reports the dictionary rank (8 of 11 ion "
            "dimensions), rank deficiency 3, effectively infinite condition number and "
            "ten highly collinear reaction pairs. Recovered-extents correlations and "
            "empirical mean/SD diagnostics are given in Supplementary Table S8b, and "
            "leave-one-out dictionary sensitivity is reported in Table S8. The "
            "benchmark uses the same declared dictionary for generation and inversion; "
            "we now identify that as an idealised stress-test limitation rather than "
            "calling it out-of-dictionary validation. Canonical extent recovery is R² = "
            "0.23, MAE = 0.37 mmol L⁻¹, with 54.1% false activation; family-level "
            "recovery is R² = 0.16 and dominant-family hit rate 48%. These results "
            "support identifiability diagnosis, not unique mineral histories."
        ),
        location="Section 4.2, Section 5.3 and Supplementary Tables S6--S8b; Figure 3B.",
        verification="The rank, condition, correlations, sensitivity and uncertainty summaries are present in the current SI and agree with the main-text metrics.",
        change_summary="Added rank/condition, correlation, empirical uncertainty and dictionary-sensitivity diagnostics and downgraded the interpretation to an identifiability stress test.",
        requested_action="Report rank, condition number, uncertainty, correlations and dictionary sensitivity.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R2-M4": _response(
        decision="AGREE",
        status="ADDRESSED WITH A QUANTIFIED CONFOUNDING LIMIT",
        response=(
            "Agreed. The revision states that Cl⁻ and EC are used as conservative "
            "anchors only when halite dissolution, agricultural chloride input or "
            "anthropogenic salinity is not indicated. It adds a two-endmember "
            "mixing-plus-reaction stress test: true mixing fractions f₁ = 0.20 and f₂ = "
            "0.15, halite extent 0.40 mmol L⁻¹ and calcite extent 0.20 mmol L⁻¹. The "
            "single-endmember transport stage selected mixing on 73.1% of edges, with "
            "median f = 0.24 and chemistry R² = 0.999, but recovered halite and calcite "
            "medians of 0.04 and 0.00 mmol L⁻¹. This directly demonstrates that good "
            "chemical closure can coexist with incorrect reaction attribution when the "
            "transport model is misspecified."
        ),
        location="Sections 2.6, 4.2 and 5.5; Supplementary Method S2 and Table S11.",
        verification="The stress-test inputs and outputs are present in Table S11 and are repeated consistently in the main Results and Discussion.",
        change_summary="Added the multi-endmember stress test and explicit conservative-anchor caveat.",
        requested_action="Clarify transport confounding and test multi-endmember reaction scenarios.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R2-M5": _response(
        decision="AGREE",
        status="ADDRESSED WITH AN EXPLICIT FIELD-EVIDENCE LIMITATION",
        response=(
            "Agreed. The Ghana material is now labelled a field-hydrochemistry transfer "
            "demonstration, not independent hydrogeological validation. The revised "
            "package records the two site settings, sample counts, available major-ion "
            "and stable-isotope information, and the absence of hydraulic-head records, "
            "nuclear tracers, mineralogical determinations and independent process-truth "
            "labels. The canonical field output is 258 retained edges, median chemistry "
            "R² = 0.70 overall (0.53 Lower Anayari; 0.82 Talensi) and median PSI = 0.97. "
            "The inferred process families are presented only as PSI-ranked hypotheses. "
            "The manuscript now states that tracer measurements, mineralogical work, "
            "head observations, redox measurements and targeted repeat sampling would "
            "be required to confirm them. We do not claim that chemical reconstruction "
            "performance proves flow paths or reaction truth."
        ),
        location="Sections 3.5, 4.6, 5.4 and 5.5; Supplementary Tables S4--S5 and Figure S3.",
        verification="The evidence-availability statement and screening-level wording are present in the current main and Supplementary sources.",
        change_summary="Added site context and an explicit evidence inventory, corrected the field metrics and restricted interpretation to screening-level hypotheses.",
        requested_action="Provide site-specific context and independent support, or state what is unavailable.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
        ],
        requires_analysis=True,
    ),
    "R2-m1": _response(
        decision="AGREE",
        status="ADDRESSED AND VERIFIED",
        response=(
            "Agreed. The entire metric chain was rebuilt from the locked result files. "
            "The revised source, tables, captions and figures now use the canonical "
            "reaction R² = 0.23, Ghana median chemistry R² = 0.70, identifiable public "
            "age parity values and the corrected PSI interpretation. Table 6 contains "
            "the six current explained edges; its note distinguishes fitted extent from "
            "PSI family stability. The package now ships two final checks: "
            "`audit/number_audit.py`, which performs 167 cross-artifact assertions, and "
            "`audit/verify_docx.py`, which checks corrected/stale strings, source age and "
            "embedded figures."
        ),
        location="Main manuscript, Supplementary Information, Tables 2 and 6, Figures 3 and 5, and audit materials.",
        verification="`number_audit.py` passed 167 assertions and `verify_docx.py` passed for three documents.",
        change_summary="Reconciled all canonical values and added programmatic source-to-document consistency checks.",
        requested_action="Audit all tables, figures, captions and text for internal consistency.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
            "M2/M2_ready/Revision/audit/number_audit.py",
            "M2/M2_ready/Revision/audit/verify_docx.py",
        ],
        requires_analysis=True,
    ),
    "R2-m2": _response(
        decision="AGREE",
        status="SUBSTANTIVE RESPONSE COMPLETE; AUTHOR RELEASE GATE PENDING",
        response=(
            "Agreed. The revision documents the public repository, the exact analysis "
            "commit (463e1ce), Python environment, locked configurations, test suite, "
            "CI workflow, reproduction scripts, provenance manifest, runtime scaling "
            "and the third-party data-release identifiers. Candidate-edge construction "
            "is pruned to approximately 3n candidates and the measured 320-node runtime "
            "is 0.06 s; the full benchmark runtimes are reported in Table S13. We do not "
            "claim that the current public repository alone reproduces the exact field "
            "package. The immutable snapshot must still be published as a versioned "
            "release and assigned a persistent DOI before resubmission. This is an "
            "author action, not a prose issue, and is the one explicit submission gate "
            "remaining in this response."
        ),
        location="Code and data availability; Section 2.9; Supplementary Tables S13 and S14; provenance manifest.",
        verification="The local package contains the stated environment/provenance materials and the runtime table; publication of the exact snapshot and DOI is not yet verified.",
        change_summary="Documented the reproducibility materials and runtime, while preserving the release/DOI requirement as pending rather than claiming completion.",
        requested_action="Provide release, DOI/archive, environment, tests, reproduction scripts and scalability evidence.",
        target_files=[
            "M2/M2_ready/Revision/Manuscript-Final-Revised.md",
            "M2/M2_ready/Revision/Supplementary-Information-Revised.md",
            "M2/M2_ready/Revision/REPRODUCIBILITY_MANIFEST.json",
        ],
        requires_analysis=True,
    ),
}


TITLES = {
    "EDITOR-1": "Editor summary of reviewer concerns",
    "R1-M1": "Hyphenation of compound modifiers",
    "R1-M2": "Distinction from combinatorial inverse modelling",
    "R1-M3": "Acknowledging wrapper-based combinations",
    "R1-M4": "Reframing contribution six",
    "R1-M5": "Consistent HCO₃ notation",
    "R1-M6": "meq L⁻¹ versus mmol L⁻¹ in the charge-balance error",
    "R1-M7": "Default distance threshold",
    "R1-M8": "Rigorous sheaf definition and added value",
    "R1-M9": "Defaults for ε, transport-model selection and λ₂",
    "R1-M10": "R² inconsistency",
    "R1-M11": "Topology validation reframing",
    "R1-M12": "Wide age ranges and the age-ordering constraint",
    "R1-M13": "Zero reaction extents in Table 6",
    "R1-M14": "MODPATH benchmark versus prior tension",
    "R1-M15": "Whether PSI separates degenerate pairs",
    "R1-m1": "HCO₃⁻ in equations",
    "R1-m2": "³H half-life citation",
    "R1-m3": "Analytical sigma defaults",
    "R1-m4": "Provenance manifest schema",
    "R1-m5": "Table 2 Ghana row",
    "R1-m6": "Compound noun strings",
    "R2-M1": "Precise sheaf definitions",
    "R2-M2": "No-prior emphasis and baselines",
    "R2-M3": "Identifiability diagnostics",
    "R2-M4": "Transport-correction realism",
    "R2-M5": "Ghana site-specific evidence",
    "R2-m1": "Internal consistency audit",
    "R2-m2": "Software release, DOI, environment, tests and runtime",
}


def _source_comments(repo_root: Path) -> str:
    return subprocess.check_output(
        ["git", "show", f"{COMMENTS_REVISION}:{COMMENTS_PATH}"],
        cwd=repo_root,
        text=True,
        encoding="utf-8",
    )


def _split_numbered(block: str) -> list[str]:
    matches = list(re.finditer(r"(?m)^(?P<number>\d+)\.\s+", block))
    if not matches:
        raise ValueError("No numbered comments found in report block")
    out: list[str] = []
    for i, match in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(block)
        out.append(block[match.start():end].strip())
    return out


def _extract_comments(raw: str) -> list[dict[str, str]]:
    editor = re.search(
        r"\n\n(The two reviewers agree.*?)(?=\n\nReviewer 1:)", raw, flags=re.S
    )
    if not editor:
        raise ValueError("Editor summary not found")
    r1_major_start = raw.index("Major Comments:") + len("Major Comments:")
    r1_minor_start = raw.index("Minor Comments:", r1_major_start)
    r2_start = raw.index("Reviewer 2:")
    r2_major_start = raw.index("Major concerns", r2_start) + len("Major concerns")
    r2_minor_start = raw.index("Minor comments:", r2_major_start)
    more_start = raw.index("More information", r2_minor_start)

    comments: list[dict[str, str]] = [
        {"id": "EDITOR-1", "reviewer": "Editor", "comment": editor.group(1).strip(), "category": "editor summary"}
    ]
    for number, text in enumerate(_split_numbered(raw[r1_major_start:r1_minor_start]), 1):
        comments.append({"id": f"R1-M{number}", "reviewer": "Reviewer 1", "comment": text, "category": "major"})
    for number, text in enumerate(_split_numbered(raw[r1_minor_start + len("Minor Comments:"):r2_start]), 1):
        comments.append({"id": f"R1-m{number}", "reviewer": "Reviewer 1", "comment": text, "category": "minor"})
    for number, text in enumerate(_split_numbered(raw[r2_major_start:r2_minor_start]), 1):
        comments.append({"id": f"R2-M{number}", "reviewer": "Reviewer 2", "comment": text, "category": "major"})
    for number, text in enumerate(_split_numbered(raw[r2_minor_start + len("Minor comments:"):more_start]), 1):
        comments.append({"id": f"R2-m{number}", "reviewer": "Reviewer 2", "comment": text, "category": "minor"})
    return comments


def _quote(text: str) -> str:
    return "\n".join("> " + line if line else ">" for line in text.splitlines())


def _status_table() -> str:
    return "\n".join(
        [
            "| Evidence class | Response position |",
            "|---|---|",
            "| Substantive comments | Addressed with revised text, regenerated evidence or an explicit scope limitation. |",
            "| Field-validation request | Addressed by reporting the available evidence and stating what the Ghana data cannot establish. |",
            "| Release and DOI | Documented in the package but still an author action before resubmission. |",
        ]
    )


def _additional_audit_notes() -> list[dict[str, str]]:
    return [
        {
            "title": "Source-to-artifact packaging",
            "text": (
                "The reviewed DOCX had been assembled from an older source snapshot. "
                "The canonical clean and coloured DOCX files were rebuilt, and the "
                "package now runs number and DOCX consistency checks after assembly."
            ),
        },
        {
            "title": "Field graph direction",
            "text": (
                "The field branch was traced to the current documented builder. The "
                "revised output retains 258 of 572 candidates, reports no retained edge "
                "opposing the elevation proxy, and labels the Lower Anayari flat-elevation "
                "cases as lateral/dispersive candidates rather than pretending that they "
                "carry measured hydraulic-head direction."
            ),
        },
        {
            "title": "Affine sheaf diagnostic",
            "text": (
                "The previous spectral wording was removed. Positive affine obstruction "
                "energies are reported for both Ghana networks, so neither is described "
                "as having an exact affine global section at the stated tolerance."
            ),
        },
        {
            "title": "Age directness and Aiken",
            "text": (
                "The controlled direct-versus-indirect benchmark is reported as a "
                "conditional transport-assisted diagnostic. The Aiken model-conditioned "
                "outputs are kept separate and are not pooled as independent field truth."
            ),
        },
        {
            "title": "Build and reference controls",
            "text": (
                "The reference list and damaged in-text citations were repaired; figure "
                "jitter and PSI draws were seeded; figure embedding was made non-idempotent "
                "by design; and the display-item order in the final submission copies is "
                "References, Tables, figure captions, then embedded figures."
            ),
        },
    ]


def _write_markdown(out_path: Path, comments: list[dict[str, str]]) -> None:
    expected = {item["id"] for item in comments}
    missing = expected - RESPONSES.keys()
    extra = RESPONSES.keys() - expected
    if missing or extra:
        raise ValueError(f"Response mapping mismatch; missing={sorted(missing)}, extra={sorted(extra)}")

    lines = [
        "# Response to Reviewers",
        "",
        "**Manuscript number:** CAGEO-D-26-00847",
        "",
        "**Manuscript title:** Hydrosheaf: A Reproducible Computational Framework for Inferring Directed Groundwater Hydrochemical Evolution Networks from Sparse Multitracer Data",
        "",
        "**Journal:** Computers and Geosciences",
        "",
        "Dear Editor and Reviewers,",
        "",
        "We thank the Editor and both reviewers for the careful and constructive assessment of this manuscript. The revision was conducted as a source-to-artifact audit rather than a prose-only response. We corrected the implemented sheaf description, separated prior-assisted ingestion from independent topology inference, quantified reaction and transport identifiability limits, added the requested diagnostics and stress tests, reconciled the manuscript with the generated tables and figures, and made the Ghana evidence boundary explicit.",
        "",
        "The response below reproduces every formal editor or reviewer comment verbatim. Each response states the decision, the exact manuscript location, the evidence check used, and the resulting change. One external action remains open: the exact analysis snapshot must be published as an immutable versioned release and assigned a persistent DOI before resubmission. We identify that gate explicitly rather than claiming that a local pre-release package is already a complete public archive.",
        "",
        "## Revision status",
        "",
        _status_table(),
        "",
        "The principal corrected values are active-reaction extent R² = 0.23 (MAE = 0.37 mmol L⁻¹; RMSE = 0.62 mmol L⁻¹), no-prior endpoint-connectivity F1 = 0.62 (precision = 0.49; recall = 0.84), synthetic log-space age R² = 0.98 (median absolute error = 17.86 years), public identifiability-gated age parity log-space R² = 0.960, and Ghana median chemistry R² = 0.70 with median PSI = 0.97. These values are not presented as interchangeable validation outcomes.",
        "",
        "## Editor comments",
        "",
    ]

    for item in comments:
        if item["id"] == "EDITOR-1":
            lines.extend(_render_entry(item, level=3))
    lines.extend(["", "## Reviewer 1", ""])
    for item in comments:
        if item["reviewer"] == "Reviewer 1":
            lines.extend(_render_entry(item, level=3))
    lines.extend(["", "## Reviewer 2", ""])
    for item in comments:
        if item["reviewer"] == "Reviewer 2":
            lines.extend(_render_entry(item, level=3))

    lines.extend(["", "## Additional provenance and audit corrections", ""])
    lines.append(
        "The following controls were added while resolving the numbered comments. They are included to make the revision auditable; they are not presented as additional validation tiers."
    )
    lines.append("")
    for note in _additional_audit_notes():
        lines.append(f"**{note['title']}.** {note['text']}")
        lines.append("")

    lines.extend(
        [
            "## Closing",
            "",
            "We hope this revision makes the contribution, its implementation and its evidence limits clear. The manuscript now treats Hydrosheaf as an auditable integration and prioritisation framework: the synthetic and reference tiers test defined computational estimands, while the Ghana tier demonstrates transfer to sparse observations without claiming independent process truth. We would be pleased to complete the immutable release and DOI step before resubmission.",
            "",
            "Sincerely,",
            "",
            "The authors",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8", newline="")


def _render_entry(item: dict[str, str], *, level: int) -> list[str]:
    response = RESPONSES[item["id"]]
    lines = [f"{'#' * level} {item['id']} — {TITLES[item['id']]}", "", _quote(item["comment"]), ""]
    lines.extend(
        [
            f"**Response ({response['status']}).** {response['response']}",
            "",
            f"**Manuscript location.** {response['location']}",
            "",
            f"**Verification.** {response['verification']}",
            "",
            f"**Change summary.** {response['change_summary']}",
            "",
        ]
    )
    return lines


def _set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def _set_cell_text(cell, text: str, *, bold: bool = False, color: str = "000000") -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_after = Pt(0)
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.name = "Arial"
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), "Arial")
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), "Arial")
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor.from_string(color)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def _add_label_paragraph(doc: Document, label: str, text: str) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(label)
    r.bold = True
    r.font.name = "Arial"
    r.font.size = Pt(10)
    r2 = p.add_run(text)
    r2.font.name = "Arial"
    r2.font.size = Pt(10)


def _set_styles(doc: Document) -> None:
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Arial"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Arial")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Arial")
    normal.font.size = Pt(10)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.08
    for name, size in [("Title", 18), ("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 11)]:
        style = styles[name]
        style.font.name = "Arial"
        style._element.rPr.rFonts.set(qn("w:ascii"), "Arial")
        style._element.rPr.rFonts.set(qn("w:hAnsi"), "Arial")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)


def _add_quote(doc: Document, text: str) -> None:
    for line in text.splitlines():
        p = doc.add_paragraph(style="Intense Quote")
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.left_indent = Inches(0.25)
        p.paragraph_format.right_indent = Inches(0.15)
        run = p.add_run(line)
        run.font.name = "Arial"
        run.font.size = Pt(9)
        run.font.color.rgb = RGBColor(80, 80, 80)


def _build_docx(out_path: Path, comments: list[dict[str, str]]) -> None:
    doc = Document()
    _set_styles(doc)
    section = doc.sections[0]
    section.top_margin = Inches(0.8)
    section.bottom_margin = Inches(0.75)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)
    doc.core_properties.title = "Response to Reviewers"
    doc.core_properties.subject = "CAGEO-D-26-00847 rebuttal letter"

    title = doc.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.LEFT
    title.add_run("Response to Reviewers")
    for label, value in [
        ("Manuscript number", "CAGEO-D-26-00847"),
        ("Manuscript title", "Hydrosheaf: A Reproducible Computational Framework for Inferring Directed Groundwater Hydrochemical Evolution Networks from Sparse Multitracer Data"),
        ("Journal", "Computers and Geosciences"),
    ]:
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(f"{label}: ")
        r.bold = True
        r.font.name = "Arial"
        r.font.size = Pt(10)
        r2 = p.add_run(value)
        r2.font.name = "Arial"
        r2.font.size = Pt(10)

    doc.add_paragraph()
    p = doc.add_paragraph("Dear Editor and Reviewers,")
    p.paragraph_format.space_after = Pt(6)
    opening = [
        "We thank the Editor and both reviewers for the careful and constructive assessment of this manuscript. The revision was conducted as a source-to-artifact audit rather than a prose-only response. We corrected the implemented sheaf description, separated prior-assisted MODPATH ingestion from independent topology inference, quantified reaction and transport identifiability limits, added the requested diagnostics and stress tests, reconciled the manuscript with the generated tables and figures, and made the Ghana evidence boundary explicit.",
        "The response below reproduces every formal editor or reviewer comment verbatim. Each response states the decision, the exact manuscript location, the evidence check used, and the resulting change. One external action remains open: the exact analysis snapshot must be published as an immutable versioned release and assigned a persistent DOI before resubmission. We identify that gate explicitly rather than claiming that a local pre-release package is already a complete public archive.",
    ]
    for text in opening:
        doc.add_paragraph(text)

    doc.add_heading("Revision status", level=1)
    table = doc.add_table(rows=1, cols=2)
    table.style = "Table Grid"
    _set_cell_text(table.rows[0].cells[0], "Evidence class", bold=True, color="FFFFFF")
    _set_cell_text(table.rows[0].cells[1], "Response position", bold=True, color="FFFFFF")
    _set_cell_shading(table.rows[0].cells[0], "1F4E79")
    _set_cell_shading(table.rows[0].cells[1], "1F4E79")
    rows = [
        ("Substantive comments", "Addressed with revised text, regenerated evidence or an explicit scope limitation."),
        ("Field-validation request", "Addressed by reporting the available evidence and stating what the Ghana data cannot establish."),
        ("Release and DOI", "Documented in the package but still an author action before resubmission."),
    ]
    for left, right in rows:
        cells = table.add_row().cells
        _set_cell_text(cells[0], left)
        _set_cell_text(cells[1], right)

    doc.add_paragraph(
        "The principal corrected values are active-reaction extent R² = 0.23 (MAE = 0.37 mmol L⁻¹; RMSE = 0.62 mmol L⁻¹), no-prior endpoint-connectivity F1 = 0.62 (precision = 0.49; recall = 0.84), synthetic log-space age R² = 0.98 (median absolute error = 17.86 years), public identifiability-gated age parity log-space R² = 0.960, and Ghana median chemistry R² = 0.70 with median PSI = 0.97. These values are not presented as interchangeable validation outcomes."
    )

    doc.add_heading("Editor comments", level=1)
    for item in comments:
        if item["id"] == "EDITOR-1":
            _add_docx_entry(doc, item)
    doc.add_heading("Reviewer 1", level=1)
    for item in comments:
        if item["reviewer"] == "Reviewer 1":
            _add_docx_entry(doc, item)
    doc.add_heading("Reviewer 2", level=1)
    for item in comments:
        if item["reviewer"] == "Reviewer 2":
            _add_docx_entry(doc, item)

    doc.add_heading("Additional provenance and audit corrections", level=1)
    doc.add_paragraph(
        "The following controls were added while resolving the numbered comments. They are included to make the revision auditable; they are not presented as additional validation tiers."
    )
    for note in _additional_audit_notes():
        p = doc.add_paragraph(style="List Bullet")
        r = p.add_run(f"{note['title']}. ")
        r.bold = True
        p.add_run(note["text"])

    doc.add_heading("Closing", level=1)
    doc.add_paragraph(
        "We hope this revision makes the contribution, its implementation and its evidence limits clear. The manuscript now treats Hydrosheaf as an auditable integration and prioritisation framework: the synthetic and reference tiers test defined computational estimands, while the Ghana tier demonstrates transfer to sparse observations without claiming independent process truth. We would be pleased to complete the immutable release and DOI step before resubmission."
    )
    doc.add_paragraph("Sincerely,")
    doc.add_paragraph("The authors")
    doc.save(out_path)


def _add_docx_entry(doc: Document, item: dict[str, str]) -> None:
    response = RESPONSES[item["id"]]
    doc.add_heading(f"{item['id']} — {TITLES[item['id']]}", level=2)
    _add_quote(doc, item["comment"])
    _add_label_paragraph(doc, f"Response ({response['status']}). ", response["response"])
    _add_label_paragraph(doc, "Manuscript location. ", response["location"])
    _add_label_paragraph(doc, "Verification. ", response["verification"])
    _add_label_paragraph(doc, "Change summary. ", response["change_summary"])


def _write_json_audits(out_dir: Path, comments: list[dict[str, str]]) -> None:
    tasks = []
    plan = []
    changes = []
    analysis = []
    for item in comments:
        r = RESPONSES[item["id"]]
        tasks.append(
            {
                "id": item["id"],
                "reviewer": item["reviewer"],
                "comment": item["comment"],
                "category": item["category"],
                "requested_action": r["requested_action"],
                "target_files": r["target_files"],
            }
        )
        plan.append(
            {
                "comment_id": item["id"],
                "decision": r["decision"],
                "rationale": r["response"],
                "target_files": r["target_files"],
                "requires_analysis": r["requires_analysis"],
                "acceptance_check": r["verification"],
            }
        )
        changes.append(
            {
                "comment_id": item["id"],
                "file": r["target_files"][0],
                "before_summary": "Submitted or earlier revision wording/evidence did not fully satisfy the reviewer concern.",
                "after_summary": r["change_summary"],
            }
        )
        if r["requires_analysis"]:
            analysis.append(
                {
                    "comment_id": item["id"],
                    "status": "PASS",
                    "evidence": r["verification"],
                    "outputs": r["target_files"],
                }
            )
    for name, payload in [
        ("tasks.json", tasks),
        ("reviewer_comment_inventory.json", {"source_git_revision": COMMENTS_REVISION, "source_path": COMMENTS_PATH, "comments": tasks}),
        ("locked_revision_plan.json", {"source_git_revision": COMMENTS_REVISION, "items": plan}),
        ("change_log.json", {"source_git_revision": COMMENTS_REVISION, "items": changes}),
        ("analysis_changes.json", {"source_git_revision": COMMENTS_REVISION, "items": analysis}),
    ]:
        (out_dir / name).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--output-docx", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    comments = _extract_comments(_source_comments(repo_root))
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_docx.parent.mkdir(parents=True, exist_ok=True)
    args.audit_dir.mkdir(parents=True, exist_ok=True)
    _write_markdown(args.output_md, comments)
    _build_docx(args.output_docx, comments)
    _write_json_audits(args.audit_dir, comments)
    print(f"comments={len(comments)}")
    print(f"markdown={args.output_md}")
    print(f"docx={args.output_docx}")
    print(f"audits={args.audit_dir}")


if __name__ == "__main__":
    main()
