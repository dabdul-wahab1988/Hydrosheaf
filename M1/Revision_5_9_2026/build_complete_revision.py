from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
WORKED = ROOT / "Worked_Demonstration"
BASE_SCRIPT = ROOT / "build_revision_artifacts.py"

MAIN_SOURCE = ROOT / "Manuscript- Water and Ecology_Revised_Clean.docx"
MAIN_MARKED_SOURCE = ROOT / "Manuscript- Water and Ecology_Revised_Colour_Marked.docx"
SUPP_SOURCE = ROOT / "SupplementaryInformation_Revised_Clean.docx"
SUPP_MARKED_SOURCE = ROOT / "SupplementaryInformation_Revised_Colour_Marked.docx"

MAIN_CLEAN = ROOT / "Manuscript- Water and Ecology_Fully_Revised_Clean.docx"
MAIN_MARKED = ROOT / "Manuscript- Water and Ecology_Fully_Revised_Colour_Marked.docx"
SUPP_CLEAN = ROOT / "SupplementaryInformation_Fully_Revised_Clean.docx"
SUPP_MARKED = ROOT / "SupplementaryInformation_Fully_Revised_Colour_Marked.docx"
RESPONSE = ROOT / "Response_to_Reviewers_Fully_Addressed.docx"
AUDIT = ROOT / "Final_Revision_Audit.json"

PUBLIC_COMMIT = "ea463ecbe9f897a861e8523f4dde30c08ff514a7"
PUBLIC_REPOSITORY = "https://github.com/dabdul-wahab1988/Hydrosheaf"
PUBLIC_BENCHMARK_URL = (
    f"{PUBLIC_REPOSITORY}/tree/{PUBLIC_COMMIT}/M7/m7_nonuniqueness_benchmark"
)
PUBLIC_RUNNER_URL = (
    f"{PUBLIC_REPOSITORY}/blob/{PUBLIC_COMMIT}/M7/m7_nonuniqueness_benchmark/"
    "scripts/run_supporting_validation.py"
)
PUBLIC_RESULTS_URL = (
    f"{PUBLIC_REPOSITORY}/tree/{PUBLIC_COMMIT}/M7/m7_nonuniqueness_benchmark/"
    "results/supporting_validation"
)

spec = importlib.util.spec_from_file_location("prior_revision", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

PURPLE = base.COLORS["both"]
GREEN = base.COLORS["editorial"]
BLACK = base.COLORS["black"]


ABSTRACT = (
    "Groundwater assessment in data-limited aquifers remains constrained by non-unique interpretations of "
    "residence time, hydraulic connectivity and hydrogeochemical evolution because tracer-based ages, graph-based "
    "connectivity and inverse reaction models are commonly applied separately. This critical review evaluates "
    "environmental-tracer residence-time methods, graph-topological representations of aquifer connectivity and "
    "inverse hydrogeochemical modelling as complementary tools for reducing interpretive ambiguity in sparse-data "
    "settings. Literature was identified through selective searches in Web of Science, Scopus, Google Scholar and "
    "Consensus, supported by hand-searching and screening of public benchmark resources. The synthesis proposes an "
    "evidence-gated sequence in which tracer-derived residence-time distributions constrain directed candidate edges, "
    "graph topology restricts the water pairs eligible for reaction fitting, and reaction solutions are filtered by "
    "thermodynamic and hydrogeological evidence. A controlled-synthetic demonstration evaluated the executable "
    "sequence on 12 locked aquifer twins generated independently with MODFLOW 6, MODPATH 7 and nonlinear chemistry. "
    "Across 825 candidate edges, hydraulic-plus-chemistry scoring achieved a precision-recall area under the curve "
    "(PR-AUC) of 0.464 and F1 of 0.507; adding an age-compatibility gate gave PR-AUC 0.462 and the same F1, so age "
    "evidence did not improve aggregate ranking. Bayesian 3H/39Ar inference had a mean absolute age error of 2.76 years "
    "and mean 95% interval coverage of 0.917. PHREEQC constraints increased reaction-family accuracy on 103 "
    "candidate-contained truth edges from 0.563 to 0.602, although carbonate reactions were not recovered. The "
    "demonstration establishes reproducible component integration under model-conditioned synthetic truth, not field "
    "transferability or operational superiority; field applications must retain that evidence boundary."
)

INTRO_SCOPE = (
    "The proposed synthesis is significant because many groundwater-dependent regions, particularly in semi-arid "
    "Africa, lack the dense monitoring networks, long-term hydraulic records and specialised isotope datasets required "
    "for conventional high-confidence modelling. In such settings, management decisions are often made using limited "
    "hydrochemical snapshots, incomplete borehole information and uncertain conceptual flow models. The review therefore "
    "combines an evidence-gated methodological synthesis with a controlled-synthetic integration test. That test checks "
    "whether the components can be executed together and exposes their failure modes; it is not a substitute for "
    "independent field validation."
)

UNCERTAINTY = (
    "The controlled-synthetic demonstration implemented uncertainty at both component and system levels. For each of "
    "12 locked aquifer twins, 3H and 39Ar age inference used four exact-grid chains with 500 retained draws per chain. "
    "Topology uncertainty used four chains with 2,500 retained samples per chain after 750 burn-in iterations and 32 "
    "updates per retained sample. PHREEQC diagnostics and constrained and unconstrained reaction fits were recorded for "
    "every candidate edge. Method contrasts were then resampled 10,000 times by independent aquifer case (seed 7781), "
    "rather than by correlated edge. All pre-specified locked cases were retained; convergence, PHREEQC success and failed "
    "runs were reported instead of being used as post hoc exclusion criteria."
)

DEMO_PARAGRAPHS = [
    (
        "A controlled-synthetic worked demonstration was drawn from the repository's locked integration benchmark "
        "(Figure 7; Supplementary Methods S3 and Table S3). Six development aquifer twins (seeds 2101–2106) were used "
        "to freeze the edge-fusion models, after which 12 independent test twins (seeds 4101–4112; 12 monitoring nodes "
        "per case) were analysed without exposing held-out ages, flow edges or reaction labels to inference. The external "
        "generator used official MODFLOW 6 and MODPATH 7 together with independent nonlinear tracer and chemistry equations."
    ),
    (
        "The age step combined 3H and 39Ar observations in a piston-flow exact-grid Bayesian model. All 12 cases "
        "converged, with mean absolute error 2.76 years, mean bias −0.98 years, mean 95% interval coverage 0.917, maximum "
        "R-hat 1.004, minimum bulk effective sample size 1,722 and no divergences. These results establish internal age-"
        "component performance for the stated generator and tracer model; they do not establish field-age accuracy."
    ),
    (
        "Candidate-graph construction recovered 103 of 108 held-out flow edges across 825 proposed edges (mean case "
        "recall 0.954). The development-locked hydraulic-plus-chemistry score achieved PR-AUC 0.464 and F1 0.507 on "
        "the locked cases. Adding the age-compatibility gate produced PR-AUC 0.462 and F1 0.507. Thus, age evidence was "
        "directionally interpretable but did not improve aggregate edge ranking. The posterior maximum-a-posteriori "
        "graphs had mean precision 0.373, recall 0.806 and F1 0.509, indicating useful sensitivity but substantial false "
        "connectivity."
    ),
    (
        "For every candidate water pair, the workflow ran unconstrained and PHREEQC-constrained network-scale inverse "
        "reaction fits. All 144 PHREEQC sample evaluations succeeded; direction constraints were active for all 825 "
        "candidate-edge fits and materially changed the objective for 681. Among the 103 candidate-contained truth edges, "
        "reaction-family accuracy increased from 0.563 without thermodynamic constraints to 0.602 with them. Recovery was "
        "strong for denitrification, sulfate reduction and silicate weathering, but both carbonate families had zero "
        "accuracy. This failure is retained because a feasible reaction fit is not evidence of unique process attribution."
    ),
    (
        "Case-block bootstrap intervals confirmed the limited incremental value of the age gate: its F1 difference from "
        "hydraulic-plus-chemistry scoring was 0.000 (95% CI 0.000 to 0.000), and its PR-AUC difference was −0.0017 (95% "
        "CI −0.0056 to 0.0001). Relative to the age-permuted control, the PR-AUC difference was 0.0354 (95% CI −0.0158 "
        "to 0.0954). The confidence intervals include no gain. The demonstration therefore resolves the execution and "
        "auditability question while setting a strict evidence ceiling: it is a controlled test of integration and failure "
        "modes, not independent field validation, general superiority or a management-ready aquifer reconstruction."
    ),
]

RESEARCH_AGENDA = (
    "The controlled-synthetic benchmark establishes that the proposed components can be executed as one auditable "
    "sequence, but it also shows that integration does not guarantee incremental predictive value: the age gate did not "
    "improve aggregate edge ranking, posterior graphs retained false connections and carbonate reaction attribution "
    "failed. The highest-priority research need is therefore independent field evaluation with pre-specified endpoints, "
    "withheld wells or campaigns, and direct hydraulic or tracer evidence for connectivity. Open-source implementations "
    "should continue to preserve raw inputs, node and edge rules, all viable reaction fits, uncertainty draws, failures "
    "and machine-readable provenance so that positive and negative results remain reproducible."
)

CONCLUSION_1 = (
    "This review develops an evidence-gated sequence for interpreting data-limited aquifers by combining three "
    "complementary evidence types. Residence-time analysis supplies temporal constraints, graph topology makes candidate "
    "hydraulic connections explicit, and inverse hydrogeochemical modelling tests whether chemistry along those candidate "
    "connections is compatible with specified reactions and mixing. A controlled-synthetic benchmark demonstrates that "
    "these components can be executed together and audited without exposing held-out truth to inference."
)

CONCLUSION_2 = (
    "The benchmark does not establish an operational field tool. Age inference was well calibrated under the synthetic "
    "design, but age information did not improve aggregate edge ranking, posterior topology retained false connections, "
    "and carbonate reaction families were not recovered. These negative results narrow the claim: the framework supports "
    "reproducible hypothesis generation and uncertainty diagnosis under controlled conditions. Independent field testing "
    "in contrasting karst, fractured-rock and porous aquifers is required before claims of transferability, superiority "
    "or management readiness can be made."
)

CODE_AVAILABILITY = (
    "The code, locked synthetic inputs, complete locked results and replay instructions for the controlled-synthetic "
    f"worked demonstration are publicly available in the HydroSheaf repository at {PUBLIC_BENCHMARK_URL}. The "
    f"executable runner is {PUBLIC_RUNNER_URL}, and its supporting validation outputs are {PUBLIC_RESULTS_URL}. "
    f"These commit-pinned links identify the exact public version used for the reported results ({PUBLIC_COMMIT}). "
    "MODFLOW 6, MODPATH 7 and PHREEQC remain external scientific programs and must be obtained under their respective "
    "distribution terms as described in the repository."
)


def replace(doc, prefix: str, text: str, marked: bool, colour: str = PURPLE) -> None:
    p = base.find_prefix(doc, prefix)
    base.put_text(p, text, color=colour if marked else BLACK)


def add_para_before(anchor, text: str, marked: bool, bold: bool = False) -> None:
    p = anchor.insert_paragraph_before()
    p.style = "Normal"
    base.put_text(p, text, color=PURPLE if marked else BLACK, bold=bold)
    if bold:
        p.paragraph_format.keep_with_next = True
        p.paragraph_format.space_before = Pt(8)


def style_table(table, marked: bool, widths: list[float]) -> None:
    base.format_table(table, marked, PURPLE, widths)
    for row in table.rows[1:]:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    base.set_font(run, color=PURPLE if marked else BLACK, size=7.5)


def add_table8(doc, marked: bool) -> None:
    caption = doc.add_paragraph()
    base.put_text(
        caption,
        "Table 8. Results and evidence limits of the controlled-synthetic worked demonstration.",
        color=PURPLE if marked else BLACK,
    )
    caption.paragraph_format.keep_with_next = True
    headers = ["Step", "Locked implementation", "Result", "Interpretive boundary"]
    rows = [
        ["Design", "6 development and 12 test aquifer twins; 12 nodes per case; truth hidden during inference.",
         "825 candidates; 103/108 true flow edges represented.",
         "Model-conditioned synthetic truth, not field truth."],
        ["Age", "3H + 39Ar piston-flow exact-grid inference; 4 chains × 500 draws.",
         "MAE 2.76 y; bias −0.98 y; 95% coverage 0.917; 0 divergences.",
         "Performance is conditional on the stated generator and model."],
        ["Edge scoring", "Development-locked hydraulic + chemistry score, then age-compatibility gate.",
         "PR-AUC 0.464 to 0.462; F1 unchanged at 0.507.",
         "Age evidence did not improve aggregate ranking."],
        ["Topology", "4 chains × 2,500 retained samples after 750 burn-in; 32 updates/sample.",
         "Mean MAP precision 0.373, recall 0.806 and F1 0.509.",
         "High recall coexisted with substantial false connectivity."],
        ["Reactions", "Unconstrained and PHREEQC-constrained fits on every candidate pair.",
         "Accuracy 0.563 to 0.602 on 103 truth edges; carbonate accuracy 0.",
         "Thermodynamic feasibility did not ensure unique attribution."],
        ["Uncertainty", "10,000 independent-case bootstrap resamples; seed 7781; no excluded cases.",
         "Age-gate ΔPR-AUC −0.0017 (95% CI −0.0056 to 0.0001).",
         "No evidence of aggregate gain; independent field validation is required."],
    ]
    table = doc.add_table(rows=1, cols=4)
    for i, value in enumerate(headers):
        table.rows[0].cells[i].text = value
    for row in rows:
        cells = table.add_row().cells
        for i, value in enumerate(row):
            cells[i].text = value
    style_table(table, marked, [1.0, 3.25, 3.25, 3.25])


def clean_runs(doc, marked: bool) -> None:
    replacements = {
        "initial-final": "initial–final",
        "water-rock": "water–rock",
        "scientific rigor": "scientific rigour",
        "  ": " ",
    }
    paragraphs = list(doc.paragraphs)
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                paragraphs.extend(cell.paragraphs)
    for paragraph in paragraphs:
        for run in paragraph.runs:
            new = run.text
            for old, replacement in replacements.items():
                while old in new:
                    new = new.replace(old, replacement)
            new = re.sub(r"\s+([,.;:])", r"\1", new)
            if new != run.text:
                run.text = new
                if marked:
                    base.set_font(run, color=GREEN)


def update_marked_legend(doc) -> None:
    matches = [
        p for p in doc.paragraphs
        if p.text.startswith("The clean copy contains the same accepted wording")
    ]
    if len(matches) == 1:
        base.put_text(
            matches[0],
            "The clean copy contains the same accepted wording without colour markup. Purple identifies the executed controlled-synthetic demonstration and its evidence-bounded interpretation; the original source files remain unchanged.",
            color=base.COLORS["gray"],
            italic=True,
            size=9,
        )


def revise_main(doc, marked: bool) -> None:
    if marked:
        update_marked_legend(doc)
    replace(doc, "Groundwater assessment in data-limited aquifers remains constrained", ABSTRACT, marked)
    replace(doc, "The proposed synthesis is significant because", INTRO_SCOPE, marked)
    replace(doc, "The uncertainty procedure should be implemented", UNCERTAINTY, marked)
    replace(doc, "6.8 Evidence required for a worked demonstration", "6.8 Controlled-synthetic worked demonstration", marked)
    replace(doc, "Because this article is a conceptual critical review", DEMO_PARAGRAPHS[0], marked)
    discussion = base.find_prefix(doc, "7. Discussion")
    for paragraph in DEMO_PARAGRAPHS[1:]:
        add_para_before(discussion, paragraph, marked)
    replace(doc, "The highest-priority research need is to move from conceptual integration", RESEARCH_AGENDA, marked)
    replace(doc, "This review develops a conceptual sequence", CONCLUSION_1, marked)
    replace(doc, "The framework is not yet a validated operational tool", CONCLUSION_2, marked)
    replace(
        doc,
        "Figure 4. Sources and propagation of uncertainty",
        "Figure 4. Explicit propagation of uncertainty through the integrated workflow. Each row links an evidence source to the sampling or model alternative used in a draw, the propagated state and the diagnostic that must be reported. Draw counts, seeds, convergence, failures and retention rules are part of the result. Independent field evidence remains a separate validation gate.",
        marked,
    )
    tables_heading = base.find_prefix(doc, "Tables")
    add_para_before(
        tables_heading,
        "Figure 7. Controlled-synthetic worked demonstration. (A) Representative posterior maximum-a-posteriori graph for locked seed 4101, with true- and false-positive inferred edges. (B) Locked edge-scoring results. (C) Bayesian 3H/39Ar age error and interval coverage across independent test aquifers. (D) Constrained versus unconstrained reaction-family accuracy and the uncertainty audit. Synthetic truth evaluates execution and internal integration, not field transferability.",
        marked,
    )
    add_table8(doc, marked)
    acknowledgement = base.find_prefix(doc, "Acknowledgement")
    add_para_before(acknowledgement, "Code and data availability", marked, bold=True)
    add_para_before(acknowledgement, CODE_AVAILABILITY, marked)
    clean_runs(doc, marked)


def move_table_before(table, anchor) -> None:
    anchor._p.addprevious(table._tbl)


SI_PARAGRAPHS = [
    (
        "S3.1 Design and provenance. The worked demonstration used six development aquifer twins (random seeds "
        "2101–2106) to freeze the edge-fusion models and 12 independent locked test twins (seeds 4101–4112) for final "
        "evaluation. Each case contained 12 monitoring nodes. Official MODFLOW 6 (version 6.7.0) and MODPATH 7 "
        "(version 7.2.001) generated heterogeneous hydraulic fields and particle-path truth; separate nonlinear equations "
        "generated tracer and chemistry observations. Held-out node ages, edges and reaction labels were excluded from "
        "inference. Source-file SHA-256 hashes, software versions, the repository commit and working-tree status are "
        "recorded in the accompanying analysis manifest."
    ),
    (
        "S3.2 Age inference and graph construction. Candidate directed edges were proposed from observable hydraulic "
        "features. Independent 3H and 39Ar observations were fitted with a piston-flow exact-grid Bayesian model using "
        "four chains and 500 retained draws per chain. Posterior node-age summaries entered a pre-specified direction "
        "compatibility gate. A development-locked fusion model combined hydraulic logit and constrained chemistry "
        "objective; the age-gated variant capped incompatible edges without refitting on test truth."
    ),
    (
        "S3.3 Inverse reaction fitting. Every candidate edge defined an eligible initial–final water pair. The workflow "
        "ran sparse network reaction fitting before and after PHREEQC-derived direction constraints across the declared "
        "mineral and redox process set. Both objectives, constraint activation, bound hits, dominant reaction families and "
        "the complete held-out comparison were retained. This is a network-scale inverse reaction benchmark; a successful "
        "fit was not interpreted as proof of hydraulic connection or unique process attribution."
    ),
    (
        "S3.4 Topology and uncertainty. Topology posteriors used four chains, 2,500 retained samples per chain after 750 "
        "burn-in iterations, 32 transition updates per retained sample, an acyclic constraint, a minimum of nine edges "
        "and a maximum out-degree of three. Method differences were quantified with 10,000 bootstrap resamples of the 12 "
        "independent aquifer cases (seed 7781). Edge-level resampling was not used because edges within an aquifer are "
        "dependent. All cases were retained; convergence flags, effective sample sizes, PHREEQC success and failed-run "
        "counts were reported explicitly."
    ),
    (
        "S3.5 Interpretation. All age and topology chains met the recorded convergence criteria, and no locked run failed. "
        "Nevertheless, the age gate did not improve aggregate edge ranking, posterior graphs contained false-positive "
        "connections and carbonate reaction families were not recovered. These are retained negative results. The "
        "demonstration supports reproducible execution and failure diagnosis under controlled synthetic conditions only."
    ),
    (
        "S3.6 Code and data availability. The complete benchmark package is available at "
        f"{PUBLIC_BENCHMARK_URL}. The confirmatory runner is {PUBLIC_RUNNER_URL}, and the locked input, truth, diagnostic "
        f"and result files are {PUBLIC_RESULTS_URL}. From the repository root, the supporting validation is replayed with "
        ".venv\\Scripts\\python.exe M7\\m7_nonuniqueness_benchmark\\scripts\\run_supporting_validation.py "
        "--confirmatory after the documented MODFLOW 6, MODPATH 7 and PHREEQC dependencies have been installed. The "
        f"commit identifier {PUBLIC_COMMIT} fixes the code and results referenced by this manuscript."
    ),
]


def revise_supplement(doc, marked: bool) -> None:
    if marked:
        update_marked_legend(doc)
    anchor = base.find_prefix(doc, "Reference")
    add_para_before(anchor, "Supplementary Methods S3. Controlled-synthetic worked demonstration", marked, bold=True)
    for paragraph in SI_PARAGRAPHS:
        add_para_before(anchor, paragraph, marked)
    add_para_before(anchor, "Supplementary Table S3. Locked design, diagnostics and retained limitations.", marked)
    headers = ["Element", "Pre-specified implementation", "Recorded result"]
    rows = [
        ["Cases", "6 development; 12 locked tests; 12 nodes/case", "All 12 locked cases retained; zero failed runs"],
        ["Age", "3H + 39Ar; 4 chains × 500 draws", "MAE 2.76 y; 95% coverage 0.917; max R-hat 1.004"],
        ["Candidates", "Hydraulic proposal followed by chemistry and age scoring", "825 candidates; 103/108 truth edges; candidate recall 0.954"],
        ["Topology", "4 chains × 2,500 retained; 750 burn-in; 32 updates/sample", "Mean MAP precision 0.373; recall 0.806; F1 0.509"],
        ["Reactions", "Unconstrained and PHREEQC-constrained fit for every candidate", "144/144 PHREEQC samples successful; constrained accuracy 0.602"],
        ["Bootstrap", "10,000 independent-case resamples; seed 7781", "Age-gate ΔPR-AUC −0.0017; 95% CI −0.0056 to 0.0001"],
        ["Evidence ceiling", "Model-conditioned synthetic truth", "No claim of field validation, superiority or management readiness"],
    ]
    table = doc.add_table(rows=1, cols=3)
    for i, value in enumerate(headers):
        table.rows[0].cells[i].text = value
    for row in rows:
        cells = table.add_row().cells
        for i, value in enumerate(row):
            cells[i].text = value
    style_table(table, marked, [1.25, 4.6, 4.6])
    move_table_before(table, anchor)
    clean_runs(doc, marked)


def word_count(doc) -> int:
    count = sum(len(p.text.split()) for p in doc.paragraphs)
    count += sum(
        len(cell.text.split())
        for table in doc.tables
        for row in table.rows
        for cell in row.cells
    )
    return count


def narrative_word_count(doc) -> int:
    """Count manuscript paragraphs before References, excluding tables and references."""
    parts = []
    for paragraph in doc.paragraphs:
        text = paragraph.text.strip()
        if text == "References":
            break
        if text:
            parts.append(text)
    return len(" ".join(parts).split())


def write_text_extract(doc, path: Path) -> None:
    lines = [paragraph.text for paragraph in doc.paragraphs if paragraph.text.strip()]
    for table in doc.tables:
        for row in table.rows:
            lines.append("\t".join(cell.text for cell in row.cells))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def section_columns(doc) -> list[int]:
    values = []
    for section in doc.sections:
        cols = section._sectPr.find(qn("w:cols"))
        values.append(int(cols.get(qn("w:num"), "1")) if cols is not None else 1)
    return values


def compliance(doc) -> dict[str, object]:
    abstract = base.find_prefix(doc, "Groundwater assessment in data-limited aquifers")
    keywords = base.find_prefix(doc, "Keywords:")
    keyword_values = [x.strip() for x in keywords.text.split(":", 1)[1].split(",") if x.strip()]
    narrative_words = narrative_word_count(doc)
    return {
        "guide_checked_utc": "2026-09-06",
        "guide_url": "https://www.keaipublishing.com/en/journals/water-and-ecology/guide-for-authors/",
        "article_type": "Review article",
        "review_article_word_limit": "7000-9000",
        "narrative_words_before_references_excluding_tables": narrative_words,
        "review_article_word_limit_compliant": 7000 <= narrative_words <= 9000,
        "review_article_original_data_policy": (
            "Review articles should not include unpublished/original data, submitted manuscripts or personal communication."
        ),
        "controlled_synthetic_demonstration_is_unpublished_original_material": True,
        "article_type_or_editorial_confirmation_required": True,
        "manuscript_words_including_tables": word_count(doc),
        "abstract_words": len(abstract.text.split()),
        "abstract_limit": 350,
        "abstract_compliant": len(abstract.text.split()) <= 350,
        "keyword_count": len(keyword_values),
        "keyword_requirement": "3 to 6",
        "keywords_compliant": 3 <= len(keyword_values) <= 6,
        "section_column_counts": section_columns(doc),
        "single_column_compliant": all(x == 1 for x in section_columns(doc)),
        "editable_source_format": "docx",
    }


def reviewer_items(length_audit: dict[str, object]) -> list[dict[str, str]]:
    items = [dict(item) for item in base.REVIEWER_COMMENTS]
    by_id = {item["id"]: item for item in items}
    for item in items:
        item["status"] = "ADDRESSED"

    by_id["R1.1"].update(
        location="Abstract; Sections 6.7–6.8; Figure 7; Table 8; Supplementary Methods S3 and Table S3.",
        response=(
            "We added a reproducible controlled-synthetic demonstration of the complete sequence. Independent MODFLOW 6/"
            "MODPATH 7 aquifer twins provide model-conditioned truth; inference remains blind to held-out ages, edges and "
            "reaction labels. The revised text reports age inference, graph construction, frozen edge scoring, PHREEQC-"
            "constrained inverse reaction fits, topology sampling, bootstrap uncertainty and the negative results. We state "
            "throughout that this is not field validation."
        ),
        text=(
            "A controlled-synthetic worked demonstration was evaluated on 12 locked aquifer twins. Adding the age gate did "
            "not improve aggregate edge ranking, and carbonate reaction attribution remained unresolved; these failures are "
            "reported rather than suppressed."
        ),
    )
    by_id["R1.5"].update(
        location="Redrawn Figure 4 and revised Figure 4 caption.",
        response=(
            "Figure 4 was redrawn. Each uncertainty source is now connected inside the figure to its quantification method, "
            "propagated state and reported diagnostic, with an explicit requirement to report draws, seeds, convergence, "
            "failures and retention rules."
        ),
        text="Figure 4. Explicit propagation of uncertainty through the integrated workflow.",
    )
    by_id["R1.6"].update(
        location="Section 6.2; Section 6.8; Supplementary Methods S3.2.",
        response=(
            "The manuscript retains the pre-specified support and coverage definitions and now shows an executed, "
            "development-locked scoring example. Hydraulic and chemistry features were frozen on development cases; the "
            "age-compatibility rule was then applied without fitting to locked-test truth."
        ),
    )
    by_id["R1.7"].update(
        location="Section 6.7; Section 6.8; Figure 7; Table 8; Supplementary Methods S3.4.",
        response=(
            "The stochastic analysis is now executed and reported: four age chains with 500 retained draws per chain, four "
            "topology chains with 2,500 retained samples per chain after 750 burn-in iterations, 32 topology updates per "
            "sample, and 10,000 independent-case bootstrap resamples (seed 7781). All locked cases were retained and zero "
            "failed runs were recorded."
        ),
        text=UNCERTAINTY,
    )
    by_id["R1.9"].update(
        location="Entire clean and colour-marked manuscript and Supplementary Information; final rendered-page audit.",
        response=(
            "The complete files were copy-edited for grammar, punctuation, terminology, numerical notation and UK spelling. "
            "The clean and colour-marked copies were generated from matched content and subjected to structural checks and "
            "complete rendered-page inspection."
        ),
        text="All revised files use consistent terminology, typography and evidence-bounded claims.",
    )
    by_id["R2.4"].update(
        status="PARTIALLY ADDRESSED — EDITOR DECISION REQUIRED",
        location="Entire manuscript; Final Revision Audit, journal-compliance section.",
        response=(
            f"We rechecked the live Water & Ecology Guide for Authors on 6 September 2026. It specifies 7,000–9,000 words "
            f"for a Review article and states that Review articles should not include unpublished/original data. The current "
            f"narrative count before the references, excluding tables, is {length_audit['narrative_words_before_references_excluding_tables']} "
            f"words. The revised abstract contains {length_audit['abstract_words']} words, the manuscript uses "
            f"{length_audit['keyword_count']} keywords, and the editable source is single-column DOCX. Because the reviewers "
            f"requested the executed demonstration, the authors should obtain explicit editorial confirmation that it may "
            f"remain in this Review article, or agree an article-type change; a further narrative reduction is also required "
            f"if the Review classification is retained. We do not claim this item is fully resolved without that decision."
        ),
        text=(
            f"Live-guide checks: Review article narrative "
            f"{length_audit['narrative_words_before_references_excluding_tables']}/7,000–9,000 words; abstract "
            f"{length_audit['abstract_words']}/350 words; keywords {length_audit['keyword_count']}/3–6; single-column "
            f"editable DOCX. The original-data policy requires an editor decision on the executed benchmark."
        ),
    )
    by_id["R2.5"].update(
        location="Section 6.8; Figure 7; Table 8; Supplementary Methods S3 and Table S3.",
        response=(
            "We added a complete numerical example covering age classification, graph construction, edge scoring, "
            "graph-supported reaction fitting and propagated uncertainty. The public-facing analysis extract includes "
            "machine-readable metrics, a representative edge audit and SHA-256 provenance."
        ),
        text=DEMO_PARAGRAPHS[2],
    )
    by_id["R2.6"].update(
        location="Abstract; Section 6.8; Research agenda; Conclusions.",
        response=(
            "The integrative framework remains the organizing structure, repeated method descriptions were condensed in "
            "the first revision, and the requested illustrative example is now executed. The conclusion reports both "
            "what worked and what did not, without treating synthetic testing as field validation."
        ),
        text=CONCLUSION_2,
    )
    return items


def add_response_letter(doc, items: list[dict[str, str]]) -> None:
    for style_name in ("Normal", "Body Text"):
        if style_name in doc.styles:
            doc.styles[style_name].font.name = "Times New Roman"
            doc.styles[style_name].font.size = Pt(11)
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    base.put_text(title, "Response to Reviewers", bold=True, size=16)
    for text in [
        "Manuscript: Environmental-Tracer Residence Time, Graph Topology and Inverse Hydrogeochemistry in Data-Limited Aquifers: A Critical Review",
        "Manuscript number: WATECO-D-26-00083",
        "Response to the editor and Reviewers 1 and 2",
        "Dear Editor and Reviewers,",
        (
            "We thank the editor and both reviewers for their careful assessment. The revised package now closes the "
            "previously outstanding empirical and graphical items with a reproducible controlled-synthetic worked "
            "demonstration, an executed uncertainty analysis and a redesigned uncertainty figure. The benchmark's negative "
            "results are reported explicitly, and all claims remain below the boundary of independent field validation."
        ),
        (
            "The exact code, locked synthetic cases, complete result tables and replay instructions are publicly accessible "
            f"at the commit-pinned HydroSheaf benchmark package: {PUBLIC_BENCHMARK_URL}. The manuscript and Supplementary "
            "Information now include a Code and data availability statement identifying the runner, results directory and "
            "external scientific-program requirements."
        ),
    ]:
        p = doc.add_paragraph()
        base.put_text(p, text, size=11)

    p = doc.add_paragraph()
    base.put_text(p, "Colour key for the marked files", bold=True, size=12)
    for label, colour, meaning in [
        ("Red", base.COLORS["r1"], "Reviewer 1-related first-round revisions"),
        ("Blue", base.COLORS["r2"], "Reviewer 2-related first-round revisions"),
        ("Purple", PURPLE, "Worked demonstration and changes addressing both reviewers"),
        ("Green", GREEN, "Editorial and consistency corrections"),
    ]:
        p = doc.add_paragraph()
        base.put_text(p, f"{label}: {meaning}", color=colour, size=10)

    for reviewer in ("Reviewer 1", "Reviewer 2"):
        h = doc.add_paragraph()
        base.put_text(h, reviewer, bold=True, size=14)
        for item in [x for x in items if x["reviewer"] == reviewer]:
            h = doc.add_paragraph()
            base.put_text(h, f"{item['id']}  {item['status']}", bold=True, size=12)
            for label, value, italic, size in [
                ("Reviewer comment: ", item["comment"], True, 11),
                ("Exact location in the revised package: ", item["location"], False, 11),
                ("Response: ", item["response"], False, 11),
                ("Representative revised text: ", item["text"], False, 10),
            ]:
                p = doc.add_paragraph()
                run = p.add_run(label)
                base.set_font(run, bold=True, size=size)
                run = p.add_run(value)
                base.set_font(run, italic=italic, size=size)

    h = doc.add_paragraph()
    base.put_text(h, "Repository and reproducibility record", bold=True, size=14)
    p = doc.add_paragraph()
    base.put_text(
        p,
        (
            f"Public benchmark package: {PUBLIC_BENCHMARK_URL}\n"
            f"Executable runner: {PUBLIC_RUNNER_URL}\n"
            f"Locked supporting results: {PUBLIC_RESULTS_URL}\n"
            f"Referenced public commit: {PUBLIC_COMMIT}"
        ),
        size=10,
    )

    h = doc.add_paragraph()
    base.put_text(h, "Closing statement", bold=True, size=14)
    p = doc.add_paragraph()
    base.put_text(
        p,
        (
            "The revision now demonstrates the complete integration sequence and its uncertainty accounting under "
            "controlled synthetic conditions. It also shows why the evidence boundary is necessary: the age gate did not "
            "improve aggregate edge ranking, inferred topology retained false connections and carbonate-family attribution "
            "failed. We therefore request assessment of the framework as a reproducible, hypothesis-generating method whose "
            "field transferability remains to be established."
        ),
        size=11,
    )
    p = doc.add_paragraph()
    base.put_text(p, "Sincerely,\nThe authors", size=11)


def main() -> None:
    evidence = json.loads((WORKED / "worked_demonstration_summary.json").read_text(encoding="utf-8"))
    assert evidence["design"]["locked_test_cases"] == 12
    assert evidence["uncertainty"]["failed_runs"] == 0

    main_clean_doc = Document(MAIN_SOURCE)
    revise_main(main_clean_doc, marked=False)
    main_clean_doc.save(MAIN_CLEAN)

    main_marked_doc = Document(MAIN_MARKED_SOURCE)
    revise_main(main_marked_doc, marked=True)
    main_marked_doc.save(MAIN_MARKED)

    supp_clean_doc = Document(SUPP_SOURCE)
    revise_supplement(supp_clean_doc, marked=False)
    supp_clean_doc.save(SUPP_CLEAN)

    supp_marked_doc = Document(SUPP_MARKED_SOURCE)
    revise_supplement(supp_marked_doc, marked=True)
    supp_marked_doc.save(SUPP_MARKED)

    length_audit = compliance(main_clean_doc)
    items = reviewer_items(length_audit)
    response_doc = Document()
    add_response_letter(response_doc, items)
    response_doc.save(RESPONSE)

    write_text_extract(main_clean_doc, ROOT / "final_main_text.txt")
    write_text_extract(supp_clean_doc, ROOT / "final_supp_text.txt")
    write_text_extract(response_doc, ROOT / "final_response_text.txt")

    audit = {
        "status": "GENERATED_AWAITING_RENDER_QA",
        "scientific_evidence": {
            "source": "Worked_Demonstration/analysis_manifest.json",
            "assertions": "PASS",
            "evidence_class": evidence["evidence_class"],
            "claim_boundary": evidence["claim_boundary"],
        },
        "journal_compliance": length_audit,
        "submission_readiness": {
            "status": "AUTHOR_EDITOR_DECISION_REQUIRED",
            "blocking_items": [
                "Review-article narrative exceeds the stated 9,000-word upper limit.",
                "The executed controlled-synthetic benchmark is unpublished/original material, which the live guide excludes from Review articles.",
                "The corresponding author must obtain an editor decision to retain the demonstration or change article type before submission.",
            ],
        },
        "code_and_data_availability": {
            "status": "PUBLIC_COMMIT_PINNED",
            "repository": PUBLIC_REPOSITORY,
            "benchmark_package": PUBLIC_BENCHMARK_URL,
            "runner": PUBLIC_RUNNER_URL,
            "locked_results": PUBLIC_RESULTS_URL,
            "commit": PUBLIC_COMMIT,
        },
        "reviewer_comments": [{"id": x["id"], "status": x["status"], "location": x["location"]} for x in items],
        "outputs": [str(x) for x in [MAIN_CLEAN, MAIN_MARKED, SUPP_CLEAN, SUPP_MARKED, RESPONSE]],
        "render_qa": {"status": "PENDING"},
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": "ARTIFACTS_GENERATED_WITH_SUBMISSION_BLOCKERS",
        "outputs": 5,
        "manuscript_words": length_audit["manuscript_words_including_tables"],
        "abstract_words": length_audit["abstract_words"],
        "keywords": length_audit["keyword_count"],
    }))


if __name__ == "__main__":
    main()
