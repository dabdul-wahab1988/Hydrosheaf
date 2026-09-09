from __future__ import annotations

import json
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(r"C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M1\Revision_5_9_2026")
MAIN = ROOT / "Manuscript- Water and Ecology.docx"
SUPP = ROOT / "SupplementaryInformation.docx"
MAIN_CLEAN = ROOT / "Manuscript- Water and Ecology_Revised_Clean.docx"
MAIN_MARKED = ROOT / "Manuscript- Water and Ecology_Revised_Colour_Marked.docx"
SUPP_CLEAN = ROOT / "SupplementaryInformation_Revised_Clean.docx"
SUPP_MARKED = ROOT / "SupplementaryInformation_Revised_Colour_Marked.docx"
RESPONSE = ROOT / "Response_to_Reviewers_Master.docx"
CHANGE_LOG = ROOT / "Revision_Change_Log.json"

COLORS = {
    "r1": "C00000",
    "r2": "0000FF",
    "both": "7030A0",
    "editorial": "008000",
    "black": "000000",
    "gray": "666666",
    "header": "1F4E78",
    "alt": "EAF2F8",
    "border": "D9D9D9",
}


def rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


def set_font(run, color=None, size=None, bold=None, italic=None):
    run.font.name = "Times New Roman"
    rpr = run._element.get_or_add_rPr()
    rpr.rFonts.set(qn("w:ascii"), "Times New Roman")
    rpr.rFonts.set(qn("w:hAnsi"), "Times New Roman")
    if color:
        run.font.color.rgb = rgb(color)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


def clear_paragraph(paragraph):
    for child in list(paragraph._p):
        if child.tag != qn("w:pPr"):
            paragraph._p.remove(child)


def clear_paragraph_borders(paragraph):
    ppr = paragraph._p.get_or_add_pPr()
    pbdr = ppr.find(qn("w:pBdr"))
    if pbdr is not None:
        ppr.remove(pbdr)


def put_text(paragraph, text, color=None, bold=None, italic=None, size=None):
    clear_paragraph(paragraph)
    run = paragraph.add_run(text)
    set_font(run, color=color, bold=bold, italic=italic, size=size)
    return run


def find_prefix(doc, prefix):
    matches = [p for p in doc.paragraphs if p.text.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one paragraph beginning {prefix!r}; found {len(matches)}")
    return matches[0]


def find_contains(doc, text):
    matches = [p for p in doc.paragraphs if text in p.text]
    if len(matches) != 1:
        raise ValueError(f"Expected one paragraph containing {text!r}; found {len(matches)}")
    return matches[0]


def replace_full(doc, prefix, new_text, marked, color, change_id, changes):
    paragraph = find_prefix(doc, prefix)
    old = paragraph.text
    put_text(paragraph, new_text, color=color if marked else COLORS["black"])
    changes.append({
        "id": change_id,
        "action": "replace",
        "location_prefix": prefix,
        "old_text": old,
        "new_text": new_text,
    })
    return paragraph


def replace_substring(doc, prefix, old_text, new_text, marked, color, change_id, changes):
    paragraph = find_prefix(doc, prefix)
    original = paragraph.text
    if old_text not in original:
        raise ValueError(f"Substring not found in paragraph {prefix!r}: {old_text!r}")
    clear_paragraph(paragraph)
    pieces = original.split(old_text)
    for i, piece in enumerate(pieces):
        if piece:
            set_font(paragraph.add_run(piece), color=COLORS["black"])
        if i < len(pieces) - 1:
            set_font(paragraph.add_run(new_text), color=color if marked else COLORS["black"])
    changes.append({
        "id": change_id,
        "action": "replace_substring",
        "location_prefix": prefix,
        "old_text": old_text,
        "new_text": new_text,
    })
    return paragraph


def remove_prefix(doc, prefix, change_id, changes):
    paragraph = find_prefix(doc, prefix)
    old = paragraph.text
    paragraph._element.getparent().remove(paragraph._element)
    changes.append({
        "id": change_id,
        "action": "delete",
        "location_prefix": prefix,
        "old_text": old,
    })


def add_before(anchor, text, marked, color, changes, change_id, heading=False):
    p = anchor.insert_paragraph_before()
    p.style = "Normal"
    put_text(p, text, color=color if marked else COLORS["black"], bold=heading)
    if heading:
        p.paragraph_format.space_before = Pt(8)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
    changes.append({"id": change_id, "action": "insert", "new_text": text})
    return p


def add_after_last(doc, text, marked, color, changes, change_id, heading=False):
    p = doc.add_paragraph()
    p.style = "Normal"
    put_text(p, text, color=color if marked else COLORS["black"], bold=heading)
    if heading:
        p.paragraph_format.space_before = Pt(8)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
    changes.append({"id": change_id, "action": "append", "new_text": text})
    return p


def set_cell_shading(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)
    shd.set(qn("w:val"), "clear")


def set_cell_borders(cell, color=COLORS["border"]):
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = qn(f"w:{edge}")
        item = borders.find(tag)
        if item is None:
            item = OxmlElement(f"w:{edge}")
            borders.append(item)
        item.set(qn("w:val"), "single")
        item.set(qn("w:sz"), "4")
        item.set(qn("w:space"), "0")
        item.set(qn("w:color"), color)


def set_cell_width(cell, width_inches):
    cell.width = Inches(width_inches)
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(int(width_inches * 1440)))
    tc_w.set(qn("w:type"), "dxa")


def repeat_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tag = OxmlElement("w:tblHeader")
    tag.set(qn("w:val"), "true")
    tr_pr.append(tag)


def format_table(table, marked, color, widths):
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    for ri, row in enumerate(table.rows):
        if ri == 0:
            repeat_header(row)
        for ci, cell in enumerate(row.cells):
            set_cell_width(cell, widths[ci])
            set_cell_borders(cell)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if ri == 0:
                set_cell_shading(cell, COLORS["header"])
            elif ri % 2 == 0:
                set_cell_shading(cell, COLORS["alt"])
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(0)
                p.paragraph_format.line_spacing = 1.0
                for run in p.runs:
                    set_font(
                        run,
                        color="FFFFFF" if ri == 0 else (color if marked else COLORS["black"]),
                        size=8.0,
                        bold=(ri == 0),
                    )


def add_table7(doc, marked, changes):
    caption = add_after_last(
        doc,
        "Table 7. Cross-method comparison for karst, fractured-rock and porous aquifers.",
        marked,
        COLORS["r1"],
        changes,
        "R1.3-caption",
    )
    caption.paragraph_format.keep_with_next = True
    headers = [
        "Aquifer setting",
        "Tracer and residence-time role",
        "Graph-topology role",
        "Inverse-modelling role",
        "Integration advantage",
        "Limitation and priority scenario",
    ]
    rows = [
        [
            "Karst",
            "Use complementary tracers and residence-time classes to identify rapid, mixed or older contributions; avoid treating one apparent age as a unique travel time.",
            "Represent conduits, recharge inlets, spring outlets and plausible conduit connections; support edges with tracer tests, geophysics and hydraulic evidence.",
            "Test reactions along graph-supported conduit or spring connections; constrain phases with mineralogy, saturation indices and redox evidence.",
            "Connect rapid-pathway evidence with reaction interpretation while retaining alternative conduit-network realisations.",
            "Conduit geometry and active pathways are incompletely observed and may change with conditions. Priority: source-to-spring connectivity and vulnerability screening where conduit or tracer evidence exists.",
        ],
        [
            "Fractured rock",
            "Use tracer suites and age classes with explicit allowance for matrix storage, fracture flow and mixed screened intervals.",
            "Represent fracture intersections or monitoring nodes, but treat the graph as a plausible realisation rather than a unique map.",
            "Test whether chemistry along supported fracture connections is compatible with the proposed reaction sequence and residence-time window.",
            "Separate geometric fracture connection from evidence for active flow and chemical evolution.",
            "Fracture connectivity is poorly observed and distributed matrix flow can be missed. Priority: rank candidate pathways and identify data needs before detailed flow modelling.",
        ],
        [
            "Porous aquifer",
            "Use lumped-parameter or numerical travel-time models where data permit, and broad age classes when tracer coverage is sparse.",
            "Use model cells, wells, recharge zones and receptor nodes; derive directed edges from heads or particle tracking when a calibrated model exists.",
            "Restrict initial-final water pairs to graph-supported relationships and report all viable reaction solutions.",
            "Link simulated travel-time structure to tracer evidence and reduce arbitrary pair selection in inverse modelling.",
            "Results depend on model structure, calibration and well-screen metadata. Priority: benchmark or regional screening where flow-model files and hydrochemical or tracer data can be linked.",
        ],
    ]
    table = doc.add_table(rows=1, cols=len(headers))
    for ci, value in enumerate(headers):
        table.rows[0].cells[ci].text = value
    for values in rows:
        cells = table.add_row().cells
        for ci, value in enumerate(values):
            cells[ci].text = value
    format_table(table, marked, COLORS["r1"], [1.05, 1.75, 1.75, 1.65, 1.75, 3.15])
    changes.append({
        "id": "R1.3-table",
        "action": "insert_table",
        "rows": len(rows) + 1,
        "columns": len(headers),
    })


def replace_table_typo(doc, marked, changes):
    old = "palegrooundwater"
    new = "palaeogroundwater"
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    if old not in p.text:
                        continue
                    original = p.text
                    clear_paragraph(p)
                    pieces = original.split(old)
                    for i, piece in enumerate(pieces):
                        if piece:
                            set_font(p.add_run(piece), color=COLORS["black"])
                        if i < len(pieces) - 1:
                            set_font(p.add_run(new), color=COLORS["editorial"] if marked else COLORS["black"])
                    changes.append({
                        "id": "R1.9-table-typo",
                        "action": "replace_cell_text",
                        "old_text": original,
                        "new_text": p.text,
                    })


def add_reference(doc, author_prefix, reference, anchor_prefix, marked, changes, change_id):
    if any(p.text.startswith(author_prefix) for p in doc.paragraphs):
        return
    try:
        anchor = find_prefix(doc, anchor_prefix)
        add_before(anchor, reference, marked, COLORS["r1"], changes, change_id)
    except ValueError:
        add_after_last(doc, reference, marked, COLORS["r1"], changes, change_id)


Suckow = (
    "Suckow, A. (2014). The age of groundwater - Definitions, models and why we do not need this term. "
    "Applied Geochemistry, 50, 222-230. https://doi.org/10.1016/j.apgeochem.2014.04.016"
)
Fentaw = (
    "Fentaw, B., Birhanu, B., Azagegn, T., & Abebe, B. (2024). Groundwater recharge sources and mechanisms "
    "in the Ethiopian central Afar rift: Insights from isotopic and hydrogeochemical tracers. Journal of African "
    "Earth Sciences, 105299. https://doi.org/10.1016/j.jafrearsci.2024.105299"
)
Kinoti = (
    "Kinoti, I., Leblanc, M., Lekula, M., Tweed, S., Kenabatho, P. K., Olioso, A., & Lubczynski, M. W. (2024). "
    "Hydrogeological conceptual model of Stampriet transboundary aquifer system in Southern Africa. Groundwater "
    "for Sustainable Development, 26, 101301. https://doi.org/10.1016/j.gsd.2024.101301"
)
Yidana = (
    "Yidana, S. M., Dzikunoo, E. A., Tetteh, J. D., & Mejida, R. A. (2024). Multiple conceptual model approach "
    "for assessing groundwater resources sustainability under multiple stresses. Water Resources Management, 38, "
    "173-191. https://doi.org/10.1007/s11269-023-03662-2"
)
Wali = (
    "Wali, S. U., Alias, N., Harun, S. B., Mohammed, I. U., Garba, M. L., & Atiku, M. (2024). Application of "
    "geochemical modelling and multiple regression analysis to reassess groundwater evolution in Kaduna Basin, NW "
    "Nigeria. Discover Water, 4, 99. https://doi.org/10.1007/s43832-024-00139-0"
)
Banda = (
    "Banda, K., Crestaz, E., Seliger, R., Mengistu, H., Sauramba, J., & Saraiva, M. (2025). Zambezi River Basin "
    "aquifer systems: Opportunities and challenges in using freely available data sources and groundwater flow "
    "modelling for spatial exploratory analysis. Groundwater for Sustainable Development, 29, 101421. "
    "https://doi.org/10.1016/j.gsd.2025.101421"
)
Gebru = (
    "Gebru, H., Gebreyohannes, T., Hagos, E., & Perilli, N. (2025). Hydrogeological assessment and steady-state "
    "groundwater flow modeling for groundwater management in the Golina River Sub-Basin, Northern Ethiopia, using "
    "MODFLOW 6. Water, 17, 949. https://doi.org/10.3390/w17070949"
)


def add_legend(doc):
    anchor = doc.paragraphs[0]
    entries = [
        ("Colour key for the marked manuscript", COLORS["black"], True, 11, False),
        ("Red: Reviewer 1-related revisions", COLORS["r1"], False, 10, False),
        ("Blue: Reviewer 2-related revisions", COLORS["r2"], False, 10, False),
        ("Purple: Revisions addressing both reviewers", COLORS["both"], False, 10, False),
        ("Green: English, reference and consistency corrections", COLORS["editorial"], False, 10, False),
        ("The clean copy contains the same accepted wording without colour markup. No numerical case-study result has been added because the supplied package contains no case-study data, code or outputs.", COLORS["gray"], False, 9, True),
    ]
    for i, (text, color, bold, size, italic) in enumerate(entries):
        p = anchor.insert_paragraph_before()
        put_text(p, text, color=color, bold=bold, italic=italic, size=size)
        p.paragraph_format.space_after = Pt(2 if i == 0 else 0)


def main_revisions(doc, marked):
    changes = []
    replace_table_typo(doc, marked, changes)

    abstract_old = (
        "Because the framework is conceptual and has not yet been tested on a field dataset, it should be regarded "
        "as a basis for future software implementation, benchmark testing and field validation rather than as a fully "
        "validated operational tool."
    )
    abstract_new = (
        "At present, the framework is conceptual: this review does not report a numerical benchmark or field validation. "
        "Its contribution is therefore methodological, and the proposed workflow should be treated as hypothesis-generating "
        "until it is implemented and tested on an auditable synthetic or public dataset and subsequently evaluated against "
        "independent field evidence."
    )
    replace_substring(
        doc,
        "Groundwater assessment in data-limited aquifers remains constrained",
        abstract_old,
        abstract_new,
        marked,
        COLORS["both"],
        "R1.1/R2.5-abstract",
        changes,
    )

    replace_full(
        doc,
        "The remainder of the manuscript is organised as follows.",
        "The remainder of the manuscript is organised as follows. Section 2 defines the critical-review design, search strategy and scope. Section 2.6 then introduces the integrative framework before the method-specific sections. Section 3 evaluates environmental-tracer residence-time methods and their uncertainty. Section 4 examines graph-topological representations of groundwater connectivity, including the distinction between physical, statistical and reaction graphs. Section 5 reviews inverse hydrogeochemical modelling and its dependence on assumed flow paths. Section 6 sets out the staged integrated workflows and their proposed uncertainty procedures. Section 7 discusses remaining gaps, transferability to African semi-arid aquifers and reporting standards. Section 8 concludes by stating the conceptual contribution and the evidence required for implementation and validation.",
        marked,
        COLORS["both"],
        "R2.2-roadmap",
        changes,
    )
    replace_full(
        doc,
        "This article is structured as a critical review rather than",
        "This article is a critical review of how environmental-tracer residence-time estimation, groundwater-connectivity graphs and inverse hydrogeochemical modelling can be combined in data-limited aquifers. It is not a systematic or scoping review and does not attempt exhaustive retrieval. Studies were selected for methodological relevance, clarity of assumptions, treatment of uncertainty, validation evidence and transferability to sparse-data settings. Foundational papers were retained where they define methods still used in current practice.",
        marked,
        COLORS["r2"],
        "R2.3-scope",
        changes,
    )
    replace_full(
        doc,
        "Literature searches were conducted across Web of Science",
        "Searches were conducted in Web of Science, Scopus and Google Scholar, with Consensus used only as a supplementary semantic-search tool. Records identified through the supplementary search were retained only after verification through publisher or DOI metadata. Searches covered three clusters: residence-time tracers and lumped-parameter models; graph representations of groundwater connectivity; and inverse hydrogeochemical modelling. The search was selective rather than exhaustive. Search sources, representative strings, screening logic and the role of each source are documented in Supplementary Table S1; the central evidence matrix is provided in Supplementary Table S2.",
        marked,
        COLORS["r2"],
        "R2.3-search",
        changes,
    )
    replace_full(
        doc,
        "This review has several acknowledged limitations.",
        "The review is limited by selective retrieval, English-language and indexing bias, and the heterogeneity of graph and inverse-modelling terminology. No formal risk-of-bias score was assigned; evidence was assessed narratively using data completeness, uncertainty treatment, validation, methodological transparency and transferability. The proposed integrated workflow is therefore a methodological synthesis, not an empirical estimate of performance.",
        marked,
        COLORS["r2"],
        "R2.3-limitations",
        changes,
    )
    remove_prefix(
        doc,
        "To standardise terminology, the following definitions are used consistently throughout this article.",
        "R1.2/R2.4-duplicate-definitions",
        changes,
    )

    anchor = find_prefix(doc, "3. Environmental-Tracer Residence-Time Methods")
    framework = [
        ("2.6 Integrative framework for the review", True),
        ("The three method domains constrain different dimensions of the same inference problem. Environmental tracers provide temporal evidence about recharge and the distribution of residence times. Graph topology provides a spatial representation of candidate connections among wells, recharge areas, discharge points and model cells. Inverse hydrogeochemical modelling tests whether the chemical difference between two connected observations is compatible with specified reactions, mixing and thermodynamic constraints. None of these outputs is sufficient on its own in a data-limited aquifer: tracer ages lack a unique spatial path, graph edges can represent hypotheses rather than active flow, and inverse models can fit a selected water pair without establishing that the pair is hydraulically connected.", False),
        ("The proposed sequence is therefore mutually constraining. Tracer evidence is first represented as an age class or residence-time distribution with explicit uncertainty. These results become node attributes and are used to evaluate the direction and plausibility of candidate edges. Graph-supported edges then define the initial-final water pairs that are eligible for inverse modelling. Reaction solutions are retained only with their phase-list, saturation-index, residence-time and mass-balance assumptions. Alternative age models, graph realisations and reaction phase lists are propagated into the final diagnostic distribution. The result is a ranked set of hypotheses and uncertainty bounds, not a claim that the aquifer has been uniquely reconstructed.", False),
        ("This framework is useful where data are incomplete because a study can begin with low-cost heads, well metadata, major ions and stable isotopes, then add age tracers, particle tracking and inverse modelling as evidence becomes available. The strength of an interpretation is determined by the agreement and coverage of independent evidence domains, while missing evidence remains visible rather than being treated as confirmation.", False),
    ]
    for i, (text, heading) in reversed(list(enumerate(framework))):
        add_before(anchor, text, marked, COLORS["both"], changes, f"R2.2-framework-{i}", heading=heading)

    replace_full(
        doc,
        "Groundwater residence-time interpretation requires clear distinction",
        "Terminology is used as follows. A groundwater age is an idealised elapsed time since recharge for a water parcel or flow contribution; it is not measured directly. An apparent age is inferred from a tracer concentration under a stated input function, correction model and lumped-parameter model. A residence-time distribution (RTD) describes the distribution of ages contributing to a sample or discharge, and the mean residence time is the mean of that distribution. When the available tracer evidence cannot support a numerical RTD, this review reports a residence-time class such as modern, mixed/intermediate or old rather than a precise age. These distinctions follow Suckow (2014) and make clear that converting a tracer concentration to an age is a modelling step, not a direct measurement.",
        marked,
        COLORS["both"],
        "R1.2/R2.4-age-terms-1",
        changes,
    )
    replace_full(
        doc,
        "Mean residence time provides a flow-system-scale interpretation",
        "These distinctions affect how temporal evidence enters the integrated framework. Numerical RTDs should be used only when tracer coverage, input functions and model assumptions support them. In sparse-data settings, age classes retain the direction of the evidence without implying a false precision. The selected class or RTD, its uncertainty and its model assumptions should be attached to graph nodes and used as a constraint on candidate edges rather than interpreted as a spatial flow map.",
        marked,
        COLORS["both"],
        "R1.2/R2.4-age-terms-2",
        changes,
    )
    replace_full(
        doc,
        "Carbon-14 remains the principal tracer",
        "Radiocarbon in dissolved inorganic carbon is useful for old groundwater but does not provide a direct travel time without correction. Carbonate dissolution, soil and geogenic carbon, redox reactions and mixing can alter the measured activity. Correction-model choice can therefore shift the inferred apparent age substantially, particularly in arid and semi-arid aquifers. The review uses 14C as one temporal constraint within a multi-tracer and geochemical framework, not as a standalone age clock. Correction assumptions, δ13C and DIC data, the selected model and the resulting uncertainty range should be reported together (Han & Plummer, 2015; Cartwright et al., 2019; Seltzer et al., 2021; Suckow, 2014).",
        marked,
        COLORS["both"],
        "R1.2/R2.4-radiocarbon-1",
        changes,
    )
    replace_full(
        doc,
        "Several approaches are available for correcting 14C groundwater ages",
        "Correction models should be compared rather than treated as interchangeable. Single-sample models are useful where data are sparse, statistical 14C-δ13C approaches require sufficient aquifer-wide observations, and full geochemical correction requires detailed chemical, isotopic and mineralogical data. For very old groundwater, 36Cl, 81Kr and 4He provide complementary constraints, but each also carries limitations related to initial ratios, specialised laboratories or uncertain subsurface production. These dependencies reinforce the need to report correction-model spread as part of age uncertainty and to combine old-water tracers with graph and reaction evidence.",
        marked,
        COLORS["both"],
        "R1.2/R2.4-radiocarbon-2",
        changes,
    )
    replace_full(
        doc,
        "Lumped parameter models provide the main link",
        "Lumped-parameter models are useful here because they translate tracer observations into a residence-time constraint without claiming to resolve detailed spatial flow. Piston-flow, exponential, dispersion, binary-mixing and Bayesian multi-tracer models imply different age structures and therefore different constraints on graph direction and reaction time. The relevant selection is the simplest model supported by the data, with model choice and input-function uncertainty retained in the reported age class or RTD.",
        marked,
        COLORS["r2"],
        "R2.4-lpm",
        changes,
    )
    replace_full(
        doc,
        "Despite these uncertainties, environmental tracers remain highly valuable",
        "For data-limited semi-arid aquifers, the minimum defensible package is one young-water indicator where feasible, one old-water indicator where feasible, stable water isotopes, major ions, field parameters and basic well-construction metadata. When this package is incomplete, tracer evidence should be reported as a residence-time class with explicit uncertainty and interpreted jointly with connectivity and hydrogeochemical evidence rather than treated as a standalone conclusion.",
        marked,
        COLORS["both"],
        "R1.9/R2.4-minimum-data",
        changes,
    )

    replace_full(
        doc,
        "Graph approaches are most physically intuitive in karst",
        "Graph approaches require aquifer-specific interpretation. In karst, nodes and edges may approximate conduits when supported by mapping, tracer tests, geophysics or hydraulic evidence. In fractured rock, edges may represent plausible fracture connections but should allow for matrix storage and mixed screened intervals. In porous aquifers, nodes more commonly represent wells, model cells, recharge zones or receptors, and directed edges are most defensible when derived from heads or particle tracking. Graph metrics are secondary to the evidence used to define an edge; the cross-method comparison is summarised in Table 7.",
        marked,
        COLORS["r1"],
        "R1.3-aquifer-comparison",
        changes,
    )
    replace_full(
        doc,
        "For porous aquifers, MODFLOW/MODPATH particle tracking",
        "For porous aquifers, MODFLOW/MODPATH particle tracking can provide graph priors when a groundwater-flow model is available. Repeated source-to-receptor pathways can be translated into directed edges whose weights represent particle-crossing frequency, flux contribution or travel-time distributions. The result is a model-conditioned prior, not an independently observed flow map, and its uncertainty should include alternative calibrated model realisations.",
        marked,
        COLORS["r2"],
        "R2.4-modflow",
        changes,
    )
    replace_full(
        doc,
        "Edge confidence is therefore essential.",
        "Edge confidence should be recorded as an evidence ledger rather than treated as an unqualified probability. For edge e, define a support score C_e = sum(w_k s_ek) / sum(w_k a_ek) and an evidence-coverage score A_e = sum(w_k a_ek) / sum(w_k), where k indexes evidence domains, w_k is a pre-specified domain weight, s_ek is the support assigned by an available domain, and a_ek indicates whether that domain is available for the edge. Support should be coded before the final interpretation using a documented scale, for example 0 for contradiction, 0.5 for ambiguous evidence and 1 for support. Missing evidence reduces A_e and is not counted as agreement. The ledger should retain the underlying head-gradient, travel-time, age, lithological, chemical, redox and model-support values. No universal threshold should be imposed across aquifers; thresholds and sensitivity analyses must be justified for the study setting. This procedure makes the graph transparent and updateable while preserving the distinction between a ranked hypothesis and a validated hydraulic connection.",
        marked,
        COLORS["r1"],
        "R1.6-edge-score",
        changes,
    )
    replace_full(
        doc,
        "Sheaf-style consistency analysis offers a useful conceptual tool",
        "Sheaf-style consistency analysis is treated here as a proposed screening concept, not as a mature groundwater-network method. It asks whether observations assigned to connected nodes can be reconciled under stated restriction rules. For a candidate downgradient edge, the rules may compare hydraulic head, tracer-derived residence-time class, stable isotopes, major ions and redox indicators. A contradiction flags an edge for review; it does not by itself disprove the edge because mixing, local recharge, pumping-induced leakage, long screens or tracer bias may produce the same pattern.",
        marked,
        COLORS["both"],
        "R1.2-sheaf-1",
        changes,
    )
    replace_full(
        doc,
        "This review treats sheaf-style analysis as a transferable screening concept",
        "The proposed consistency analysis is not presented as a formal theory or probability model. Formal groundwater implementations, calibration data and independent validation remain necessary before a consistency score can be interpreted as a validated probability of hydraulic connection. Its immediate role is to flag edges where age, chemistry or redox evolution contradicts the assumed hydraulic connection.",
        marked,
        COLORS["both"],
        "R1.2-sheaf-2",
        changes,
    )
    replace_full(
        doc,
        "A major source of confusion in groundwater studies is the use",
        "Three graph objects are distinguished. A physical water-map graph represents candidate hydraulic connectivity among wells, springs, recharge zones, discharge points, conduits, fractures or model cells. Its edges require hydrogeological support such as head gradients, tracer tests, particle tracks, geological continuity, geophysics or mapped conduit and fracture networks. A hydrogeochemical statistical graph connects samples because their compositions are similar under a specified distance, correlation, clustering or dimensional-reduction rule. It describes chemical similarity or community structure, not water movement. A reaction graph represents transformations among aqueous species, minerals, gases or redox processes within a geochemical model. It describes chemical process structure and does not establish hydraulic connection between samples.",
        marked,
        COLORS["r1"],
        "R1.4-graph-types-1",
        changes,
    )
    replace_full(
        doc,
        "These graph types must not be conflated.",
        "These objects answer different questions and cannot substitute for one another. A statistical edge may group samples or generate a candidate hypothesis, but it must not be used as proof of a flow path. A reaction solution may satisfy mass balance for an initial and final water pair, but it must not be used as evidence that the pair is hydraulically connected. A physical graph constructed from distance or map proximity alone must not be reported as an active-flow map. In the proposed workflow, only physically or hydrogeologically constrained edges define candidate flow paths; statistical graphs support grouping and hypothesis generation, and reaction graphs support process diagnosis after a candidate edge has been established.",
        marked,
        COLORS["r1"],
        "R1.4-graph-types-2",
        changes,
    )
    add_before(
        find_prefix(doc, "5. Inverse Hydrogeochemistry and Reaction-Pathway Diagnosis"),
        "Typical misuse includes treating compositional similarity between two wells as direct hydraulic connectivity without hydraulic, tracer, geological or particle-tracking support; treating a chemically feasible PHREEQC or NETPATH solution as proof that the selected samples are an upgradient-downgradient pair; treating distance-based edges as active pathways without considering barriers, vertical separation, pumping or transient flow; and treating a statistical edge and a reaction solution as independent confirmation when both were derived from the same hydrochemical variables.",
        marked,
        COLORS["r1"],
        changes,
        "R1.4-misuse-cases",
    )

    replace_full(
        doc,
        "PHREEQC and NETPATH are the most widely used tools",
        "PHREEQC and NETPATH are used here as examples of inverse mass-balance tools, not as the contribution of the paper. They estimate reaction and mixing terms between a defined initial and final water pair, but they do not independently derive the hydraulic path connecting those waters. Their relevance to the proposed framework is conditional: graph evidence selects defensible candidate pairs, mineralogy and saturation indices constrain the phase list, residence-time evidence constrains reaction-time plausibility, and all viable solutions are retained.",
        marked,
        COLORS["r2"],
        "R2.4-inverse-tools-1",
        changes,
    )
    replace_full(
        doc,
        "The strength of these models is that they translate qualitative hydrochemical interpretation",
        "The strength of inverse models is that they translate qualitative hydrochemical interpretation into explicit mass-transfer estimates. Their limitation is that the modeller must define the initial water, final water, candidate phases and uncertainty tolerances before the model is run. Inverse modelling therefore tests whether a selected pair and phase list can satisfy mass balance; it does not discover the flow path. A chemically valid result remains hydrogeologically ambiguous if the pair is not independently supported.",
        marked,
        COLORS["r2"],
        "R2.4-inverse-tools-2",
        changes,
    )

    replace_full(
        doc,
        "Environmental-tracer residence time, graph topology and inverse hydrogeochemical modelling each address",
        "Section 2.6 introduced the integration logic. This section specifies how the three evidence domains can be implemented as a mutually constraining sequence: residence-time evidence supplies temporal constraints, graph topology makes candidate connections explicit, and inverse modelling tests chemical evolution only along graph-supported pairs. The sequence does not eliminate non-uniqueness; it makes its sources visible and testable.",
        marked,
        COLORS["both"],
        "R2.2-late-integration",
        changes,
    )
    replace_full(
        doc,
        "The proposed integrated framework addresses this problem",
        "The workflow is mutually constraining. Residence-time classes are assigned to graph nodes and used to test whether inferred flow directions are plausible. Graph topology restricts inverse modelling to hydraulically defensible initial-final sample pairs. Inverse models then test whether chemistry along each graph-supported edge is consistent with mineralogy, saturation indices, redox state and residence-time evidence. Alternative age models, graph realisations and phase lists are retained so that the final output is a distribution of ranked hypotheses rather than a single preferred reconstruction.",
        marked,
        COLORS["both"],
        "R2.2-integration-sequence",
        changes,
    )
    replace_full(
        doc,
        "The final workflow propagates uncertainty across residence-time estimation",
        "The uncertainty procedure should be implemented as a joint ensemble rather than as separate error bars. Each draw samples tracer observations within analytical error and detection-limit models, selects an input-function history and age-model class according to pre-specified alternatives, generates a graph realisation from the edge-evidence ledger, and selects a phase list using the documented mineralogical, saturation-index and redox filters. Inverse modelling is then run for every retained graph-supported initial-final pair. The output of each draw records edge inclusion, age class or RTD parameters, mass-balance residuals, viable reaction solutions and consistency scores. The analysis should report the number of draws, random seed, distributions or discrete alternatives, convergence or stability checks, failed runs and the rule for retaining solutions. These procedures quantify uncertainty only when the required observations and assumptions are available; they remain proposed and not empirically calibrated in this conceptual review.",
        marked,
        COLORS["r1"],
        "R1.7-uncertainty-implementation",
        changes,
    )

    discussion_anchor = find_prefix(doc, "7. Discussion")
    add_before(
        discussion_anchor,
        "6.8 Evidence required for a worked demonstration",
        marked,
        COLORS["both"],
        changes,
        "R1.1/R2.5-case-study-heading",
        heading=True,
    )
    add_before(
        discussion_anchor,
        "Because this article is a conceptual critical review, it does not report a public-data or synthetic case study. A future demonstration should identify the dataset and access route, describe node attributes and preprocessing, assign residence-time classes, construct the candidate directed graph, calculate edge support and coverage, select graph-supported initial-final pairs, run the inverse model across pre-specified phase lists, and propagate measurement, structural and model uncertainty. It should report edge-inclusion frequencies, age-class distributions, viable reaction-solution counts and cross-method consistency. Such a demonstration would test workflow logic or computational reproducibility; it would not by itself establish field transferability or independent validation.",
        marked,
        COLORS["both"],
        changes,
        "R1.1/R2.5-case-study-scope",
    )

    add_before(
        find_prefix(doc, "Future integrated studies should meet minimum reporting standards"),
        "Recent African studies reinforce the need for aquifer-specific conceptualisation and explicit model-structure uncertainty. Kinoti et al. (2024) developed a three-dimensional hydrostratigraphic and flow conceptual model for the Stampriet Transboundary Aquifer System and identified regional flow divides and structural barriers. Yidana et al. (2024) compared six transient groundwater-flow models derived from alternative conceptualisations in southern Ghana, illustrating how boundary and vertical-structure uncertainty affects sustainability assessment. Fentaw et al. (2024) combined stable isotopes, hydrochemistry and geological structures to propose a groundwater-flow model for the central Afar Rift. Wali et al. (2024) used geochemical modelling and multiple regression to examine groundwater evolution in the Kaduna Basin, providing a recent example of process interpretation in a Nigerian basement setting. At basin scale, Banda et al. (2025) used freely available datasets and a FEFLOW model to examine the opportunities and limitations of exploratory groundwater modelling in the Zambezi River Basin, while Gebru et al. (2025) applied MODFLOW 6 in the Golina River sub-basin of northern Ethiopia. These studies do not validate the present integrated framework and do not all use graph theory; their value here is to show that regional flow structure, model alternatives, geochemical evidence and data scarcity must be treated explicitly before a connectivity diagnosis is transferred across African aquifer settings.",
        marked,
        COLORS["r1"],
        changes,
        "R1.8-africa-literature",
    )
    replace_full(
        doc,
        "This review argues that data-limited aquifer interpretation remains non-unique",
        "This review develops a conceptual sequence for interpreting data-limited aquifers by combining three complementary evidence types. Residence-time analysis supplies temporal constraints, graph topology makes candidate hydraulic connections explicit, and inverse hydrogeochemical modelling tests whether chemistry along those candidate connections is compatible with specified reactions and mixing. The sequence is useful because each domain constrains a different source of non-uniqueness; it does not make any one domain a substitute for the others.",
        marked,
        COLORS["both"],
        "R2.6-conclusion-1",
        changes,
    )
    replace_full(
        doc,
        "The framework is not yet a validated field tool.",
        "The framework is not yet a validated operational tool. Its proposed edge scores, ensemble procedures and sheaf-style consistency checks require implementation, benchmark testing and independent field evaluation. Until those tests are completed, outputs should be reported as hypotheses or uncertainty-aware decision support, not as management-ready reconstructions. The next empirical step is a reproducible public-data or synthetic demonstration followed by testing in contrasting karst, fractured-rock and porous aquifers, including data-limited African settings.",
        marked,
        COLORS["both"],
        "R1.1/R2.5/R2.6-conclusion-2",
        changes,
    )

    replace_full(
        doc,
        "Figure 1. Conceptual integration",
        "Figure 1. Conceptual integration of environmental-tracer residence time, directed graph topology and inverse hydrogeochemical modelling. Residence-time evidence constrains timing and recharge source, graph topology constrains plausible hydraulic connectivity, and inverse hydrogeochemical modelling tests reaction and mixing processes along graph-supported edges. Uncertainty propagation and independent validation are required before the integrated diagnosis can be treated as management-ready.",
        marked,
        COLORS["editorial"],
        "R1.9-figure1-caption",
        changes,
    )
    replace_full(
        doc,
        "Figure 3. Point-based versus graph-constrained groundwater interpretation",
        "Figure 3. Point-based versus graph-constrained groundwater interpretation. The paired schematic shows how isolated sample interpretation differs from network edge-based process inference.",
        marked,
        COLORS["editorial"],
        "R1.9-figure3-caption",
        changes,
    )
    replace_full(
        doc,
        "Figure 4. Sources of uncertainty",
        "Figure 4. Sources and propagation of uncertainty in data-limited aquifer diagnosis. Feasible quantification routes include analytical-error and detection-limit models for tracer measurements, alternative recharge and atmospheric-input histories for age modelling, stochastic graph or particle-track ensembles for connectivity, phase-list and saturation-index sensitivity for inverse reactions, and multi-model or parameter-sensitivity analysis for model selection and scale. These methods quantify uncertainty only when the required observations and model assumptions are available. Independent validation is required before an integrated diagnosis is treated as management-ready.",
        marked,
        COLORS["r1"],
        "R1.5-figure4-caption",
        changes,
    )

    for author, ref, anchor_prefix, change_id in [
        ("Banda, K.", Banda, "Banerjee, A.", "R1.8-ref-Banda"),
        ("Fentaw, B.", Fentaw, "Ferrer, N.", "R1.8-ref-Fentaw"),
        ("Gebru, H.", Gebru, "Gerber, C.", "R1.8-ref-Gebru"),
        ("Kinoti, I.", Kinoti, "Keesari, T.", "R1.8-ref-Kinoti"),
        ("Suckow, A.", Suckow, "Taccari, M.", "R1.2-ref-Suckow"),
        ("Wali, S. U.", Wali, "Wang, J.", "R1.8-ref-Wali"),
        ("Yidana, S. M.", Yidana, "Yu, X.", "R1.8-ref-Yidana"),
    ]:
        add_reference(doc, author, ref, anchor_prefix, marked, changes, change_id)

    add_table7(doc, marked, changes)
    if marked:
        add_legend(doc)
    return changes


def revise_supplement(doc, marked):
    changes = []
    replace_full(
        doc,
        "Note: Supplementary Table S1 documents the search strategy",
        "Note: Supplementary Table S1 documents the selective search strategy used for this critical review. It does not report PRISMA-style hit counts, deduplication counts or final inclusion numbers because the review did not aim to exhaustively identify every eligible study. Instead, the search was used to identify representative, methodologically relevant and high-quality studies capable of supporting critical comparison across environmental-tracer residence-time estimation, graph-topological connectivity and inverse hydrogeochemical modelling.",
        marked,
        COLORS["r2"],
        "SI-R2.3-S1-note",
        changes,
    )
    if any("Cid-Escobar" in p.text for p in doc.paragraphs):
        replace_full(
            doc,
            "Cid-Escobar, D.",
            "Cid-Escobar, D., Folch, A., Ferrer, N., Katuva, J., & Sanchez-Vila, X. (2024). An assessment tool to improve rural groundwater access: Integrating hydrogeological modelling with socio-technical factors. Science of the Total Environment, 912, 168864. https://doi.org/10.1016/j.scitotenv.2023.168864",
            marked,
            COLORS["editorial"],
            "SI-cross-document-Cid-Escobar",
            changes,
        )
    for author, ref, anchor_prefix, change_id in [
        ("Banda, K.", Banda, "Banerjee, A.", "SI-R1.8-ref-Banda"),
        ("Fentaw, B.", Fentaw, "Ferrer, N.", "SI-R1.8-ref-Fentaw"),
        ("Gebru, H.", Gebru, "Gerber, C.", "SI-R1.8-ref-Gebru"),
        ("Kinoti, I.", Kinoti, "Keesari, T.", "SI-R1.8-ref-Kinoti"),
        ("Suckow, A.", Suckow, "Thiros, N.", "SI-R1.2-ref-Suckow"),
        ("Yidana, S. M.", Yidana, "Yu, X.", "SI-R1.8-ref-Yidana"),
    ]:
        add_reference(doc, author, ref, anchor_prefix, marked, changes, change_id)
    if marked:
        add_legend(doc)
    return changes


REVIEWER_COMMENTS = [
    {
        "id": "R1.1",
        "reviewer": "Reviewer 1",
        "comment": "It remains purely conceptual. Provide case studies and numerical tests using public datasets (e.g., USGS), demonstrating the full workflow from age classification, graph construction, edge scoring, inverse modelling and uncertainty propagation.",
        "status": "PARTIALLY ADDRESSED",
        "location": "Abstract, final sentence beginning \"At present, the framework is conceptual...\"; new Section 6.8, \"Evidence required for a worked demonstration\"; Conclusion, final paragraph beginning \"The framework is not yet a validated operational tool...\"",
        "response": "We agree that the original manuscript did not contain a worked public-data or synthetic demonstration. The revision now states this limitation plainly, adds the exact minimum contents and audit trail required for a future demonstration, and prevents Figures 3 and 6 from being interpreted as case studies. This comment is not claimed as fully resolved because the supplied package contains no case-study data, code or outputs. The revised manuscript therefore does not invent numerical results or field-validation claims.",
        "text": "At present, the framework is conceptual: this review does not report a numerical benchmark or field validation. Its contribution is therefore methodological, and the proposed workflow should be treated as hypothesis-generating until it is implemented and tested on an auditable synthetic or public dataset and subsequently evaluated against independent field evidence.",
    },
    {
        "id": "R1.2",
        "reviewer": "Reviewer 1",
        "comment": "It should unify professional terms, such as groundwater age, apparent age, mean residence time and lumped parameter models, etc. Objectively explain the application maturity of the beam consistency theory and clarify that it is only used as a screening approach and is not a mature analytical tool.",
        "status": "ADDRESSED IN TEXT; EMPIRICAL MATURITY REMAINS LIMITED",
        "location": "Section 3.1, first two paragraphs; Section 4.4, paragraphs beginning \"Sheaf-style consistency analysis...\" and \"The phrase beam consistency theory...\"; added Suckow (2014) reference.",
        "response": "The terminology is now consolidated in Section 3.1. The revision distinguishes observations from model-derived apparent ages and RTDs, and reports residence-time classes when numerical RTDs are not defensible. The reviewer’s phrase \"beam consistency theory\" does not occur in the manuscript; the nearest manuscript term is \"sheaf-style consistency analysis.\" The revision makes that interpretation explicit and limits the method to screening until formal implementation, calibration and independent validation exist.",
        "text": "Terminology is used as follows. A groundwater age is an idealised elapsed time since recharge for a water parcel or flow contribution; it is not measured directly. An apparent age is inferred from a tracer concentration under a stated input function, correction model and lumped-parameter model. A residence-time distribution (RTD) describes the distribution of ages contributing to a sample or discharge, and the mean residence time is the mean of that distribution. When the available tracer evidence cannot support a numerical RTD, this review reports a residence-time class such as modern, mixed/intermediate or old rather than a precise age. These distinctions follow Suckow (2014) and make clear that converting a tracer concentration to an age is a modelling step, not a direct measurement.",
    },
    {
        "id": "R1.3",
        "reviewer": "Reviewer 1",
        "comment": "The manuscript should systematically compare graph topology, tracer dating, and inverse modelling across karst, fractured, and porous aquifers. For each type, discuss adaptability, advantages, limitations, and priority application scenarios. A summary table is recommended.",
        "status": "ADDRESSED IN TEXT AND TABLE 7",
        "location": "Section 4.1, paragraph beginning \"Graph approaches require aquifer-specific interpretation...\"; new Table 7, \"Cross-method comparison for karst, fractured-rock and porous aquifers.\"",
        "response": "We added a setting-specific comparison that separates the role of tracers, graph topology and inverse modelling, and states the integration advantage, limitation and priority scenario for karst, fractured-rock and porous aquifers. The text also clarifies that graph edges have different physical interpretations across these settings.",
        "text": "Graph approaches require aquifer-specific interpretation. In karst, nodes and edges may approximate conduits when supported by mapping, tracer tests, geophysics or hydraulic evidence. In fractured rock, edges may represent plausible fracture connections but should allow for matrix storage and mixed screened intervals. In porous aquifers, nodes more commonly represent wells, model cells, recharge zones or receptors, and directed edges are most defensible when derived from heads or particle tracking. Graph metrics are secondary to the evidence used to define an edge; the cross-method comparison is summarised in Table 7.",
    },
    {
        "id": "R1.4",
        "reviewer": "Reviewer 1",
        "comment": "Strictly distinguish the essential differences among the three types of graph models, including physical water maps, hydrogeochemical statistical graphs, and reaction graphs, and supplement typical cases of misuse in the field of groundwater.",
        "status": "ADDRESSED IN TEXT",
        "location": "Section 4.5, revised graph-definition paragraphs; added paragraph listing four misuse cases.",
        "response": "The revised Section 4.5 defines the object represented by nodes and edges, the evidence required, and the question answered by each graph type. It then lists typical misuse patterns without attributing them to any particular cited study.",
        "text": "A statistical edge may group samples or generate a candidate hypothesis, but it must not be used as proof of a flow path. A reaction solution may satisfy mass balance for an initial and final water pair, but it must not be used as evidence that the pair is hydraulically connected. A physical graph constructed from distance or map proximity alone must not be reported as an active-flow map.",
    },
    {
        "id": "R1.5",
        "reviewer": "Reviewer 1",
        "comment": "Fig. 4 should indicate feasible quantification methods for each stage.",
        "status": "PARTIALLY ADDRESSED",
        "location": "Figure 4 caption; Section 6.7 uncertainty procedure.",
        "response": "The Figure 4 caption now identifies feasible quantification routes for tracer measurement, age modelling, connectivity, inverse reactions and model selection, and states the data conditions required for their use. The accompanying uncertainty procedure is also specified in Section 6.7. The supplied Figure 4 artwork is a separate PNG and has not yet been redrawn to place those routes inside each graphical stage; therefore this comment remains partial until the updated figure asset is produced and checked.",
        "text": "Feasible quantification routes include analytical-error and detection-limit models for tracer measurements, alternative recharge and atmospheric-input histories for age modelling, stochastic graph or particle-track ensembles for connectivity, phase-list and saturation-index sensitivity for inverse reactions, and multi-model or parameter-sensitivity analysis for model selection and scale.",
    },
    {
        "id": "R1.6",
        "reviewer": "Reviewer 1",
        "comment": "Section 6.3, edge confidence scoring is too vague.",
        "status": "ADDRESSED IN TEXT AS A PROPOSED, UNCALIBRATED PROCEDURE",
        "location": "Section 4.3, revised paragraph beginning \"Edge confidence should be recorded as an evidence ledger...\"; Section 6.3 uses the resulting staged rule.",
        "response": "The revision replaces the qualitative description with a support score, a separate evidence-coverage score, definitions for available evidence, a documented support scale, treatment of missing evidence, and a requirement for setting-specific thresholds and sensitivity analysis. It is explicitly presented as a proposed procedure rather than a validated probability.",
        "text": "For edge e, define a support score C_e = sum(w_k s_ek) / sum(w_k a_ek) and an evidence-coverage score A_e = sum(w_k a_ek) / sum(w_k), where k indexes evidence domains, w_k is a pre-specified domain weight, s_ek is the support assigned by an available domain, and a_ek indicates whether that domain is available for the edge.",
    },
    {
        "id": "R1.7",
        "reviewer": "Reviewer 1",
        "comment": "Section 6.7 mentions Monte Carlo and stochastic graph ensembles but lacks implementation details.",
        "status": "ADDRESSED AS A REPRODUCIBLE PROTOCOL; NOT EXECUTED",
        "location": "Section 6.7, revised paragraph beginning \"The uncertainty procedure should be implemented as a joint ensemble...\"",
        "response": "The revision specifies the uncertainty draw, the alternative input-function and age-model choices, graph realisations, phase-list filters, inverse-model execution, output ledger, random seed, convergence or stability checks, failed runs and solution-retention rules. Because no implementation or analysis package was supplied, these details are a proposed reproducibility protocol and not reported empirical results.",
        "text": "Each draw samples tracer observations within analytical error and detection-limit models, selects an input-function history and age-model class according to pre-specified alternatives, generates a graph realisation from the edge-evidence ledger, and selects a phase list using the documented mineralogical, saturation-index and redox filters.",
    },
    {
        "id": "R1.8",
        "reviewer": "Reviewer 1",
        "comment": "There is a shortage of literature on the modeling of the African groundwater network from 2024 to 2025 in the references. It is recommended to supplement it.",
        "status": "ADDRESSED IN SECTION 7.5 AND REFERENCES",
        "location": "Section 7.5, new paragraph beginning \"Recent African studies reinforce...\"; main and supplementary reference lists.",
        "response": "We added a focused 2024-2025 synthesis covering conceptual flow modelling, multiple conceptual models, isotope and hydrochemical interpretation, freely available data with FEFLOW, and MODFLOW 6. The paragraph states that these studies do not validate the present integrated framework and do not all use graph theory. The same verified bibliographic records were added to the main and supplementary reference lists where absent.",
        "text": "Recent African studies reinforce the need for aquifer-specific conceptualisation and explicit model-structure uncertainty. Kinoti et al. (2024) developed a three-dimensional hydrostratigraphic and flow conceptual model for the Stampriet Transboundary Aquifer System and identified regional flow divides and structural barriers. Yidana et al. (2024) compared six transient groundwater-flow models derived from alternative conceptualisations in southern Ghana, illustrating how boundary and vertical-structure uncertainty affects sustainability assessment. Fentaw et al. (2024) combined stable isotopes, hydrochemistry and geological structures to propose a groundwater-flow model for the central Afar Rift. Wali et al. (2024) used geochemical modelling and multiple regression to examine groundwater evolution in the Kaduna Basin, providing a recent example of process interpretation in a Nigerian basement setting. At basin scale, Banda et al. (2025) used freely available datasets and a FEFLOW model to examine the opportunities and limitations of exploratory groundwater modelling in the Zambezi River Basin, while Gebru et al. (2025) applied MODFLOW 6 in the Golina River sub-basin of northern Ethiopia. These studies do not validate the present integrated framework and do not all use graph theory; their value here is to show that regional flow structure, model alternatives, geochemical evidence and data scarcity must be treated explicitly before a connectivity diagnosis is transferred across African aquifer settings.",
    },
    {
        "id": "R1.9",
        "reviewer": "Reviewer 1",
        "comment": "Polish the English.",
        "status": "ADDRESSED IN THE REVISED PASSAGES; FULL COPY-EDIT SHOULD STILL BE COMPLETED BEFORE SUBMISSION",
        "location": "Targeted revisions in Sections 2.1-2.5, 3.1, 3.3, 3.5, 4.1, 4.4-4.5, 5.2, 6.7, 7.5 and 8; Figure 1, Figure 3 and Figure 4 captions; Table 1 carbon-14 row.",
        "response": "The revision corrects the specifically identified grammatical errors, removes duplicated minimum-data wording, corrects \"palegrooundwater\" to \"palaeogroundwater\", standardises selected hyphenation and improves sentence economy in the revised sections. Because the original manuscript was not subjected to a complete line-by-line copy-edit in this pass, a final English edit remains a submission requirement.",
        "text": "Figure 1. Conceptual integration of environmental-tracer residence time, directed graph topology and inverse hydrogeochemical modelling. Residence-time evidence constrains timing and recharge source, graph topology constrains plausible hydraulic connectivity, and inverse hydrogeochemical modelling tests reaction and mixing processes along graph-supported edges. Uncertainty propagation and independent validation are required before the integrated diagnosis can be treated as management-ready.",
    },
    {
        "id": "R2.1",
        "reviewer": "Reviewer 2",
        "comment": "The manuscript is described as a critical review, but it is not clear what distinguishes it as such or what new insight the reader gains from it.",
        "status": "ADDRESSED IN TEXT",
        "location": "Section 2.1, revised scope paragraph; new Section 2.6; Abstract and Conclusion.",
        "response": "The revised scope paragraph defines the article as a critical review by its evaluative purpose: it compares assumptions, uncertainty, validation evidence, transferability and integration potential rather than simply cataloguing methods. Section 2.6 now states the central added insight: temporal, spatial and reaction evidence should be mutually constraining, and missing evidence should remain visible rather than being treated as confirmation.",
        "text": "The strength of an interpretation is determined by the agreement and coverage of independent evidence domains, while missing evidence remains visible rather than being treated as confirmation.",
    },
    {
        "id": "R2.2",
        "reviewer": "Reviewer 2",
        "comment": "The components are introduced sequentially. The manuscript should begin by introducing the conceptual framework that explains why the approaches should be combined, how they complement each other, and what additional insights are gained from integration.",
        "status": "ADDRESSED IN STRUCTURE AND TEXT",
        "location": "New Section 2.6, \"Integrative framework for the review\"; revised roadmap in the Introduction; condensed Section 6.1.",
        "response": "The integrative framework has been moved ahead of the method-specific sections. The roadmap now identifies Section 2.6 as the conceptual entry point, and Sections 3-5 are framed as answers to three linked questions: what temporal constraint is available, what spatial relation can be defended, and what reaction inference remains possible after those constraints are applied.",
        "text": "The proposed sequence is therefore mutually constraining. Tracer evidence is first represented as an age class or residence-time distribution with explicit uncertainty. These results become node attributes and are used to evaluate the direction and plausibility of candidate edges. Graph-supported edges then define the initial-final water pairs that are eligible for inverse modelling.",
    },
    {
        "id": "R2.3",
        "reviewer": "Reviewer 2",
        "comment": "The manuscript is considerably longer than necessary and devotes too much space to meta-level discussion, especially in Section 2. The search-strategy wording should be more direct and transparent.",
        "status": "ADDRESSED IN SECTION 2",
        "location": "Sections 2.1, 2.2 and 2.5; the former long search-strategy paragraph beginning \"Literature searches were conducted...\" has been replaced.",
        "response": "Section 2 has been reduced to the review type, selection logic, search sources, scope and limitations. Detailed search strings and the evidence matrix remain in Supplementary Tables S1 and S2. The revised main-text search statement is direct and no longer presents a selective critical-review search as a PRISMA corpus.",
        "text": "Searches covered three clusters: residence-time tracers and lumped-parameter models; graph representations of groundwater connectivity; and inverse hydrogeochemical modelling. The search was selective rather than exhaustive. Search sources, representative strings, screening logic and the role of each source are documented in Supplementary Table S1; the central evidence matrix is provided in Supplementary Table S2.",
    },
    {
        "id": "R2.4",
        "reviewer": "Reviewer 2",
        "comment": "Sections 3-5 should be substantially shortened and should focus on integration rather than general descriptions of established topics and packages such as MODFLOW, NETPATH and PHREEQC. The requested groundwater-age reference should also be added.",
        "status": "ADDRESSED IN TARGETED SECTIONS; FINAL LENGTH CHECK REQUIRED",
        "location": "Sections 3.1, 3.3-3.5, 4.1, 4.3-4.5 and 5.2; Suckow (2014) added to both reference lists.",
        "response": "The revised passages retain only details that affect integration, uncertainty or misuse. General package descriptions were shortened, while the dependence of inverse modelling on selected water pairs and phase lists was retained. The age terminology and radiocarbon passages now focus on the distinction between observation and model output and on how uncertainty enters the integrated workflow.",
        "text": "PHREEQC and NETPATH are used here as examples of inverse mass-balance tools, not as the contribution of the paper. They estimate reaction and mixing terms between a defined initial and final water pair, but they do not independently derive the hydraulic path connecting those waters.",
    },
    {
        "id": "R2.5",
        "reviewer": "Reviewer 2",
        "comment": "The graph-topology section would benefit from a practical example, and the manuscript lacks illustrative case studies demonstrating the proposed workflow.",
        "status": "PARTIALLY ADDRESSED",
        "location": "Section 4.5 graph distinctions; new Section 6.8, \"Evidence required for a worked demonstration\"; Conclusion.",
        "response": "The revision explains exactly how a worked public-data or synthetic demonstration would be constructed and makes clear that the existing figures are conceptual schematics. The practical example itself is not present because no underlying data, code or outputs were supplied. This comment therefore remains partially addressed and must be completed with an auditable demonstration before submission if the journal requires it.",
        "text": "Because this article is a conceptual critical review, it does not report a public-data or synthetic case study. A future demonstration should identify the dataset and access route, describe node attributes and preprocessing, assign residence-time classes, construct the candidate directed graph, calculate edge support and coverage, select graph-supported initial-final pairs, run the inverse model across pre-specified phase lists, and propagate measurement, structural and model uncertainty.",
    },
    {
        "id": "R2.6",
        "reviewer": "Reviewer 2",
        "comment": "The manuscript should be streamlined, several sections shortened, and, if possible, an illustrative synthetic case study should demonstrate how the workflow can be applied.",
        "status": "PARTIALLY ADDRESSED",
        "location": "New Section 2.6; targeted shortening in Sections 2-5; revised Section 6.7; revised Conclusion; Section 6.8 records the remaining evidence-dependent requirement.",
        "response": "The manuscript has been reorganised around the integration problem, repetitive or general material has been condensed in the identified sections, and the conclusion now states the contribution and evidence ceiling directly. The illustrative case study remains outstanding for the same reason stated under R1.1 and R2.5: the supplied package contains no data or analysis outputs from which a responsible numerical example could be written.",
        "text": "This review develops a conceptual sequence for interpreting data-limited aquifers by combining three complementary evidence types. Residence-time analysis supplies temporal constraints, graph topology makes candidate hydraulic connections explicit, and inverse hydrogeochemical modelling tests whether chemistry along those candidate connections is compatible with specified reactions and mixing.",
    },
]


def add_response_letter(doc):
    for style_name in ("Normal", "Body Text"):
        try:
            style = doc.styles[style_name]
            style.font.name = "Times New Roman"
            style.font.size = Pt(11)
        except KeyError:
            pass

    # The bundled Word Title style carries a decorative bottom rule. Remove
    # it so the response letter uses the same restrained manuscript style.
    title_style = doc.styles["Title"]
    title_style.font.color.rgb = rgb(COLORS["black"])
    title_ppr = title_style._element.pPr
    if title_ppr is not None:
        title_pbdr = title_ppr.find(qn("w:pBdr"))
        if title_pbdr is not None:
            title_ppr.remove(title_pbdr)

    title = doc.add_paragraph(style="Title")
    clear_paragraph_borders(title)
    put_text(title, "Response to Reviewers", color=COLORS["black"], size=16)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for text in [
        "Manuscript: Environmental-Tracer Residence Time, Graph Topology and Inverse Hydrogeochemistry in Data-Limited Aquifers: A Critical Review",
        "Manuscript number: WATECO-D-26-00083",
        "Response to the editor and Reviewers 1 and 2",
    ]:
        p = doc.add_paragraph()
        put_text(p, text, size=11)

    for text in [
        "Dear Editor and Reviewers,",
        "We thank the editor and both reviewers for their careful assessment. The comments identified a clear distinction between the manuscript's conceptual contribution and the evidence required for an operational groundwater workflow. We have reorganised the manuscript around the integration problem, shortened targeted method descriptions, clarified terminology and the screening status of sheaf-style consistency analysis, specified the proposed edge-confidence and ensemble procedures, added a cross-aquifer comparison, and expanded the recent African literature.",
        "We have not inserted numerical case-study results, software outputs or field-validation claims because the supplied revision package contains no case-study dataset, implementation or analysis outputs. The response below identifies this limitation explicitly. Comments that can be resolved by text, structure or reference changes are described as addressed; requests that require new computation or a revised graphical asset are identified as partial and are not presented as complete.",
    ]:
        p = doc.add_paragraph()
        put_text(p, text, size=11)

    p = doc.add_paragraph()
    put_text(p, "Colour key for the marked manuscript", bold=True, size=12)
    for label, color, description in [
        ("Red", COLORS["r1"], "Reviewer 1-related revisions"),
        ("Blue", COLORS["r2"], "Reviewer 2-related revisions"),
        ("Purple", COLORS["both"], "Revisions addressing both reviewers"),
        ("Green", COLORS["editorial"], "English, reference and consistency corrections"),
    ]:
        p = doc.add_paragraph()
        put_text(p, f"{label}: {description}", color=color, size=10)
    p = doc.add_paragraph()
    put_text(p, "The clean revised copies contain the same accepted wording without colour markup. Colour markup identifies changed or newly added text; it is not a substitute for the original-versus-revised comparison.", color=COLORS["gray"], italic=True, size=9)

    for reviewer in ("Reviewer 1", "Reviewer 2"):
        h = doc.add_paragraph()
        put_text(h, reviewer, bold=True, size=14)
        for item in [x for x in REVIEWER_COMMENTS if x["reviewer"] == reviewer]:
            h = doc.add_paragraph()
            put_text(h, f"{item['id']}  {item['status']}", bold=True, size=12)
            p = doc.add_paragraph()
            r = p.add_run("Reviewer comment: ")
            set_font(r, bold=True, size=11)
            r = p.add_run(item["comment"])
            set_font(r, italic=True, size=11)
            p = doc.add_paragraph()
            r = p.add_run("Exact location in the revised manuscript: ")
            set_font(r, bold=True, size=11)
            r = p.add_run(item["location"])
            set_font(r, size=11)
            p = doc.add_paragraph()
            r = p.add_run("Response: ")
            set_font(r, bold=True, size=11)
            r = p.add_run(item["response"])
            set_font(r, size=11)
            p = doc.add_paragraph()
            r = p.add_run("Revised or added text: ")
            set_font(r, bold=True, size=11)
            r = p.add_run(item["text"])
            set_font(r, size=10)

    h = doc.add_paragraph()
    put_text(h, "Cross-document checks", bold=True, size=14)
    for text in [
        "The Cid-Escobar reference was reconciled between the main manuscript and Supplementary Information.",
        "Suckow (2014) and the recent African references used in Section 7.5 were added consistently where absent.",
        "Supplementary Tables S1 and S2 remain the locations for search-detail and evidence-matrix information removed from the main text.",
        "The marked and clean copies were generated from the same source revision; only the colour markup and marked-copy legend differ.",
        "Figure 4 still requires a graphical redraw if the reviewer expects each quantification route to appear inside the figure rather than in its caption. This is recorded as partial, not silently claimed as complete.",
    ]:
        p = doc.add_paragraph(style="List Bullet")
        put_text(p, text, size=11)

    h = doc.add_paragraph()
    put_text(h, "Closing statement", bold=True, size=14)
    for text in [
        "We appreciate the reviewers' identification of the manuscript's main weakness: the integration rationale was present but appeared too late, while the practical evidence boundary was not operationally explicit. The revision brings the framework forward, makes the evidence dependencies visible, and provides a reproducible specification for the analyses that remain to be implemented.",
        "We respectfully ask that the revised manuscript be assessed with the distinction between a critical methodological review and an empirical validation study in mind. We nevertheless recognise that the worked public-data or synthetic demonstration requested by both reviewers would strengthen the paper materially and remains the principal outstanding revision if the journal requires an applied demonstration.",
        "Sincerely,\nThe authors",
    ]:
        p = doc.add_paragraph()
        put_text(p, text, size=11)


def main():
    main_clean_doc = Document(MAIN)
    main_clean_changes = main_revisions(main_clean_doc, marked=False)
    main_clean_doc.save(MAIN_CLEAN)

    main_marked_doc = Document(MAIN)
    main_marked_changes = main_revisions(main_marked_doc, marked=True)
    main_marked_doc.save(MAIN_MARKED)

    supp_clean_doc = Document(SUPP)
    supp_clean_changes = revise_supplement(supp_clean_doc, marked=False)
    supp_clean_doc.save(SUPP_CLEAN)

    supp_marked_doc = Document(SUPP)
    supp_marked_changes = revise_supplement(supp_marked_doc, marked=True)
    supp_marked_doc.save(SUPP_MARKED)

    response_doc = Document()
    add_response_letter(response_doc)
    response_doc.save(RESPONSE)

    log = {
        "source_main": str(MAIN),
        "source_supplement": str(SUPP),
        "outputs": [str(MAIN_CLEAN), str(MAIN_MARKED), str(SUPP_CLEAN), str(SUPP_MARKED), str(RESPONSE)],
        "colour_key": {
            "red": "Reviewer 1-related revisions",
            "blue": "Reviewer 2-related revisions",
            "purple": "Revisions addressing both reviewers",
            "green": "English, reference and consistency corrections",
        },
        "comment_status": [
            {"id": x["id"], "status": x["status"], "location": x["location"]}
            for x in REVIEWER_COMMENTS
        ],
        "main_change_count": len(main_clean_changes),
        "supplement_change_count": len(supp_clean_changes),
        "note": "No case-study numerical results or software outputs were added because none were present in the supplied revision package.",
    }
    CHANGE_LOG.write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "outputs": [p.name for p in [MAIN_CLEAN, MAIN_MARKED, SUPP_CLEAN, SUPP_MARKED, RESPONSE]],
        "main_changes": len(main_clean_changes),
        "supplement_changes": len(supp_clean_changes),
        "comment_count": len(REVIEWER_COMMENTS),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
