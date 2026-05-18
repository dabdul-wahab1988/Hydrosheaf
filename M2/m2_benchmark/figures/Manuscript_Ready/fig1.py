from pathlib import Path
from textwrap import wrap

from PIL import Image, ImageDraw, ImageFont


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_PATH = SCRIPT_DIR / "Manuscript_Fig1_Architecture.png"

CANVAS_W = 9000
CANVAS_H = 6200

TEXT = "#0F172A"
MUTED = "#111827"
NODE_FILL = "#FFFFFF"
NODE_STROKE = "#334155"
EDGE = "#64748B"
EDGE_DARK = "#475569"

GROUPS = {
    "entry": ("Entry Layer", (140, 120, 2600, 980), "#EFF6FF", "#1E6091"),
    "input": ("Input Data Layer", (2860, 120, 6000, 980), "#F0FDF4", "#2F7D4E"),
    "conditioning": ("Data Conditioning", (140, 1200, 2300, 1600), "#FFF7E8", "#A15C00"),
    "topology": ("Topology and Prior Engines", (2560, 1200, 3400, 1600), "#F4ECFF", "#6D4E9E"),
    "geochem": ("Geochemical Context and Constraints", (6080, 1200, 2780, 1600), "#FFF1F1", "#B13A3A"),
    "solver": ("Core Inverse Solver", (140, 2920, 8520, 1700), "#ECFEFF", "#267587"),
    "optional": ("Optional Inference and Validation", (140, 4780, 4250, 1300), "#FAF0FF", "#8A4D99"),
    "output": ("Output Layer", (4600, 4780, 4060, 1300), "#F1F5F9", "#52708A"),
}


def load_font(file_name: str, size: int) -> ImageFont.FreeTypeFont:
    try:
        return ImageFont.truetype(file_name, size)
    except OSError:
        return ImageFont.load_default()


F_GROUP = load_font("DejaVuSans-Bold.ttf", 72)
F_TITLE = load_font("DejaVuSans-Bold.ttf", 58)
F_SUB = load_font("DejaVuSans.ttf", 50)
F_SMALL = load_font("DejaVuSans.ttf", 44)

img = Image.new("RGB", (CANVAS_W, CANVAS_H), "#FFFFFF")
draw = ImageDraw.Draw(img)
nodes: dict[str, tuple[int, int, int, int, str, str]] = {}


def text_width(text: str, font: ImageFont.ImageFont) -> float:
    return draw.textlength(text, font=font)


def fit_lines(text: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    lines: list[str] = []
    for raw_line in str(text).splitlines():
        if not raw_line:
            lines.append("")
            continue
        words = raw_line.split()
        current = words[0]
        for word in words[1:]:
            candidate = f"{current} {word}"
            if text_width(candidate, font) <= max_width:
                current = candidate
            else:
                lines.append(current)
                current = word
        lines.append(current)
    return lines


def font_height(font: ImageFont.ImageFont) -> int:
    bbox = draw.textbbox((0, 0), "Ag", font=font)
    return bbox[3] - bbox[1]


def add_node(key: str, x: int, y: int, w: int, h: int, title: str, subtitle: str) -> None:
    nodes[key] = (x, y, w, h, title, subtitle)


def anchor(key: str, side: str) -> tuple[float, float]:
    x, y, w, h, _, _ = nodes[key]
    if side == "left":
        return x, y + h / 2
    if side == "right":
        return x + w, y + h / 2
    if side == "top":
        return x + w / 2, y
    if side == "bottom":
        return x + w / 2, y + h
    if side == "center":
        return x + w / 2, y + h / 2
    raise ValueError(f"Unsupported anchor side: {side}")


def draw_group(title: str, box: tuple[int, int, int, int], fill: str, outline: str) -> None:
    x, y, w, h = box
    draw.rounded_rectangle(
        [x, y, x + w, y + h],
        radius=18,
        fill=fill,
        outline=outline,
        width=4,
    )
    tw = text_width(title, F_GROUP)
    draw.text((x + (w - tw) / 2, y + 32), title, font=F_GROUP, fill=TEXT)


def draw_node(key: str) -> None:
    x, y, w, h, title, subtitle = nodes[key]
    draw.rounded_rectangle(
        [x, y, x + w, y + h],
        radius=14,
        fill=NODE_FILL,
        outline=NODE_STROKE,
        width=3,
    )

    max_width = w - 80
    title_lines = fit_lines(title, F_TITLE, max_width)
    sub_font = F_SMALL if w < 1150 or h < 260 else F_SUB
    subtitle_lines = fit_lines(subtitle, sub_font, max_width)

    title_step = font_height(F_TITLE) + 10
    sub_step = font_height(sub_font) + 8
    gap = 10 if subtitle_lines else 0
    total_h = len(title_lines) * title_step + gap + len(subtitle_lines) * sub_step
    y_cursor = y + (h - total_h) / 2

    for line in title_lines:
        lw = text_width(line, F_TITLE)
        draw.text((x + (w - lw) / 2, y_cursor), line, font=F_TITLE, fill=TEXT)
        y_cursor += title_step

    y_cursor += gap
    for line in subtitle_lines:
        lw = text_width(line, sub_font)
        draw.text((x + (w - lw) / 2, y_cursor), line, font=sub_font, fill=MUTED)
        y_cursor += sub_step


def draw_dashed_line(points: list[tuple[float, float]], fill: str, width: int, dash_len: int = 18, gap_len: int = 14) -> None:
    for start, end in zip(points, points[1:]):
        x1, y1 = start
        x2, y2 = end
        dx = x2 - x1
        dy = y2 - y1
        length = (dx * dx + dy * dy) ** 0.5
        if length == 0:
            continue
        step_x = dx / length
        step_y = dy / length
        pos = 0.0
        while pos < length:
            seg_end = min(pos + dash_len, length)
            draw.line(
                [
                    (x1 + step_x * pos, y1 + step_y * pos),
                    (x1 + step_x * seg_end, y1 + step_y * seg_end),
                ],
                fill=fill,
                width=width,
            )
            pos += dash_len + gap_len


def draw_arrow_head(start: tuple[float, float], end: tuple[float, float], fill: str, scale: int = 18) -> None:
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    if length == 0:
        return
    ux = dx / length
    uy = dy / length
    px = -uy
    py = ux
    base_x = x2 - ux * scale
    base_y = y2 - uy * scale
    half_w = scale * 0.55
    draw.polygon(
        [
            (x2, y2),
            (base_x + px * half_w, base_y + py * half_w),
            (base_x - px * half_w, base_y - py * half_w),
        ],
        fill=fill,
    )


def draw_arrow(
    points: list[tuple[float, float]],
    fill: str = EDGE,
    width: int = 4,
    dashed: bool = False,
    both: bool = False,
) -> None:
    if dashed:
        draw_dashed_line(points, fill, width)
    else:
        draw.line(points, fill=fill, width=width, joint="curve")
    draw_arrow_head(points[-2], points[-1], fill)
    if both:
        draw_arrow_head(points[1], points[0], fill)


def elbow(
    start_key: str,
    end_key: str,
    start_side: str = "right",
    end_side: str = "left",
    fill: str = EDGE,
    width: int = 4,
    both: bool = False,
    dashed: bool = False,
    mid: float | None = None,
) -> None:
    start = anchor(start_key, start_side)
    end = anchor(end_key, end_side)
    points = [start]
    if start_side in {"left", "right"}:
        mid_x = mid if mid is not None else start[0] + (end[0] - start[0]) * 0.5
        points.extend([(mid_x, start[1]), (mid_x, end[1])])
    else:
        mid_y = mid if mid is not None else start[1] + (end[1] - start[1]) * 0.5
        points.extend([(start[0], mid_y), (end[0], mid_y)])
    points.append(end)
    draw_arrow(points, fill=fill, width=width, both=both, dashed=dashed)


def route(points: list[tuple[float, float]], fill: str = EDGE, width: int = 4, dashed: bool = False, both: bool = False) -> None:
    draw_arrow(points, fill=fill, width=width, dashed=dashed, both=both)


# Node placement mirrors the reference Figure 1 layout.
add_node("CLI", 250, 320, 700, 240, "Command line", "hydrosheaf.cli:main")
add_node("API", 1040, 320, 760, 240, "Python API", "fit_network_pipeline()\nfit_network_with_priors()")
add_node("AUTO", 1890, 320, 700, 240, "Auto workflow", "workflows.auto.analyze_dataset()")
add_node("CFG", 320, 730, 2160, 280, "Config dataclass", "ion order, weights, feature flags, PHREEQC, isotopes, UQ, 3D, temporal")

x0, y0 = 3020, 320
bw, bh = 1800, 260
gapx, gapy = 130, 130
add_node("SAMPLES", x0, y0, bw, bh, "Sample chemistry", "site_id, major ions, pH, EC/TDS, isotopes, coordinates, heads")
add_node("EDGES", x0 + bw + gapx, y0, bw, bh, "User edges or candidate edges", "u, v, edge_id, attrs")
add_node("ENDM", x0 + 2 * (bw + gapx), y0, bw, bh, "Endmembers", "manual, JSON, sample-derived, latent virtual endmembers")
add_node("TS", x0, y0 + bh + gapy, bw, bh, "Time-series data", "TemporalNode, TimeSeriesSample")
add_node("VADOSE", x0 + bw + gapx, y0 + bh + gapy, bw, bh, "Vadose inputs", "profile, forcing, water table, vadose links")
add_node("PRIORS", x0 + 2 * (bw + gapx), y0 + bh + gapy, bw, bh, "External physics priors", "MODPATH / CSV / JSON\np_uv, tau mean/std/quantiles")

cx, cy = 320, 1370
cw, ch = 1940, 280
conditioning_nodes = [
    ("UNITS", "Unit handling", "mg/L, meq/L to mmol/L"),
    ("SCHEMA", "Schema parsing and validation", "parse_numeric(), vector_from_sample()"),
    ("QC", "QC checks", "missing values, detection limits, charge balance, non-negative chemistry"),
    ("AUTODISABLE", "Auto-disable missing modules", "PHREEQC, isotope, nitrate-source flags"),
]
for idx, (key, title, subtitle) in enumerate(conditioning_nodes):
    add_node(key, cx, cy + idx * (ch + 90), cw, ch, title, subtitle)

px, py = 2730, 1390
tw, th = 1000, 280
add_node("HEADS", px, py, tw, th, "Head inference", "heuristic, Bayesian linear, Bayesian MCMC")
add_node("VADOSEPRI", px + 1130, py, tw, th, "Vadose Richards-1D priors", "recharge, TTD, tau summaries")
add_node("TEMPPRI", px + 2260, py, tw, th, "Temporal residence-time priors", "cross-correlation, Bayesian lag, TTD convolution, tracer decay")
add_node("GRAPH2D", px, py + 430, tw, th, "2D graph inference", "infer_edges_probabilistic()\ndistance, heads, p_uv")
add_node("GRAPH3D", px + 1130, py + 430, tw, th, "3D graph inference", "Node3D, Edge3D, layers, screen overlap, anisotropy")
add_node("PHYSAPPLY", px + 2260, py + 430, tw, th, "apply_physics_priors()", "override, merge, or only")
add_node("SHEAFTOPO", px + 565, py + 970, 1600, 300, "Sheaf topology refinement", "isotope cost, Cl cost, age cost, global section residual")

ex, ey = 6250, 1390
ew, eh = 1180, 280
add_node("REDOX", ex, ey, ew, eh, "Redox constraints", "get_redox_constraints()")
add_node("PHREEQC", ex + 1290, ey, ew, eh, "PHREEQC speciation", "run_phreeqc()\nSI, ionic strength, charge balance")
add_node("RXNDICT", ex, ey + 430, ew, eh, "Reaction dictionary", "mineral reactions, nitrate input,\ndenitrification, cation exchange")
add_node("PENALTIES", ex + 1290, ey + 430, ew, eh, "Auxiliary penalties", "EC/TDS, Gibbs, isotope consistency, kinetic feasibility")
add_node("BOUNDS", ex + 645, ey + 910, ew, eh, "Thermodynamic bounds", "build_edge_bounds()\ndissolution / precipitation / free")

fx, fy = 280, 3180
fw, fh = 1320, 540
fgap = 92
add_node("FITNET", fx, fy, fw, fh, "fit_network()", "iterate primary edges")
add_node("EDGELOOP", fx + (fw + fgap), fy, fw, fh, "Per-edge preparation", "u/v sample lookup, ion vectors, lateral endmembers, tau estimate")
add_node("TRANSPORT", fx + 2 * (fw + fgap), fy, fw, fh, "Transport candidates", "evaporation gamma\nmixing fraction f")
add_node("LASSO", fx + 3 * (fw + fgap), fy, fw, fh, "Sparse reaction fitting", "fit_reactions()\ncoordinate descent + L1 penalty")
add_node("OBJECTIVE", fx + 4 * (fw + fgap), fy, fw, fh, "Candidate objective", "residual + L1 + PHREEQC/redox bounds + EC/TDS + isotope + Gibbs + kinetic penalties")
add_node("EDGERES", fx + 5 * (fw + fgap), fy, fw, fh, "EdgeResult", "transport model, probabilities, reaction extents, residuals, constraints, QC flags, edge metadata")

ox, oy = 300, 5020
ow, oh = 1240, 310
og = 120
add_node("CAL", ox, oy, ow, oh, "Calibration sidecar", "hydrosheaf-cal / PESTGLM\ntransport, vadose, kinetic, isotope adapters")
add_node("UQ", ox + ow + og, oy, ow, oh, "Uncertainty", "bootstrap, Monte Carlo, Bayesian PyMC")
add_node("NITRATE", ox + 2 * (ow + og), oy, ow, oh, "Nitrate source inference", "dual isotope / MCMC priority or hydrochemical evidence fusion")
add_node("NUCLEAR", ox + ow / 2, oy + 560, ow, oh, "Network age inference", "tritium / nuclide age model\nBayesian network aging")
add_node("RTVAL", ox + ow / 2 + ow + og, oy + 560, ow, oh, "Reactive transport validation", "PHREEQC kinetic forward check\nRMSE, NSE, PBIAS, consistency")

hx, hy = 4800, 5020
add_node("TABLES", hx, hy, 1340, 320, "Tables and summaries", "edge_results_table(), summarize_network()")
add_node("PLOTS", hx + 1500, hy, 1180, 320, "Plots", "Gibbs, ILR, TTD, breakthrough, posterior ridges, 3D network")
add_node("REPORT", hx + 2860, hy, 1180, 320, "Interpretation reports", "process maps, diagnostics, science-facing results")
add_node("EXPORT", hx + 220, hy + 570, 900, 250, "Exports", "JSON, CSV")


for group in GROUPS.values():
    draw_group(*group)

# Entry and input routes.
for source in ["CLI", "API", "AUTO"]:
    elbow(source, "CFG", "bottom", "top", fill=EDGE, width=3)
elbow("AUTO", "SAMPLES", "right", "left", fill=EDGE, width=3, mid=2750)
elbow("SAMPLES", "TS", "bottom", "top", fill=EDGE, width=3)
elbow("EDGES", "VADOSE", "bottom", "top", fill=EDGE, width=3)
elbow("ENDM", "PRIORS", "bottom", "top", fill=EDGE, width=3)

# Data-conditioning flow.
elbow("UNITS", "SCHEMA", "bottom", "top", fill=EDGE_DARK, width=4)
elbow("SCHEMA", "QC", "bottom", "top", fill=EDGE_DARK, width=4)
elbow("QC", "AUTODISABLE", "bottom", "top", fill=EDGE_DARK, width=4)
route([anchor("SAMPLES", "bottom"), (3900, 1180), (2440, 1180), (2440, 1510), anchor("UNITS", "right")], fill=EDGE, width=3)

# Topology and prior-engine routes.
elbow("TS", "HEADS", "bottom", "top", fill=EDGE, width=3)
elbow("VADOSE", "VADOSEPRI", "bottom", "top", fill=EDGE, width=3)
elbow("PRIORS", "PHYSAPPLY", "bottom", "top", fill=EDGE, width=3)
elbow("HEADS", "GRAPH2D", "bottom", "top", fill=EDGE_DARK, width=4)
elbow("HEADS", "GRAPH3D", "right", "left", fill=EDGE, width=3, mid=3790)
elbow("VADOSEPRI", "GRAPH3D", "bottom", "top", fill=EDGE, width=3)
elbow("TEMPPRI", "PHYSAPPLY", "bottom", "top", fill=EDGE, width=3)
elbow("GRAPH2D", "GRAPH3D", "right", "left", fill=EDGE_DARK, width=4, both=True)
elbow("GRAPH3D", "PHYSAPPLY", "right", "left", fill=EDGE_DARK, width=4)
elbow("GRAPH2D", "SHEAFTOPO", "bottom", "top", fill=EDGE, width=3)
elbow("GRAPH3D", "SHEAFTOPO", "bottom", "top", fill=EDGE, width=3)
elbow("PHYSAPPLY", "SHEAFTOPO", "bottom", "top", fill=EDGE, width=3)
route([anchor("SHEAFTOPO", "bottom"), (4095, 2910), (940, 2910), anchor("FITNET", "top")], fill=EDGE, width=3)
route([anchor("EDGES", "bottom"), (5850, 1180), (3230, 1180), anchor("GRAPH2D", "top")], fill=EDGE, width=3)

# Geochemical-context routes.
elbow("REDOX", "PHREEQC", "right", "left", fill=EDGE_DARK, width=4)
elbow("REDOX", "RXNDICT", "bottom", "top", fill=EDGE, width=3)
elbow("PHREEQC", "PENALTIES", "bottom", "top", fill=EDGE, width=3)
elbow("RXNDICT", "PENALTIES", "right", "left", fill=EDGE, width=3)
elbow("RXNDICT", "BOUNDS", "bottom", "top", fill=EDGE, width=3)
elbow("PHREEQC", "BOUNDS", "bottom", "top", fill=EDGE, width=3)
route([anchor("BOUNDS", "bottom"), (7485, 2910), (2350, 2910), anchor("EDGELOOP", "top")], fill=EDGE, width=3)
route([anchor("PENALTIES", "bottom"), (8130, 2910), (6590, 2910), anchor("OBJECTIVE", "top")], fill=EDGE, width=3)

# Core inverse-solver route.
elbow("FITNET", "EDGELOOP", "right", "left", fill=EDGE_DARK, width=5, both=True)
elbow("EDGELOOP", "TRANSPORT", "right", "left", fill=EDGE_DARK, width=5, both=True)
elbow("TRANSPORT", "LASSO", "right", "left", fill=EDGE_DARK, width=5)
elbow("LASSO", "OBJECTIVE", "right", "left", fill=EDGE_DARK, width=5)
elbow("OBJECTIVE", "EDGERES", "right", "left", fill=EDGE_DARK, width=5, both=True)

# Optional validation and output routes.
route([anchor("CAL", "top"), (920, 4680), (920, 1010), anchor("CFG", "bottom")], fill=EDGE, width=3, dashed=True)
elbow("EDGERES", "UQ", "bottom", "top", fill=EDGE, width=3)
elbow("EDGERES", "NITRATE", "bottom", "top", fill=EDGE, width=3)
elbow("NITRATE", "RTVAL", "bottom", "top", fill=EDGE, width=3)
elbow("RTVAL", "NUCLEAR", "left", "right", fill=EDGE, width=3)
route([anchor("EDGERES", "bottom"), (8000, 4680), (5500, 4680), anchor("TABLES", "top")], fill=EDGE, width=3)
elbow("NITRATE", "TABLES", "right", "left", fill=EDGE, width=3)
elbow("RTVAL", "TABLES", "right", "left", fill=EDGE, width=3)
elbow("UQ", "TABLES", "right", "left", fill=EDGE, width=3)
elbow("TABLES", "EXPORT", "bottom", "top", fill=EDGE_DARK, width=4)
elbow("TABLES", "PLOTS", "right", "left", fill=EDGE, width=3)
elbow("PLOTS", "REPORT", "right", "left", fill=EDGE, width=3)

for key in nodes:
    draw_node(key)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
img.save(OUTPUT_PATH, dpi=(300, 300))
