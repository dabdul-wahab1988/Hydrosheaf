from PIL import Image, ImageDraw, ImageFont
import math

def orthogonal_points(a, b, start_side='right', end_side='left', mid_shift=0.5):
    x1, y1, w1, h1 = a
    x2, y2, w2, h2 = b
    def get_pt(side, x, y, w, h):
        if side == 'right': return x + w, y + h / 2
        if side == 'left': return x, y + h / 2
        if side == 'top': return x + w / 2, y
        if side == 'bottom': return x + w / 2, y + h
        return x, y
    p1 = get_pt(start_side, x1, y1, w1, h1)
    p2 = get_pt(end_side, x2, y2, w2, h2)
    if start_side in ['left', 'right']:
        xm = p1[0] + (p2[0] - p1[0]) * mid_shift
        return [p1, (xm, p1[1]), (xm, p2[1]), p2]
    else:
        ym = p1[1] + (p2[1] - p1[1]) * mid_shift
        return [p1, (p1[0], ym), (p2[0], ym), p2]

nodes = {}
def add_node(key, x, y, w, h, title, subtitle):
    nodes[key] = (x, y, w, h, title, subtitle)

def add_edge(u, v, color, width, dash=False, start_side='right', end_side='left', mid_shift=0.5, label=None, label_offset=(0,0)):
    a = nodes[u][:4]
    b = nodes[v][:4]
    pts = orthogonal_points(a, b, start_side=start_side, end_side=end_side, mid_shift=mid_shift)
    return dict(a=a, b=b, points=pts, color=color, width=width, dash=dash, label=label, label_offset=label_offset)

EDGE = '#1E293B'
EDGE_LIGHT = '#94A3B8'

# --- Node placement ---
add_node('CLI', 250, 320, 700, 240, 'Command line', 'hydrosheaf.cli:main')
add_node('API', 1040, 320, 760, 240, 'Python API', 'fit_network_pipeline()\nfit_network_with_priors()')
add_node('AUTO', 1890, 320, 700, 240, 'Auto workflow', 'workflows.auto.analyze_dataset()')
add_node('CFG', 320, 730, 2160, 280, 'Config dataclass', 'ion order, weights, feature flags, PHREEQC, isotopes, UQ, 3D, temporal')

# Input Data
x0, y0 = 3020, 320
bw, bh = 1800, 260
gapx, gapy = 130, 130
add_node('SAMPLES', x0, y0, bw, bh, 'Sample chemistry', 'site_id, major ions, pH, EC/TDS, isotopes, coordinates, heads')
add_node('EDGES', x0 + bw + gapx, y0, bw, bh, 'User edges or candidate edges', 'u, v, edge_id, attrs')
add_node('ENDM', x0 + 2 * (bw + gapx), y0, bw, bh, 'Endmembers', 'manual, JSON, sample-derived, latent virtual endmembers')
add_node('TS', x0, y0 + bh + gapy, bw, bh, 'Time-series data', 'TemporalNode, TimeSeriesSample')
add_node('VADOSE', x0 + bw + gapx, y0 + bh + gapy, bw, bh, 'Vadose inputs', 'profile, forcing, water table, vadose links')
add_node('PRIORS', x0 + 2 * (bw + gapx), y0 + bh + gapy, bw, bh, 'External physics priors', 'MODPATH / CSV / JSON\np_uv, tau mean/std/quantiles')

# Data Conditioning
cx, cy = 320, 1370
cw, ch = 1940, 280
for i, (k, t, b) in enumerate([
    ('UNITS', 'Unit handling', 'mg/L, meq/L to mmol/L'),
    ('SCHEMA', 'Schema parsing and validation', 'parse_numeric(), vector_from_sample()'),
    ('QC', 'QC checks', 'missing values, detection limits, charge balance, non-negative chemistry'),
    ('AUTODISABLE', 'Auto-disable missing modules', 'PHREEQC, isotope, nitrate-source flags'),
]):
    add_node(k, cx, cy + i * (ch + 90), cw, ch, t, b)

# Topology and prior engines
px, py = 2730, 1390
tw, th = 1000, 280
add_node('HEADS', px, py, tw, th, 'Head inference', 'heuristic, Bayesian linear, Bayesian MCMC')
add_node('VADOSEPRI', px + 1130, py, tw, th, 'Vadose Richards-1D priors', 'recharge, TTD, tau summaries')
add_node('TEMPPRI', px + 2260, py, tw, th, 'Temporal residence-time priors', 'cross-correlation, Bayesian lag, TTD convolution, tracer decay')
add_node('GRAPH2D', px, py + 430, tw, th, '2D graph inference', 'infer_edges_probabilistic()\ndistance, heads, p_uv')
add_node('GRAPH3D', px + 1130, py + 430, tw, th, '3D graph inference', 'Node3D, Edge3D, layers, screen overlap, anisotropy')
add_node('PHYSAPPLY', px + 2260, py + 430, tw, th, 'apply_physics_priors()', 'override, merge, or only')
add_node('SHEAFTOPO', px + 565, py + 970, 1600, 300, 'Sheaf topology refinement', 'isotope cost, Cl cost, age cost, global section residual')

# Geochemical context
ex, ey = 6250, 1390
ew, eh = 1180, 280
add_node('REDOX', ex, ey, ew, eh, 'Redox constraints', 'get_redox_constraints()')
add_node('PHREEQC', ex + 1290, ey, ew, eh, 'PHREEQC speciation', 'run_phreeqc()\nSI, ionic strength, charge balance')
add_node('RXNDICT', ex, ey + 430, ew, eh, 'Reaction dictionary', 'minerals, NO3src, denit, cation exchange, geologic bias, logic-gate pruning')
add_node('PENALTIES', ex + 1290, ey + 430, ew, eh, 'Auxiliary penalties', 'EC/TDS, Gibbs, isotope consistency, kinetic feasibility')
add_node('BOUNDS', ex + 645, ey + 910, ew, eh, 'Thermodynamic bounds', 'build_edge_bounds()\ndissolution / precipitation / free')

# Core solver
fx, fy = 280, 3180
fw, fh = 1320, 540
fgap = 92
add_node('FITNET', fx, fy, fw, fh, 'fit_network()', 'iterate primary edges')
add_node('EDGELOOP', fx + (fw + fgap), fy, fw, fh, 'Per-edge preparation', 'u/v sample lookup, ion vectors, lateral endmembers, tau estimate')
add_node('TRANSPORT', fx + 2 * (fw + fgap), fy, fw, fh, 'Transport candidates', 'evaporation gamma\nmixing fraction f')
add_node('LASSO', fx + 3 * (fw + fgap), fy, fw, fh, 'Sparse reaction fitting', 'fit_reactions()\ncoordinate descent + L1 penalty')
add_node('OBJECTIVE', fx + 4 * (fw + fgap), fy, fw, fh, 'Candidate objective', 'residual + L1 + PHREEQC/redox bounds + EC/TDS + isotope + Gibbs + kinetic penalties')
add_node('EDGERES', fx + 5 * (fw + fgap), fy, fw, fh, 'EdgeResult', 'transport model, reaction extents, residuals, metadata')

# Optional
ox, oy = 300, 4620
ow, oh = 1240, 310
og = 120
add_node('CAL', ox, oy, ow, oh, 'Calibration sidecar', 'hydrosheaf-cal / PESTGLM')
add_node('UQ', ox + ow + og, oy, ow, oh, 'Uncertainty', 'bootstrap, Monte Carlo, Bayesian PyMC')
add_node('NITRATE', ox + 2 * (ow + og), oy, ow, oh, 'Nitrate source inference', 'dual isotope evidence fusion')
add_node('NUCLEAR', ox + ow / 2, oy + 560, ow, oh, 'Network age inference', 'tritium / nuclide age model')
add_node('RTVAL', ox + ow / 2 + ow + og, oy + 560, ow, oh, 'Reactive transport validation', 'PHREEQC kinetic forward check')

# Output
hx, hy = 4800, 4620
add_node('TABLES', hx, hy, 1340, 320, 'Tables and summaries', 'edge_results_table(), summarize_network()')
add_node('EXPORT', hx + 220, hy + 570, 900, 250, 'Exports', 'JSON, CSV')
add_node('PLOTS', hx + 1500, hy, 1180, 320, 'Plots', 'Gibbs, ILR, breakthrough, posterior ridges')
add_node('REPORT', hx + 2860, hy, 1180, 320, 'Interpretation reports', 'process maps, diagnostics')

# --- Drawing ---
CANVAS_W, CANVAS_H = 8800, 6000
img = Image.new('RGB', (CANVAS_W, CANVAS_H), '#FFFFFF')
d = ImageDraw.Draw(img)
try:
    F_TITLE = ImageFont.truetype("DejaVuSans-Bold.ttf", 100)
    F_SUB = ImageFont.truetype("DejaVuSans.ttf", 70)
except:
    F_TITLE = ImageFont.load_default()
    F_SUB = ImageFont.load_default()

def draw_node(k):
    x, y, w, h, t, s = nodes[k]
    d.rectangle([x, y, x + w, y + h], fill='#F8FAFC', outline='#1E293B', width=8)
    d.text((x + 40, y + 40), t, font=F_TITLE, fill='#0F172A')
    d.text((x + 40, y + 160), s, font=F_SUB, fill='#475569')

edges = []
for a in ['CLI', 'API', 'AUTO']: edges.append(add_edge(a, 'CFG', EDGE_LIGHT, 8, start_side='bottom', end_side='top'))
for a, b in [('SAMPLES', 'UNITS'), ('UNITS', 'SCHEMA'), ('SCHEMA', 'QC'), ('QC', 'AUTODISABLE')]: edges.append(add_edge(a, b, EDGE_LIGHT, 8, start_side='left', end_side='left', mid_shift=-0.2))
for a, b in [('SCHEMA', 'FITNET'), ('AUTODISABLE', 'FITNET'), ('SHEAFTOPO', 'FITNET')]: edges.append(add_edge(a, b, EDGE_LIGHT, 8, start_side='bottom', end_side='top'))
for a, b in [('FITNET', 'EDGELOOP'), ('EDGELOOP', 'TRANSPORT'), ('TRANSPORT', 'LASSO'), ('LASSO', 'OBJECTIVE'), ('OBJECTIVE', 'EDGERES')]: edges.append(add_edge(a, b, EDGE, 12))

for e in edges:
    d.line(e['points'], fill=e['color'], width=e['width'])
for k in nodes: draw_node(k)

img.save('Manuscript_Fig1_Architecture.png')
