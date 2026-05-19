path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# Tighten edge annotation further: only show 0.60 <= psi < 0.95
old = '''            if is_top:
                # Add PSI text label only for high-confidence edges to avoid clutter
                if d["psi"] >= 0.60:
                    mid_x = (pos_dict[u][0] + pos_dict[v][0]) / 2
                    mid_y = (pos_dict[u][1] + pos_dict[v][1]) / 2
                    ax.annotate(f"{d['psi']*100:.0f}%", (mid_x, mid_y),
                                color=FAMILY_COLORS[d["family"]], fontsize=6, fontweight="bold",
                                alpha=0.9,
                                path_effects=[path_effects.withStroke(linewidth=2, foreground="white", alpha=0.9)])'''

new = '''            if is_top:
                # Add PSI text label only for mid-high confidence edges to avoid clutter
                if 0.60 <= d["psi"] < 0.95:
                    mid_x = (pos_dict[u][0] + pos_dict[v][0]) / 2
                    mid_y = (pos_dict[u][1] + pos_dict[v][1]) / 2
                    ax.annotate(f"{d['psi']*100:.0f}%", (mid_x, mid_y),
                                color=FAMILY_COLORS[d["family"]], fontsize=6, fontweight="bold",
                                alpha=0.9,
                                path_effects=[path_effects.withStroke(linewidth=2, foreground="white", alpha=0.9)])'''

content = content.replace(old, new)

# Add white bounding box to node labels for better overlap tolerance
old_node = '''            for node, (x, y) in pos_dict.items():
                if node in G and node in high_psi_nodes:
                    ax.annotate(node, (x, y), xytext=(5, 5), textcoords='offset points',
                               fontsize=6, fontweight="bold",
                               path_effects=[path_effects.withStroke(linewidth=2, foreground="white")])'''

new_node = '''            for node, (x, y) in pos_dict.items():
                if node in G and node in high_psi_nodes:
                    ax.annotate(node, (x, y), xytext=(4, 4), textcoords='offset points',
                               fontsize=5.5, fontweight="bold",
                               bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.75))'''

content = content.replace(old_node, new_node)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print('Patched Fig4 overlap further')
