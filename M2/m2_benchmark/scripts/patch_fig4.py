path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# Replace edge annotation block
old_edge = '''            if is_top:
                # Add PSI text label to the edge
                mid_x = (pos_dict[u][0] + pos_dict[v][0]) / 2
                mid_y = (pos_dict[u][1] + pos_dict[v][1]) / 2
                # Fading: Use alpha=0.3 for low-probability labels
                label_alpha = 1.0 if d["psi"] >= 0.50 else 0.3
                ax.annotate(f"{d['psi']*100:.0f}%", (mid_x, mid_y),
                            color=FAMILY_COLORS[d["family"]], fontsize=6.5, fontweight="bold",
                            alpha=label_alpha,
                            path_effects=[path_effects.withStroke(linewidth=2, foreground="white", alpha=label_alpha)])'''

new_edge = '''            if is_top:
                # Add PSI text label only for high-confidence edges to avoid clutter
                if d["psi"] >= 0.60:
                    mid_x = (pos_dict[u][0] + pos_dict[v][0]) / 2
                    mid_y = (pos_dict[u][1] + pos_dict[v][1]) / 2
                    ax.annotate(f"{d['psi']*100:.0f}%", (mid_x, mid_y),
                                color=FAMILY_COLORS[d["family"]], fontsize=6, fontweight="bold",
                                alpha=0.9,
                                path_effects=[path_effects.withStroke(linewidth=2, foreground="white", alpha=0.9)])'''

content = content.replace(old_edge, new_edge)

# Replace node annotation block
old_node = '''        if is_top:
            for node, (x, y) in pos_dict.items():
                if node in G:
                    ax.annotate(node, (x, y), xytext=(6, 6), textcoords='offset points',
                               fontsize=6.5, fontweight="bold",
                               path_effects=[path_effects.withStroke(linewidth=2, foreground="white")])'''

new_node = '''        if is_top:
            # Only label nodes connected to at least one high-PSI edge to reduce overlap
            high_psi_nodes = set()
            for u, v, d in G.edges(data=True):
                if d["psi"] >= 0.60:
                    high_psi_nodes.add(u)
                    high_psi_nodes.add(v)
            for node, (x, y) in pos_dict.items():
                if node in G and node in high_psi_nodes:
                    ax.annotate(node, (x, y), xytext=(5, 5), textcoords='offset points',
                               fontsize=6, fontweight="bold",
                               path_effects=[path_effects.withStroke(linewidth=2, foreground="white")])'''

content = content.replace(old_node, new_node)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print('Patched Fig4 annotation density')
