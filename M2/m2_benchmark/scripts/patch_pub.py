import sys
path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Fix Figure 4: reduce annotation font sizes
for i, line in enumerate(lines):
    if 'fontsize=8, fontweight="bold",' in line and 'path_effects' in lines[i+1]:
        lines[i] = line.replace('fontsize=8,', 'fontsize=6.5,')
    if 'fontsize=7.5, fontweight="bold",' in line and 'path_effects' in lines[i+1]:
        lines[i] = line.replace('fontsize=7.5,', 'fontsize=6.5,')

# Fix Figure 6: smarter annotation placement
for i, line in enumerate(lines):
    if 'ax_right.annotate(fr"Optimal' in line:
        # Replace the annotation block
        new_block = """        # Smart placement: shift text left if lambda is on the right half of axis
        x_min, x_max = ax.get_xlim()
        mid_geom = (x_min * x_max) ** 0.5
        text_x = best_lambda * 2.0 if best_lambda < mid_geom else best_lambda / 2.5
        text_ha = "left" if best_lambda < mid_geom else "right"
        aicc_min = df["aicc"].min()
        y_range = ax_right.get_ylim()[1] - ax_right.get_ylim()[0]
        ax_right.annotate(fr"Optimal $\lambda$ = {best_lambda:.4f}", xy=(best_lambda, aicc_min), xycoords="data",
                          xytext=(text_x, aicc_min + y_range * 0.12), fontsize=10, fontweight="bold", color="#1f2937",
                          ha=text_ha, va="bottom", bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#d1d5db", alpha=0.9),
                          arrowprops=dict(arrowstyle="->", color="#6b7280", lw=1))\n"""
        lines[i] = new_block
        j = i + 1
        while j < len(lines) and ('bbox=dict' in lines[j] or 'alpha=0.9' in lines[j] or 'xycoords=' in lines[j] or 'rotation=90' in lines[j]):
            lines[j] = ''
            j += 1
        break

# Fix Figure 7: increase width and wrap x labels, add xlabel
for i, line in enumerate(lines):
    if 'fig, ax = plt.subplots(figsize=(10, 11))' in line:
        lines[i] = line.replace('(10, 11)', '(12, 10)')
    if 'Lower Anayari (Basement)' in line and 'set_xticklabels' in line:
        lines[i] = '    ax.set_xticklabels(["Lower Anayari\\n(Basement)", "Talensi\\n(Mining/Agriculture)"], rotation=0, fontweight="bold", fontsize=FONT_LABEL)\n'
    if 'ax.set_ylabel("")' in line and i+1 < len(lines) and 'ax.set_xlabel("")' in lines[i+1]:
        lines[i+1] = '    ax.set_xlabel("Geologic Province", fontsize=FONT_LABEL, fontweight="bold")\n'

with open(path, 'w', encoding='utf-8') as f:
    f.writelines(lines)
print('Patched make_publication_figures.py')
