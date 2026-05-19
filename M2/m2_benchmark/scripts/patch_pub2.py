path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# Fix the edge PSI label font size specifically
old = 'ax.annotate(f"{d[\'psi\']*100:.0f}%", (mid_x, mid_y),\n                            color=FAMILY_COLORS[d["family"]], fontsize=8, fontweight="bold",'
new = 'ax.annotate(f"{d[\'psi\']*100:.0f}%", (mid_x, mid_y),\n                            color=FAMILY_COLORS[d["family"]], fontsize=6.5, fontweight="bold",'
content = content.replace(old, new)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print('Fixed edge label font size')
