path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# PANEL B replacements
lines[247] = '    df_b = pd.concat(plot_data_b).reset_index(drop=True)\n'
lines.insert(248, '    df_b_active = df_b[df_b[\"True\"] > 0.01].copy()\n')
lines[250] = '    sns.scatterplot(data=df_b_active, x=\"True\", y=\"Inferred\", hue=\"Mineral\", palette=\"viridis\",\n'
lines[253] = '    lim = max(df_b_active[\"True\"].max(), df_b_active[\"Inferred\"].max(), 0.1) * 1.15\n'
lines[254] = '    ax_b_main.plot([0, lim], [0, lim], \"k--\", alpha=0.4, lw=1.5, zorder=0)\n'
lines[256] = '    sns.kdeplot(data=df_b_active, x=\"True\", ax=ax_b_histx, fill=True, alpha=0.3, hue=\"Mineral\", legend=False)\n'
lines[257] = '    sns.kdeplot(data=df_b_active, y=\"Inferred\", ax=ax_b_histy, fill=True, alpha=0.3, hue=\"Mineral\", legend=False)\n'

metrics_block = '''    # Metrics annotation
    r2_b = np.corrcoef(df_b_active[\"True\"], df_b_active[\"Inferred\"])[0, 1] ** 2
    mae_b = np.mean(np.abs(df_b_active[\"True\"] - df_b_active[\"Inferred\"]))
    n_b = len(df_b_active)
    ax_b_main.text(0.05, 0.95, f\"R2 = {r2_b:.2f}\\nMAE = {mae_b:.2f} mmol/L\\nn = {n_b}\",
                   transform=ax_b_main.transAxes, fontsize=FONT_ANNOTATE, fontweight="bold",
                   verticalalignment='top', bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db", alpha=0.9))
'''
lines.insert(261, metrics_block)
lines[267] = '    ax_b_histx.set_title("B. Reaction Recovery", fontsize=FONT_TITLE, fontweight="bold", pad=20)\n'
lines[265] = '    ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")\n'

with open(path, 'w', encoding='utf-8') as f:
    f.writelines(lines)
print('Patched panel B')
