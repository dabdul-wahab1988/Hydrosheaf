path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
i = 0
while i < len(lines):
    line = lines[i]

    # === PANEL B FIX ===
    if 'df_b = pd.concat(plot_data_b).reset_index(drop=True)' in line:
        new_lines.append(line)
        new_lines.append('    df_b_active = df_b[df_b["True"] > 0.01].copy()\n')
        i += 1
        continue

    if 'sns.scatterplot(data=df_b, x="True", y="Inferred", hue="Mineral", palette="viridis",' in line:
        # Skip the next line (continuation) and replace both
        new_lines.append('    sns.scatterplot(data=df_b_active, x="True", y="Inferred", hue="Mineral", palette="viridis",\n')
        new_lines.append('                    s=90, alpha=0.7, edgecolor="white", linewidth=0.5, ax=ax_b_main)\n')
        i += 2  # skip original continuation line
        continue

    if 'lim = max(df_b["True"].max(), df_b["Inferred"].max(), 0.1) * 1.1' in line:
        new_lines.append('    lim = max(df_b_active["True"].max(), df_b_active["Inferred"].max(), 0.1) * 1.15\n')
        i += 1
        continue

    if 'ax_b_main.plot([-0.05, lim], [-0.05, lim], "k--", alpha=0.3, lw=1.5, zorder=0)' in line:
        new_lines.append('    ax_b_main.plot([0, lim], [0, lim], "k--", alpha=0.4, lw=1.5, zorder=0)\n')
        i += 1
        continue

    if 'sns.kdeplot(data=df_b, x="True", ax=ax_b_histx, fill=True, alpha=0.3, hue="Mineral", legend=False)' in line:
        new_lines.append('    sns.kdeplot(data=df_b_active, x="True", ax=ax_b_histx, fill=True, alpha=0.3, hue="Mineral", legend=False)\n')
        i += 1
        continue

    if 'sns.kdeplot(data=df_b, y="Inferred", ax=ax_b_histy, fill=True, alpha=0.3, hue="Mineral", legend=False)' in line:
        new_lines.append('    sns.kdeplot(data=df_b_active, y="Inferred", ax=ax_b_histy, fill=True, alpha=0.3, hue="Mineral", legend=False)\n')
        i += 1
        continue

    if 'ax_b_histy.axis("off")' in line:
        new_lines.append(line)
        # Insert metrics block after this
        new_lines.append('    # Metrics annotation\n')
        new_lines.append('    r2_b = np.corrcoef(df_b_active["True"], df_b_active["Inferred"])[0, 1] ** 2\n')
        new_lines.append('    mae_b = np.mean(np.abs(df_b_active["True"] - df_b_active["Inferred"]))\n')
        new_lines.append('    n_b = len(df_b_active)\n')
        new_lines.append('    ax_b_main.text(0.05, 0.95, f"R² = {r2_b:.2f}\\nMAE = {mae_b:.2f} mmol/L\\nn = {n_b}",\n')
        new_lines.append('                   transform=ax_b_main.transAxes, fontsize=FONT_ANNOTATE, fontweight="bold",\n')
        new_lines.append("                   verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#d1d5db', alpha=0.9))\n")
        i += 1
        continue

    if 'ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="upper left")' in line:
        new_lines.append('    ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")\n')
        i += 1
        continue

    # === PANEL D FIX ===
    if 'df_d = pd.read_csv(iso_path)' in line:
        # Replace the entire panel D block
        new_lines.append('        df_d = pd.read_csv(iso_path)\n')
        new_lines.append('        df_d["category"] = df_d["true_shift"].apply(lambda x: "Active shift" if abs(x) > 0.01 else "No shift")\n')
        new_lines.append('        sns.scatterplot(data=df_d, x="true_shift", y="inf_shift", hue="category",\n')
        new_lines.append('                        palette={"Active shift": "#db2777", "No shift": "#94a3b8"},\n')
        new_lines.append('                        s=70, alpha=0.6, edgecolor="white", linewidth=0.5, ax=ax_d)\n')
        new_lines.append('\n')
        new_lines.append('        active = df_d[df_d["true_shift"].abs() > 0.01]\n')
        new_lines.append('        if len(active) > 1:\n')
        new_lines.append('            sns.regplot(data=active, x="true_shift", y="inf_shift", ax=ax_d, color="#db2777",\n')
        new_lines.append('                        scatter=False,\n')
        new_lines.append("                        line_kws={'color': 'black', 'ls': '--', 'lw': 1.5})\n")
        new_lines.append('\n')
        new_lines.append('        d_min = min(df_d["true_shift"].min(), df_d["inf_shift"].min()) - 0.3\n')
        new_lines.append('        d_max = max(df_d["true_shift"].max(), df_d["inf_shift"].max()) + 0.3\n')
        new_lines.append('        ax_d.plot([d_min, d_max], [d_min, d_max], "k:", alpha=0.3, lw=1.5)\n')
        new_lines.append('\n')
        new_lines.append('        ax_d.set_title("D. Isotopic Forensic Consistency Check", fontsize=FONT_TITLE, fontweight="bold", pad=20)\n')
        new_lines.append('        ax_d.set_xlabel(r"True $\\Delta \\delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")\n')
        new_lines.append('        ax_d.set_ylabel(r"Inferred $\\Delta \\delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")\n')
        new_lines.append('        ax_d.tick_params(labelsize=FONT_TICK)\n')
        new_lines.append('\n')
        new_lines.append('        r2_iso_all = np.corrcoef(df_d["true_shift"], df_d["inf_shift"])[0,1]**2\n')
        new_lines.append('        mae_iso = np.mean(np.abs(df_d["true_shift"] - df_d["inf_shift"]))\n')
        new_lines.append('        n_active = (df_d["true_shift"].abs() > 0.01).sum()\n')
        new_lines.append('        if len(active) > 1:\n')
        new_lines.append('            r2_iso_pos = np.corrcoef(active["true_shift"], active["inf_shift"])[0,1]**2\n')
        new_lines.append('            stats_text = f"R² (all) = {r2_iso_all:.2f}\\nR² (active) = {r2_iso_pos:.2f}\\nMAE = {mae_iso:.2f} permil\\nn_active = {n_active}"\n')
        new_lines.append('        else:\n')
        new_lines.append('            stats_text = f"R² = {r2_iso_all:.2f}\\nMAE = {mae_iso:.2f} permil"\n')
        new_lines.append('        ax_d.text(0.05, 0.95, stats_text, transform=ax_d.transAxes,\n')
        new_lines.append('                  fontsize=FONT_ANNOTATE, fontweight="bold", verticalalignment="top",\n')
        new_lines.append("                  bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#d1d5db', alpha=0.9))\n")
        new_lines.append('        ax_d.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")\n')
        # Skip all original panel D lines until we hit 'else:'
        i += 1
        while i < len(lines) and 'ax_d.text(0.5, 0.5, "Isotope results not found"' not in lines[i] and 'else:' not in lines[i]:
            i += 1
        continue

    # Also need to skip old panel D r2 annotation and bbox lines
    if '# Add R2 annotation' in line:
        # Skip the next 3 lines (r2_iso calculation + text annotation)
        i += 4
        continue

    new_lines.append(line)
    i += 1

with open(path, 'w', encoding='utf-8') as f:
    f.writelines(new_lines)
print('Patched make_publication_figures.py cleanly')
