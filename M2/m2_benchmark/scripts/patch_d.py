path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Replace panel D block lines 301-319 (0-indexed: 300-318)
new_d = '''    if iso_path.exists():
        df_d = pd.read_csv(iso_path)
        # Categorise points for clarity
        df_d["category"] = df_d["true_shift"].apply(lambda x: "Active shift" if abs(x) > 0.01 else "No shift")
        sns.scatterplot(data=df_d, x="true_shift", y="inf_shift", hue="category",
                        palette={"Active shift": "#db2777", "No shift": "#94a3b8"},
                        s=70, alpha=0.6, edgecolor="white", linewidth=0.5, ax=ax_d)

        # Regression line for active shifts only
        active = df_d[df_d["true_shift"].abs() > 0.01]
        if len(active) > 1:
            sns.regplot(data=active, x="true_shift", y="inf_shift", ax=ax_d, color="#db2777",
                        scatter=False,
                        line_kws={'color': 'black', 'ls': '--', 'lw': 1.5})

        # 1:1 Line
        d_min = min(df_d["true_shift"].min(), df_d["inf_shift"].min()) - 0.3
        d_max = max(df_d["true_shift"].max(), df_d["inf_shift"].max()) + 0.3
        ax_d.plot([d_min, d_max], [d_min, d_max], "k:", alpha=0.3, lw=1.5)

        ax_d.set_title("D. Isotopic Forensic Consistency Check", fontsize=FONT_TITLE, fontweight="bold", pad=20)
        ax_d.set_xlabel(r"True $\\Delta \\delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
        ax_d.set_ylabel(r"Inferred $\\Delta \\delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
        ax_d.tick_params(labelsize=FONT_TICK)

        # Metrics annotation
        r2_iso_all = np.corrcoef(df_d["true_shift"], df_d["inf_shift"])[0,1]**2
        mae_iso = np.mean(np.abs(df_d["true_shift"] - df_d["inf_shift"]))
        n_active = (df_d["true_shift"].abs() > 0.01).sum()
        if len(active) > 1:
            r2_iso_pos = np.corrcoef(active["true_shift"], active["inf_shift"])[0,1]**2
            stats_text = f"R2 all = {r2_iso_all:.2f}\\nR2 active = {r2_iso_pos:.2f}\\nMAE = {mae_iso:.2f} permil\\nn_active = {n_active}"
        else:
            stats_text = f"R2 = {r2_iso_all:.2f}\\nMAE = {mae_iso:.2f} permil"
        ax_d.text(0.05, 0.95, stats_text, transform=ax_d.transAxes,
                  fontsize=FONT_ANNOTATE, fontweight="bold", verticalalignment='top',
                  bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db", alpha=0.9))
        ax_d.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")
'''

# Replace lines 300-318 (inclusive) with new block
lines[300:319] = [new_d]

with open(path, 'w', encoding='utf-8') as f:
    f.writelines(lines)
print('Patched panel D')
