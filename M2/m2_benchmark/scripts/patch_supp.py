path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_supplementary_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

old_block = '''    # 2x2 Grid, one empty panel
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ((ax1, ax2), (ax3, ax4)) = axes

    # (a) RMSE Histogram
    ax1.hist(df["rmse_mmolL"], bins=20, color=COLOR_PRIMARY, alpha=0.7, edgecolor="white")
    ax1.axvline(df["rmse_mmolL"].median(), color=COLOR_ACCENT, ls="--", lw=2,
                label=f"Median: {df['rmse_mmolL'].median():.3f}")
    ax1.set_title("(a) Geochemical Forward RMSE", fontsize=FONT_TITLE, fontweight="bold")
    ax1.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax1.set_ylabel("Frequency", fontsize=FONT_LABEL, fontweight="bold")
    ax1.tick_params(labelsize=FONT_TICK)
    ax1.legend(fontsize=FONT_LEGEND)

    # (b) Model Efficiency (NSE)
    df_clean = df.dropna(subset=['nse'])
    ax2.scatter(df_clean["rmse_mmolL"], df_clean["nse"], s=40, alpha=0.4, color="#0f766e")
    ax2.axhline(0.5, color=COLOR_ACCENT, linestyle="--", lw=1.5, label="Threshold (NSE=0.5)")
    ax2.set_title("(b) Model Efficiency (Nash-Sutcliffe)", fontsize=FONT_TITLE, fontweight="bold")
    ax2.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylabel("NSE", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylim(-0.5, 1.05)
    ax2.tick_params(labelsize=FONT_TICK)
    ax2.legend(loc='lower right', fontsize=FONT_LEGEND)

    # (c) Mineral Saturation Indices
    si_cols = [c for c in df.columns if c.startswith("si_")]
    if si_cols:
        labels = [c.replace("si_", "").capitalize() for c in si_cols]
        ax3.boxplot([df[c].dropna() for c in si_cols], patch_artist=True,
                    boxprops=dict(facecolor="#f1f5f9", color="#334155"),
                    medianprops=dict(color=COLOR_ACCENT))
        ax3.set_xticklabels(labels, rotation=45, fontsize=FONT_TICK)
        ax3.axhline(0, color="black", ls="-", lw=1)
        ax3.axhline(0.2, color="gray", ls=":", alpha=0.5)
        ax3.axhline(-0.2, color="gray", ls=":", alpha=0.5)
        ax3.set_title("(c) Mineral Saturation Indices", fontsize=FONT_TITLE, fontweight="bold")
        ax3.set_ylabel("Saturation Index (SI)", fontsize=FONT_LABEL, fontweight="bold")
        ax3.tick_params(axis='y', labelsize=FONT_TICK)

    # Hide the 4th panel
    ax4.axis("off")

    for ax in [ax1, ax2, ax3]:
        ax.grid(True, which="major", ls="-", alpha=GRID_ALPHA)'''

new_block = '''    # 1x3 Grid for clean supplementary layout
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    ax1, ax2, ax3 = axes

    # (a) RMSE Histogram
    ax1.hist(df["rmse_mmolL"], bins=20, color=COLOR_PRIMARY, alpha=0.7, edgecolor="white")
    ax1.axvline(df["rmse_mmolL"].median(), color=COLOR_ACCENT, ls="--", lw=2,
                label=f"Median: {df['rmse_mmolL'].median():.3f}")
    ax1.set_title("(a) Geochemical Forward RMSE", fontsize=FONT_TITLE, fontweight="bold")
    ax1.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax1.set_ylabel("Frequency", fontsize=FONT_LABEL, fontweight="bold")
    ax1.tick_params(labelsize=FONT_TICK)
    ax1.legend(fontsize=FONT_LEGEND)

    # (b) Model Efficiency (NSE)
    df_clean = df.dropna(subset=['nse'])
    ax2.scatter(df_clean["rmse_mmolL"], df_clean["nse"], s=40, alpha=0.4, color="#0f766e")
    ax2.axhline(0.5, color=COLOR_ACCENT, linestyle="--", lw=1.5, label="Threshold (NSE=0.5)")
    ax2.set_title("(b) Model Efficiency (Nash-Sutcliffe)", fontsize=FONT_TITLE, fontweight="bold")
    ax2.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylabel("NSE", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylim(-0.5, 1.05)
    ax2.tick_params(labelsize=FONT_TICK)
    ax2.legend(loc='lower right', fontsize=FONT_LEGEND)

    # (c) Percent Bias (PBIAS) Distribution
    pbias_clean = df["pbias_percent"].dropna()
    ax3.hist(pbias_clean, bins=20, color=COLOR_WARNING, alpha=0.7, edgecolor="white")
    ax3.axvline(pbias_clean.median(), color=COLOR_ACCENT, ls="--", lw=2,
                label=f"Median: {pbias_clean.median():.1f}%")
    ax3.axvline(0, color="black", ls="-", lw=1)
    ax3.set_title("(c) Percent Bias (PBIAS)", fontsize=FONT_TITLE, fontweight="bold")
    ax3.set_xlabel("PBIAS (%)", fontsize=FONT_LABEL, fontweight="bold")
    ax3.set_ylabel("Frequency", fontsize=FONT_LABEL, fontweight="bold")
    ax3.tick_params(labelsize=FONT_TICK)
    ax3.legend(fontsize=FONT_LEGEND)

    for ax in [ax1, ax2, ax3]:
        ax.grid(True, which="major", ls="-", alpha=GRID_ALPHA)'''

content = content.replace(old_block, new_block)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print('Patched make_supplementary_figures.py')
