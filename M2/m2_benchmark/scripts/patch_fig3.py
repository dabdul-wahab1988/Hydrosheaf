import re

path = r'C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\scripts\make_publication_figures.py'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# ===== PANEL B FIX =====
old_b = '''      minerals = ["calcite", "gypsum", "halite"]
      plot_data_b = []
      for m in minerals:
          # Construct truth mapping for each edge explicitly (including false edges as 0)
          truth_map = {e["edge_id"]: e["reactions"].get(m, 0.0) for e in truth["generation_edges"]}
          temp = baseline[["edge_id", f"reaction_{m}"]].copy()
          temp["True"] = temp["edge_id"].apply(lambda x: truth_map.get(x, 0.0))
          temp = temp.rename(columns={f"reaction_{m}": "Inferred"})
          temp["Mineral"] = PROCESS_LABEL_MAP.get(m, m.capitalize())
          plot_data_b.append(temp[["Inferred", "True", "Mineral"]])

      df_b = pd.concat(plot_data_b).reset_index(drop=True)
      sns.scatterplot(data=df_b, x="True", y="Inferred", hue="Mineral", palette="viridis",
                      s=80, alpha=0.6, edgecolor="white", linewidth=0.5, ax=ax_b_main)

      lim = max(df_b["True"].max(), df_b["Inferred"].max(), 0.1) * 1.1
      ax_b_main.plot([-0.05, lim], [-0.05, lim], "k--", alpha=0.3, lw=1.5, zorder=0)

      # Marginals
      sns.kdeplot(data=df_b, x="True", ax=ax_b_histx, fill=True, alpha=0.3, hue="Mineral", legend=False)
      sns.kdeplot(data=df_b, y="Inferred", ax=ax_b_histy, fill=True, alpha=0.3, hue="Mineral", legend=False)

      ax_b_histx.axis("off")
      ax_b_histy.axis("off")
      ax_b_main.set_xlabel("True Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
      ax_b_main.set_ylabel("Inferred Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
      ax_b_main.tick_params(labelsize=FONT_TICK)
      ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="upper left")

      # Title on the top axis (histogram) to avoid overlap
      ax_b_histx.set_title("B. Reaction Recovery & Signal Sparsity", fontsize=FONT_TITLE, fontweight="bold", pad=20)'''

new_b = '''      minerals = ["calcite", "gypsum", "halite"]
      plot_data_b = []
      for m in minerals:
          truth_map = {e["edge_id"]: e["reactions"].get(m, 0.0) for e in truth["generation_edges"]}
          temp = baseline[["edge_id", f"reaction_{m}"]].copy()
          temp["True"] = temp["edge_id"].apply(lambda x: truth_map.get(x, 0.0))
          temp = temp.rename(columns={f"reaction_{m}": "Inferred"})
          temp["Mineral"] = PROCESS_LABEL_MAP.get(m, m.capitalize())
          plot_data_b.append(temp[["Inferred", "True", "Mineral"]])

      df_b = pd.concat(plot_data_b).reset_index(drop=True)
      # Only plot edges where the mineral is truth-active (True > 0.01) to avoid visual clutter from false positives
      df_b_active = df_b[df_b["True"] > 0.01].copy()

      sns.scatterplot(data=df_b_active, x="True", y="Inferred", hue="Mineral", palette="viridis",
                      s=90, alpha=0.7, edgecolor="white", linewidth=0.5, ax=ax_b_main)

      lim = max(df_b_active["True"].max(), df_b_active["Inferred"].max(), 0.1) * 1.15
      ax_b_main.plot([0, lim], [0, lim], "k--", alpha=0.4, lw=1.5, zorder=0)

      # Marginals (only for active reactions)
      sns.kdeplot(data=df_b_active, x="True", ax=ax_b_histx, fill=True, alpha=0.3, hue="Mineral", legend=False)
      sns.kdeplot(data=df_b_active, y="Inferred", ax=ax_b_histy, fill=True, alpha=0.3, hue="Mineral", legend=False)

      # Metrics annotation
      r2_b = np.corrcoef(df_b_active["True"], df_b_active["Inferred"])[0, 1] ** 2
      mae_b = np.mean(np.abs(df_b_active["True"] - df_b_active["Inferred"]))
      n_b = len(df_b_active)
      ax_b_main.text(0.05, 0.95, fr"$R^2$ = {r2_b:.2f}" + "\n" + fr"MAE = {mae_b:.2f} mmol/L" + "\n" + fr"$n$ = {n_b}",
                     transform=ax_b_main.transAxes, fontsize=FONT_ANNOTATE, fontweight="bold",
                     verticalalignment='top', bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db", alpha=0.9))

      ax_b_histx.axis("off")
      ax_b_histy.axis("off")
      ax_b_main.set_xlabel("True Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
      ax_b_main.set_ylabel("Inferred Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
      ax_b_main.tick_params(labelsize=FONT_TICK)
      ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")

      # Title on the top axis (histogram) to avoid overlap
      ax_b_histx.set_title("B. Reaction Recovery (Truth-Active Only)", fontsize=FONT_TITLE, fontweight="bold", pad=20)'''

content = content.replace(old_b, new_b)

# ===== PANEL D FIX =====
old_d = '''      if iso_path.exists():
          df_d = pd.read_csv(iso_path)
          # REVIEWER FIX: Plot the actual isotopic SHIFT (Delta d18O) instead of absolute states
          sns.regplot(data=df_d, x="true_shift", y="inf_shift", ax=ax_d, color="#db2777",
                      scatter_kws={'s': 80, 'alpha': 0.6, 'edgecolor': 'white'},
                      line_kws={'color': 'black', 'ls': '--', 'lw': 1.5})

          # 1:1 Line
          d_lim = [min(df_d["true_shift"].min(), df_d["inf_shift"].min()) - 0.2,
                   max(df_d["true_shift"].max(), df_d["inf_shift"].max()) + 0.2]
          ax_d.plot(d_lim, d_lim, "k:", alpha=0.3)

          ax_d.set_title("D. Isotopic Forensic Consistency Check", fontsize=FONT_TITLE, fontweight="bold", pad=20)
          ax_d.set_xlabel(r"True $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
          ax_d.set_ylabel(r"Inferred $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
          ax_d.tick_params(labelsize=FONT_TICK)

          # Add R2 annotation
          r2_iso = np.corrcoef(df_d["true_shift"], df_d["inf_shift"])[0,1]**2
          ax_d.text(0.05, 0.95, fr"$R^2 = {r2_iso:.3f}$", transform=ax_d.transAxes,
                    fontsize=FONT_LEGEND, fontweight="bold", verticalalignment='top',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#d1d5db"))'''

new_d = '''      if iso_path.exists():
          df_d = pd.read_csv(iso_path)
          # Plot all points, but color by whether true shift is non-zero
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
          ax_d.set_xlabel(r"True $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
          ax_d.set_ylabel(r"Inferred $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
          ax_d.tick_params(labelsize=FONT_TICK)

          # Add metrics annotation
          r2_iso_all = np.corrcoef(df_d["true_shift"], df_d["inf_shift"])[0,1]**2
          mae_iso = np.mean(np.abs(df_d["true_shift"] - df_d["inf_shift"]))
          n_active = (df_d["true_shift"].abs() > 0.01).sum()
          if len(active) > 1:
              r2_iso_pos = np.corrcoef(active["true_shift"], active["inf_shift"])[0,1]**2
              stats_text = fr"$R^2$ (all) = {r2_iso_all:.2f}" + "\n" + fr"$R^2$ (active) = {r2_iso_pos:.2f}" + "\n" + fr"MAE = {mae_iso:.2f} permil" + "\n" + fr"$n_{{active}}$ = {n_active}"
          else:
              stats_text = fr"$R^2$ = {r2_iso_all:.2f}" + "\n" + fr"MAE = {mae_iso:.2f} permil"
          ax_d.text(0.05, 0.95, stats_text, transform=ax_d.transAxes,
                    fontsize=FONT_ANNOTATE, fontweight="bold", verticalalignment='top',
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db", alpha=0.9))
          ax_d.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")'''

content = content.replace(old_d, new_d)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)
print('Patched Fig3 panels B and D')
