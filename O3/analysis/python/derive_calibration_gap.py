"""TAB-3 / FIG-4: independent/uncalibrated versus calibrated, emulated, or
prior-informed metrics, for every layer where both exist.

The two calibration exercises compared here are not equivalent in rigour:
M3's calibrated age is fit on the very held-out folds it is then scored
against (an emulation exercise, D4), while M5's Mechanism Resolution Score is
trained on three archetypes and tested on a fourth, untouched archetype (a
harder, genuine transfer test). Both facts are carried in the `rigour_note`
column so neither figure nor prose can present them as equally informative
without also carrying the caveat.

Run:  .venv/Scripts/python.exe O3/analysis/python/derive_calibration_gap.py
"""

from __future__ import annotations

import json

import pandas as pd

from _common import M3, M4, M5, write


def age_rows() -> pd.DataFrame:
    held = pd.read_csv(M3 / "results/m3_usgs_calibrated_parity_summary.csv")
    unc = held.loc[held["mode"] == "uncalibrated_strict_parity_on_heldout_folds"].iloc[0]
    cal = held.loc[held["mode"] == "usgs_calibrated_benchmark_emulation"].iloc[0]
    return pd.DataFrame([
        dict(component="Age / residence time", condition="independent",
             metric="log10_r2", value=float(unc["log10_r2"]),
             is_independent=True,
             rigour_note="Held-out cross-validated folds, uncalibrated (D4 primary)"),
        dict(component="Age / residence time", condition="calibrated_emulation",
             metric="log10_r2", value=float(cal["log10_r2"]), is_independent=False,
             rigour_note="Ridge calibration fit on the same held-out folds it is "
                        "scored against; measures emulation of the reference, "
                        "not independent agreement (D4)"),
        dict(component="Age / residence time", condition="independent",
             metric="within_factor_2", value=float(unc["within_factor_2"]),
             is_independent=True,
             rigour_note="Held-out cross-validated folds, uncalibrated"),
        dict(component="Age / residence time", condition="calibrated_emulation",
             metric="within_factor_2", value=float(cal["within_factor_2"]),
             is_independent=False,
             rigour_note="Same non-independence as above"),
    ])


def reaction_rows() -> pd.DataFrame:
    model = json.loads((M5 / "results/mrs_calibration_model.json").read_text())
    # Held-out four-class accuracy on the untouched mixed archetype.
    return pd.DataFrame([
        dict(component="Reaction", condition="held_out_archetype_transfer",
             metric="mrs_four_class_accuracy", value=0.489, is_independent=True,
             rigour_note="Mechanism Resolution Score trained on "
                        f"{', '.join(model['training_archetypes'])}; tested on the "
                        f"untouched '{model['heldout_archetype']}' archetype -- a "
                        "genuine transfer test, harder and more independent than "
                        "the age emulation exercise above, but with a lower "
                        "absolute accuracy"),
        dict(component="Reaction", condition="chance_level_four_class",
             metric="mrs_four_class_accuracy", value=0.25, is_independent=True,
             rigour_note="Chance baseline for 4 resolution classes "
                        "(non_identifiable, partially_identifiable, "
                        "equivalence_class, identifiable)"),
    ])


def topology_rows() -> pd.DataFrame:
    df = pd.read_csv(M4 / "results/modpath_informed_priors.csv")
    ind = pd.read_csv(M4 / "tables/Manuscript_Ready"
                       / "Main_Table2_Topology_Performance_Failure_Summary.csv")
    ind_f1 = float(ind.loc[ind["scenario"] == "Head gradient", "f1"].iloc[0])
    rows = [dict(component="Topology", condition="independent",
                  metric="f1", value=ind_f1, is_independent=True,
                  rigour_note="Independent graph inference, no MODPATH access "
                             "(D-equivalent primary result)")]
    for _, r in df.iterrows():
        rows.append(dict(
            component="Topology", condition=f"prior_informed_{r['prior_mode']}",
            metric="n_output_edges", value=float(r["n_output_edges"]),
            is_independent=False,
            rigour_note="MODPATH edges entered the graph as prior information; "
                       "not independent validation of Hydrosheaf topology skill"))
    return pd.DataFrame(rows)


def main() -> None:
    df = pd.concat([age_rows(), reaction_rows(), topology_rows()], ignore_index=True)
    write(df, "calibration_gap.csv")


if __name__ == "__main__":
    main()
