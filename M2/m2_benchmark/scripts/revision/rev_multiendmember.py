"""Revision analysis 6 (CAGEO-D-26-00847): multi-endmember mixing + reaction scenario.

Answers Reviewer 2 Major 4: how the framework behaves when mixing is driven by
TWO endmembers while reactions co-occur (halite dissolution), testing whether
halite dissolution, evapoconcentration, and mixing are distinguishable.

Scenario: downstream = (1 - f1 - f2) * upstream + f1*E1 + f2*E2 + reactions
(halite + calcite extents), with E1/E2 the benchmark's own endmembers. The
framework's transport stage offers a single-endmember mix plus evaporation, so
this is a deliberate model-misspecification stress test. We report transport
model selection, mixing-fraction recovery, and reaction-family recovery.

Output: M2/m2_benchmark/results/revision/multiendmember_recovery.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
BENCHMARK_ROOT = PROJECT_ROOT / "M2" / "m2_benchmark"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(BENCHMARK_ROOT / "scripts"))

from run_m2_benchmark import (  # noqa: E402
    build_noiseless_nodes,
    df_to_samples,
    edge_inputs_for_variant,
    load_truth,
    make_config,
    make_realisations,
    node_vectors,
    reaction_maps,
)
from hydrosheaf import fit_network  # noqa: E402

OUT_DIR = BENCHMARK_ROOT / "results" / "revision"

FAMILY = {
    "calcite": "carbonate", "dolomite": "carbonate", "gypsum": "evaporite",
    "halite": "evaporite", "albite": "silicate", "pyrite_net": "redox",
    "sulfate_reduction": "redox", "iron_reduction": "redox", "denit": "redox",
    "NO3src": "nitrate", "CaNa_exch": "exchange", "NaCa_exch": "exchange",
    "MgNa_exch": "exchange", "NaMg_exch": "exchange",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    truth = load_truth(BENCHMARK_ROOT / "config" / "ground_truth.yaml")
    config = make_config(truth)
    nodes = build_noiseless_nodes(truth, config)
    config = make_config(truth, node_vectors(nodes, truth["ion_order"]))
    vectors = node_vectors(nodes, truth["ion_order"])
    rmap, labels = reaction_maps(config)
    stoich = {lbl: np.array(vec, dtype=float) for lbl, vec in rmap.items()}

    E1 = np.array(vectors["LAT_SAL"], dtype=float)
    E2 = np.array(vectors["LAT_CARB"], dtype=float)
    f1_true, f2_true = 0.20, 0.15
    ext_halite, ext_calcite = 0.40, 0.20

    gen_edges = truth["generation_edges"]
    n_edges = len(gen_edges)
    n_realisations = 20
    rng = np.random.default_rng(20260819)

    rows = []
    realisations_df = make_realisations(truth, nodes, n_realisations)
    for realisation in range(n_realisations):
        realised = realisations_df[realisations_df["realisation"] == realisation].copy()
        samples = df_to_samples(realised, truth["ion_order"])
        # rebuild a scenario frame where every generation edge follows the
        # two-endmember + reaction model
        scenario = []
        for edge in gen_edges:
            u, v = str(edge["u"]), str(edge["v"])
            row = realised[realised["site_id"] == u].iloc[0].to_dict()
            row = dict(row)
            row["site_id"] = v
            xu = np.array([float(row.get(ion, 0.0)) for ion in truth["ion_order"]], dtype=float)
            noise = 1.0 + rng.normal(0.0, 0.04, size=len(xu))
            xv = ((1.0 - f1_true - f2_true) * xu + f1_true * E1 + f2_true * E2
                  + ext_halite * stoich["halite"] + ext_calcite * stoich["calcite"]) * noise
            for ion, val in zip(truth["ion_order"], xv):
                row[ion] = max(0.0, float(val))
            scenario.append(row)
        scenario_df = pd.DataFrame(scenario)
        edges = edge_inputs_for_variant(truth, "full", realisation, scenario_df)
        results = fit_network(df_to_samples(scenario_df, truth["ion_order"]), edges, config)
        for res in results:
            zmap = dict(zip(res.z_labels, res.z_extents))
            fam_ext = {}
            for lbl, val in zmap.items():
                fam = FAMILY.get(lbl, lbl)
                fam_ext[fam] = fam_ext.get(fam, 0.0) + abs(float(val))
            rows.append(
                {
                    "realisation": realisation,
                    "edge_id": res.edge_id,
                    "transport_model": res.transport_model,
                    "gamma": res.gamma,
                    "f": res.f,
                    "chemistry_r2": res.chemistry_r2,
                    "objective_score": res.objective_score,
                    "extent_halite": float(zmap.get("halite", 0.0)),
                    "extent_calcite": float(zmap.get("calcite", 0.0)),
                    "family_evaporite_extent": fam_ext.get("evaporite", 0.0),
                    "family_carbonate_extent": fam_ext.get("carbonate", 0.0),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "multiendmember_recovery.csv", index=False)
    summary = pd.DataFrame(
        [
            {"metric": "n_edges", "value": n_edges},
            {"metric": "n_realisations", "value": n_realisations},
            {"metric": "true_mixing_fraction_total", "value": f1_true + f2_true},
            {"metric": "true_halite_extent_mmolL", "value": ext_halite},
            {"metric": "true_calcite_extent_mmolL", "value": ext_calcite},
            {"metric": "mix_model_selected_fraction", "value": float((out.transport_model == "mix").mean())},
            {"metric": "evap_model_selected_fraction", "value": float((out.transport_model == "evap").mean())},
            {"metric": "median_recovered_f", "value": float(out["f"].median())},
            {"metric": "median_chemistry_r2", "value": float(out["chemistry_r2"].median())},
            {"metric": "median_extent_halite_mmolL", "value": float(out["extent_halite"].median())},
            {"metric": "median_extent_calcite_mmolL", "value": float(out["extent_calcite"].median())},
            {"metric": "median_family_evaporite_extent", "value": float(out["family_evaporite_extent"].median())},
            {"metric": "median_family_carbonate_extent", "value": float(out["family_carbonate_extent"].median())},
        ]
    )
    summary.to_csv(OUT_DIR / "multiendmember_recovery_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote {OUT_DIR / 'multiendmember_recovery*.csv'}")


if __name__ == "__main__":
    main()
