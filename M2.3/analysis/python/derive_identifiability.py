"""Structural identifiability of the reaction-extent inverse problem.

Added in response to peer review, which correctly objected that a failed
inversion by one solver on one scenario does not establish that the target is
non-identifiable. This script supplies the structural argument: it builds the
stoichiometry matrix from the declared reaction dictionary and reports its rank,
null-space dimension and conditioning against the observable ion set.

It also reports the leakage fraction as a function of threshold together with a
propagated-measurement-noise baseline, so that the reported leakage can be
compared with what analytical noise alone would produce.

Run:  .venv/Scripts/python.exe M2.3/analysis/python/derive_identifiability.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "M2.3" / "manuscript" / "artifacts" / "data"
BENCH = ROOT / "M2/m2_benchmark"


def stoichiometry_matrix() -> tuple[np.ndarray, list[str], list[str]]:
    """Assemble the ion-by-reaction stoichiometry matrix A."""
    dic = pd.read_csv(BENCH / "data/hydrosheaf_reaction_dictionary.csv")
    truth = yaml.safe_load((BENCH / "config/ground_truth.yaml").read_text(encoding="utf-8"))
    ions = list(truth["ion_order"])
    reactions = list(dic["reaction_label"])
    A = np.zeros((len(ions), len(reactions)))
    for j, raw in enumerate(dic["stoichiometry_mmolL_per_extent"]):
        for ion, coeff in json.loads(raw).items():
            if ion in ions:
                A[ions.index(ion), j] = float(coeff)
    return A, ions, reactions


def identifiability_report() -> pd.DataFrame:
    A, ions, reactions = stoichiometry_matrix()
    n_ion, n_rx = A.shape
    rank = int(np.linalg.matrix_rank(A))
    sv = np.linalg.svd(A, compute_uv=False)
    nonzero = sv[sv > 1e-12]
    cond = float(nonzero.max() / nonzero.min()) if nonzero.size else np.inf

    # Ions that actually vary in the benchmark chemistry constrain the system;
    # an ion absent from every reaction contributes no constraint.
    active_ions = [ions[i] for i in range(n_ion) if np.any(np.abs(A[i]) > 0)]
    A_active = A[[ions.index(i) for i in active_ions], :]
    rank_active = int(np.linalg.matrix_rank(A_active))

    rows = [dict(
        quantity="Reaction dictionary terms", value=n_rx,
        note="Unknown extents to be recovered"),
        dict(quantity="Ions in the declared order", value=n_ion,
             note="Maximum number of mass-balance constraints"),
        dict(quantity="Ions appearing in at least one reaction", value=len(active_ions),
             note="Constraints that actually bind: " + ", ".join(active_ions)),
        dict(quantity="Rank of the stoichiometry matrix", value=rank,
             note="Independent constraints available"),
        dict(quantity="Rank over binding ions only", value=rank_active,
             note="Unchanged by dropping non-participating ions"),
        dict(quantity="Null-space dimension", value=n_rx - rank,
             note="Directions in extent space that leave the chemistry unchanged"),
        dict(quantity="Condition number of the stoichiometry matrix",
             value=round(cond, 2),
             note="Ratio of largest to smallest non-zero singular value"),
    ]
    return pd.DataFrame(rows)


def null_space_examples() -> pd.DataFrame:
    """Concrete combinations of reactions that are chemically indistinguishable."""
    A, ions, reactions = stoichiometry_matrix()
    _, _, vt = np.linalg.svd(A)
    rank = int(np.linalg.matrix_rank(A))
    null = vt[rank:]
    rows = []
    for k, vec in enumerate(null):
        # Report the terms carrying most of the null direction's weight.
        idx = np.argsort(-np.abs(vec))
        keep = [i for i in idx if abs(vec[i]) > 0.15][:5]
        combo = "; ".join(f"{reactions[i]} {vec[i]:+.2f}" for i in keep)
        rows.append(dict(null_direction=k + 1, dominant_terms=combo,
                         max_residual_ion_change=float(np.abs(A @ vec).max())))
    return pd.DataFrame(rows)


def leakage_curve() -> pd.DataFrame:
    """Leakage against threshold, with a propagated-noise baseline."""
    truth = yaml.safe_load((BENCH / "config/ground_truth.yaml").read_text(encoding="utf-8"))
    rel_sigma = float(truth["noise"]["major_ion_rel_sigma"])
    rx = pd.read_csv(BENCH / "results/reaction_recovery.csv")
    inactive = rx[rx["true_extent_mmolL"].abs() <= 0.01]
    active = rx[rx["true_extent_mmolL"].abs() > 0.01]

    # Noise floor: the extent change a single reaction must produce for its ion
    # signature to exceed analytical noise on a typical major-ion concentration.
    A, ions, reactions = stoichiometry_matrix()
    typical_conc = 2.0  # mmol/L, order of magnitude of the benchmark major ions
    noise_mmol = rel_sigma * typical_conc
    per_rx_gain = np.abs(A).max(axis=0)
    detection_floor = float(np.median(noise_mmol / per_rx_gain[per_rx_gain > 0]))

    rows = []
    for thr in [0.01, 0.02, 0.05, 0.10, 0.20, 0.50]:
        rows.append(dict(
            threshold_mmolL=thr,
            leakage_fraction=float((inactive["recovered_extent_mmolL"].abs() > thr).mean()),
            active_terms_detected_fraction=float(
                (active["recovered_extent_mmolL"].abs() > thr).mean()),
            n_inactive=int(len(inactive)), n_active=int(len(active))))
    df = pd.DataFrame(rows)
    df["analytical_noise_floor_mmolL"] = round(detection_floor, 4)
    df["major_ion_relative_sigma"] = rel_sigma
    return df


def replication_structure() -> pd.DataFrame:
    """Distinguish scored rows from independent units of information."""
    rx = pd.read_csv(BENCH / "results/reaction_recovery.csv")
    tr = pd.read_csv(BENCH / "results/transport_recovery.csv")
    age = pd.read_csv(BENCH / "results/age_inference_validation.csv")
    rows = [
        dict(target="Reaction extents", scored_rows=len(rx),
             realisations=int(rx["realisation"].nunique()),
             distinct_units=int(rx.groupby(["edge_id", "reaction_label"]).ngroups),
             distinct_truth_values=int(rx["true_extent_mmolL"].round(6).nunique()),
             truth_varies_across_realisations=bool(
                 rx.groupby(["edge_id", "reaction_label"])
                   ["true_extent_mmolL"].nunique().max() > 1)),
        dict(target="Transport parameters", scored_rows=len(tr),
             realisations=int(tr["realisation"].nunique()),
             distinct_units=int(tr.groupby(["edge_id", "parameter"]).ngroups),
             distinct_truth_values=int(tr["true_value"].round(6).nunique()),
             truth_varies_across_realisations=bool(
                 tr.groupby(["edge_id", "parameter"])["true_value"].nunique().max() > 1)),
        dict(target="Residence time", scored_rows=len(age),
             realisations=int(age["realisation"].nunique()),
             distinct_units=int(age.groupby("site_id").ngroups),
             distinct_truth_values=int(age["true_mrt_years"].round(6).nunique()),
             truth_varies_across_realisations=bool(
                 age.groupby("site_id")["true_mrt_years"].nunique().max() > 1)),
    ]
    return pd.DataFrame(rows)


def clustered_recovery() -> pd.DataFrame:
    """Recovery statistics computed per distinct truth unit, then aggregated.

    Averaging within each (edge, reaction) or (edge, parameter) unit before
    summarising prevents the 100 noise replicates from being counted as
    independent information.
    """
    rx = pd.read_csv(BENCH / "results/reaction_recovery.csv")
    tr = pd.read_csv(BENCH / "results/transport_recovery.csv")

    act = rx[rx["true_extent_mmolL"].abs() > 0.01]
    per_unit = (act.groupby(["edge_id", "reaction_label"])
                   .agg(truth=("true_extent_mmolL", "first"),
                        mean_recovered=("recovered_extent_mmolL", "mean"))
                   .reset_index())
    err = per_unit["mean_recovered"] - per_unit["truth"]
    r2_unit = float(1 - np.sum(err ** 2)
                    / np.sum((per_unit["truth"] - per_unit["truth"].mean()) ** 2))

    t_ok = tr[tr["absolute_error"].notna()]
    t_unit = (t_ok.groupby(["edge_id", "parameter"])
                  .agg(truth=("true_value", "first"),
                       mean_recovered=("recovered_value", "mean"))
                  .reset_index())
    t_err = (t_unit["mean_recovered"] - t_unit["truth"]).abs()

    return pd.DataFrame([
        dict(target="Reaction extents (active)", n_units=len(per_unit),
             statistic="R-squared over unit means", value=round(r2_unit, 4),
             row_level_equivalent=-0.0231),
        dict(target="Transport parameters", n_units=len(t_unit),
             statistic="Median absolute error of unit means",
             value=round(float(t_err.median()), 4),
             row_level_equivalent=0.0577),
    ])


def main() -> None:
    identifiability_report().to_csv(OUT / "reaction_identifiability.csv", index=False)
    null_space_examples().to_csv(OUT / "reaction_null_space.csv", index=False)
    leakage_curve().to_csv(OUT / "reaction_leakage_curve.csv", index=False)
    replication_structure().to_csv(OUT / "benchmark_replication_structure.csv", index=False)
    clustered_recovery().to_csv(OUT / "clustered_recovery.csv", index=False)

    print("== identifiability ==")
    print(identifiability_report().to_string(index=False))
    print("\n== null directions ==")
    print(null_space_examples().to_string(index=False))
    print("\n== leakage vs threshold ==")
    print(leakage_curve().to_string(index=False))
    print("\n== replication structure ==")
    print(replication_structure().to_string(index=False))
    print("\n== clustered recovery ==")
    print(clustered_recovery().to_string(index=False))


if __name__ == "__main__":
    main()
