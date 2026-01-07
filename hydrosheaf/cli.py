"""Command-line interface for hydrosheaf."""

import argparse
import csv
from typing import Dict, Iterable, List, Mapping, Tuple

from .config import Config, DEFAULT_ION_ORDER
from .data.endmembers import load_endmembers_json
from .data.minerals import MINERAL_LIBRARY
from .data.schema import build_endmember_vectors, normalize_samples
from .data.units import mgL_to_mmolL, meqL_to_mmolL
from .inference.network_fit import fit_network, infer_edges
from .models.ec_tds import calibrate_ec_tds
from .isotopes import fit_lmwl
from .outputs.export import export_edge_results_csv, export_edge_results_json
from .outputs.interpret import print_interpretation_report
import json


def _parse_weights(value: str) -> List[float]:
    weights = [float(item) for item in value.split(",")]
    if len(weights) != 10:
        raise ValueError("weights must have 10 comma-separated values.")
    return weights


def _parse_ion_order(value: str) -> List[str]:
    order = [item.strip() for item in value.split(",") if item.strip()]
    if len(order) != 10:
        raise ValueError("ion-order must have 10 comma-separated values.")
    return order


def _parse_endmembers(values: Iterable[str]) -> Dict[str, List[float]]:
    endmembers: Dict[str, List[float]] = {}
    for entry in values:
        if ":" not in entry:
            raise ValueError("endmember format must be name:v1,v2,...")
        name, raw = entry.split(":", 1)
        vector = [float(item) for item in raw.split(",")]
        if len(vector) != 10:
            raise ValueError(f"endmember '{name}' must have 10 values.")
        endmembers[name] = vector
    return endmembers


def _read_samples(path: str) -> List[Dict[str, object]]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: List[Dict[str, object]] = []
        for row in reader:
            parsed: Dict[str, object] = {}
            for key, value in row.items():
                if value is None:
                    continue
                parsed[key] = value if value != "" else None
            rows.append(parsed)
    return rows


def _read_edges(path: str) -> List[Tuple[str, str, str]]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        edges: List[Tuple[str, str, str]] = []
        for row in reader:
            if "u" not in row or "v" not in row:
                raise ValueError("edges CSV must include 'u' and 'v' columns.")
            edge_id = row.get("edge_id") or f"{row['u']}->{row['v']}"
            edges.append((edge_id, row["u"], row["v"]))
    return edges


def _required_columns(ion_order: List[str]) -> List[str]:
    return ["sample_id", "site_id", *ion_order, "EC", "TDS", "pH"]


def _validate_required(samples: List[Mapping[str, object]], ion_order: List[str]) -> None:
    required = _required_columns(ion_order)
    for entry in samples:
        missing = [col for col in required if col not in entry]
        if missing:
            raise ValueError(f"Sample {entry.get('sample_id', '?')} missing columns: {missing}")


def _convert_samples(
    samples: List[Dict[str, object]],
    ion_order: List[str],
    unit: str,
) -> List[Dict[str, object]]:
    if unit == "mmol/L":
        return samples
    converted: List[Dict[str, object]] = []
    for entry in samples:
        updated = dict(entry)
        for ion in ion_order:
            if ion not in updated:
                continue
            value = updated[ion]
            if value is None:
                continue
            value = float(value)
            if unit == "mg/L":
                updated[ion] = mgL_to_mmolL(value, ion)
            elif unit == "meq/L":
                updated[ion] = meqL_to_mmolL(value, ion)
        converted.append(updated)
    return converted


def _convert_endmembers(
    endmembers: Dict[str, List[float]],
    ion_order: List[str],
    unit: str,
) -> Dict[str, List[float]]:
    if unit == "mmol/L":
        return endmembers
    converted: Dict[str, List[float]] = {}
    for name, vector in endmembers.items():
        new_vector: List[float] = []
        for value, ion in zip(vector, ion_order):
            if unit == "mg/L":
                new_vector.append(mgL_to_mmolL(value, ion))
            elif unit == "meq/L":
                new_vector.append(meqL_to_mmolL(value, ion))
        converted[name] = new_vector
    return converted


def main() -> None:
    parser = argparse.ArgumentParser(description="Run hydrosheaf edge inference.")
    parser.add_argument("--samples", help="Path to samples CSV.")
    parser.add_argument("--edges", help="Path to edges CSV.")
    parser.add_argument("--output", help="Output file path.")
    parser.add_argument("--format", choices=["json", "csv"], default="json")
    parser.add_argument("--interpret", action="store_true", help="Print a geochemical interpretation report after inference.")
    parser.add_argument("--report-only", type=str, help="Print interpretation report from an existing JSON results file and exit.")
    parser.add_argument("--lambda-sparse", type=float, default=0.0)
    parser.add_argument("--lambda-l1", type=float, default=0.0)
    parser.add_argument("--allow-signed", action="store_true")
    parser.add_argument("--reaction-max-iter", type=int, default=300)
    parser.add_argument("--reaction-tol", type=float, default=1e-6)
    parser.add_argument("--charge-balance-limit", type=float, default=0.1)
    parser.add_argument("--ec-tds-penalty-limit", type=float, default=0.0)
    parser.add_argument("--ec-tds-penalty-enabled", action="store_true")
    parser.add_argument(
        "--missing-policy",
        choices=["skip", "impute_zero"],
        default="skip",
        help="How to handle missing ions.",
    )
    parser.add_argument(
        "--detection-policy",
        choices=["half", "zero", "value", "drop"],
        default="half",
        help="How to handle detection limits like '<0.01'.",
    )
    parser.add_argument("--eta-ec", type=float, default=0.0)
    parser.add_argument("--eta-tds", type=float, default=0.0)
    parser.add_argument("--weights", type=str)
    parser.add_argument(
        "--signed-reaction",
        action="append",
        default=[],
        help="Allow signed extent for a reaction label (repeatable).",
    )
    parser.add_argument(
        "--ion-order",
        type=str,
        help="Comma-separated ion order (8 entries).",
    )
    parser.add_argument(
        "--unit",
        choices=["mmol/L", "mg/L", "meq/L"],
        default="mmol/L",
        help="Input sample units for ions.",
    )
    parser.add_argument(
        "--infer-edges",
        action="store_true",
        help="Infer edges from coordinates when edges CSV is not supplied.",
    )
    parser.add_argument(
        "--infer-edges-method",
        choices=["probabilistic", "simple"],
        default="probabilistic",
        help="Edge inference strategy.",
    )
    parser.add_argument("--max-neighbors", type=int, default=1)
    parser.add_argument("--allow-uphill", action="store_true")
    parser.add_argument("--head-key", type=str, default="hydraulic_head")
    parser.add_argument("--elevation-key", type=str, default="elevation")
    parser.add_argument("--edge-p-min", type=float, default=0.75)
    parser.add_argument("--edge-radius-km", type=float, default=5.0)
    parser.add_argument("--edge-max-neighbors", type=int, default=3)
    parser.add_argument("--edge-sigma-meas", type=float, default=0.5)
    parser.add_argument("--edge-sigma-dtw", type=float, default=1.0)
    parser.add_argument("--edge-sigma-elev", type=float, default=1.0)
    parser.add_argument("--edge-sigma-topo", type=float, default=10.0)
    parser.add_argument("--edge-gradient-min", type=float, default=1e-4)
    parser.add_argument("--edge-head-key", type=str, default="head_meas")
    parser.add_argument("--edge-dtw-key", type=str, default="dtw")
    parser.add_argument("--edge-elevation-key", type=str, default="elevation")
    parser.add_argument("--edge-aquifer-key", type=str, default="aquifer_unit")
    parser.add_argument("--edge-screen-depth-key", type=str, default="screen_depth")
    parser.add_argument("--edge-well-depth-key", type=str, default="well_depth")
    parser.add_argument("--edge-depth-mismatch", type=float, default=20.0)
    parser.add_argument("--isotope-enabled", action="store_true")
    parser.add_argument("--isotope-weight", type=float, default=1.0)
    parser.add_argument("--isotope-d-excess-weight", type=float, default=0.0)
    parser.add_argument("--isotope-d18o-key", type=str, default="18O")
    parser.add_argument("--isotope-d2h-key", default="2H", help="Column name for delta-2H.")
    parser.add_argument("--auto-lmwl", action="store_true", help="Automatically calibrate LMWL from samples.")
    parser.add_argument("--iso-consistency-weight", type=float, default=0.0, help="Weight for isotope-chemistry consistency penalty.")
    
    # Mineral Library Options
    parser.add_argument("--lmwl-a", type=float)
    parser.add_argument("--lmwl-b", type=float)
    parser.add_argument("--fit-lmwl", action="store_true")
    parser.add_argument(
        "--endmember",
        action="append",
        default=[],
        help="Mixing endmember name:v1,v2,... (8 values).",
    )
    parser.add_argument(
        "--endmember-id",
        action="append",
        default=[],
        help="Mixing endmember site_id from samples (repeatable).",
    )
    parser.add_argument(
        "--endmembers-json",
        type=str,
        help="Path to endmembers.json (see proposal).",
    )
    parser.add_argument(
        "--calibrate-ec-tds",
        action="store_true",
        help="Fit EC/TDS linear models from samples.",
    )
    parser.add_argument("--phreeqc-enabled", action="store_true")
    parser.add_argument("--phreeqc-mode", choices=["phreeqpython", "subprocess"], default="phreeqpython")
    parser.add_argument("--phreeqc-database", type=str, default=Config().phreeqc_database)
    parser.add_argument("--phreeqc-executable", type=str, default="")
    parser.add_argument("--temp-default-c", type=float, default=25.0)
    parser.add_argument("--si-threshold", type=float, default=0.2)
    parser.add_argument("--constraints-hard", action="store_true")
    
    # Nitrate Source Discrim
    parser.add_argument("--nitrate-source-enabled", action="store_true", help="Enable nitrate source discrimination (manure vs fertilizer).")
    
    # Mineral Library Options
    parser.add_argument(
        "--minerals",
        type=str,
        help="Comma-separated list of active minerals (e.g., 'calcite,albite,pyrite_oxidation_aerobic'). by default it uses a standard set.",
    )
    parser.add_argument(
        "--list-minerals",
        action="store_true",
        help="Print available minerals in the library and exit.",
    )

    args = parser.parse_args()

    if args.list_minerals:
        print("Available Minerals in Library:")
        for name, stoich in sorted(MINERAL_LIBRARY.items()):
            formula = ", ".join(f"{ion}:{coeff}" for ion, coeff in stoich.items())
            print(f"  - {name}: {formula}")
        return

    if args.report_only:
        with open(args.report_only, 'r') as f:
            data = json.load(f)
        print_interpretation_report(data)
        return

    if not args.samples or not args.output:
        parser.error("--samples and --output are required unless --report-only or --list-minerals is used.")

    ion_order = _parse_ion_order(args.ion_order) if args.ion_order else DEFAULT_ION_ORDER.copy()
    samples = _read_samples(args.samples)
    _validate_required(samples, ion_order)

    config = Config(
        ion_order=ion_order,
        lambda_sparse=args.lambda_sparse,
        lambda_l1=args.lambda_l1,
        allow_signed_reactions=args.allow_signed,
        reaction_max_iter=args.reaction_max_iter,
        reaction_tol=args.reaction_tol,
        charge_balance_limit=args.charge_balance_limit,
        ec_tds_penalty_limit=args.ec_tds_penalty_limit,
        ec_tds_penalty_enabled=args.ec_tds_penalty_enabled,
        missing_policy=args.missing_policy,
        detection_limit_policy=args.detection_policy,
        eta_ec=args.eta_ec,
        eta_tds=args.eta_tds,
        unit="mmol/L",
        phreeqc_enabled=args.phreeqc_enabled,
        phreeqc_mode=args.phreeqc_mode,
        phreeqc_database=args.phreeqc_database,
        phreeqc_executable=args.phreeqc_executable,
        temp_default_c=args.temp_default_c,
        si_threshold_tau=args.si_threshold,
        constraints_hard=args.constraints_hard,
        edge_p_min=args.edge_p_min,
        edge_radius_km=args.edge_radius_km,
        edge_max_neighbors=args.edge_max_neighbors,
        edge_sigma_meas=args.edge_sigma_meas,
        edge_sigma_dtw=args.edge_sigma_dtw,
        edge_sigma_elev=args.edge_sigma_elev,
        edge_sigma_topo=args.edge_sigma_topo,
        edge_gradient_min=args.edge_gradient_min,
        edge_head_key=args.edge_head_key,
        edge_dtw_key=args.edge_dtw_key,
        edge_elevation_key=args.edge_elevation_key,
        edge_aquifer_key=args.edge_aquifer_key,
        edge_screen_depth_key=args.edge_screen_depth_key,
        edge_well_depth_key=args.edge_well_depth_key,
        edge_depth_mismatch=args.edge_depth_mismatch,
        isotope_enabled=args.isotope_enabled,
        isotope_weight=args.isotope_weight,
        isotope_d_excess_weight=args.isotope_d_excess_weight,
        isotope_d18o_key=args.isotope_d18o_key,
        nitrate_source_enabled=args.nitrate_source_enabled,
    )
    config.isotope_d2h_key = args.isotope_d2h_key
    config.auto_lmwl = args.auto_lmwl
    config.isotope_consistency_weight = args.iso_consistency_weight
    
    # Gibbs Diagram settings
    if args.edges:
        edges = _read_edges(args.edges)
    elif args.infer_edges:
        edges = infer_edges(
            samples,
            max_neighbors=args.max_neighbors,
            allow_uphill=args.allow_uphill,
            head_key=args.head_key,
            elevation_key=args.elevation_key,
            method=args.infer_edges_method,
            config=config,
        )
    else:
        raise ValueError("Provide --edges or use --infer-edges to generate edges.")
    if args.weights:
        config.weights = _parse_weights(args.weights)
    if args.eta_ec or args.eta_tds:
        config.ec_tds_penalty_enabled = True
    if args.isotope_enabled:
        config.isotope_enabled = True
    if args.endmember:
        config.mixing_endmembers = _parse_endmembers(args.endmember)
    samples = normalize_samples(samples, ion_order, config.detection_limit_policy)
    if args.endmember_id:
        endmember_map = build_endmember_vectors(samples, args.endmember_id, ion_order)
        missing_ids = [end_id for end_id in args.endmember_id if end_id not in endmember_map]
        if missing_ids:
            raise ValueError(f"Endmember IDs not found in samples: {missing_ids}")
        for end_id, vector in endmember_map.items():
            config.mixing_endmembers.setdefault(end_id, vector)
    if args.endmembers_json:
        endmember_map, meta = load_endmembers_json(args.endmembers_json)
        if meta.get("ion_order") and meta.get("ion_order") != ion_order:
            raise ValueError("endmembers.json ion_order does not match CLI ion_order")
        for end_id, vector in endmember_map.items():
            config.mixing_endmembers.setdefault(end_id, vector)
    if args.signed_reaction:
        config.signed_reaction_labels = args.signed_reaction
    
    if args.minerals:
        config.active_minerals = [m.strip() for m in args.minerals.split(",") if m.strip()]


    if args.unit != "mmol/L":
        samples = _convert_samples(samples, ion_order, args.unit)
        if config.mixing_endmembers:
            config.mixing_endmembers = _convert_endmembers(
                config.mixing_endmembers,
                ion_order,
                args.unit,
            )

    if args.calibrate_ec_tds:
        calibrate_ec_tds(samples, config)
        config.ec_tds_penalty_enabled = True

    if args.lmwl_a is not None and args.lmwl_b is not None:
        config.lmwl_a = args.lmwl_a
        config.lmwl_b = args.lmwl_b
        config.lmwl_defined = True
        config.isotope_enabled = True
    elif args.fit_lmwl:
        a, b = fit_lmwl(samples, d18o_key=config.isotope_d18o_key, d2h_key=config.isotope_d2h_key)
        config.lmwl_a = a
        config.lmwl_b = b
        config.lmwl_defined = True
        config.isotope_enabled = True

    if config.isotope_enabled and config.auto_lmwl:
        try:
            intercept, slope = fit_lmwl(
                samples,
                d18o_key=config.isotope_d18o_key,
                d2h_key=config.isotope_d2h_key
            )
            config.lmwl_a = intercept
            config.lmwl_b = slope
            config.lmwl_defined = True
            print(f"Auto-calibrated LMWL: d2H = {slope:.2f} * d18O + {intercept:.2f}")
        except Exception as e:
            print(f"Warning: Failed to auto-calibrate LMWL: {e}")

    results = fit_network(samples, edges, config)

    if args.format == "csv":
        export_edge_results_csv(results, args.output)
    else:
        export_edge_results_json(results, args.output)

    if args.interpret:
        # Convert results to list of dicts for the interpreter if it's not already
        from .outputs.tables import edge_results_table
        table_data = edge_results_table(results)
        print_interpretation_report(table_data)


if __name__ == "__main__":
    main()
