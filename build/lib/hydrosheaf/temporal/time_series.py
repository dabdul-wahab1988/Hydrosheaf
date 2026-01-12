"""
Time-series data loading and management.
"""

import csv
from datetime import datetime
from typing import Dict, List

from . import TemporalNode, TimeSeriesSample


def load_time_series_csv(
    path: str,
    ion_order: List[str],
    time_column: str = "timestamp",
    time_format: str = "%Y-%m-%d",
    node_id_column: str = "sample_id",
) -> Dict[str, TemporalNode]:
    """
    Load time-series data from CSV.

    Parameters
    ----------
    path : str
        Path to CSV file with columns: node_id, timestamp, Ca, Mg, Na, ...
    ion_order : List[str]
        Ion order matching Hydrosheaf config (default 10 ions)
    time_column : str
        Name of timestamp column
    time_format : str
        strptime format for parsing timestamps
    node_id_column : str
        Column identifying which node each sample belongs to

    Returns
    -------
    Dict[str, TemporalNode]
        Mapping from node_id to TemporalNode with sorted samples

    Implementation Notes
    --------------------
    1. Parse CSV row by row
    2. Group samples by node_id
    3. Parse timestamps using datetime.strptime(row[time_column], time_format)
    4. Extract concentrations in ion_order
    5. Sort each node's samples by timestamp
    6. Return dict of TemporalNode objects
    """
    # Read CSV
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Group by node_id
    node_samples: Dict[str, List[TimeSeriesSample]] = {}

    for row in rows:
        node_id = row.get(node_id_column, "")
        if not node_id:
            continue

        # Parse timestamp
        timestamp_str = row.get(time_column, "")
        if not timestamp_str:
            continue

        try:
            timestamp = datetime.strptime(timestamp_str, time_format)
        except ValueError:
            # Try common alternative formats
            for fmt in ["%Y-%m-%d %H:%M:%S", "%Y/%m/%d", "%m/%d/%Y"]:
                try:
                    timestamp = datetime.strptime(timestamp_str, fmt)
                    break
                except ValueError:
                    continue
            else:
                print(f"Warning: Could not parse timestamp '{timestamp_str}', skipping")
                continue

        # Extract concentrations
        concentrations = []
        for ion in ion_order:
            value_str = row.get(ion, "")
            if value_str == "" or value_str is None:
                concentrations.append(0.0)  # Missing value
            else:
                try:
                    concentrations.append(float(value_str))
                except ValueError:
                    concentrations.append(0.0)

        # Optional metadata
        temperature_c = None
        ph = None
        ec_uS_cm = None
        isotopes = {}

        if "temperature" in row or "temp_c" in row:
            temp_str = row.get("temperature") or row.get("temp_c", "")
            try:
                temperature_c = float(temp_str)
            except ValueError:
                pass

        if "pH" in row or "ph" in row:
            ph_str = row.get("pH") or row.get("ph", "")
            try:
                ph = float(ph_str)
            except ValueError:
                pass

        if "EC" in row or "ec" in row:
            ec_str = row.get("EC") or row.get("ec", "")
            try:
                ec_uS_cm = float(ec_str)
            except ValueError:
                pass

        # Isotopes
        for iso_key in ["18O", "2H", "d18O", "d2H"]:
            if iso_key in row:
                try:
                    isotopes[iso_key] = float(row[iso_key])
                except ValueError:
                    pass

        # Create sample
        sample = TimeSeriesSample(
            sample_id=row.get("sample_id", f"{node_id}_{timestamp_str}"),
            node_id=node_id,
            timestamp=timestamp,
            concentrations=concentrations,
            temperature_c=temperature_c,
            ph=ph,
            ec_uS_cm=ec_uS_cm,
            isotopes=isotopes if isotopes else None,
        )

        if node_id not in node_samples:
            node_samples[node_id] = []
        node_samples[node_id].append(sample)

    # Sort samples by timestamp and create TemporalNode objects
    temporal_nodes = {}
    for node_id, samples in node_samples.items():
        samples.sort(key=lambda s: s.timestamp)
        temporal_nodes[node_id] = TemporalNode(node_id=node_id, samples=samples)

    return temporal_nodes


def compute_node_statistics(node: TemporalNode, ion_order: List[str]) -> None:
    """
    Compute mean, std, and trend for each ion in a temporal node.

    Modifies node in-place by setting:
    - mean_concentration
    - std_concentration
    - trend_coefficients

    Parameters
    ----------
    node : TemporalNode
        Node with time-series samples
    ion_order : List[str]
        Ion order (for reference)
    """
    if not node.samples:
        return

    n_samples = len(node.samples)
    n_ions = len(node.samples[0].concentrations)

    # Extract concentration matrix
    conc_matrix = [sample.concentrations for sample in node.samples]

    # Compute mean and std
    mean_conc = [0.0] * n_ions
    std_conc = [0.0] * n_ions

    for j in range(n_ions):
        values = [conc_matrix[i][j] for i in range(n_samples)]
        mean_conc[j] = sum(values) / n_samples
        if n_samples > 1:
            variance = sum((v - mean_conc[j]) ** 2 for v in values) / (n_samples - 1)
            std_conc[j] = variance**0.5
        else:
            std_conc[j] = 0.0

    node.mean_concentration = mean_conc
    node.std_concentration = std_conc

    # Compute linear trend (simple least squares)
    if n_samples >= 2:
        # Convert timestamps to days since first sample
        t0 = node.samples[0].timestamp
        times = [(s.timestamp - t0).total_seconds() / 86400.0 for s in node.samples]

        trend_coeff = [0.0] * n_ions
        for j in range(n_ions):
            values = [conc_matrix[i][j] for i in range(n_samples)]

            # y = a + b*t, solve for b (trend)
            mean_t = sum(times) / n_samples
            mean_y = mean_conc[j]

            numerator = sum((times[i] - mean_t) * (values[i] - mean_y) for i in range(n_samples))
            denominator = sum((times[i] - mean_t) ** 2 for i in range(n_samples))

            if denominator > 1e-10:
                trend_coeff[j] = numerator / denominator
            else:
                trend_coeff[j] = 0.0

        node.trend_coefficients = trend_coeff


def compute_temporal_gradients(node: TemporalNode) -> List[List[float]]:
    """
    Compute rate of change for each ion at each time step.

    Parameters
    ----------
    node : TemporalNode
        Node with at least 2 time points

    Returns
    -------
    List[List[float]]
        Shape (N_times - 1, N_ions). Gradient in mmol/L/day.

    Mathematical Implementation
    ---------------------------
    For each consecutive pair (t_k, t_{k+1}):
        dC_j/dt ≈ (C_j(t_{k+1}) - C_j(t_k)) / (t_{k+1} - t_k)

    Convert time difference to days.
    """
    if len(node.samples) < 2:
        return []

    gradients = []

    for i in range(len(node.samples) - 1):
        sample_k = node.samples[i]
        sample_k1 = node.samples[i + 1]

        dt_seconds = (sample_k1.timestamp - sample_k.timestamp).total_seconds()
        dt_days = dt_seconds / 86400.0

        if dt_days < 1e-6:
            # Avoid division by zero
            gradients.append([0.0] * len(sample_k.concentrations))
            continue

        grad = [
            (sample_k1.concentrations[j] - sample_k.concentrations[j]) / dt_days
            for j in range(len(sample_k.concentrations))
        ]

        gradients.append(grad)

    return gradients
