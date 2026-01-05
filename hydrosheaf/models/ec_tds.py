"""EC/TDS prediction utilities."""

from typing import Iterable, List, Mapping, Sequence, Tuple

from ..config import Config
from ..data.schema import vector_from_sample


def _sum(values: Iterable[float]) -> float:
    return float(sum(values))


def predict_ec(values: Iterable[float], config: Config) -> float:
    a, b = config.ec_model
    return a * _sum(values) + b


def predict_tds(values: Iterable[float], config: Config) -> float:
    a, b = config.tds_model
    return a * _sum(values) + b


def ec_tds_penalty(values: Iterable[float], obs: Mapping[str, float], config: Config) -> float:
    penalty = 0.0
    if config.eta_ec and "EC" in obs:
        est = predict_ec(values, config)
        penalty += config.eta_ec * (float(obs["EC"]) - est) ** 2
    if config.eta_tds and "TDS" in obs:
        est = predict_tds(values, config)
        penalty += config.eta_tds * (float(obs["TDS"]) - est) ** 2
    return penalty


def predict_ec_tds(values: Iterable[float], config: Config) -> Tuple[float, float]:
    return predict_ec(values, config), predict_tds(values, config)


def fit_linear_model(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float]:
    if not x or not y or len(x) != len(y):
        return 0.0, 0.0
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denom = sum((xi - x_mean) ** 2 for xi in x)
    if denom == 0:
        return 0.0, y_mean
    slope = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y)) / denom
    intercept = y_mean - slope * x_mean
    return slope, intercept


def calibrate_ec_tds(samples: Sequence[Mapping[str, float]], config: Config) -> Config:
    x_values_ec: List[float] = []
    x_values_tds: List[float] = []
    ec_values: List[float] = []
    tds_values: List[float] = []
    for sample in samples:
        values, normalized = vector_from_sample(
            sample,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if values is None:
            continue
        total = _sum(values)
        if normalized.get("EC") is not None:
            x_values_ec.append(total)
            ec_values.append(float(normalized["EC"]))
        if normalized.get("TDS") is not None:
            x_values_tds.append(total)
            tds_values.append(float(normalized["TDS"]))

    if x_values_ec and ec_values:
        config.ec_model = fit_linear_model(x_values_ec, ec_values)
    if x_values_tds and tds_values:
        config.tds_model = fit_linear_model(x_values_tds, tds_values)
    return config
