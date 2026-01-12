"""Export helpers."""

import csv
import json
from typing import Iterable, List

from .tables import edge_results_table
from ..inference.edge_fit import EdgeResult


def export_edge_results_csv(results: List[EdgeResult], path: str) -> None:
    rows = edge_results_table(results)
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def export_edge_results_json(results: List[EdgeResult], path: str) -> None:
    rows = edge_results_table(results)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)
