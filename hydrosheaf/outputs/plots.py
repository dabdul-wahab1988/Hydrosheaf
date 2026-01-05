"""Optional plotting helpers."""

from typing import List, Optional

from ..inference.edge_fit import EdgeResult


def plot_edge_anomalies(results: List[EdgeResult], path: Optional[str] = None) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for plotting.") from exc

    edge_ids = [result.edge_id for result in results]
    values = [result.anomaly_norm for result in results]

    fig, ax = plt.subplots(figsize=(max(6.0, len(edge_ids) * 0.6), 4.0))
    ax.bar(edge_ids, values, color="#4C72B0")
    ax.set_xlabel("Edge")
    ax.set_ylabel("Anomaly norm")
    ax.set_title("Edge anomaly norms")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
    else:
        plt.show()
