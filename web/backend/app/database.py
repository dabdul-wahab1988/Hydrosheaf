"""
SQLite persistence layer for Hydrosheaf web application.
Provides persistent storage for analysis jobs, sample datasets, and networks.
"""

import sqlite3
import json
import os
from datetime import datetime
from typing import Optional, List, Dict, Any
from contextlib import contextmanager
import logging

logger = logging.getLogger(__name__)

# Database file location
DB_PATH = os.path.join(os.path.dirname(__file__), "..", "hydrosheaf_data.db")


@contextmanager
def get_db_connection():
    """Context manager for database connections."""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    try:
        yield conn
    finally:
        conn.close()


def init_db():
    """Initialize database tables."""
    with get_db_connection() as conn:
        cursor = conn.cursor()

        # Analysis jobs table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS analysis_jobs (
                job_id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'pending',
                progress INTEGER DEFAULT 0,
                current_step TEXT,
                dataset_id TEXT,
                config TEXT,
                results TEXT,
                error TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """
        )

        # Sample datasets table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS sample_datasets (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                samples TEXT NOT NULL,
                sample_count INTEGER NOT NULL,
                created_at TEXT NOT NULL
            )
        """
        )

        # Networks table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS networks (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                nodes TEXT NOT NULL,
                edges TEXT NOT NULL,
                inference_results TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """
        )

        conn.commit()
        logger.info("Database initialized successfully")


# ===== Analysis Jobs CRUD =====


def create_job(
    job_id: str, name: str, dataset_id: str, config: Dict[str, Any]
) -> Dict[str, Any]:
    """Create a new analysis job."""
    now = datetime.utcnow().isoformat()
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO analysis_jobs (job_id, name, status, progress, dataset_id, config, created_at, updated_at)
            VALUES (?, ?, 'pending', 0, ?, ?, ?, ?)
        """,
            (job_id, name, dataset_id, json.dumps(config), now, now),
        )
        conn.commit()

    return {
        "job_id": job_id,
        "name": name,
        "status": "pending",
        "progress": 0,
        "dataset_id": dataset_id,
        "config": config,
        "created_at": now,
        "updated_at": now,
    }


def get_job(job_id: str) -> Optional[Dict[str, Any]]:
    """Get an analysis job by ID."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM analysis_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()
        if row:
            return _row_to_job(row)
    return None


def get_all_jobs() -> List[Dict[str, Any]]:
    """Get all analysis jobs."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM analysis_jobs ORDER BY created_at DESC")
        return [_row_to_job(row) for row in cursor.fetchall()]


def update_job_status(
    job_id: str,
    status: str,
    progress: int = None,
    current_step: str = None,
    error: str = None,
) -> bool:
    """Update job status and progress."""
    now = datetime.utcnow().isoformat()
    with get_db_connection() as conn:
        cursor = conn.cursor()

        updates = ["status = ?", "updated_at = ?"]
        params = [status, now]

        if progress is not None:
            updates.append("progress = ?")
            params.append(progress)

        if current_step is not None:
            updates.append("current_step = ?")
            params.append(current_step)

        if error is not None:
            updates.append("error = ?")
            params.append(error)

        params.append(job_id)

        cursor.execute(
            f"""
            UPDATE analysis_jobs SET {', '.join(updates)} WHERE job_id = ?
        """,
            params,
        )
        conn.commit()
        return cursor.rowcount > 0


def update_job_results(job_id: str, results: Dict[str, Any]) -> bool:
    """Update job with analysis results."""
    now = datetime.utcnow().isoformat()
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            UPDATE analysis_jobs SET results = ?, status = 'completed', progress = 100, updated_at = ?
            WHERE job_id = ?
        """,
            (json.dumps(results), now, job_id),
        )
        conn.commit()
        return cursor.rowcount > 0


def delete_job(job_id: str) -> bool:
    """Delete an analysis job."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM analysis_jobs WHERE job_id = ?", (job_id,))
        conn.commit()
        return cursor.rowcount > 0


def _row_to_job(row: sqlite3.Row) -> Dict[str, Any]:
    """Convert a database row to a job dictionary."""
    return {
        "job_id": row["job_id"],
        "name": row["name"],
        "status": row["status"],
        "progress": row["progress"],
        "current_step": row["current_step"],
        "dataset_id": row["dataset_id"],
        "config": json.loads(row["config"]) if row["config"] else {},
        "results": json.loads(row["results"]) if row["results"] else None,
        "error": row["error"],
        "created_at": row["created_at"],
        "updated_at": row["updated_at"],
    }


# ===== Sample Datasets CRUD =====


def create_dataset(
    dataset_id: str, name: str, samples: List[Dict[str, Any]], description: str = None
) -> Dict[str, Any]:
    """Create a new sample dataset."""
    now = datetime.utcnow().isoformat()
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO sample_datasets (id, name, description, samples, sample_count, created_at)
            VALUES (?, ?, ?, ?, ?, ?)
        """,
            (dataset_id, name, description, json.dumps(samples), len(samples), now),
        )
        conn.commit()

    return {
        "id": dataset_id,
        "name": name,
        "description": description,
        "samples": samples,
        "sample_count": len(samples),
        "created_at": now,
    }


def get_dataset(dataset_id: str) -> Optional[Dict[str, Any]]:
    """Get a dataset by ID."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM sample_datasets WHERE id = ?", (dataset_id,))
        row = cursor.fetchone()
        if row:
            return _row_to_dataset(row)
    return None


def get_all_datasets() -> List[Dict[str, Any]]:
    """Get all datasets (without full sample data for listing)."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT id, name, description, sample_count, created_at FROM sample_datasets ORDER BY created_at DESC"
        )
        return [
            {
                "id": row["id"],
                "name": row["name"],
                "description": row["description"],
                "sample_count": row["sample_count"],
                "created_at": row["created_at"],
            }
            for row in cursor.fetchall()
        ]


def delete_dataset(dataset_id: str) -> bool:
    """Delete a dataset."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM sample_datasets WHERE id = ?", (dataset_id,))
        conn.commit()
        return cursor.rowcount > 0


def _row_to_dataset(row: sqlite3.Row) -> Dict[str, Any]:
    """Convert a database row to a dataset dictionary."""
    return {
        "id": row["id"],
        "name": row["name"],
        "description": row["description"],
        "samples": json.loads(row["samples"]) if row["samples"] else [],
        "sample_count": row["sample_count"],
        "created_at": row["created_at"],
    }


# ===== Networks CRUD =====


def create_network(
    network_id: str,
    name: str,
    nodes: List[Dict[str, Any]],
    edges: List[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Create a new network."""
    now = datetime.utcnow().isoformat()
    edges = edges or []
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO networks (id, name, nodes, edges, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?)
        """,
            (network_id, name, json.dumps(nodes), json.dumps(edges), now, now),
        )
        conn.commit()

    return {
        "id": network_id,
        "name": name,
        "nodes": nodes,
        "edges": edges,
        "inference_results": None,
        "created_at": now,
        "updated_at": now,
    }


def get_network(network_id: str) -> Optional[Dict[str, Any]]:
    """Get a network by ID."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM networks WHERE id = ?", (network_id,))
        row = cursor.fetchone()
        if row:
            return _row_to_network(row)
    return None


def get_all_networks() -> List[Dict[str, Any]]:
    """Get all networks."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM networks ORDER BY created_at DESC")
        return [_row_to_network(row) for row in cursor.fetchall()]


def update_network(
    network_id: str,
    nodes: List[Dict[str, Any]] = None,
    edges: List[Dict[str, Any]] = None,
    inference_results: Dict[str, Any] = None,
) -> bool:
    """Update network nodes, edges, or inference results."""
    now = datetime.utcnow().isoformat()
    with get_db_connection() as conn:
        cursor = conn.cursor()

        updates = ["updated_at = ?"]
        params = [now]

        if nodes is not None:
            updates.append("nodes = ?")
            params.append(json.dumps(nodes))

        if edges is not None:
            updates.append("edges = ?")
            params.append(json.dumps(edges))

        if inference_results is not None:
            updates.append("inference_results = ?")
            params.append(json.dumps(inference_results))

        params.append(network_id)

        cursor.execute(
            f"""
            UPDATE networks SET {', '.join(updates)} WHERE id = ?
        """,
            params,
        )
        conn.commit()
        return cursor.rowcount > 0


def delete_network(network_id: str) -> bool:
    """Delete a network."""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM networks WHERE id = ?", (network_id,))
        conn.commit()
        return cursor.rowcount > 0


def _row_to_network(row: sqlite3.Row) -> Dict[str, Any]:
    """Convert a database row to a network dictionary."""
    return {
        "id": row["id"],
        "name": row["name"],
        "nodes": json.loads(row["nodes"]) if row["nodes"] else [],
        "edges": json.loads(row["edges"]) if row["edges"] else [],
        "inference_results": (
            json.loads(row["inference_results"]) if row["inference_results"] else None
        ),
        "created_at": row["created_at"],
        "updated_at": row["updated_at"],
    }
