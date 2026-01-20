"""
WebSocket connection manager for real-time analysis progress updates.
"""

import asyncio
from typing import Dict, Set
from fastapi import WebSocket
import json
import logging

logger = logging.getLogger(__name__)


class ConnectionManager:
    """
    Manages WebSocket connections for real-time progress updates.
    Clients connect with a job_id to receive updates for that specific job.
    """

    def __init__(self):
        # Map of job_id -> set of connected WebSockets
        self.active_connections: Dict[str, Set[WebSocket]] = {}
        # Lock for thread-safe operations
        self._lock = asyncio.Lock()

    async def connect(self, websocket: WebSocket, job_id: str):
        """Accept a WebSocket connection and register it for a job."""
        await websocket.accept()
        async with self._lock:
            if job_id not in self.active_connections:
                self.active_connections[job_id] = set()
            self.active_connections[job_id].add(websocket)
        logger.info(f"WebSocket connected for job {job_id}")

    async def disconnect(self, websocket: WebSocket, job_id: str):
        """Remove a WebSocket connection."""
        async with self._lock:
            if job_id in self.active_connections:
                self.active_connections[job_id].discard(websocket)
                if not self.active_connections[job_id]:
                    del self.active_connections[job_id]
        logger.info(f"WebSocket disconnected for job {job_id}")

    async def broadcast_progress(
        self, job_id: str, progress: int, step: str, status: str = "running"
    ):
        """
        Broadcast progress update to all connected clients for a job.

        Args:
            job_id: The analysis job ID
            progress: Progress percentage (0-100)
            step: Current step description
            status: Job status (pending, running, completed, failed)
        """
        async with self._lock:
            connections = self.active_connections.get(job_id, set()).copy()

        if not connections:
            return

        message = json.dumps(
            {
                "type": "progress",
                "job_id": job_id,
                "progress": progress,
                "step": step,
                "status": status,
            }
        )

        # Send to all connected clients
        disconnected = []
        for websocket in connections:
            try:
                await websocket.send_text(message)
            except Exception as e:
                logger.warning(f"Failed to send WebSocket message: {e}")
                disconnected.append(websocket)

        # Clean up disconnected clients
        if disconnected:
            async with self._lock:
                for ws in disconnected:
                    if job_id in self.active_connections:
                        self.active_connections[job_id].discard(ws)

    async def send_completion(self, job_id: str, status: str, error: str = None):
        """
        Send completion message to all connected clients for a job.

        Args:
            job_id: The analysis job ID
            status: Final status (completed or failed)
            error: Error message if failed
        """
        async with self._lock:
            connections = self.active_connections.get(job_id, set()).copy()

        if not connections:
            return

        message = {
            "type": "completion",
            "job_id": job_id,
            "status": status,
            "progress": 100 if status == "completed" else 0,
        }

        if error:
            message["error"] = error

        message_str = json.dumps(message)

        for websocket in connections:
            try:
                await websocket.send_text(message_str)
            except Exception as e:
                logger.warning(f"Failed to send completion message: {e}")


# Global instance for use across the application
manager = ConnectionManager()


def broadcast_progress_sync(
    job_id: str, progress: int, step: str, status: str = "running"
):
    """
    Synchronous wrapper for broadcasting progress.
    Creates a new event loop if needed for background tasks.
    """
    try:
        loop = asyncio.get_event_loop()
        if loop.is_running():
            # We're in an async context, schedule the coroutine
            asyncio.create_task(
                manager.broadcast_progress(job_id, progress, step, status)
            )
        else:
            loop.run_until_complete(
                manager.broadcast_progress(job_id, progress, step, status)
            )
    except RuntimeError:
        # No event loop, create one
        asyncio.run(manager.broadcast_progress(job_id, progress, step, status))


def send_completion_sync(job_id: str, status: str, error: str = None):
    """
    Synchronous wrapper for sending completion message.
    """
    try:
        loop = asyncio.get_event_loop()
        if loop.is_running():
            asyncio.create_task(manager.send_completion(job_id, status, error))
        else:
            loop.run_until_complete(manager.send_completion(job_id, status, error))
    except RuntimeError:
        asyncio.run(manager.send_completion(job_id, status, error))
