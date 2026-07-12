"""FastAPI app entry: mount API routers for loop control plane.

Run:
  uvicorn server.app:app --host 127.0.0.1 --port 8765
"""

from __future__ import annotations

from fastapi import FastAPI

from server.routes_events import router as events_router
from server.routes_loop import router as loop_router
from server.routes_name import router as name_router
from server.routes_sessions import router as sessions_router

app = FastAPI(title="ChemAgent Control Plane", version="0.1.0")
app.include_router(loop_router)
app.include_router(sessions_router)
app.include_router(name_router)
app.include_router(events_router)


@app.get("/health")
def health() -> dict[str, str]:
    """Liveness probe."""
    return {"status": "ok"}
