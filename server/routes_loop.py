"""Loop control routes: state, start, pause, stop."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Body
from pydantic import BaseModel, Field

from server.deps import get_controller

router = APIRouter(prefix="/api/v1/loop", tags=["loop"])


class StopBody(BaseModel):
    """Optional force stop."""

    force: bool = Field(default=False)


@router.get("/state")
def loop_state() -> dict[str, Any]:
    """Return STATE.json fields plus controller status."""
    return get_controller().state_payload()


@router.post("/start")
def loop_start() -> dict[str, str]:
    """Start or resume the agent loop in a background thread."""
    return get_controller().start()


@router.post("/pause")
def loop_pause() -> dict[str, str]:
    """Soft-pause between cycles."""
    return get_controller().pause()


@router.post("/stop")
def loop_stop(body: StopBody | None = Body(default=None)) -> dict[str, str]:
    """Stop loop; write STOP file; optional force join."""
    force = bool(body.force) if body else False
    return get_controller().stop(force=force)
