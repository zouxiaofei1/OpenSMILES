"""SSE event stream from EventBus via per-client queue."""

from __future__ import annotations

import asyncio
import json
import queue
from collections.abc import AsyncIterator
from typing import Any

from fastapi import APIRouter, Request
from fastapi.responses import StreamingResponse

from server.deps import event_queue

router = APIRouter(prefix="/api/v1", tags=["events"])


def _format_sse(event: dict[str, Any]) -> str:
    etype = event.get("type", "message")
    data = json.dumps(event, ensure_ascii=False)
    return f"event: {etype}\ndata: {data}\n\n"


async def _sse_gen(q: queue.Queue, request: Request) -> AsyncIterator[str]:
    while True:
        if await request.is_disconnected():
            break
        try:
            item = q.get_nowait()
            yield _format_sse(item)
        except queue.Empty:
            yield ": keepalive\n\n"
            await asyncio.sleep(0.5)


@router.get("/events")
async def events_stream(request: Request) -> StreamingResponse:
    """Subscribe to EventBus and stream text/event-stream."""
    q = event_queue()
    return StreamingResponse(
        _sse_gen(q, request),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )
