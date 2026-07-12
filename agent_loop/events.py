"""In-process pub/sub event bus for loop → console fan-out."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

Handler = Callable[[str, dict[str, Any]], None]


class EventBus:
    """Simple subscribe/publish bus (same process)."""

    def __init__(self) -> None:
        self._handlers: list[Handler] = []

    def subscribe(self, fn: Handler) -> None:
        """Register a handler(type, payload)."""
        self._handlers.append(fn)

    def publish(self, event_type: str, payload: dict[str, Any]) -> None:
        """Deliver event to all subscribers."""
        for fn in list(self._handlers):
            fn(event_type, payload)
