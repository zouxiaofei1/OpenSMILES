"""Shared dependencies: EventBus, StateStore, LoopController, namer."""

from __future__ import annotations

import queue
import threading
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

from agent_loop.config import LoopConfig
from agent_loop.events import EventBus
from agent_loop.loop import AgentLoop
from agent_loop.state import StateStore
from namepredict.namer import SMILESNNamer
from namepredict.types import NameResult

MEMORY_ROOT = Path("agent_loop/memory")
SESSIONS_DIR = MEMORY_ROOT / "sessions"
STOP_FILE = MEMORY_ROOT / "STOP"


def get_bus() -> EventBus:
    """Return process-wide EventBus."""
    return _BUS


def get_store() -> StateStore:
    """Return process-wide StateStore."""
    return _STORE


def get_namer() -> SMILESNNamer:
    """Return process-wide SMILESNNamer."""
    return _NAMER


def get_controller() -> "LoopController":
    """Return process-wide LoopController."""
    return _CONTROLLER


def name_result_dict(result: NameResult) -> dict[str, Any]:
    """Serialize NameResult for JSON responses."""
    return asdict(result)


class LoopController:
    """Owns outer loop thread: pause sleeps; stop via STOP file + flag."""

    def __init__(
        self,
        store: StateStore,
        bus: EventBus,
        memory: Path = MEMORY_ROOT,
    ) -> None:
        self.store = store
        self.bus = bus
        self.memory = Path(memory)
        self._lock = threading.Lock()
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._paused = False
        self._status = "stopped"

    def status(self) -> str:
        """Return running | paused | stopped."""
        with self._lock:
            return self._status

    def state_payload(self) -> dict[str, Any]:
        """STATE.json fields plus controller status."""
        data = dict(self.store.load())
        data["status"] = self.status()
        return data

    def start(self) -> dict[str, str]:
        """Start background loop if not already running."""
        with self._lock:
            if self._thread is not None and self._thread.is_alive():
                self._paused = False
                self._status = "running"
                return {"status": self._status}
            self._begin_thread()
            return {"status": self._status}

    def pause(self) -> dict[str, str]:
        """Soft-pause between cycles."""
        with self._lock:
            self._paused = True
            if self._status == "running":
                self._status = "paused"
            return {"status": self._status}

    def stop(self, force: bool = False) -> dict[str, str]:
        """Write STOP file, signal thread, optional join."""
        self._write_stop()
        with self._lock:
            self._stop.set()
            self._paused = False
            self._status = "stopped"
            thread = self._thread
        self._join_thread(thread, force)
        with self._lock:
            self._thread = None
            return {"status": self._status}

    def _begin_thread(self) -> None:
        self._stop.clear()
        self._paused = False
        self._status = "running"
        self._clear_stop()
        self._thread = threading.Thread(target=self._run_loop, daemon=True)
        self._thread.start()

    def _run_loop(self) -> None:
        agent = self._make_agent()
        while not self._stop.is_set():
            if not self._cycle_step(agent):
                break
        self._mark_stopped()

    def _cycle_step(self, agent: AgentLoop) -> bool:
        if self._paused:
            time.sleep(0.2)
            return True
        if self._should_halt(agent):
            return False
        agent.run_once()
        return not self._should_halt(agent)

    def _mark_stopped(self) -> None:
        with self._lock:
            self._status = "stopped"
            self._thread = None

    def _make_agent(self) -> AgentLoop:
        cfg = LoopConfig(cwd=Path(".").resolve(), mock_pi=True)
        return AgentLoop(cfg, bus=self.bus, store=self.store)

    def _should_halt(self, agent: AgentLoop) -> bool:
        if self._stop.is_set():
            return True
        from agent_loop.loop import _should_stop

        return _should_stop(self.store.load(), agent.config, self.memory)

    def _write_stop(self) -> None:
        self.memory.mkdir(parents=True, exist_ok=True)
        STOP_FILE.write_text("1\n", encoding="utf-8")

    def _clear_stop(self) -> None:
        if STOP_FILE.is_file():
            STOP_FILE.unlink()

    def _join_thread(self, thread: threading.Thread | None, force: bool) -> None:
        if thread is None or not thread.is_alive():
            return
        timeout = 0.5 if force else 5.0
        thread.join(timeout=timeout)


_BUS = EventBus()
_STORE = StateStore(MEMORY_ROOT)
_NAMER = SMILESNNamer()
_CONTROLLER = LoopController(_STORE, _BUS, MEMORY_ROOT)


def event_queue() -> queue.Queue:
    """Subscribe bus to a new queue for SSE consumers."""
    q: queue.Queue = queue.Queue(maxsize=256)

    def _handler(event_type: str, payload: dict[str, Any]) -> None:
        try:
            q.put_nowait({"type": event_type, "payload": payload})
        except queue.Full:
            pass

    _BUS.subscribe(_handler)
    return q
