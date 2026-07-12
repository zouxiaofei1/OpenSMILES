"""Unit tests for LoopController lifecycle and EventBus.unsubscribe."""

from __future__ import annotations

import threading
import time
from pathlib import Path

from agent_loop.events import EventBus
from agent_loop.state import StateStore
from server.deps import LoopController


def _stuck_worker(gate: threading.Event) -> None:
    gate.wait(timeout=30.0)


def _attach_live_thread(ctrl: LoopController, gate: threading.Event) -> threading.Thread:
    t = threading.Thread(target=_stuck_worker, args=(gate,), daemon=True)
    t.start()
    with ctrl._lock:
        ctrl._thread = t
        ctrl._status = "running"
        ctrl._paused = False
    return t


def test_start_refused_while_thread_alive(tmp_path: Path):
    bus = EventBus()
    ctrl = LoopController(StateStore(tmp_path), bus, tmp_path)
    gate = threading.Event()
    t = _attach_live_thread(ctrl, gate)
    try:
        out = ctrl.start()
        assert out.get("error") == "prior thread still alive"
        assert ctrl._thread is t
        assert t.is_alive()
    finally:
        gate.set()
        t.join(timeout=2.0)


def test_force_stop_keeps_live_thread_and_blocks_start(tmp_path: Path):
    bus = EventBus()
    ctrl = LoopController(StateStore(tmp_path), bus, tmp_path)
    gate = threading.Event()
    t = _attach_live_thread(ctrl, gate)
    try:
        out = ctrl.stop(force=True)
        assert out["status"] == "stopping"
        assert ctrl._thread is t
        assert t.is_alive()
        again = ctrl.start()
        assert again.get("error") == "prior thread still alive"
        assert (tmp_path / "STOP").is_file()
    finally:
        gate.set()
        t.join(timeout=2.0)
        time.sleep(0.05)
        assert ctrl.status() == "stopped"


def test_event_bus_unsubscribe_removes_handler():
    bus = EventBus()
    seen: list[str] = []

    def handler(event_type: str, payload: dict) -> None:
        seen.append(event_type)

    bus.subscribe(handler)
    bus.publish("a", {})
    bus.unsubscribe(handler)
    bus.publish("b", {})
    assert seen == ["a"]


def test_stop_writes_stop_under_controller_memory(tmp_path: Path):
    bus = EventBus()
    ctrl = LoopController(StateStore(tmp_path), bus, tmp_path)
    out = ctrl.stop(force=True)
    assert out["status"] == "stopped"
    assert (tmp_path / "STOP").is_file()
    assert (tmp_path / "STOP").read_text(encoding="utf-8").strip() == "1"
