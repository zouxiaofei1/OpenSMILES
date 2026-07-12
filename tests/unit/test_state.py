"""Unit tests for agent_loop state store and event bus."""

from __future__ import annotations

from pathlib import Path

from agent_loop.config import LoopConfig
from agent_loop.events import EventBus
from agent_loop.state import StateStore, default_state


def test_loop_config_defaults():
    cfg = LoopConfig(data_path=Path("data/merged.json"), cwd=Path("."))
    assert cfg.k == 5
    assert cfg.target_dual == 0.99
    assert cfg.max_iters == 100
    assert cfg.drop_tol == 0.005
    assert cfg.n_rules == 3
    assert cfg.mock_pi is False


def test_state_store_load_missing_returns_defaults(tmp_path: Path):
    store = StateStore(tmp_path)
    state = store.load()
    assert state["iter"] == 0
    assert state["best_dual"] == 0.0
    assert state["last_dual"] == 0.0
    assert state["no_improve"] == 0
    assert state["last_cluster"] is None
    assert "target" in state
    assert "K" in state


def test_state_store_save_and_load_roundtrip(tmp_path: Path):
    store = StateStore(tmp_path)
    payload = default_state()
    payload["iter"] = 3
    payload["best_dual"] = 0.42
    payload["last_dual"] = 0.40
    payload["no_improve"] = 1
    payload["last_cluster"] = "alcohol"
    store.save(payload)
    loaded = store.load()
    assert loaded["iter"] == 3
    assert loaded["best_dual"] == 0.42
    assert loaded["last_dual"] == 0.40
    assert loaded["no_improve"] == 1
    assert loaded["last_cluster"] == "alcohol"
    assert (tmp_path / "STATE.json").is_file()


def test_state_store_append_progress(tmp_path: Path):
    store = StateStore(tmp_path)
    store.append_progress("first line")
    store.append_progress("second line")
    text = (tmp_path / "progress.md").read_text(encoding="utf-8")
    assert "first line" in text
    assert "second line" in text
    assert text.index("first line") < text.index("second line")


def test_event_bus_subscribe_and_publish():
    bus = EventBus()
    seen: list[tuple[str, dict]] = []

    def handler(event_type: str, payload: dict) -> None:
        seen.append((event_type, payload))

    bus.subscribe(handler)
    bus.publish("cycle_start", {"iter": 1})
    bus.publish("cycle_end", {"iter": 1, "decision": "commit"})
    assert seen == [
        ("cycle_start", {"iter": 1}),
        ("cycle_end", {"iter": 1, "decision": "commit"}),
    ]


def test_event_bus_multiple_subscribers():
    bus = EventBus()
    a: list[str] = []
    b: list[str] = []
    bus.subscribe(lambda t, p: a.append(t))
    bus.subscribe(lambda t, p: b.append(t))
    bus.publish("tick", {})
    assert a == ["tick"]
    assert b == ["tick"]
