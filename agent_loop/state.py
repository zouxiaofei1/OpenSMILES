"""Persistent loop state (STATE.json) and progress.md append log."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def default_state() -> dict[str, Any]:
    """Return a fresh STATE payload with design-doc keys."""
    return {
        "iter": 0,
        "best_dual": 0.0,
        "last_dual": 0.0,
        "no_improve": 0,
        "last_cluster": None,
        "target": 0.99,
        "K": 5,
    }


class StateStore:
    """Load/save STATE.json and append lines to progress.md under root."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.state_path = self.root / "STATE.json"
        self.progress_path = self.root / "progress.md"

    def load(self) -> dict[str, Any]:
        """Load STATE.json or return defaults if missing."""
        if not self.state_path.is_file():
            return default_state()
        return self._read_json()

    def save(self, state: dict[str, Any]) -> None:
        """Write STATE.json (creates parent dirs)."""
        self.root.mkdir(parents=True, exist_ok=True)
        text = json.dumps(state, indent=2, ensure_ascii=False) + "\n"
        self.state_path.write_text(text, encoding="utf-8")

    def append_progress(self, line: str) -> None:
        """Append one line (plus newline) to progress.md."""
        self.root.mkdir(parents=True, exist_ok=True)
        with self.progress_path.open("a", encoding="utf-8") as fh:
            fh.write(line.rstrip("\n") + "\n")

    def _read_json(self) -> dict[str, Any]:
        raw = self.state_path.read_text(encoding="utf-8")
        data = json.loads(raw)
        if not isinstance(data, dict):
            raise ValueError("STATE.json must be an object")
        return data
