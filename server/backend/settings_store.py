"""Persistent app settings (LLM config + theme + concurrency).

Stored as JSON at server/settings.json (gitignored; contains the API key).
The full API key is only ever read/written to disk — never exposed through
public() or the API responses (only a last-4 hint).
"""
from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SETTINGS_PATH = ROOT / "server" / "settings.json"

DEFAULTS: dict[str, Any] = {
    "model": "deepseek-v4-pro",
    "base_url": "https://api.deepseek.com",
    "api_key": "",
    "concurrency": 1,
    "theme": "dark",
}


class SettingsStore:
    def __init__(self, path: Path | None = None) -> None:
        self.path = Path(path or SETTINGS_PATH)
        self._lock = threading.RLock()

    def load(self) -> dict[str, Any]:
        with self._lock:
            merged = dict(DEFAULTS)
            if not self.path.is_file():
                return merged
            try:
                data = json.loads(self.path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                return merged
            if isinstance(data, dict):
                for k in DEFAULTS:
                    if data.get(k) is not None:
                        merged[k] = data[k]
            return merged

    def get_all(self) -> dict[str, Any]:
        return self.load()

    def save(self, updates: dict[str, Any]) -> dict[str, Any]:
        with self._lock:
            current = self.load()
            for k, v in updates.items():
                if v is not None and k in DEFAULTS:
                    current[k] = v
            self.path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.path.with_suffix(".json.tmp")
            tmp.write_text(
                json.dumps(current, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            tmp.replace(self.path)
            return current

    def public(self, data: dict[str, Any] | None = None) -> dict[str, Any]:
        """Redacted view: api_key only as configured flag + last-4 hint."""
        data = data or self.load()
        key = str(data.get("api_key") or "")
        return {
            "model": data.get("model"),
            "base_url": data.get("base_url"),
            "api_key_configured": bool(key),
            "api_key_hint": key[-4:] if len(key) >= 4 else None,
            "concurrency": int(data.get("concurrency", 1)),
            "theme": data.get("theme", "dark"),
        }


_store = SettingsStore()


def get_store() -> SettingsStore:
    return _store
