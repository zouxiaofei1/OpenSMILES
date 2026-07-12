"""Local Anthropic API key store (never log full key)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DEFAULT_SECRETS_PATH = Path("agent_loop/memory/secrets.json")
_KEY = "anthropic_api_key"


class SecretsStore:
    def __init__(self, path: Path | None = None) -> None:
        self.path = Path(path or DEFAULT_SECRETS_PATH)

    def is_configured(self) -> bool:
        return self.get_key() is not None

    def get_key(self) -> str | None:
        data = self._read()
        raw = data.get(_KEY)
        if not isinstance(raw, str):
            return None
        key = raw.strip()
        return key or None

    def set_key(self, key: str) -> None:
        text = (key or "").strip()
        if not text:
            self.clear_key()
            return
        self._write({_KEY: text})

    def clear_key(self) -> None:
        self._write({})

    def public_status(self) -> dict[str, Any]:
        key = self.get_key()
        if key is None:
            return {"configured": False, "hint": None}
        return {"configured": True, "hint": key[-4:] if len(key) >= 4 else key}

    def _read(self) -> dict[str, Any]:
        if not self.path.is_file():
            return {}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {}
        return data if isinstance(data, dict) else {}

    def _write(self, data: dict[str, Any]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
        self.path.write_text(text, encoding="utf-8")
