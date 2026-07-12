from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any


@dataclass
class NameResult:
    en: str
    zh: str
    success: bool
    source: str = "iupac"
    time_ms: float = 0.0
    meta: dict[str, Any] = field(default_factory=dict)
