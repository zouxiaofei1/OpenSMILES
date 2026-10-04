"""核心类型定义：命名结果等共享数据结构。"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any


@dataclass
class NameResult:
    """一次命名流水线的输出结果（双语名称、成功标志与元数据）。"""
    en: str
    zh: str
    success: bool
    source: str = "iupac"
    time_ms: float = 0.0
    meta: dict[str, Any] = field(default_factory=dict)
