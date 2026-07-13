from __future__ import annotations
import re

_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",").replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    return (name or "").strip()
