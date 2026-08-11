"""Shared functional-group eligibility predicates + aliphatic FG filters (L2)."""
from __future__ import annotations


def _c_idxs(entries, n: int) -> list[int] | None:
    xs = [e["c_idx"] for e in entries or [] if "c_idx" in e]
    return xs if len(xs) == n and len(set(xs)) == n else None


def _no_fgs(info: dict, keys: tuple) -> bool:
    return not any(info.get(k) for k in keys)
