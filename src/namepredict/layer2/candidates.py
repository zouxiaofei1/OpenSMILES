"""Layer2 母体候选收集（P-44），环母体来自规则驱动候选。"""
from __future__ import annotations
from namepredict.layer2.principal_parent import rule_driven_parent_candidates


def _candidate_key(candidate: dict) -> tuple:
    """构造候选去重键：(kind, chain)。"""
    return candidate.get("kind"), tuple(candidate.get("chain") or [])

def _collect_candidates(info: dict) -> list[dict]:
    """收集并去重规则驱动的母体候选（剔除 None）。"""
    raw = rule_driven_parent_candidates(info)
    return [c for c in raw if c is not None]
