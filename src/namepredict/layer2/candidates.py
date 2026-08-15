"""Layer2 母体候选收集（由 layer2.scoring 评分，P-44）；环母体来自 `rule_driven_parent_candidates`；kind_registry 中原 `ring_producers`/`unsat_producers` 引导链已作为死代码移除。"""
from __future__ import annotations

from namepredict.layer2.parent_candidate import with_principal_group_contract
from namepredict.layer2.principal_parent import rule_driven_parent_candidates


def _candidate_key(candidate: dict) -> tuple:
    """构造候选去重键：(kind, chain)。"""
    return candidate.get("kind"), tuple(candidate.get("chain") or [])


def _dedupe_parents(cands: list[dict]) -> list[dict]:
    """按 (kind, chain) 键去重，并应用主官能团契约。"""
    seen: set[tuple] = set()
    out: list[dict] = []
    for raw in cands:
        candidate = with_principal_group_contract(raw)
        key = _candidate_key(candidate)
        if key not in seen:
            seen.add(key)
            out.append(candidate)
    return out


def _candidate_result(raw: list[dict]) -> list[dict]:
    """候选过滤后透传（当前无额外处理）。"""
    return raw


def _collect_candidates(info: dict) -> list[dict]:
    """收集并去重规则驱动的母体候选（剔除 None）。"""
    raw = rule_driven_parent_candidates(info)
    return _dedupe_parents(_candidate_result([c for c in raw if c is not None]))
