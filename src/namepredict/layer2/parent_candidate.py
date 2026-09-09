"""类型化 P-44 母体 facts 与旧版生产者边界适配器。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG
from namepredict.layer2.principal import legacy_rank
_FIXED_MULTI: dict[str, int] = {}
_DYNAMIC_IDS: dict[str, str] = {}


def _as_fg(kind: str) -> FG | None:
    """kind 字符串转 FG 枚举，非法返回 None。"""
    try:
        return FG(kind)
    except ValueError:
        return None


def _kind_rank(kind: str) -> int:
    """kind 的主官能团兼容等级（旧式 parent 无 pef 时兜底）。"""
    return legacy_rank(_as_fg(kind))


@dataclass(frozen=True, order=True)
class P44Facts:
    """P-44 打分事实：主官能团类等级与个数。"""
    principal_group_class: int
    principal_group_count: int


@dataclass(frozen=True)
class ParentCandidate:
    """候选母体 dict 及其 P-44 打分事实。"""
    parent: dict
    facts: P44Facts


def principal_contract_kind(kind: str) -> str:
    """返回 kind 的主官能团契约模式（dynamic/fixed/single/none）。"""
    if kind in _DYNAMIC_IDS:
        return "dynamic"
    if kind in _FIXED_MULTI:
        return "fixed"
    return "single" if _kind_rank(kind) else "none"


def _legacy_count(parent: dict, kind: str) -> int:
    """按契约模式计算主官能团个数。"""
    mode = principal_contract_kind(kind)
    if mode == "dynamic":
        return len(parent.get(_DYNAMIC_IDS[kind]) or ())
    if mode == "fixed":
        return _FIXED_MULTI[kind]
    return 1 if mode == "single" else 0


def with_principal_group_contract(parent: dict) -> dict:
    """保证候选 dict 带 principal_group_count 字段。"""
    if "principal_group_count" in parent:
        return parent
    kind = parent.get("kind") or ""
    return {**parent, "principal_group_count": _legacy_count(parent, kind)}


def from_parent_dict(parent: dict) -> ParentCandidate:
    """由母体 dict 构造 ParentCandidate 与打分 facts。"""
    kind = parent.get("kind") or ""
    if "principal_group_count" not in parent:
        raise ValueError(f"principal_group_count missing for {kind}")
    pef = parent.get("principal_expression_facts")  # 打分收敛：主官能团等级直接取 actual FG class rank（kind 只兜底，无 facts 时）。
    rank = legacy_rank(pef.group_class) if pef else _kind_rank(kind)
    facts = P44Facts(rank, int(parent["principal_group_count"]))
    return ParentCandidate(parent, facts)


def principal_key(parent: dict) -> P44Facts:
    """提取母体的主官能团打分键 facts。"""
    return from_parent_dict(parent).facts
