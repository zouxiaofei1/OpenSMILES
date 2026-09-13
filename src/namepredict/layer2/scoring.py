"""可组合的母体候选评分（简化的 IUPAC P-44 优先规则）：每个候选母体 dict 被评分为元组"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2 import kind_registry as _kr



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


def principal_key(parent: dict) -> P44Facts:
    """提取母体的主官能团打分键 facts。"""
    return ParentCandidate(parent, P44Facts(None, int(parent["principal_group_count"]))).facts


def _p44_1_1(parent: dict) -> tuple[int, int]:
    """P-44.1.1 键：(FG 等级, 基团个数)。"""
    facts = principal_key(parent)
    return facts.principal_group_class, facts.principal_group_count
