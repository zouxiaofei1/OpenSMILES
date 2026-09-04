"""环骨架上类型化主基团表达式的能力策略。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG
from namepredict.layer2.ring_scaffold import ScaffoldIdentity


@dataclass(frozen=True)
class RingExpressionPolicy:
    naming_classes: frozenset[str]
    group_class: FG
    relations: frozenset[str]


_POLICIES = (
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    # 未注册稠环(fused_hetero/fused): kind 已正交化, 主 FG 走 L5 fused_tree 组装, 需放开 typed 表达(否则环酮被 _unsupported_typed_ring 拦截)。
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.NITRILE, frozenset({"exocyclic"})),
)


def supports_ring_expression(scaffold: ScaffoldIdentity, facts) -> bool:
    """判断 scaffold 是否支持该主基团的 typed 环表达。"""
    return any(scaffold.naming_class in policy.naming_classes
               and facts.group_class is policy.group_class
               and facts.relation.value in policy.relations for policy in _POLICIES)
