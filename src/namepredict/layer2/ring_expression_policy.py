"""环骨架上类型化主基团表达式的能力策略。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG
from namepredict.layer2.ring_scaffold import ScaffoldIdentity


@dataclass(frozen=True)
class RingExpressionPolicy:
    """一条能力策略：命名类集合 + 主基团类 + 允许的关系集合。"""
    naming_classes: frozenset[str]
    group_class: FG
    relations: frozenset[str]


_CYCLIC_CORE = frozenset({"fused_hetero", "fused", "bridged",
                          "mono_spiro", "fused_bridged_spiro"})  # 多环骨架共用同一套环表达策略

_POLICIES = (
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.KETONE, frozenset({"in_skeleton"})), 
    RingExpressionPolicy(frozenset({"monohetero", "fused56", "purine", "carbazole",  
                                    "acridine", "phenothiazine", "benzodioxole",
                                    "anthra", "phenanthrene", "pyrene"}),
                         FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"xanthene", "steroid"}), FG.KETONE, frozenset({"in_skeleton"})),  #
    RingExpressionPolicy(frozenset({"xanthene", "steroid"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_CYCLIC_CORE, FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_CYCLIC_CORE, FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_CYCLIC_CORE, FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_CYCLIC_CORE, FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(_CYCLIC_CORE, FG.NITRILE, frozenset({"exocyclic"})),
)


def supports_ring_expression(scaffold: ScaffoldIdentity, facts) -> bool:
    """判断 scaffold 是否支持该主基团的 typed 环表达。"""
    return any(scaffold.naming_class in policy.naming_classes
               and facts.group_class is policy.group_class
               and facts.relation.value in policy.relations for policy in _POLICIES)
