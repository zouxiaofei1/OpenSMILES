"""环骨架上类型化主基团表达式的能力策略。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.constants import HW_CLASS
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
_HW = frozenset({HW_CLASS})  # 生成式 Hantzsch-Widman 杂单环（P-22.2.2）

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
    # 生成式 HW 杂单环（P-22.2.2）：环内嵌 FG 同环系内核，环外 FG 走后缀
    RingExpressionPolicy(_HW, FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_HW, FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_HW, FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_HW, FG.AMIDE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_HW, FG.ESTER, frozenset({"in_skeleton"})),
    RingExpressionPolicy(_HW, FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(_HW, FG.NITRILE, frozenset({"exocyclic"})),
    RingExpressionPolicy(_HW, FG.ALDEHYDE, frozenset({"exocyclic"})),
    RingExpressionPolicy(_HW, FG.SULFONAMIDE, frozenset({"exocyclic"})),
)


_POLICY_INDEX = frozenset(
    (nc, p.group_class, rel) for p in _POLICIES for nc in p.naming_classes for rel in p.relations
)  # 展平为 O(1) 查询集，避免每次全表扫描


def supports_ring_expression(scaffold: ScaffoldIdentity, facts) -> bool:
    """判断 scaffold 是否支持该主基团的 typed 环表达。"""
    return (scaffold.naming_class, facts.group_class, facts.relation.value) in _POLICY_INDEX
