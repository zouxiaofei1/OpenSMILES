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


_POLICIES = (
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ACID, frozenset({"exocyclic"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.AMINE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.KETONE, frozenset({"in_skeleton"})),  # 保留稠环母体的环内酮（色烯-4-酮/喹啉酮型）：缺此条则环酮被 _unsupported_typed_ring 整体拦截，候选数为 0、直接空输出（chromone/coumaranone 等）。
    RingExpressionPolicy(frozenset({"monohetero", "fused56", "purine", "carbazole",  # 单杂环/稠杂环（吡啶酮、吡唑酮、乙内酰脲、嘌呤二酮、吖啶酮型）与纯碳稠环的环内酮：L1 已把环内 C=O 归酮，缺此条则这些命名类的环酮同样被 _unsupported_typed_ring 拦截成空输出（单杂环不在 naph_family，苯并噁唑/嘌呤等是各自命名类）。
                                    "acridine", "phenothiazine", "benzodioxole",
                                    "anthra", "phenanthrene", "pyrene"}),
                         FG.KETONE, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"xanthene", "steroid"}), FG.KETONE, frozenset({"in_skeleton"})),  # 传统编号保留母体（呫吨/甾体）：缺此条则环内酮、环内醇被 _unsupported_typed_ring 拦截。
    RingExpressionPolicy(frozenset({"xanthene", "steroid"}), FG.ALCOHOL, frozenset({"in_skeleton"})),
    RingExpressionPolicy(frozenset({"fused_hetero", "fused"}), FG.ALCOHOL, frozenset({"in_skeleton"})),  # 未注册稠环(fused_hetero/fused): kind 已正交化, 主 FG 走 L5 fused_tree 组装, 需放开 typed 表达(否则环酮被 _unsupported_typed_ring 拦截)。
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
