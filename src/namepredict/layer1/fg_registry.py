"""FG 跨层注册元数据（L1 key、L2 优先级、L4/L5）唯一事实来源。
下游表由 ``FG_SPECS`` 派生，避免重复登记。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FgSpec:
    fg: str  # FunctionalGroupClass 值（权威枚举字符串，如 "alcohol"）；同作 L1 analyzer 列表 key
    p41: int = 0  # P-41 主官能团等级（0 = 非主官能团）
    path: tuple[int, ...] = ()  # P-43 优先级路径
    expr: str = "suffix"  # 表达类型：suffix / prefix_only / legacy_compat
    anchors: tuple[str, ...] = ()  # occurrence payload 锚点 key（空 = 不收集锚点）
    parent_anchor_fields: tuple[str, str] | None = None  # parent 锚点字段（单, 复）
    locant_source: str = "attachment"  # 位次原子来源：attachment 取 principal_expression_facts 骨架内附着原子；attachment_exocyclic 仅环外表达时取；anchor_field 取 parent_anchor_fields 首位语义字段（固定 locant 1 锚点）


FG_SPECS: tuple[FgSpec, ...] = (  # 全部 L1 检测列表对应的 FG 类别
    FgSpec("radical", p41=1, anchors=("center_idx",), parent_anchor_fields=("radical_c_idx")),
    FgSpec("acyl", p41=1, anchors=("center_idx",), parent_anchor_fields=("acyl_c_idx")),
    FgSpec("acid", p41=7, path=(1,), anchors=("center_idx",)),
    FgSpec("phosphate", p41=9, path=(1,), anchors=("p_idx",)),
    FgSpec("ester", p41=9, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("acyl_halide", p41=10, anchors=("center_idx",)),
    FgSpec("amide", p41=11, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("nitrile", p41=14, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("aldehyde", p41=15, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("ketone", p41=16, anchors=("center_idx",)),
    FgSpec("alcohol", p41=17, path=(1,), anchors=("surr_idx",)),
    FgSpec("thiol", p41=17, path=(2,), anchors=("surr_idx",)),
    FgSpec("amine", p41=19, anchors=("surr_idx",)),
)

