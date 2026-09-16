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
    FgSpec("oxoacid", p41=9, path=(0,), anchors=("center_idx",)),  # 含氧酸中心 P/S 合一类：oxo_kind 由 L1 payload 归一
    FgSpec("sulfonamide", p41=11, path=(1,), anchors=("center_idx",)),  # P-41 类 11：磺酰胺与酰胺同组，排在酰胺之后
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


# ── P-41 表 4.1 类 7：含氧酸的"酸式"判据 ──
# 中心带 O⁻（阴离子）或中心为碳锚定酸式（膦酸/磺酸）的含氧酸按类 7 参与主基团竞争，
# 排在类 9 酯与类 11 酰胺之前；中性磷酸/硫酸酯仍属类 9（按相应酸排羧酸酯之后）。
OXO_ACID_P41 = 8  # 类 7 内排在羧酸（acid=7）之后的等级
OXO_ACID_KINDS = frozenset({"sulfonic"})  # 碳锚定且自带酸式氢的 oxo_kind；phosphonate 另按 n_oh 判
OXO_ACID_KIND_BY_H = {"phosphonate": 1}  # 碳锚定 P 酸：n_oh ≥ 该值时按酸式（膦酸/膦酸氢酯）


def oxoacid_is_acid(payload: dict) -> bool:
    """含氧酸 occurrence 是否按 P-41 类 7 的酸参与主基团竞争。"""
    if int(payload.get("n_om") or 0) > 0:  # 中心带 O⁻：按酸根处理
        return True
    kind = payload.get("oxo_kind")
    if kind in OXO_ACID_KINDS:
        return True
    need = OXO_ACID_KIND_BY_H.get(kind)
    return need is not None and int(payload.get("n_oh") or 0) >= need
