"""FG 跨层注册元数据（L1 key、L2 优先级、L4/L5）唯一事实来源。
下游表由 ``FG_SPECS`` 派生，避免重复登记。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FgSpec:
    fg: str  # FunctionalGroupClass 枚举值；同作 L1 列表 key
    p41: int = 0  # P-41 主官能团等级（0 = 非主官能团）
    path: tuple[int, ...] = ()  # P-43 优先级路径
    expr: str = "suffix"  # 表达类型：suffix/prefix_only/legacy_compat
    anchors: tuple[str, ...] = ()  # occurrence payload 锚点 key（空 = 不收集锚点）
    parent_anchor_fields: tuple[str, str] | None = None  # parent 锚点字段（单, 复）
    locant_source: str = "attachment"  # 位次原子来源（attachment/anchor_field 等）


FG_SPECS: tuple[FgSpec, ...] = (  # 全部 L1 检测列表对应的 FG 类别
    FgSpec("radical", p41=1, anchors=("center_idx",), parent_anchor_fields=("radical_c_idx")),
    FgSpec("acyl", p41=1, anchors=("center_idx",), parent_anchor_fields=("acyl_c_idx")),
    FgSpec("azanide", p41=4, anchors=("center_idx",)),
    FgSpec("cation", p41=6, anchors=("center_idx",)),
    FgSpec("acid", p41=7, path=(1,), anchors=("center_idx",)),
    FgSpec("oxoacid", p41=9, path=(0,), anchors=("center_idx",)),  # 含氧酸中心 P/S 合一类：oxo_kind 由 L1 payload 归一
    FgSpec("sulfonamide", p41=11, path=(1,), anchors=("center_idx",)),  # P-41 类 11：磺酰胺与酰胺同组，排在酰胺之后
    FgSpec("ester", p41=9, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("acyl_halide", p41=10, anchors=("center_idx",)),
    FgSpec("amide", p41=11, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("nitrile", p41=14, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("aldehyde", p41=15, anchors=("center_idx",), locant_source="attachment_exocyclic"),
    FgSpec("ketone", p41=16, anchors=("center_idx",)),
    FgSpec("thione", p41=16, path=(1,), anchors=("center_idx",)),
    FgSpec("alcohol", p41=17, path=(1,), anchors=("surr_idx",)),
    FgSpec("thiol", p41=17, path=(2,), anchors=("surr_idx",)),
    FgSpec("amine", p41=19, anchors=("surr_idx",)),
    FgSpec("heterane", p41=36, anchors=("center_idx",)),
)


OXO_ACID_P41 = 8  # 类 7 内排在羧酸（acid=7）之后的等级
OXO_ACID_KINDS = frozenset({"sulfonic"})  # 碳锚定且带酸式氢的 oxo_kind
OXO_ACID_KIND_BY_H = {"phosphonate": 1, "boronic": 1}  # 碳锚定 P/B 酸：n_oh ≥ 该值时按酸式（膦酸/膦酸氢酯/硼酸）


def oxoacid_is_acid(payload: dict) -> bool:
    """含氧酸 occurrence 是否按 P-41 类 7 的酸参与主基团竞争。"""
    if int(payload.get("n_om") or 0) > 0:  # 中心带 O⁻：按酸根处理
        return True
    kind = payload.get("oxo_kind")
    if kind in OXO_ACID_KINDS:
        return True
    need = OXO_ACID_KIND_BY_H.get(kind)
    return need is not None and int(payload.get("n_oh") or 0) >= need
