"""FG 跨层注册元数据（L1 key、L2 优先级、L4/L5）唯一事实来源。
下游表由 ``FG_SPECS`` 派生，避免重复登记。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FgSpec:
    """一个官能团类别在命名管线各层的注册元数据（字段纯数据，不依赖下游类型）。"""

    fg: str  # FunctionalGroupClass 值（权威枚举字符串，如 "alcohol"）
    list_key: str  # L1 analyzer 列表复数 key（如 "hydroxyls"）
    p41: int = 0  # P-41 主官能团等级（0 = 非主官能团）
    path: tuple[int, ...] = ()  # P-43 优先级路径
    compat: int = 0  # 兼容等级（legacy_rank 消费）
    expr: str = "suffix"  # 表达类型：suffix / prefix_only / legacy_compat
    anchors: tuple[str, ...] = ()  # occurrence payload 锚点 key（空 = 不收集锚点）
    parent_anchor_fields: tuple[str, str] | None = None  # parent 锚点字段（单, 复）
    chain: bool = False  # 可作链式主官能团（_CHAIN_FG）
    multi: bool = False  # 支持数量后缀（_MULTI_FG / L5 mult_ok）
    rs: bool = False  # 支持 R/S（stereo._RS_KINDS）
    keep_locant: bool = False  # 取代基位次保留（assembler_prefixes._KEEP_LOCANT_KINDS）
    locant_kind: str | None = None  # fg_locants 记录 kind（"oh"/"amine"/…）；None = 不产位次记录
    locant_source: str = "attachment"  # 位次原子来源：attachment 取 principal_expression_facts 骨架内附着原子；attachment_exocyclic 仅环外表达时取；anchor_field 取 parent_anchor_fields 首位语义字段（固定 locant 1 锚点）
    oh_parent: bool = False  # 醇类母体（抑制羟基取代基提取）
    nh2_parent: bool = False  # 胺类母体
    oxo_parent: bool = False  # 酮类母体


FG_SPECS: tuple[FgSpec, ...] = (  # 全部 L1 检测列表对应的 FG 类别
    FgSpec("radical", "radicals", p41=1,  anchors=("center_idx",),
           parent_anchor_fields=("radical_c_idx", "radical_c_idxs"),
           locant_kind="radical", locant_source="anchor_field"),
    FgSpec("acyl", "acyls", p41=1,  anchors=("center_idx",),  # 酰基残基
           parent_anchor_fields=("acyl_c_idx", "acyl_c_idxs"),
           chain=True, rs=True, keep_locant=True, locant_kind="acyl"),
    FgSpec("acid", "carboxyls", p41=7, path=(1,),  anchors=("center_idx",),
           chain=True, multi=True, rs=True, keep_locant=True, locant_kind="acid"),
    FgSpec("phosphate", "phosphates", p41=9, path=(1,),  anchors=("p_idx",),  chain=True),
    FgSpec("ester", "esters", p41=9,  anchors=("center_idx",),
           chain=True, multi=True,rs=True, locant_kind="ester", locant_source="attachment_exocyclic"),
    FgSpec("acyl_halide", "acyl_chlorides", p41=10,  anchors=("center_idx",), chain=True),
    FgSpec("amide", "amides", p41=11, anchors=("center_idx",),
           chain=True, multi=True,rs=True, locant_kind="amide", locant_source="attachment_exocyclic"),
    FgSpec("nitrile", "nitriles", p41=14,anchors=("center_idx",),
           chain=True, rs=True, locant_kind="nitrile", locant_source="attachment_exocyclic"),
    FgSpec("aldehyde", "aldehydes", p41=15,  anchors=("center_idx",),
           chain=True, rs=True, locant_kind="aldehyde", locant_source="attachment_exocyclic"),
    FgSpec("ketone", "ketones", p41=16,  anchors=("center_idx",),
           chain=True, multi=True, rs=True, keep_locant=True, locant_kind="ketone",
           oxo_parent=True),
    FgSpec("alcohol", "hydroxyls", p41=17, path=(1,), anchors=("surr_idx",),  
           chain=True, multi=True, rs=True, locant_kind="oh", oh_parent=True),
    FgSpec("thiol", "thiols", p41=17, path=(2,),  anchors=("surr_idx",), 
           chain=True, multi=True, rs=True, locant_kind="sh"),
    FgSpec("amine", "amines", p41=19,anchors=("surr_idx",),  
           chain=True, multi=True, rs=True, locant_kind="amine", nh2_parent=True),
)

def chain_fgs() -> frozenset[str]:
    """可作链式主官能团的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.chain)


def multi_fgs() -> frozenset[str]:
    """支持数量后缀的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.multi)


def srs_fgs() -> frozenset[str]:
    """支持 R/S 的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.rs)


def keep_locant_fgs() -> frozenset[str]:
    """取代基位次保留的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.keep_locant)
