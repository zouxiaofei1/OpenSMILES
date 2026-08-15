"""FG 主官能团元数据单一事实来源。

集中每个官能团类别的跨层注册元数据（L1 检测列表 key、L2 优先级与锚点字段、
L4 位次记录 kind、L5 命名/立体/位次集合），下游表由 ``FG_SPECS`` 派生，
避免同一 FG 在多个文件各登记一遍（多重注册收敛点）。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FgSpec:
    """一个官能团类别在命名管线各层的注册元数据。

    字段全部为纯数据，不依赖下游类型，consumer 自行转换成所需形态。
    """

    fg: str  # FunctionalGroupClass 值（权威枚举字符串，如 "alcohol"）
    list_key: str  # L1 analyzer 列表复数 key（如 "hydroxyls"）
    p41: int = 0  # P-41 主官能团等级（0 = 非主官能团）
    path: tuple[int, ...] = ()  # P-43 优先级路径
    expr: str = "suffix"  # 表达类型：suffix / prefix_only / legacy_compat
    compat: int = 0  # 兼容等级（legacy_rank 消费）
    anchors: tuple[str, ...] = ()  # occurrence payload 锚点 key（空 = 不收集锚点）
    parent_anchor_fields: tuple[str, str] | None = None  # parent 锚点字段（单, 复）
    chain: bool = False  # 可作链式主官能团（_CHAIN_FG）
    multi: bool = False  # 支持数量后缀（_MULTI_FG / L5 mult_ok）
    rs: bool = False  # 支持 R/S（stereo._RS_KINDS）
    keep_locant: bool = False  # 取代基位次保留（assembler_prefixes._KEEP_LOCANT_KINDS）
    locant_kind: str | None = None  # fg_locants 记录 kind（"oh"/"amine"/…）
    oh_parent: bool = False  # 醇类母体（抑制羟基取代基提取）
    nh2_parent: bool = False  # 胺类母体
    oxo_parent: bool = False  # 酮类母体


_LEGACY = "legacy_compat"
_PREFIX = "prefix_only"

# 全部 L1 检测列表对应的 FG 类别（含非主官能团：醚/硫醚/季铵/异氰酸酯/硝基…）。
# 新增一个官能团类别：在本表加一条 + L1 analyzer 检测 + L5 chain_engine._KIND_TABLE 命名规格。
FG_SPECS: tuple[FgSpec, ...] = (
    FgSpec("radical", "radicals", p41=1, compat=1, anchors=("c_idx",),
           parent_anchor_fields=("radical_c_idx", "radical_c_idxs"),
           locant_kind="radical"),
    FgSpec("acid", "carboxyls", p41=7, path=(1,), compat=14, anchors=("c_idx",),
           parent_anchor_fields=("cooh_c_idx", "cooh_c_idxs"),
           chain=True, multi=True, rs=True, keep_locant=True, locant_kind="acid"),
    FgSpec("anhydride", "anhydrides", p41=8, compat=12),
    FgSpec("ester", "esters", p41=9, compat=11, anchors=("c_idx",),
           parent_anchor_fields=("ester_c_idx", "ester_c_idxs"),
           chain=True, rs=True),
    FgSpec("acyl_halide", "acyl_chlorides", p41=10, compat=10, anchors=("c_idx",)),
    FgSpec("amide", "amides", p41=11, compat=9, anchors=("c_idx",),
           parent_anchor_fields=("amide_c_idx", "amide_c_idxs"),
           chain=True, rs=True, locant_kind="amide"),
    FgSpec("nitrile", "nitriles", p41=14, compat=8, anchors=("c_idx",),
           parent_anchor_fields=("nitrile_c_idx", "nitrile_c_idxs"),
           chain=True, rs=True),
    FgSpec("aldehyde", "aldehydes", p41=15, compat=7, anchors=("c_idx",),
           parent_anchor_fields=("aldehyde_c_idx", "aldehyde_c_idxs"),
           chain=True, rs=True),
    FgSpec("ketone", "ketones", p41=16, compat=6, anchors=("c_idx",),
           parent_anchor_fields=("ketone_c_idx", "ketone_c_idxs"),
           chain=True, multi=True, rs=True, keep_locant=True, locant_kind="ketone",
           oxo_parent=True),
    FgSpec("alcohol", "hydroxyls", p41=17, path=(1,), compat=5, anchors=("c_idx",),
           parent_anchor_fields=("oh_c_idx", "oh_c_idxs"),
           chain=True, multi=True, rs=True, locant_kind="oh", oh_parent=True),
    FgSpec("thiol", "thiols", p41=17, path=(2,), compat=4, anchors=("c_idx",),
           parent_anchor_fields=("sh_c_idx", "sh_c_idxs"),
           chain=True, multi=True, rs=True, locant_kind="sh"),
    FgSpec("amine", "amines", p41=19, compat=3, anchors=("c_idx", "c_idxs"),
           parent_anchor_fields=("amine_c_idx", "amine_c_idxs"),
           chain=True, multi=True, rs=True, locant_kind="amine", nh2_parent=True),
    FgSpec("quaternary_ammonium", "quaternary_ammoniums"),
    FgSpec("isocyanate", "isocyanates", p41=41, expr=_LEGACY, compat=8),
    FgSpec("isothiocyanate", "isothiocyanates", p41=41, expr=_LEGACY, compat=8),
    FgSpec("sulfide", "sulfides", p41=41, path=(2,), expr=_LEGACY, compat=2),
    FgSpec("ether", "ethers", p41=41, path=(1,), expr=_PREFIX),
)

def chain_fgs() -> frozenset[str]:
    """可作链式主官能团的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.chain)


def multi_fgs() -> frozenset[str]:
    """支持数量后缀的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.multi)


def rs_fgs() -> frozenset[str]:
    """支持 R/S 的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.rs)


def keep_locant_fgs() -> frozenset[str]:
    """取代基位次保留的 FG 值集合。"""
    return frozenset(sp.fg for sp in FG_SPECS if sp.keep_locant)
