"""由旧条目列表构建的带类型 Layer 1 官能团清单。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from namepredict.layer1.fg_registry import FG_SPECS


class FunctionalGroupClass(str, Enum):
    """官能团类别枚举（对应 fg_registry 的 p41 优先级体系）。"""
    RADICAL = "radical"
    ACYL = "acyl"
    AZANIDE = "azanide"  # 氮负离子母体（P-72.2.2.2）：N⁻ 自任阴离子母体，余者作前缀
    CATION = "cation"
    ACID = "acid"
    OXOACID = "oxoacid"
    SULFONAMIDE = "sulfonamide"
    ESTER = "ester"
    ACYL_HALIDE = "acyl_halide"
    AMIDE = "amide"
    NITRILE = "nitrile"
    ALDEHYDE = "aldehyde"
    KETONE = "ketone"
    THIONE = "thione"  # 酮的硫族类似物 C=S（P-64.6.1）：同类 16，让位于 C=O
    ALCOHOL = "alcohol"
    THIOL = "thiol"
    AMINE = "amine"
    HETERANE = "heterane"  # 非碳母体氢化物（P-21 类 21–39）：杂原子自任母体
    NONE = 'alkane'


OXO_FG_CLASSES = frozenset({FunctionalGroupClass.OXOACID, FunctionalGroupClass.SULFONAMIDE})


@dataclass(frozen=True)
class FunctionalGroupOccurrence:
    """单个官能团出现：id、类别、特征原子、母体锚点、附加负载与 P-41 仲裁状态。"""
    id: str
    group_class: FunctionalGroupClass
    characteristic_atoms: frozenset[int]
    parent_anchors: frozenset[int]
    payload: dict
    demoted: bool = False  # 被 P-41 仲裁降级为前缀叶（P-61.1.3）


@dataclass(frozen=True)
class FunctionalGroupInventory:
    """官能团清单：承载全部出现并提供按类查询/统计。"""
    entries: tuple[FunctionalGroupOccurrence, ...]
    has_anion: bool = False  # 分子内存在负形式电荷原子（P-41 表 4.1 类 4）：阴离子 > 阳离子

    def occurrences(self, group_class: FunctionalGroupClass) -> tuple[FunctionalGroupOccurrence, ...]:
        """返回给定官能团类的全部出现（不含被 P-41 仲裁降级的叶条目）。"""
        return tuple(e for e in self.entries if e.group_class == group_class and not e.demoted)

    def demoted_entries(self) -> tuple[FunctionalGroupOccurrence, ...]:
        """返回被 P-41 仲裁降级为前缀叶的全部条目。"""
        return tuple(e for e in self.entries if e.demoted)

_FG_KEYS = tuple(sp.fg for sp in FG_SPECS)  # 类别注册来源见 fg_registry.FG_SPECS

_ANCHOR_KEYS = {FunctionalGroupClass(sp.fg): sp.anchors for sp in FG_SPECS if sp.anchors}  # 只有声明 anchors 的 FG 才收集锚点

from namepredict.constants import C, N, O


def _idx(payload: dict, key: str) -> int | None:
    """取 payload 中单个索引键（缺失返回 None）。"""
    v = payload.get(key)
    return int(v) if v is not None else None

def _oxoacid_atoms(mol, payload: dict) -> set[int]:
    """含氧酸：锚点碳 + 中心 P/S + 中心的非碳邻居（碳臂留给链/取代基侧）。"""
    z = _idx(payload, "oxo_z")
    if z is None:
        return set()
    out = {z} | {n.GetIdx() for n in mol.GetAtomWithIdx(z).GetNeighbors()
                 if n.GetAtomicNum() not in (1, C)}
    anchor = _idx(payload, "center_idx")
    return out | ({anchor} if anchor is not None else set())


def center_surr_atoms(payload: dict) -> frozenset[int]:
    """通用 FG 特征原子：中心原子与全部周边原子之并。"""
    out = {int(i) for i in (payload.get("surr_idx") or ())}
    center = _idx(payload, "center_idx")
    if center is not None:
        out.add(center)
    return frozenset(out)


def _cation_atoms(_mol, payload: dict) -> set[int]:
    """阳离子/杂原子烃：特征原子只有中心本身，周边原子全留给取代基侧。"""
    center = _idx(payload, "center_idx")
    return {center} if center is not None else set()


def _nitrile_atoms(mol, payload: dict) -> set[int]:
    """腈：特征原子只含腈碳与三键氮（P-66.1.5）。"""
    from rdkit.Chem import BondType

    center = _idx(payload, "center_idx")
    if center is None or mol is None:
        return set()
    out = {center}
    for nb in mol.GetAtomWithIdx(center).GetNeighbors():
        if nb.GetAtomicNum() != N:
            continue
        if mol.GetBondBetweenAtoms(center, nb.GetIdx()).GetBondType() == BondType.TRIPLE:
            out.add(nb.GetIdx())
    return out


FG_ATOM_FNS = {  # 不走通用规则的例外类别
    "azanide": _cation_atoms,  # 氮负离子母体只占 N⁻ 本身（P-72.2.2.2）：全部臂退为前缀
    "oxoacid": _oxoacid_atoms,
    "sulfonamide": _oxoacid_atoms,  # 含氧酸合一类的 P-41 酰胺分组，特征原子同一判据
    "cation": _cation_atoms,  # 单核母体阳离子作母体时只占一个原子（P-73.1.1）
    "heterane": _cation_atoms,  # 杂原子烃母体同理只占杂原子（P-21）：周边臂一律退为取代基
    "nitrile": _nitrile_atoms,  # 腈的 R 侧连接原子不并入所有权（P-66.1.5）
}

def _indices(payload: dict, keys: tuple[str, ...]) -> frozenset[int]:
    """从 payload 按键收集全部整数值作为索引集合。"""
    values = (payload.get(key) for key in keys)
    flat = [x for value in values for x in (value if isinstance(value, (list, tuple, set, frozenset)) else [value])]
    return frozenset(x for x in flat if isinstance(x, int))


def _characteristic_atoms(group_class: FunctionalGroupClass, mol, payload: dict) -> frozenset[int]:
    """按 FG 类别取特征原子集（未登记时退回通用并集）。"""
    fn = FG_ATOM_FNS.get(group_class.value)
    if fn is not None and mol is not None:
        return frozenset(fn(mol, payload))
    return center_surr_atoms(payload)


def _one(key: str, index: int, payload: dict, mol, demoted: bool = False) -> FunctionalGroupOccurrence:
    """将单条官能团 dict 组装为带类型的出现。"""
    group_class = FunctionalGroupClass(key)
    anchors = _indices(payload, _ANCHOR_KEYS.get(group_class, ()))
    return FunctionalGroupOccurrence(f"{key}:{index}", group_class,
                                     _characteristic_atoms(group_class, mol, payload), anchors, payload, demoted)


def build_inventory(lists: dict, mol=None, demoted: frozenset[str] = frozenset(),
                    has_anion: bool | None = None) -> FunctionalGroupInventory:
    """由官能团列表构建带类型的清单；demoted 为降级 id 集。"""
    entries = tuple(_one(key, i, item, mol, f"{key}:{i}" in demoted)
                    for key in _FG_KEYS for i, item in enumerate(lists.get(key) or ()))
    if has_anion is None:
        has_anion = mol is not None and any(a.GetFormalCharge() < 0 for a in mol.GetAtoms())
    return FunctionalGroupInventory(entries, has_anion)

def inventory_from_info(info: dict) -> FunctionalGroupInventory:
    """从分析信息中取出清单；缺失即显式失败。"""
    inventory = info.get("fg_inventory")
    if not isinstance(inventory, FunctionalGroupInventory):
        raise KeyError("info 缺少 fg_inventory：L1 出口只以 FunctionalGroupInventory 承载官能团事实")
    return inventory
