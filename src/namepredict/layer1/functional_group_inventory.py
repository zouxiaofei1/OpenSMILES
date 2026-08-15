"""由旧条目列表构建的带类型 Layer 1 官能团清单。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class FunctionalGroupClass(str, Enum):
    RADICAL = "radical"
    ACID = "acid"
    ANHYDRIDE = "anhydride"
    ESTER = "ester"
    ACYL_HALIDE = "acyl_halide"
    AMIDE = "amide"
    NITRILE = "nitrile"
    ALDEHYDE = "aldehyde"
    KETONE = "ketone"
    ALCOHOL = "alcohol"
    THIOL = "thiol"
    AMINE = "amine"
    QUATERNARY_AMMONIUM = "quaternary_ammonium"
    ISOCYANATE = "isocyanate"
    ISOTHIOCYANATE = "isothiocyanate"
    ETHER = "ether"
    SULFIDE = "sulfide"
    NONE = 'alkane'


@dataclass(frozen=True)
class FunctionalGroupOccurrence:
    id: str
    group_class: FunctionalGroupClass
    characteristic_atoms: frozenset[int]
    parent_anchors: frozenset[int]
    payload: dict

    @property
    def atoms(self) -> frozenset[int]:
        """返回特征原子集合。"""
        return self.characteristic_atoms


@dataclass(frozen=True)
class FunctionalGroupInventory:
    entries: tuple[FunctionalGroupOccurrence, ...]

    def occurrences(self, group_class: FunctionalGroupClass) -> tuple[FunctionalGroupOccurrence, ...]:
        """返回给定官能团类的全部出现。"""
        return tuple(e for e in self.entries if e.group_class == group_class)

    def has(self, group_class: FunctionalGroupClass) -> bool:
        """判断是否存在给定官能团类。"""
        return bool(self.occurrences(group_class))

    def count(self, group_class: FunctionalGroupClass) -> int:
        """统计给定官能团类的出现次数。"""
        return len(self.occurrences(group_class))


_LIST_CLASSES = {
    "radicals": FunctionalGroupClass.RADICAL,
    "carboxyls": FunctionalGroupClass.ACID,
    "anhydrides": FunctionalGroupClass.ANHYDRIDE,
    "esters": FunctionalGroupClass.ESTER,
    "acyl_chlorides": FunctionalGroupClass.ACYL_HALIDE,
    "amides": FunctionalGroupClass.AMIDE,
    "nitriles": FunctionalGroupClass.NITRILE,
    "aldehydes": FunctionalGroupClass.ALDEHYDE,
    "ketones": FunctionalGroupClass.KETONE,
    "hydroxyls": FunctionalGroupClass.ALCOHOL,
    "thiols": FunctionalGroupClass.THIOL,
    "amines": FunctionalGroupClass.AMINE,
    "quaternary_ammoniums": FunctionalGroupClass.QUATERNARY_AMMONIUM,
    "isocyanates": FunctionalGroupClass.ISOCYANATE,
    "isothiocyanates": FunctionalGroupClass.ISOTHIOCYANATE,
    "ethers": FunctionalGroupClass.ETHER,
    "sulfides": FunctionalGroupClass.SULFIDE,
}


_ANCHOR_KEYS = {
    FunctionalGroupClass.RADICAL: ("c_idx",),
    FunctionalGroupClass.ACID: ("c_idx",),
    FunctionalGroupClass.ESTER: ("c_idx",),
    FunctionalGroupClass.ACYL_HALIDE: ("c_idx",),
    FunctionalGroupClass.AMIDE: ("c_idx",),
    FunctionalGroupClass.NITRILE: ("c_idx",),
    FunctionalGroupClass.ALDEHYDE: ("c_idx",),
    FunctionalGroupClass.KETONE: ("c_idx",),
    FunctionalGroupClass.ALCOHOL: ("c_idx",),
    FunctionalGroupClass.THIOL: ("c_idx",),
    # 胺：2°/3° 的 N 上所有碳邻居都是候选母体臂锚点（P-62.2 多臂选优）。
    FunctionalGroupClass.AMINE: ("c_idx", "c_idxs"),
}


def _indices(payload: dict, keys: tuple[str, ...]) -> frozenset[int]:
    """从 payload 按键收集全部整数值作为索引集合。"""
    values = (payload.get(key) for key in keys)
    flat = [x for value in values for x in (value if isinstance(value, (list, tuple, set, frozenset)) else [value])]
    return frozenset(x for x in flat if isinstance(x, int))


def _atom_ids(payload: dict) -> frozenset[int]:
    """收集 payload 中所有 idx/idxs 键的原子索引集合。"""
    keys = tuple(k for k in payload if k.endswith("idx") or k.endswith("idxs"))
    return _indices(payload, keys)


def _one(key: str, index: int, payload: dict) -> FunctionalGroupOccurrence:
    """将单条官能团 dict 组装为带类型的出现。"""
    group_class = _LIST_CLASSES[key]
    anchors = _indices(payload, _ANCHOR_KEYS.get(group_class, ()))
    return FunctionalGroupOccurrence(f"{key}:{index}", group_class, _atom_ids(payload), anchors, payload)


def build_inventory(lists: dict) -> FunctionalGroupInventory:
    """由官能团列表构建带类型的 FunctionalGroupInventory。"""
    entries = tuple(_one(key, i, item) for key in _LIST_CLASSES for i, item in enumerate(lists.get(key) or ()))
    return FunctionalGroupInventory(entries)


def inventory_from_info(info: dict) -> FunctionalGroupInventory:
    """从分析信息中取出清单，缺失时回退构建。"""
    inventory = info.get("fg_inventory")
    return inventory if isinstance(inventory, FunctionalGroupInventory) else build_inventory(info)
