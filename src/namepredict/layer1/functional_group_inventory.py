"""由旧条目列表构建的带类型 Layer 1 官能团清单。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from namepredict.layer1.fg_registry import FG_SPECS


class FunctionalGroupClass(str, Enum):
    RADICAL = "radical"
    ACYL = "acyl"
    ACID = "acid"
    PHOSPHATE = "phosphate"
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


# FG 类别注册唯一事实来源在 fg_registry.FG_SPECS；此处派生，不再逐条手写。
_LIST_CLASSES = {sp.list_key: FunctionalGroupClass(sp.fg) for sp in FG_SPECS}

# 锚点 key（occurrence payload）：只有声明 anchors 的 FG 才收集；胺含多臂锚点（P-62.2）。
_ANCHOR_KEYS = {FunctionalGroupClass(sp.fg): sp.anchors for sp in FG_SPECS if sp.anchors}


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
