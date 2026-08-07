"""Typed Layer 1 functional-group inventory built from legacy entry lists."""
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
    HYDRAZINE = "hydrazine"
    AMINE = "amine"
    QUATERNARY_AMMONIUM = "quaternary_ammonium"
    BORONIC = "boronic"
    CARBAMATE = "carbamate"
    CARBONATE = "carbonate"
    GUANIDINE = "guanidine"
    UREA = "urea"
    ISOCYANATE = "isocyanate"
    ISOTHIOCYANATE = "isothiocyanate"
    PHOSPHATE = "phosphate"
    PHOSPHONIC = "phosphonic"
    ETHER = "ether"
    SULFIDE = "sulfide"
    SULFOXIDE = "sulfoxide"
    SULFONE = "sulfone"
    SULFONAMIDE = "sulfonamide"
    SULFONATE = "sulfonate"
    SULFONYL_HALIDE = "sulfonyl_halide"
    SULFONIC_ACID = "sulfonic_acid"


@dataclass(frozen=True)
class FunctionalGroupOccurrence:
    id: str
    group_class: FunctionalGroupClass
    characteristic_atoms: frozenset[int]
    parent_anchors: frozenset[int]
    payload: dict

    @property
    def atoms(self) -> frozenset[int]:
        return self.characteristic_atoms


@dataclass(frozen=True)
class FunctionalGroupInventory:
    entries: tuple[FunctionalGroupOccurrence, ...]

    def occurrences(self, group_class: FunctionalGroupClass) -> tuple[FunctionalGroupOccurrence, ...]:
        return tuple(e for e in self.entries if e.group_class == group_class)

    def has(self, group_class: FunctionalGroupClass) -> bool:
        return bool(self.occurrences(group_class))

    def count(self, group_class: FunctionalGroupClass) -> int:
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
    "hydrazines": FunctionalGroupClass.HYDRAZINE,
    "amines": FunctionalGroupClass.AMINE,
    "quaternary_ammoniums": FunctionalGroupClass.QUATERNARY_AMMONIUM,
    "boronics": FunctionalGroupClass.BORONIC,
    "carbamates": FunctionalGroupClass.CARBAMATE,
    "carbonates": FunctionalGroupClass.CARBONATE,
    "guanidines": FunctionalGroupClass.GUANIDINE,
    "ureas": FunctionalGroupClass.UREA,
    "isocyanates": FunctionalGroupClass.ISOCYANATE,
    "isothiocyanates": FunctionalGroupClass.ISOTHIOCYANATE,
    "phosphates": FunctionalGroupClass.PHOSPHATE,
    "phosphonics": FunctionalGroupClass.PHOSPHONIC,
    "ethers": FunctionalGroupClass.ETHER,
    "sulfides": FunctionalGroupClass.SULFIDE,
    "sulfoxides": FunctionalGroupClass.SULFOXIDE,
    "sulfones": FunctionalGroupClass.SULFONE,
    "sulfonamides": FunctionalGroupClass.SULFONAMIDE,
    "sulfonates": FunctionalGroupClass.SULFONATE,
    "sulfonyl_chlorides": FunctionalGroupClass.SULFONYL_HALIDE,
    "sulfonic_acids": FunctionalGroupClass.SULFONIC_ACID,
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
    FunctionalGroupClass.AMINE: ("c_idx",),
}


def _indices(payload: dict, keys: tuple[str, ...]) -> frozenset[int]:
    values = (payload.get(key) for key in keys)
    flat = [x for value in values for x in (value if isinstance(value, (list, tuple, set, frozenset)) else [value])]
    return frozenset(x for x in flat if isinstance(x, int))


def _atom_ids(payload: dict) -> frozenset[int]:
    keys = tuple(k for k in payload if k.endswith("idx") or k.endswith("idxs"))
    return _indices(payload, keys)


def _one(key: str, index: int, payload: dict) -> FunctionalGroupOccurrence:
    group_class = _LIST_CLASSES[key]
    anchors = _indices(payload, _ANCHOR_KEYS.get(group_class, ()))
    return FunctionalGroupOccurrence(f"{key}:{index}", group_class, _atom_ids(payload), anchors, payload)


def build_inventory(lists: dict) -> FunctionalGroupInventory:
    entries = tuple(_one(key, i, item) for key in _LIST_CLASSES for i, item in enumerate(lists.get(key) or ()))
    return FunctionalGroupInventory(entries)


def inventory_from_info(info: dict) -> FunctionalGroupInventory:
    inventory = info.get("fg_inventory")
    return inventory if isinstance(inventory, FunctionalGroupInventory) else build_inventory(info)
