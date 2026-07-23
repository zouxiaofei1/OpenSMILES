"""Typed P-44 parent facts and the legacy producer boundary adapter."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2 import kind_registry as _kr

_FIXED_MULTI = {
    "diacid": 2, "diester": 2, "dione": 2, "diol": 2, "diamine": 2,
    "triol": 3, "triamine": 3, "tetraamine": 4, "benzenediol": 2,
    "benzenediamine": 2, "cycloalkanediol": 2, "cycloalkanedione": 2,
    "anthraquinone": 2, "benzoquinone": 2, "ortho_benzoquinone": 2,
    "naphthalenediol": 2, "quinolinediol": 2,
}
_DYNAMIC_IDS = {
    "polycarboxylic": "cooh_c_idxs", "benzene_polycarboxylic": "cooh_c_idxs",
    "cycloalkane_polycarboxylic": "cooh_c_idxs",
}


@dataclass(frozen=True, order=True)
class P44Facts:
    principal_group_class: int
    principal_group_count: int


@dataclass(frozen=True)
class ParentCandidate:
    parent: dict
    facts: P44Facts


def principal_contract_kind(kind: str) -> str:
    if kind in _DYNAMIC_IDS:
        return "dynamic"
    if kind in _FIXED_MULTI:
        return "fixed"
    return "single" if _kr.fg_rank(kind) else "none"


def _legacy_count(parent: dict, kind: str) -> int:
    mode = principal_contract_kind(kind)
    if mode == "dynamic":
        return len(parent.get(_DYNAMIC_IDS[kind]) or ())
    if mode == "fixed":
        return _FIXED_MULTI[kind]
    return 1 if mode == "single" else 0


def with_principal_group_contract(parent: dict) -> dict:
    if "principal_group_count" in parent:
        return parent
    kind = parent.get("kind") or ""
    return {**parent, "principal_group_count": _legacy_count(parent, kind)}


def from_parent_dict(parent: dict) -> ParentCandidate:
    kind = parent.get("kind") or ""
    if "principal_group_count" not in parent:
        raise ValueError(f"principal_group_count missing for {kind}")
    facts = P44Facts(_kr.fg_rank(kind), int(parent["principal_group_count"]))
    return ParentCandidate(parent, facts)


def principal_key(parent: dict) -> P44Facts:
    return from_parent_dict(parent).facts


def principal_phases(parents: list[dict]) -> list[list[dict]]:
    keys = sorted({principal_key(parent) for parent in parents}, reverse=True)
    return [[parent for parent in parents if principal_key(parent) == key] for key in keys]
