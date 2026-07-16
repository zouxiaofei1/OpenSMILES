"""Shared aryl/cyclo side stem fields for special FG parents (L2 pack).

Precompute EN/ZH into side dicts so L5 never walks via L2 private APIs.
"""
from __future__ import annotations


def aryl_en_zh(mol, c_idx: int, parent: int) -> tuple[str, str]:
    """Phenyl stem EN/ZH from attach carbon and parent heteroatom."""
    from namepredict.layer2.aryl_sub import _phenyl_at, _phenyl_name
    ph = _phenyl_at(mol, c_idx, parent)
    if ph is None:
        return "phenyl", "苯基"
    en, zh, _ = _phenyl_name(mol, ph, c_idx)
    return en, zh


def with_aryl_names(side: dict, mol, c_idx: int, parent: int) -> dict:
    """Copy side and attach en/zh from aryl stem (keep ph/c for debug)."""
    en, zh = aryl_en_zh(mol, c_idx, parent)
    return {**side, "en": en, "zh": zh}


def with_cyclo_names(side: dict) -> dict:
    """Attach cycloalkyl en/zh from size field."""
    from namepredict.layer2.side_cycloalkyl import _cycloalkyl_names
    names = _cycloalkyl_names(int(side.get("size") or 0))
    if names is None:
        return side
    return {**side, "en": names[0], "zh": names[1]}
