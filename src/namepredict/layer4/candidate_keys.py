"""L4 并列候选母体比较键：P-44.1.1 principal 特征基团位次集合、P-45.2.2 前缀取代基位次集合。"""
from __future__ import annotations

from namepredict.layer4.locant_key import locant_key


def suffix_locant_set(numbered: dict) -> tuple:
    """P-44.1.1 键：principal 特征基团（后缀）附着原子的位次集合（升序、保留重复，P-14.3.5）。"""
    from namepredict.layer4.numbering_engine import _principal_atoms

    parent = numbered.get("parent") or {}
    chain = parent.get("chain") or []
    labels = (parent.get("numbering_scaffold") or {}).get("labels") or []
    locs = []
    for atom in _principal_atoms(parent):
        if atom in chain:
            i = chain.index(atom)
            locs.append(labels[i] if len(labels) == len(chain) else i + 1)
    return tuple(sorted(locant_key(x) for x in locs))


def prefix_locant_set(numbered: dict) -> tuple:
    """P-45.2.2 键：以前缀引用的取代基位次集合（升序、保留重复，P-14.3.5）；无位次前缀不入键。"""
    locs = [s["locant"] for s in (numbered.get("substituents") or []) if s.get("locant") is not None]
    return tuple(sorted(locant_key(x) for x in locs))
