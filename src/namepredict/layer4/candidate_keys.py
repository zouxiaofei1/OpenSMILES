"""L4 并列候选母体比较键：P-44.1.1 后缀集、P-45.2.2 前缀集。"""
from __future__ import annotations

from namepredict.layer4.locant_calc import atom_locant, locant_key


def suffix_locant_set(numbered: dict) -> tuple:
    """P-44.1.1：principal 特征基团位次集合（P-14.3.5）。"""
    from namepredict.layer4.numbering_engine import _principal_atoms

    parent = numbered.get("parent") or {}
    chain = parent.get("chain") or []
    labels = (parent.get("numbering_scaffold") or {}).get("labels") or []
    facts = {"labels": labels}
    locs = [loc for atom in _principal_atoms(parent)
            if (loc := atom_locant(chain, atom, facts)) is not None]
    return tuple(sorted(locant_key(x) for x in locs))


def prefix_locant_set(numbered: dict) -> tuple:
    """P-45.2.2：前缀取代基位次集合（P-14.3.5）；无位次前缀不入键。"""
    locs = [s["locant"] for s in (numbered.get("substituents") or []) if s.get("locant") is not None]
    return tuple(sorted(locant_key(x) for x in locs))
