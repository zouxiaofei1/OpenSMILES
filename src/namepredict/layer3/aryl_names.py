"""基于 L2 芳基拓扑事实组装侧链名称。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.side_facts import ArylArmFact, ArylArmKind


_STEMS = {
    ArylArmKind.DIRECT_C: ("phenyl", "苯基"),
    ArylArmKind.DIRECT_O: ("phenoxy", "苯氧基"),
    ArylArmKind.METHYLENE_C: ("benzyl", "苄基"),
    ArylArmKind.O_METHYLENE_C: ("benzyloxy", "苄氧基"),
}


def _ring_attachment(mol: Mol, fact: ArylArmFact) -> int:
    if fact.kind in (ArylArmKind.DIRECT_C, ArylArmKind.DIRECT_O):
        return fact.outer
    for idx in fact.ring:
        if any(n.GetIdx() == fact.outer for n in mol.GetAtomWithIdx(idx).GetNeighbors()):
            return idx
    return min(fact.ring)


def _stem(en: str, zh: str, kind: ArylArmKind) -> tuple[str, str]:
    target_en, target_zh = _STEMS[kind]
    if en.endswith("phenyl"):
        en = en[:-6] + target_en
    if zh.endswith("苯基"):
        zh = zh[:-2] + target_zh
    return en, zh


def _ring_parent(fact: ArylArmFact) -> int:
    if fact.kind == ArylArmKind.DIRECT_C:
        return fact.attachment
    bridged = (ArylArmKind.METHYLENE_C, ArylArmKind.O_METHYLENE_C)
    return fact.outer if fact.kind in bridged else fact.bridge or -1


def aryl_arm_name(mol: Mol, fact: ArylArmFact, owned_atoms) -> tuple[str, str, bool]:
    """Name an aryl arm via the generic cut→free-name→yl pipeline.

    The benzene ring is the free parent; leaves (halo/Me/alkoxy/nitro/
    OH/NH2/CN/nested Ph) are collected topologically so the cut submol
    carries the full arm. name_as_substituent re-names the anchored ring as a
    phenyl radical (P-41), numbering the attach carbon locant 1 (P-29.6.2).

    owned_atoms is the parent boundary (parent.owned_atoms). The cut atom
    set is the full non-parent connected component (side_atoms) minus the
    kind-encoded arm atoms (fact.atoms - fact.ring, e.g. phenoxy's O or
    benzyl's CH2), which the _stem suffix carries instead.
    """
    from namepredict.layer3.as_substituent import name_as_substituent
    from namepredict.tools.block_cut import side_atoms

    attach = _ring_attachment(mol, fact)
    seed = frozenset(fact.atoms)
    arm = seed - frozenset(fact.ring)
    parent = _ring_parent(fact)
    atoms = side_atoms(mol, frozenset(owned_atoms), attach, seed) - arm - {parent}
    hit = name_as_substituent(mol, attach, atoms)
    if hit is None:
        return "", "", False
    en, zh, paren = hit
    en, zh = _stem(en, zh, fact.kind)
    return en, zh, paren
