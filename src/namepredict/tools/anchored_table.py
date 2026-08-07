"""Anchored canonical-SMILES lookup table for simple substituents.

build_anchor_submol marks the attach site with a dummy atom (*), so the
canonical SMILES encodes both the substituent shape and the attachment site:
  isopropyl *C(C)C vs n-propyl *CCC; 4-Cl *c1ccc(Cl)cc1 vs 3-Cl *c1cccc(Cl)c1
  vs 2-Cl *c1ccccc1Cl.

Each key maps to a (registry_key, en, zh, paren) tuple:
  - registry_key is a retained_substituents key when the name is name_mode
    sensitive (isopropyl → propan-2-yl under pin), else None.
  - en/zh are the general names (also the pin names when registry_key is None).
A hit proves the anchored key uniquely identifies a simple substituent.  The
full naming path remains the fallback when the key is absent.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.submol_build import build_anchor_submol


# anchor SMILES → (registry_key | None, en, zh, requires_parentheses)
ANCHOR_TABLE: dict[str, tuple[str | None, str, str, bool]] = {
    # linear n-alkyl (C1–C12; not name_mode sensitive)
    "*C": (None, "methyl", "甲基", False),
    "*CC": (None, "ethyl", "乙基", False),
    "*CCC": (None, "propyl", "丙基", False),
    "*CCCC": (None, "butyl", "丁基", False),
    "*CCCCC": (None, "pentyl", "戊基", False),
    "*CCCCCC": (None, "hexyl", "己基", False),
    "*CCCCCCC": (None, "heptyl", "庚基", False),
    "*CCCCCCCC": (None, "octyl", "辛基", False),
    "*CCCCCCCCC": (None, "nonyl", "壬基", False),
    "*CCCCCCCCCC": (None, "decyl", "癸基", False),
    "*CCCCCCCCCCC": (None, "undecyl", "十一烷基", False),
    # branched retained (registry keys, name_mode sensitive)
    "*C(C)C": ("isopropyl", "isopropyl", "异丙基", False),
    "*C(C)(C)C": ("tert-butyl", "tert-butyl", "叔丁基", False),
    "*CC(C)C": ("isobutyl", "isobutyl", "异丁基", False),
    "*C(C)CC": ("sec-butyl", "sec-butyl", "仲丁基", False),
    "*CC(C)(C)C": ("neopentyl", "neopentyl", "新戊基", False),
    "*CCC(C)C": ("isopentyl", "isopentyl", "异戊基", False),
    "*C(C)(C)CC": ("2-methylbutan-2-yl", "2-methylbutan-2-yl", "2-甲基丁-2-基", False),
    # alkenyl retained
    "*C=C": ("vinyl", "vinyl", "乙烯基", False),
    "*C=C-C": ("allyl", "allyl", "烯丙基", False),
    "*C=CC": ("allyl", "allyl", "烯丙基", False),
    "*CC=C": ("allyl", "allyl", "烯丙基", False),
    "*C(=C)C": ("isopropenyl", "isopropenyl", "异丙烯基", False),
    # haloalkyl (compound prefixes need parentheses)
    "*CCl": (None, "chloromethyl", "氯甲基", True),
    "*CBr": (None, "bromomethyl", "溴甲基", True),
    "*CCCl": (None, "2-chloroethyl", "2-氯乙基", True),
    "*CCCCl": (None, "3-chloropropyl", "3-氯丙基", True),
    "*CCCCCl": (None, "4-chlorobutyl", "4-氯丁基", True),
    "*CCCBr": (None, "3-bromopropyl", "3-溴丙基", True),
    "*CCCCBr": (None, "4-bromobutyl", "4-溴丁基", True),
    # 1-cycloalkylethyl (C1(ring)-C(C)H-)
    "*C(C)C1CCCCC1": (None, "1-cyclohexylethyl", "1-环己基乙基", True),
    # CF3
    "*C(F)(F)F": (None, "trifluoromethyl", "三氟甲基", False),
    # piperidinyl
    "*C1CCNCC1": (None, "piperidin-4-yl", "哌啶-4-基", True),
    "*C1CCCNC1": (None, "piperidin-3-yl", "哌啶-3-基", True),
    "*[C@@H]1CCCCN1": (None, "piperidin-2-yl", "哌啶-2-基", True),
    # alkenyl-branched retained
    "*C=C(C)C": ("isobutyl", "isobutyl", "异丁基", False),
    "*C=CC(C)C": ("isopentyl", "isopentyl", "异戊基", False),
    "*CC=C(C)C": ("3-methylbut-2-enyl", "3-methylbut-2-enyl", "3-甲基丁-2-烯基", False),
    # cycloalkyl (not name_mode sensitive)
    "*C1CC1": (None, "cyclopropyl", "环丙基", False),
    "*C1CCC1": (None, "cyclobutyl", "环丁基", False),
    "*C1CCCC1": (None, "cyclopentyl", "环戊基", False),
    "*C1CCCCC1": (None, "cyclohexyl", "环己基", False),
    "*C1CCCCCC1": (None, "cycloheptyl", "环庚基", False),
    "*C1CCCCCCC1": (None, "cyclooctyl", "环辛基", False),
    # aryl
    "*c1ccccc1": (None, "phenyl", "苯基", False),
    "*c1ccc(Cl)cc1": (None, "4-chlorophenyl", "4-氯苯基", True),
    "*c1cccc(Cl)c1": (None, "3-chlorophenyl", "3-氯苯基", True),
    "*c1ccccc1Cl": (None, "2-chlorophenyl", "2-氯苯基", True),
    # single-atom halogens (always single-bonded; no bond-type ambiguity)
    "*F": (None, "fluoro", "氟", False),
    "*Cl": (None, "chloro", "氯", False),
    "*Br": (None, "bromo", "溴", False),
    "*I": (None, "iodo", "碘", False),
    # multi-atom FG leaves with unique bond topology
    "*[N+](=O)[O-]": (None, "nitro", "硝基", False),
    "*N=C=O": (None, "isocyanato", "异氰酸根合", False),
    "*N=C=S": (None, "isothiocyanato", "异硫氰酸根合", False),
}


def pick_root(mol: Mol, atoms: frozenset[int]) -> int:
    """Substituent-side attach atom: the atom in `atoms` bonded to an outside
    heavy atom (the parent).  Falls back to the lowest index."""
    for a in atoms:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in atoms:
                return a
    return min(atoms)


def anchored_key(mol: Mol, atoms: frozenset[int], attach_old: int | None = None) -> str | None:
    """Anchored canonical SMILES for a substituent atom set.

    attach_old is the substituent-side attach atom; auto-picked via pick_root
    when None (B-path dicts carry the parent-side attach_idx, not this one).
    """
    from rdkit.Chem import MolToSmiles

    if attach_old is None:
        attach_old = pick_root(mol, atoms)
    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
    *, name_mode: str = "general",
) -> tuple[str, str, bool] | None:
    """Look up a substituent atom set; returns (en, zh, paren) under name_mode.

    registry-keyed entries resolve through retained_substituents.resolve_name
    so pin mode yields e.g. propan-2-yl for isopropyl.  None when no hit.
    """
    from namepredict.layer3.retained_substituents import resolve_name

    key = anchored_key(mol, atoms, attach_old)
    if key is None:
        return None
    hit = ANCHOR_TABLE.get(key)
    if hit is None:
        return None
    reg_key, en, zh, paren = hit
    if reg_key is not None:
        en, zh = resolve_name(reg_key, name_mode=name_mode)
    return en, zh, paren
