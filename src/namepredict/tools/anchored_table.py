"""Anchored canonical-SMILES lookup table for simple substituents.

build_anchor_submol marks the attach site with a dummy atom (*), so the
canonical SMILES encodes both the substituent shape and the attachment site:
  isopropyl *C(C)C vs n-propyl *CCC; 4-Cl *c1ccc(Cl)cc1 vs 3-Cl *c1cccc(Cl)c1
  vs 2-Cl *c1ccccc1Cl.

Each key maps to a (registry_key, en, zh, paren, kind) tuple:
  - registry_key is a retained_substituents key when the name is name_mode
    sensitive (isopropyl → propan-2-yl under pin), else None.
  - en/zh are the general names (also the pin names when registry_key is None).
  - paren: compound prefix needs parentheses in assembled names.
  - kind classifies the shape: "alkyl" (pure-carbon side chains, the only kind
    the side-alkyl extractor may claim), "aryl", "halo", or "leaf" (hetero
    atom prefixes/functions).  The L3 namer accepts all kinds; the extractor
    filters to "alkyl" so cyano/nitroso etc. are not misclaimed as alkyl sides.
A hit proves the anchored key uniquely identifies a simple substituent.  The
full naming path remains the fallback when the key is absent.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.submol_build import build_anchor_submol


# anchor SMILES → (registry_key | None, en, zh, requires_parentheses, kind)
ANCHOR_TABLE: dict[str, tuple[str | None, str, str, bool, str]] = {
    # linear n-alkyl (C1–C12; not name_mode sensitive)
    "*C": (None, "methyl", "甲基", False, "alkyl"),
    "*CC": (None, "ethyl", "乙基", False, "alkyl"),
    "*CCC": (None, "propyl", "丙基", False, "alkyl"),
    "*CCCC": (None, "butyl", "丁基", False, "alkyl"),
    "*CCCCC": (None, "pentyl", "戊基", False, "alkyl"),
    "*CCCCCC": (None, "hexyl", "己基", False, "alkyl"),
    "*CCCCCCC": (None, "heptyl", "庚基", False, "alkyl"),
    "*CCCCCCCC": (None, "octyl", "辛基", False, "alkyl"),
    "*CCCCCCCCC": (None, "nonyl", "壬基", False, "alkyl"),
    "*CCCCCCCCCC": (None, "decyl", "癸基", False, "alkyl"),
    "*CCCCCCCCCCC": (None, "undecyl", "十一烷基", False, "alkyl"),
    # branched retained (registry keys, name_mode sensitive)
    "*C(C)C": ("isopropyl", "isopropyl", "异丙基", False, "alkyl"),
    "*C(C)(C)C": ("tert-butyl", "tert-butyl", "叔丁基", False, "alkyl"),
    "*CC(C)C": ("isobutyl", "isobutyl", "异丁基", False, "alkyl"),
    "*C(C)CC": ("sec-butyl", "sec-butyl", "仲丁基", False, "alkyl"),
    "*CC(C)(C)C": ("neopentyl", "neopentyl", "新戊基", False, "alkyl"),
    "*CCC(C)C": ("isopentyl", "isopentyl", "异戊基", False, "alkyl"),
    "*C(C)(C)CC": ("2-methylbutan-2-yl", "2-methylbutan-2-yl", "2-甲基丁-2-基", False, "alkyl"),
    # alkenyl retained
    "*C=C": ("vinyl", "vinyl", "乙烯基", False, "alkyl"),
    "*C=C-C": ("allyl", "allyl", "烯丙基", False, "alkyl"),
    "*C=CC": ("allyl", "allyl", "烯丙基", False, "alkyl"),
    "*CC=C": ("allyl", "allyl", "烯丙基", False, "alkyl"),
    "*C(=C)C": ("isopropenyl", "isopropenyl", "异丙烯基", False, "alkyl"),
    # haloalkyl (compound prefixes need parentheses)
    "*CCl": (None, "chloromethyl", "氯甲基", True, "alkyl"),
    "*CBr": (None, "bromomethyl", "溴甲基", True, "alkyl"),
    "*CCCl": (None, "2-chloroethyl", "2-氯乙基", True, "alkyl"),
    "*CCCCl": (None, "3-chloropropyl", "3-氯丙基", True, "alkyl"),
    "*CCCCCl": (None, "4-chlorobutyl", "4-氯丁基", True, "alkyl"),
    "*CCCBr": (None, "3-bromopropyl", "3-溴丙基", True, "alkyl"),
    "*CCCCBr": (None, "4-bromobutyl", "4-溴丁基", True, "alkyl"),
    # 1-cycloalkylethyl (C1(ring)-C(C)H-)
    "*C(C)C1CCCCC1": (None, "1-cyclohexylethyl", "1-环己基乙基", True, "alkyl"),
    # CF3
    "*C(F)(F)F": (None, "trifluoromethyl", "三氟甲基", False, "alkyl"),
    # heteroatom prefixes — alkoxy / thio / sulfinyl / sulfonyl (registry-keyed)
    "*OC": ("methoxy", "methoxy", "甲氧基", False, "leaf"),
    "*SC": ("methylsulfanyl", "methylsulfanyl", "甲硫基", False, "leaf"),
    "*S(C)=O": ("methylsulfinyl", "methylsulfinyl", "甲亚磺酰基", False, "leaf"),
    "*S(C)(=O)=O": ("methylsulfonyl", "methylsulfonyl", "甲磺酰基", False, "leaf"),
    "*S(=O)(=O)O": ("sulfo", "sulfo", "磺基", False, "leaf"),
    "*S(=O)(=O)c1ccc(C)cc1": ("tosyl", "tosyl", "对甲苯磺酰基", False, "leaf"),
    "*S(=O)(=O)C(F)(F)F": ("triflyl", "triflyl", "三氟甲磺酰基", False, "leaf"),
    # nitrogen leaves (registry-keyed)
    "*N=O": ("nitroso", "nitroso", "亚硝基", False, "leaf"),
    "*N=[N+]=[N-]": ("azido", "azido", "叠氮基", False, "leaf"),
    "*[N+]#[C-]": ("isocyano", "isocyano", "异氰基", False, "leaf"),
    "*C#N": ("cyano", "cyano", "氰基", False, "leaf"),
    # piperidinyl
    "*C1CCNCC1": (None, "piperidin-4-yl", "哌啶-4-基", True, "alkyl"),
    "*C1CCCNC1": (None, "piperidin-3-yl", "哌啶-3-基", True, "alkyl"),
    "*[C@@H]1CCCCN1": (None, "piperidin-2-yl", "哌啶-2-基", True, "alkyl"),
    # alkenyl-branched retained
    "*C=C(C)C": ("isobutyl", "isobutyl", "异丁基", False, "alkyl"),
    "*C=CC(C)C": ("isopentyl", "isopentyl", "异戊基", False, "alkyl"),
    "*CC=C(C)C": ("3-methylbut-2-enyl", "3-methylbut-2-enyl", "3-甲基丁-2-烯基", False, "alkyl"),
    # cycloalkyl (not name_mode sensitive)
    "*C1CC1": (None, "cyclopropyl", "环丙基", False, "alkyl"),
    "*C1CCC1": (None, "cyclobutyl", "环丁基", False, "alkyl"),
    "*C1CCCC1": (None, "cyclopentyl", "环戊基", False, "alkyl"),
    "*C1CCCCC1": (None, "cyclohexyl", "环己基", False, "alkyl"),
    "*C1CCCCCC1": (None, "cycloheptyl", "环庚基", False, "alkyl"),
    "*C1CCCCCCC1": (None, "cyclooctyl", "环辛基", False, "alkyl"),
    # aryl
    "*c1ccccc1": (None, "phenyl", "苯基", False, "aryl"),
    "*c1ccc(Cl)cc1": (None, "4-chlorophenyl", "4-氯苯基", True, "aryl"),
    "*c1cccc(Cl)c1": (None, "3-chlorophenyl", "3-氯苯基", True, "aryl"),
    "*c1ccccc1Cl": (None, "2-chlorophenyl", "2-氯苯基", True, "aryl"),
    # single-atom halogens (always single-bonded; no bond-type ambiguity)
    "*F": (None, "fluoro", "氟", False, "halo"),
    "*Cl": (None, "chloro", "氯", False, "halo"),
    "*Br": (None, "bromo", "溴", False, "halo"),
    "*I": (None, "iodo", "碘", False, "halo"),
    # multi-atom FG leaves with unique bond topology
    "*[N+](=O)[O-]": (None, "nitro", "硝基", False, "leaf"),
    "*N=C=O": (None, "isocyanato", "异氰酸根合", False, "leaf"),
    "*N=C=S": (None, "isothiocyanato", "异硫氰酸根合", False, "leaf"),
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


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> tuple[str, tuple[str | None, str, str, bool, str]] | None:
    key = anchored_key(mol, atoms, attach_old)
    if key is None:
        return None
    hit = ANCHOR_TABLE.get(key)
    return None if hit is None else (key, hit)


def _resolve(
    hit: tuple[str | None, str, str, bool, str], *, name_mode: str,
) -> tuple[str, str, bool, str]:
    from namepredict.layer3.retained_substituents import resolve_name

    reg_key, en, zh, paren, kind = hit
    if reg_key is not None:
        en, zh = resolve_name(reg_key, name_mode=name_mode)
    return en, zh, paren, kind


def anchored_entry(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None, *, name_mode: str = "general",
) -> tuple[str, str, bool, str] | None:
    """Resolve an atom set to (en, zh, paren, kind) under name_mode, or None."""
    got = _table_hit(mol, atoms, attach_old)
    return None if got is None else _resolve(got[1], name_mode=name_mode)


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
    *, name_mode: str = "general",
) -> tuple[str, str, bool] | None:
    """Look up a substituent atom set; returns (en, zh, paren) under name_mode.

    registry-keyed entries resolve through retained_substituents.resolve_name
    so pin mode yields e.g. propan-2-yl for isopropyl.  None when no hit.
    """
    entry = anchored_entry(mol, atoms, attach_old, name_mode=name_mode)
    return None if entry is None else (entry[0], entry[1], entry[2])
