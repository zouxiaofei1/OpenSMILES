"""Alkenoic / alkenedioic acid name assembly and E/Z stereo prefix."""
from __future__ import annotations

from rdkit.Chem import BondStereo, Mol

from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH, zh_stem


def _stereo_tag(st) -> str:
    if st == BondStereo.STEREOE:
        return "(E)-"
    if st == BondStereo.STEREOZ:
        return "(Z)-"
    return ""


def _bond_stereo(mol: Mol | None, double_bond) -> str:
    if mol is None or not double_bond:
        return ""
    c1, c2 = double_bond
    bond = mol.GetBondBetweenAtoms(int(c1), int(c2))
    return _stereo_tag(bond.GetStereo()) if bond is not None else ""


def _ez_prefix(numbered: dict) -> str:
    parent = numbered.get("parent") or {}
    return _bond_stereo(parent.get("mol"), parent.get("double_bond"))


def _bond_min_loc(chain: list[int], pair) -> int | None:
    if not pair or pair[0] not in chain or pair[1] not in chain:
        return None
    return min(chain.index(pair[0]) + 1, chain.index(pair[1]) + 1)


def _ez_letter(tag: str) -> str:
    """'(E)-' → 'E'; empty → ''."""
    return tag[1] if len(tag) >= 3 and tag[0] == "(" else ""


def _ez_bond_part(mol, chain: list[int], bond) -> tuple[int, str] | None:
    loc = _bond_min_loc(chain, bond)
    letter = _ez_letter(_bond_stereo(mol, bond))
    return (loc, letter) if loc is not None and letter else None


def _ez_parts(mol, chain: list[int], bonds) -> list[tuple[int, str]] | None:
    parts = [_ez_bond_part(mol, chain, b) for b in bonds]
    if any(p is None for p in parts):
        return None
    return sorted(parts, key=lambda x: x[0])


def _ez_multi_prefix(numbered: dict) -> str:
    """Multi-ene prefix: (2E,6Z)- when all bonds have stereo; else ''."""
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), parent.get("chain") or []
    bonds = list(parent.get("double_bonds") or [])
    if mol is None or not bonds or not chain:
        return ""
    parts = _ez_parts(mol, chain, bonds)
    if not parts:
        return ""
    return f"({','.join(f'{loc}{let}' for loc, let in parts)})-"


def _ez_for_alkenol(numbered: dict) -> str:
    parent = numbered.get("parent") or {}
    if parent.get("double_bonds"):
        return _ez_multi_prefix(numbered)
    return _ez_prefix(numbered)


def _unsat_acid_pair(n, locant, ez, en_sfx, zh_sfx, min_n=2) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or locant is None or n < min_n:
        return None
    return f"{ez}{en[:-3]}-{locant}-{en_sfx}", f"{ez}{zh_stem(zh)}-{locant}-{zh_sfx}"


def _alkenedioic_names(n: int, locant: int | None, ez: str) -> tuple[str, str] | None:
    return _unsat_acid_pair(n, locant, ez, "enedioic acid", "烯二酸", 3)


def alkenedioic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    return _alkenedioic_names(n, numbered.get("ene_locant"), _ez_prefix(numbered))


def _alkenoic_acid_names(n: int, locant: int | None, ez: str) -> tuple[str, str] | None:
    return _unsat_acid_pair(n, locant, ez, "enoic acid", "烯酸")


def _ene_mult_acid(k: int) -> tuple[str, str]:
    en = {
        2: "dienoic acid", 3: "trienoic acid", 4: "tetraenoic acid",
        5: "pentaenoic acid", 6: "hexaenoic acid",
    }.get(k, "")
    zh = {
        2: "二烯酸", 3: "三烯酸", 4: "四烯酸", 5: "五烯酸", 6: "六烯酸",
    }.get(k, "")
    return en, zh


def _polyalkenoic_acid_names(n: int, locs, ez: str) -> tuple[str, str] | None:
    """octadeca-9,11-dienoic acid / 十八-9,11-二烯酸."""
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    me, mz = _ene_mult_acid(len(locs or []))
    if not en or not zh or not me or not locs or len(locs) < 2 or n < 4:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"{ez}{en[:-3]}a-{loc}-{me}", f"{ez}{zh_stem(zh)}-{loc}-{mz}"


def alkenoic_acid_names(n: int, numbered: dict) -> tuple[str, str] | None:
    locs = numbered.get("ene_locants")
    if locs and len(locs) >= 2:
        return _polyalkenoic_acid_names(n, locs, _ez_for_alkenol(numbered))
    return _alkenoic_acid_names(
        n, numbered.get("ene_locant"), _ez_prefix(numbered),
    )


def alkenamide_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """prop-2-enamide / (E)-but-2-enamide with E/Z when stereo defined."""
    return _unsat_acid_pair(
        n, numbered.get("ene_locant"), _ez_prefix(numbered), "enamide", "烯酰胺",
    )


def unsat_carbonyl_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Dispatch alkenoic / alkenedioic / alkenamide parent stems."""
    if kind == "alkenedioic":
        return alkenedioic_names(n, numbered)
    if kind == "alkenoic_acid":
        return alkenoic_acid_names(n, numbered)
    return alkenamide_names(n, numbered) if kind == "alkenamide" else None
