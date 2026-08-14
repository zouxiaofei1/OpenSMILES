"""Layer5 stereodescriptor prefixes: E/Z (P-91.2/P-93.4) + CIP R/S (P-92/P-93).

Merged from stereo_ez.py, stereo_rs.py and _stereo_common.py (same stereo
prefix cluster); _split_stereo_lead is the shared stereo-block splitter.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import BondStereo, Mol


def _split_stereo_lead(name: str) -> tuple[str, str]:
    """Split leading '(…)-' stereo block from a name/stem."""
    if not name.startswith("("):
        return "", name
    close = name.find(")-")
    if close < 0:
        return "", name
    return name[: close + 2], name[close + 2 :]


# --- E/Z stereodescriptors -------------------------------------------------

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


def _ez_parts(mol, chain: list[int], bonds) -> list[tuple[int, str]]:
    """Collect (loc, letter) only for bonds with defined stereo (partial OK)."""
    parts = [_ez_bond_part(mol, chain, b) for b in bonds]
    return sorted((p for p in parts if p is not None), key=lambda x: x[0])


def _ez_multi_prefix(numbered: dict) -> str:
    """Multi-ene prefix from bonds that have stereo: (2E,6Z)- or (14Z)-."""
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), parent.get("chain") or []
    bonds = list(parent.get("double_bonds") or [])
    if mol is None or not bonds or not chain:
        return ""
    parts = _ez_parts(mol, chain, bonds)
    if not parts:
        return ""
    return f"({','.join(f'{loc}{let}' for loc, let in parts)})-"


def ez_for_parent(numbered: dict) -> str:
    """E/Z prefix for parent: multi-ene when double_bonds, else single bond."""
    parent = numbered.get("parent") or {}
    if parent.get("double_bonds"):
        return _ez_multi_prefix(numbered)
    return _ez_prefix(numbered)


# --- CIP R/S stereodescriptors ----------------------------------------------

_RS_KINDS = frozenset({
    "acid", "alcohol",
    "diol", "triol", "amine", "diamine",
    "ketone", "ester",
    "diacid", "amide", "nitrile", "aldehyde", "thiol",
    "piperidine", "pyrrolidine", "piperazine", "morpholine",
    "oxolane", "oxane",
})

# Single-center sat-hetero parents omit locant: (R)- not (3R)-.
_RS_OMIT_LOC = frozenset({
    "piperidine", "pyrrolidine", "piperazine", "morpholine",
    "oxolane", "oxane",
})


def _assign_cip(mol: Mol) -> None:
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)


def _cip_code(atom) -> str | None:
    if not atom.HasProp("_CIPCode"):
        return None
    code = atom.GetProp("_CIPCode")
    return code if code in ("R", "S") else None


def _cip_on_chain(mol: Mol, chain: list[int]) -> list[tuple[int, str]]:
    """Return (locant, R/S) for chiral centers on parent chain."""
    _assign_cip(mol)
    out: list[tuple[int, str]] = []
    for loc, idx in enumerate(chain, 1):
        code = _cip_code(mol.GetAtomWithIdx(int(idx)))
        if code:
            out.append((loc, code))
    return out


def _collapsed_parent(parent: dict) -> bool:
    """Skip RS when parent C1 collapses from a larger/ring molecule."""
    n = int(parent.get("n_carbons") or 0)
    if n > 1:
        return False
    mol = parent.get("mol")
    if mol is None:
        return False
    if mol.GetRingInfo().NumRings() > 0:
        return True
    return mol.GetNumHeavyAtoms() > n + 2


def _rs_parts(numbered: dict) -> list[tuple[int, str]]:
    parent = numbered.get("parent") or {}
    if parent.get("kind") not in _RS_KINDS or _collapsed_parent(parent):
        return []
    mol, chain = parent.get("mol"), parent.get("chain") or []
    if mol is None or not chain:
        return []
    return _cip_on_chain(mol, chain)


def _parse_token(tok: str) -> tuple[int | None, str] | None:
    """'E'→(None,'E'); '8R'→(8,'R'); '2E'→(2,'E')."""
    if tok in ("E", "Z", "R", "S"):
        return None, tok
    i = 0
    while i < len(tok) and tok[i].isdigit():
        i += 1
    if i and tok[i:] in ("E", "Z", "R", "S"):
        return int(tok[:i]), tok[i:]
    return None


def _parse_stereo(tag: str) -> list[tuple[int | None, str]]:
    """Parse '(E)-' / '(2E,6Z)-' / '(E,8R)-' into (loc, letter) list."""
    if not tag.startswith("(") or not tag.endswith(")-"):
        return []
    raw = tag[1:-2]
    out: list[tuple[int | None, str]] = []
    for tok in raw.split(","):
        p = _parse_token(tok.strip())
        if p:
            out.append(p)
    return out


def _fmt_part(loc: int | None, letter: str) -> str:
    return letter if loc is None else f"{loc}{letter}"


def _format_stereo(parts: list[tuple[int | None, str]]) -> str:
    if not parts:
        return ""
    ordered = sorted(parts, key=lambda x: (x[0] is not None, x[0] or 0))
    body = ",".join(_fmt_part(loc, let) for loc, let in ordered)
    return f"({body})-"


def _merge_parts(
    old: list[tuple[int | None, str]], rs: list[tuple[int | None, str]],
) -> list[tuple[int | None, str]]:
    """Keep non-RS stereo; add R/S by locant."""
    keep = [(loc, let) for loc, let in old if let not in ("R", "S")]
    return keep + list(rs)


def _omit_locants(rs: list[tuple[int, str]]) -> list[tuple[int | None, str]]:
    """Drop locants when a single center needs bare (R)/(S)."""
    if len(rs) != 1:
        return [(loc, let) for loc, let in rs]
    return [(None, rs[0][1])]


def _display_rs(
    kind: str | None, rs: list[tuple[int, str]],
) -> list[tuple[int | None, str]]:
    if kind in _RS_OMIT_LOC:
        return _omit_locants(rs)
    return [(loc, let) for loc, let in rs]


def _with_rs(name: str, rs: list[tuple[int | None, str]]) -> str:
    if not rs:
        return name
    tag, stem = _split_stereo_lead(name)
    parts = _merge_parts(_parse_stereo(tag), rs)
    return f"{_format_stereo(parts)}{stem}"


def _ester_en_rs(en: str, rs: list[tuple[int | None, str]]) -> str:
    """Insert RS after alkyl word: 'methyl X' → 'methyl (2S)-X'."""
    if not rs or " " not in en:
        return _with_rs(en, rs)
    alkyl, acyl = en.split(" ", 1)
    return f"{alkyl} {_with_rs(acyl, rs)}"


def apply_rs_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    rs_raw = _rs_parts(numbered)
    if not rs_raw:
        return en, zh
    kind = (numbered.get("parent") or {}).get("kind")
    rs = _display_rs(kind, rs_raw)
    if kind == "ester":
        return _ester_en_rs(en, rs), _with_rs(zh, rs)
    return _with_rs(en, rs), _with_rs(zh, rs)
