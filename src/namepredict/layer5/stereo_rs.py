"""CIP R/S stereodescriptor prefixes (IUPAC P-92 / P-93)."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

_RS_KINDS = frozenset({
    "acid", "alkenoic_acid", "alcohol", "alkenol",
    "diol", "triol", "amine", "diamine",
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


def _rs_parts(numbered: dict) -> list[tuple[int, str]]:
    parent = numbered.get("parent") or {}
    if parent.get("kind") not in _RS_KINDS:
        return []
    mol, chain = parent.get("mol"), parent.get("chain") or []
    if mol is None or not chain:
        return []
    return _cip_on_chain(mol, chain)


def _strip_stereo(name: str) -> tuple[str, str]:
    """Split leading '(…)-' stereo block from name."""
    if not name.startswith("("):
        return "", name
    close = name.find(")-")
    if close < 0:
        return "", name
    return name[: close + 2], name[close + 2 :]


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
    old: list[tuple[int | None, str]], rs: list[tuple[int, str]],
) -> list[tuple[int | None, str]]:
    """Keep non-RS stereo; add R/S by locant."""
    keep = [(loc, let) for loc, let in old if let not in ("R", "S")]
    return keep + [(loc, let) for loc, let in rs]


def _with_rs(name: str, rs: list[tuple[int, str]]) -> str:
    if not rs:
        return name
    tag, stem = _strip_stereo(name)
    parts = _merge_parts(_parse_stereo(tag), rs)
    return f"{_format_stereo(parts)}{stem}"


def apply_rs_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    rs = _rs_parts(numbered)
    if not rs:
        return en, zh
    return _with_rs(en, rs), _with_rs(zh, rs)
