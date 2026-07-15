"""CIP R/S stereodescriptor prefixes (IUPAC P-92 / P-93)."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

_RS_KINDS = frozenset({
    "acid", "alcohol",
    "diol", "triol", "amine", "diamine",
    "ketone", "ester",
    "piperidine", "pyrrolidine", "piperazine", "morpholine",
    "oxolane", "oxane",
    "piperidinecarboxylic", "pyrrolidinecarboxylic",
    "piperazinecarboxylic", "morpholinecarboxylic",
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
    tag, stem = _strip_stereo(name)
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
