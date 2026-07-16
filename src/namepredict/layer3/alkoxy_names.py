"""Alkoxy substituent names: linear C1–C4 and PEG tails (L3)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_facts import outer_alkoxy

_ALKOXY_BASE_EN = {1: "methoxy", 2: "ethoxy", 3: "propoxy", 4: "butoxy"}
_ALKOXY_BASE_ZH = {1: "甲氧基", 2: "乙氧基", 3: "丙氧基", 4: "丁氧基"}
_BRANCHED_ALKOXY = {31: ("isopropoxy", "异丙氧基"), 41: ("isobutoxy", "异丁氧基")}


def _peg_k1_en(n_term: int) -> str:
    return "2-methoxyethoxy" if n_term == 1 else "2-ethoxyethoxy"


def _peg_alkoxy_en(peg_k: int, n_term: int) -> str:
    if peg_k == 0:
        return _ALKOXY_BASE_EN[n_term]
    if peg_k == 1:
        return _peg_k1_en(n_term)
    return f"2-({_peg_alkoxy_en(peg_k - 1, n_term)})ethoxy"


def _peg_alkoxy_zh(peg_k: int, n_term: int) -> str:
    base = _ALKOXY_BASE_ZH[n_term]
    if peg_k == 0:
        return base
    if peg_k == 1:
        return f"2-{base}乙氧基"
    return f"2-({_peg_alkoxy_zh(peg_k - 1, n_term)})乙氧基"


def _peg_code_ok(s: str, n_term: int, peg_k: int) -> bool:
    if n_term not in (1, 2) or not (1 <= peg_k <= 3):
        return False
    return peg_k == 1 or set(s[:-2]) == {"1"}


def _decode_alkoxy_code(code: int) -> tuple[int, int] | None:
    """Map code → (peg_k, n_terminal); simple 1–4 or PEG 12/22/112/122/…"""
    if code in _ALKOXY_BASE_EN:
        return 0, code
    s = str(code)
    if len(s) < 2 or s[-1] != "2" or s[-2] not in "12":
        return None
    n_term, peg_k = int(s[-2]), len(s) - 1
    return (peg_k, n_term) if _peg_code_ok(s, n_term, peg_k) else None


def _alkoxy_names(code: int) -> tuple[str, str] | None:
    if code in _BRANCHED_ALKOXY:
        return _BRANCHED_ALKOXY[code]
    dec = _decode_alkoxy_code(code)
    if dec is None:
        return None
    peg_k, n_term = dec
    return _peg_alkoxy_en(peg_k, n_term), _peg_alkoxy_zh(peg_k, n_term)


def _make_alkoxy(attach: int, o_idx: int, atoms: list[int], n: int) -> dict | None:
    names = _alkoxy_names(n)
    if names is None:
        return None
    en, zh = names
    n_c = {31: 3, 41: 4}.get(n, n)
    return {"kind": "alkoxy", "attach_idx": attach, "atoms": [o_idx] + atoms,
            "n_carbons": n_c, "en": en, "zh": zh, "paren": False}


def _alkoxy_ends(e: dict, chain_set: set[int]) -> tuple[int, int, int] | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in chain_set) == (c2 in chain_set):
        return None
    ring_c, outer = (c1, c2) if c1 in chain_set else (c2, c1)
    return o, ring_c, outer


def _one_ring_alkoxy(mol: Mol, e: dict, chain_set: set[int]) -> dict | None:
    ends = _alkoxy_ends(e, chain_set)
    if ends is None:
        return None
    o, ring_c, outer = ends
    if not mol.GetAtomWithIdx(ring_c).IsInRing() or mol.GetAtomWithIdx(outer).GetIsAromatic():
        return None
    fact = outer_alkoxy(mol, outer, o)
    return _make_alkoxy(ring_c, o, list(fact.atoms), fact.code) if fact else None


def _extract_alkoxys(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _one_ring_alkoxy(mol, e, chain_set)
        if one is not None:
            out.append(one)
    return out
