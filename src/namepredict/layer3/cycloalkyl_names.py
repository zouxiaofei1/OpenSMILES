"""Cycloalkyl substituent naming (L3; P-29.6)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_facts import (
    cycloalkyl_side,
    cycloalkylethyl_side,
    saturated_heterocycle_side,
)

_CYCLO_EN = {3: "cyclopropyl", 4: "cyclobutyl", 5: "cyclopentyl", 6: "cyclohexyl", 7: "cycloheptyl", 8: "cyclooctyl"}
_CYCLO_ZH = {3: "环丙基", 4: "环丁基", 5: "环戊基", 6: "环己基", 7: "环庚基", 8: "环辛基"}


def _cyclo_names(size: int) -> tuple[str, str] | None:
    return (_CYCLO_EN[size], _CYCLO_ZH[size]) if size in _CYCLO_EN else None


def _make_cycloalkyl_branch(
    attach: int, atoms: list[int], en: str, zh: str, paren: bool = False,
) -> dict:
    return {
        "kind": "alkyl", "n_carbons": len(atoms), "attach_idx": attach,
        "atoms": atoms, "en": en, "zh": zh, "paren": paren,
    }


def _one_cycloalkylethyl(mol: Mol, attach: int, start: int,
                         chain_set: set[int]) -> dict | None:
    fact = cycloalkylethyl_side(mol, start, chain_set)
    atoms = list(fact.atoms) if fact else []
    names = _cyclo_names(len(atoms) - 2)
    if not names:
        return None
    return _make_cycloalkyl_branch(attach, atoms, f"1-{names[0]}ethyl", f"1-{names[1]}乙基", True)


def _one_monocycloalkyl(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    fact = cycloalkyl_side(mol, start, chain_set)
    atoms = list(fact.atoms) if fact else []
    names = _cyclo_names(len(atoms))
    return _make_cycloalkyl_branch(attach, atoms, *names) if names else None


_SAT_HETERO = {
    (5, frozenset(((7, 1),))): ("pyrrolidin", "吡咯烷", 7),
    (6, frozenset(((7, 1),))): ("piperidin", "哌啶", 7),
    (6, frozenset(((8, 1),))): ("oxan", "噁烷", 8),
    (5, frozenset(((8, 1),))): ("oxolan", "氧杂环戊烷", 8),
    (6, frozenset(((8, 1), (7, 1)))): ("morpholin", "吗啉", 7),
}


def _ring_distance(mol: Mol, ring: set[int], start: int, target: int) -> int:
    frontier, seen, distance = {start}, {start}, 0
    while target not in frontier:
        frontier = {n.GetIdx() for i in frontier for n in mol.GetAtomWithIdx(i).GetNeighbors()
                    if n.GetIdx() in ring and n.GetIdx() not in seen}
        seen |= frontier
        distance += 1
    return distance


def _sat_hetero_names(mol: Mol, fact) -> tuple[str, str]:
    en, zh, z = _SAT_HETERO[fact.signature]
    hetero = next(i for i in fact.ring if mol.GetAtomWithIdx(i).GetAtomicNum() == z)
    locant = _ring_distance(mol, set(fact.ring), hetero, fact.root) + 1
    return f"{en}-{locant}-yl", f"{zh}-{locant}-基"


def _one_sat_hetero_side(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    fact = saturated_heterocycle_side(mol, start, chain_set)
    if fact is None:
        return None
    en, zh = _sat_hetero_names(mol, fact)
    return _make_cycloalkyl_branch(attach, list(fact.ring), en, zh, True)


def _one_cycloalkyl_side(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    return (
        _one_cycloalkylethyl(mol, attach, start, chain_set)
        or _one_monocycloalkyl(mol, attach, start, chain_set)
        or _one_sat_hetero_side(mol, attach, start, chain_set)
    )
