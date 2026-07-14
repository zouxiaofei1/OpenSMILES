"""Cycloalkyl substituent naming (L3; P-29.6)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_cycloalkyl import (
    _cycloalkyl_names,
    _is_1_cycloalkylethyl,
    _is_monocycloalkyl,
)
from namepredict.layer2.side_sat_hetero import is_sat_hetero_side


def _make_cycloalkyl_branch(
    attach: int, atoms: list[int], en: str, zh: str, paren: bool = False,
) -> dict:
    return {
        "kind": "alkyl", "n_carbons": len(atoms), "attach_idx": attach,
        "atoms": atoms, "en": en, "zh": zh, "paren": paren,
    }


def _one_cycloalkylethyl(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    eth = _is_1_cycloalkylethyl(mol, start, chain_set)
    if eth is None:
        return None
    names = _cycloalkyl_names(len(eth) - 2)
    if not names:
        return None
    en, zh = f"1-{names[0]}ethyl", f"1-{names[1]}乙基"
    return _make_cycloalkyl_branch(attach, eth, en, zh, True)


def _one_monocycloalkyl(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    cyc = _is_monocycloalkyl(mol, start, chain_set)
    if cyc is None:
        return None
    names = _cycloalkyl_names(len(cyc))
    if not names:
        return None
    return _make_cycloalkyl_branch(attach, cyc, names[0], names[1])


def _one_sat_hetero_side(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    got = is_sat_hetero_side(mol, start, chain_set)
    if got is None:
        return None
    atoms, en, zh = got
    return _make_cycloalkyl_branch(attach, atoms, en, zh, True)


def _one_cycloalkyl_side(
    mol: Mol, attach: int, start: int, chain_set: set[int],
) -> dict | None:
    return (
        _one_cycloalkylethyl(mol, attach, start, chain_set)
        or _one_monocycloalkyl(mol, attach, start, chain_set)
        or _one_sat_hetero_side(mol, attach, start, chain_set)
    )
