"""Typed L3 substituent namer: ordered retained → rooted_tree → recursive."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence

from namepredict.layer2.claimable_block import ClaimedBlock
from namepredict.layer3.as_substituent import name_as_substituent


@dataclass(frozen=True)
class SubstituentName:
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained", "rooted_tree", or "recursive"


class SubstituentBackend(Protocol):
    name: str

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None: ...


def _retained_hit(claim: ClaimedBlock, en: str, zh: str, paren: bool) -> SubstituentName:
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="retained",
    )


def _try_cycloalkyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.side_cycloalkyl import _cycloalkyl_names, _is_monocycloalkyl

    parent = {claim.attach_parent}
    atoms = _is_monocycloalkyl(mol, claim.root, parent)
    if atoms is None or set(atoms) != set(claim.atoms):
        return None
    names = _cycloalkyl_names(len(atoms))
    return None if names is None else _retained_hit(claim, names[0], names[1], False)


def _try_sub_phenyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.aryl_sub import _phenyl_at
    from namepredict.layer3.ring_namer import recursive_ph_name

    ph = _phenyl_at(mol, claim.root, claim.attach_parent)
    if ph is None:
        return None
    en, zh, paren, atoms = recursive_ph_name(mol, ph, claim.root, claim.attach_parent)
    return None if set(atoms) != set(claim.atoms) or not en else _retained_hit(claim, en, zh, paren)


def _methoxy_o(mol, claim: ClaimedBlock) -> int | None:
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() == 8:
        return claim.root
    return next((i for i in claim.atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 8), None)


def _is_plain_methyl_on_o(mol, o_idx: int, c_idx: int) -> bool:
    if mol.GetAtomWithIdx(c_idx).GetAtomicNum() != 6:
        return False
    heavies = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == o_idx


def _try_methoxy(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """O–CH3 side attached at claim.attach_parent (ether O or chain C via O root)."""
    if len(claim.atoms) != 2:
        return None
    o_idx = _methoxy_o(mol, claim)
    if o_idx is None:
        return None
    c_idxs = [i for i in claim.atoms if i != o_idx]
    if len(c_idxs) != 1 or not _is_plain_methyl_on_o(mol, o_idx, c_idxs[0]):
        return None
    return _retained_hit(claim, "methoxy", "甲氧基", False)


def _methylsulfanyl_s(mol, claim: ClaimedBlock) -> int | None:
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() == 16:
        return claim.root
    return next((i for i in claim.atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 16), None)


def _is_plain_methyl_on_s(mol, s_idx: int, c_idx: int) -> bool:
    if mol.GetAtomWithIdx(c_idx).GetAtomicNum() != 6:
        return False
    heavies = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _try_methylsulfanyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """CH3–S side attached at claim.attach_parent (thioether S or chain C via S root)."""
    if len(claim.atoms) != 2:
        return None
    s_idx = _methylsulfanyl_s(mol, claim)
    if s_idx is None:
        return None
    c_idxs = [i for i in claim.atoms if i != s_idx]
    if len(c_idxs) != 1 or not _is_plain_methyl_on_s(mol, s_idx, c_idxs[0]):
        return None
    return _retained_hit(claim, "methylsulfanyl", "甲硫基", False)


def _retained_name(mol, claim: ClaimedBlock) -> SubstituentName | None:
    return (
        _try_cycloalkyl(mol, claim)
        or _try_sub_phenyl(mol, claim)
        or _try_methoxy(mol, claim)
        or _try_methylsulfanyl(mol, claim)
    )


class RetainedBackend:
    """Retained simple leaves: cycloalkyl, substituted phenyl, methoxy."""

    name = "retained"

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return _retained_name(mol, claim)


def _rooted_tree_name(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.side_alkyl_sys import build_rooted_alkyl_tree
    from namepredict.layer3.alkyl_sys_names import name_rooted_alkyl

    tree = build_rooted_alkyl_tree(mol, root=claim.root, atoms=claim.atoms)
    hit = None if tree is None else name_rooted_alkyl(tree)
    if hit is None or not hit[0] or not hit[1]:
        return None
    en, zh, paren = hit
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="rooted_tree",
    )


class RootedTreeBackend:
    """Pure saturated-carbon rooted tree (max_atoms=12, max_depth=3)."""

    name = "rooted_tree"

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return _rooted_tree_name(mol, claim)


class RecursiveBackend:
    """Bounded recursive cut → free-name → yl_form."""

    name = "recursive"

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        hit = name_as_substituent(mol, claim.root, claim.atoms, depth=depth)
        return None if hit is None else _from_yl(claim, hit)


def _from_yl(claim: ClaimedBlock, hit: tuple[str, str, bool]) -> SubstituentName | None:
    en, zh, paren = hit
    if not en or not zh:
        return None
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="recursive",
    )


def _default_backends() -> list[SubstituentBackend]:
    return [RetainedBackend(), RootedTreeBackend(), RecursiveBackend()]


class SubstituentNamer:
    def __init__(self, backends: Sequence[SubstituentBackend] | None = None) -> None:
        self._backends = list(backends) if backends is not None else _default_backends()

    def name(self, mol, claim: ClaimedBlock, *, depth: int = 0) -> SubstituentName | None:
        for backend in self._backends:
            hit = backend.try_name(mol, claim, depth=depth)
            if hit is not None:
                return hit
        return None
