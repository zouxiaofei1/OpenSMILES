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


class RetainedBackend:
    """True simple leaves; Task 4 default is no-op (order tests use fakes)."""

    name = "retained"

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return None


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
