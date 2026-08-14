"""带类型的 L3 取代基命名器：ordered retained → rooted_tree → recursive。"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.claimable_block import ClaimedBlock
from namepredict.layer3.as_substituent import name_as_substituent
from namepredict.tools.anchored_table import anchored_lookup


@dataclass(frozen=True)
class SubstituentName:
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained"、"rooted_tree" 或 "recursive"


class SubstituentBackend(Protocol):
    name: str

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None: ...


def _retained_hit(claim: ClaimedBlock, en: str, zh: str, paren: bool) -> SubstituentName:
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="retained",
    )


def _try_anchored_lookup(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """通过锚定 canonical-SMILES 表解析简单的 retained 叶子。

    锚定键唯一标识形状 + 连接位点，因此字典查找即可覆盖 alkyl/cycloalkyl/
    aryl/halo/alkoxy/sulfinyl/sulfonyl/N-leaves。未命中则落入 rooted_tree/
    recursive 链。
    """
    hit = anchored_lookup(mol, claim.atoms, claim.root, name_mode=name_mode)
    if hit is None:
        return None
    en, zh, paren = hit
    return _retained_hit(claim, en, zh, paren)


def _retained_name(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """简单 retained 叶子：alkyl/cycloalkyl/phenyl/halo/alkoxy/sulfinyl/
    sulfonyl/N-leaves，全部通过锚定 canonical-SMILES 表解析。"""
    return _try_anchored_lookup(mol, claim, name_mode=name_mode)


class RetainedBackend:
    """锚定表 retained 叶子：alkyl/aryl/halo/alkoxy/sulfinyl/..."""

    name = "retained"

    def __init__(self, *, name_mode: str = "general") -> None:
        self._name_mode = name_mode

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return _retained_name(mol, claim, name_mode=self._name_mode)

class RecursiveBackend:
    """有界递归 cut → free-name → yl_form。"""

    name = "recursive"

    def __init__(self, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> None:
        self._name_mode = name_mode
        self._cache = cache

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        hit = name_as_substituent(mol, claim.root, claim.atoms, depth=depth, name_mode=self._name_mode, cache=self._cache)
        return None if hit is None else _from_yl(claim, hit)


def _from_yl(claim: ClaimedBlock, hit: tuple[str, str, bool]) -> SubstituentName | None:
    en, zh, paren = hit
    if not en or not zh:
        return None
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="recursive",
    )


def _default_backends(name_mode: str = "general", cache: CommonNameCache | None = None) -> list[SubstituentBackend]:
    return [RetainedBackend(name_mode=name_mode),  RecursiveBackend(name_mode=name_mode, cache=cache)]


class SubstituentNamer:
    def __init__(self, backends: Sequence[SubstituentBackend] | None = None, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> None:
        self._backends = list(backends) if backends is not None else _default_backends(name_mode, cache=cache)

    def name(self, mol, claim: ClaimedBlock, *, depth: int = 0) -> SubstituentName | None:
        for backend in self._backends:
            hit = backend.try_name(mol, claim, depth=depth)
            if hit is not None:
                return hit
        return None
