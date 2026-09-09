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
    """取代基命名结果：claim、双语名、是否需括号与后端来源。"""
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained"、"rooted_tree" 或 "recursive"


class SubstituentBackend(Protocol):
    """取代基命名后端协议：各后端实现 try_name 尝试为 claim 命名。"""
    name: str

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        """各后端尝试为 claim 命名，失败返回 None。"""
        ...


def _retained_hit(claim: ClaimedBlock, en: str, zh: str, paren: bool) -> SubstituentName:
    """以 retained 后端类型封装 SubstituentName。"""
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="retained",
    )


def _try_anchored_lookup(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """尝试通过锚定 canonical-SMILES 表解析 retained 叶子。"""
    hit = anchored_lookup(mol, claim.atoms, claim.root, name_mode=name_mode)
    if hit is None:
        return None
    en, zh, paren = hit
    return _retained_hit(claim, en, zh, paren)


def _retained_name(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """简单 retained 叶子（alkyl/cycloalkyl/phenyl/halo/alkoxy/sulfinyl/sulfonyl/N-leaves），经锚定 canonical-SMILES 表解析。"""
    return _try_anchored_lookup(mol, claim, name_mode=name_mode)


class RetainedBackend:
    """锚定表 retained 叶子：alkyl/aryl/halo/alkoxy/sulfinyl/..."""

    name = "retained"

    def __init__(self, *, name_mode: str = "general") -> None:
        """记录命名模式。"""
        self._name_mode = name_mode

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        """通过锚定表 retained 叶子尝试命名。"""
        return _retained_name(mol, claim, name_mode=self._name_mode)

class RecursiveBackend:
    """有界递归 cut → free-name → yl_form。"""

    name = "recursive"

    def __init__(self, *, name_mode: str = "general", cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> None:
        """记录命名模式、共享缓存与根分子上下文。"""
        self._name_mode = name_mode
        self._cache = cache
        self._root_ctx = root_ctx

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        """通过递归 cut → free-name → yl_form 尝试命名。"""
        hit = name_as_substituent(mol, claim.root, claim.atoms, depth=depth, name_mode=self._name_mode, cache=self._cache, root_ctx=self._root_ctx)
        return None if hit is None else _from_yl(claim, hit)


def _from_yl(claim: ClaimedBlock, hit: tuple[str, str, bool]) -> SubstituentName | None:
    """将 -yl 双语结果封装为递归后端的 SubstituentName。"""
    en, zh, paren = hit
    if not en or not zh:
        return None
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="recursive",
    )


def _default_backends(name_mode: str = "general", cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> list[SubstituentBackend]:
    """构造默认命名后端列表（retained → recursive）。"""
    return [RetainedBackend(name_mode=name_mode),  RecursiveBackend(name_mode=name_mode, cache=cache, root_ctx=root_ctx)]


class SubstituentNamer:
    """按序尝试各后端为 claim 命名，返回首个命中的命名器。"""

    def __init__(self, backends: Sequence[SubstituentBackend] | None = None, *, name_mode: str = "general", cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> None:
        """初始化后端列表，缺省时用默认后端；root_ctx 供递归取代基命名回根分子校正 R/S。"""
        self._backends = list(backends) if backends is not None else _default_backends(name_mode, cache=cache, root_ctx=root_ctx)

    def name(self, mol, claim: ClaimedBlock, *, depth: int = 0) -> SubstituentName | None:
        """按序尝试各后端命名 claim，返回首个命中。"""
        for backend in self._backends:
            hit = backend.try_name(mol, claim, depth=depth)
            if hit is not None:
                return hit
        return None
