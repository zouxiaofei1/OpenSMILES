"""L3 取代基命名器：ordered retained → recursive。"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence

from namepredict.tools.common_names import CommonNameCache
from namepredict.layer3.claimable_block import ClaimedBlock
from namepredict.layer3.as_substituent import name_as_substituent
from namepredict.tools.anchored_table import anchored_lookup


@dataclass(frozen=True)
class SubstituentName:
    """取代基命名结果：claim、双语名与是否需括号。"""
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool


class SubstituentBackend(Protocol):
    """取代基命名后端协议：各后端实现 try_name 尝试为 claim 命名。"""
    name: str




def _retained_hit(claim: ClaimedBlock, en: str, zh: str, paren: bool) -> SubstituentName:
    """封装 retained 叶子命中结果为 SubstituentName。"""
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren,
    )


def _try_anchored_lookup(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """尝试通过锚定 canonical-SMILES 表解析 retained 叶子。"""
    hit = anchored_lookup(mol, claim.atoms, claim.root)
    if hit is None:
        return None
    en, zh, paren = hit
    return _retained_hit(claim, en, zh, paren)


def _retained_name(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """retained 叶子（alkyl/aryl/halo 等）经锚定表解析。"""
    return _try_anchored_lookup(mol, claim)


class RetainedBackend:
    """锚定表 retained 叶子：alkyl/aryl/halo 等。"""

    name = "retained"

    def try_name(self, mol, claim: ClaimedBlock) -> SubstituentName | None:
        """通过锚定表 retained 叶子尝试命名。"""
        return _retained_name(mol, claim)

class RecursiveBackend:
    """有界递归 cut → free-name → yl_form。"""

    name = "recursive"

    def __init__(self, *, cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> None:
        """记录共享缓存与根分子上下文。"""
        self._cache = cache
        self._root_ctx = root_ctx

    def try_name(self, mol, claim: ClaimedBlock) -> SubstituentName | None:
        """通过递归 cut → free-name → yl_form 尝试命名。"""
        hit = name_as_substituent(mol, claim.root, claim.atoms, cache=self._cache, root_ctx=self._root_ctx)
        return None if hit is None else _from_yl(claim, hit)


def _from_yl(claim: ClaimedBlock, hit: tuple[str, str, bool]) -> SubstituentName | None:
    """将 -yl 双语结果封装为递归后端的 SubstituentName。"""
    en, zh, paren = hit
    if not en or not zh:
        return None
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren,
    )


def _default_backends(cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> list[SubstituentBackend]:
    """构造默认命名后端列表（retained → recursive）。"""
    return [RetainedBackend(),  RecursiveBackend(cache=cache, root_ctx=root_ctx)]


class SubstituentNamer:
    """按序尝试各后端为 claim 命名，返回首个命中的命名器。"""

    def __init__(self, backends: Sequence[SubstituentBackend] | None = None, *, cache: CommonNameCache | None = None, root_ctx: tuple | None = None) -> None:
        """初始化后端列表，缺省用默认后端；root_ctx 供 R/S 校正。"""
        self._backends = list(backends) if backends is not None else _default_backends(cache=cache, root_ctx=root_ctx)

    def name(self, mol, claim: ClaimedBlock) -> SubstituentName | None:
        """按序尝试各后端命名 claim，返回首个命中。"""
        for backend in self._backends:
            hit = backend.try_name(mol, claim)
            if hit is not None:
                return hit
        return None
