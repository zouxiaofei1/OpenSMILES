"""L3 dual systematic -yl names for rooted saturated-carbon trees."""
from __future__ import annotations

from namepredict.layer2.side_alkyl_sys import RootedAlkylTree

_CHAIN_EN = {
    1: "meth", 2: "eth", 3: "prop", 4: "but", 5: "pent", 6: "hex",
    7: "hept", 8: "oct", 9: "non", 10: "dec", 11: "undec", 12: "dodec",
}
_CHAIN_ZH = {
    1: "甲", 2: "乙", 3: "丙", 4: "丁", 5: "戊", 6: "己",
    7: "庚", 8: "辛", 9: "壬", 10: "癸", 11: "十一烷", 12: "十二烷",
}
_BRANCH_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_BRANCH_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
from namepredict.constants import MULT_EN as _MULT_EN, MULT_ZH as _MULT_ZH


def _paths_from(tree: RootedAlkylTree, node: int) -> list[list[int]]:
    kids = tree.children.get(node) or ()
    if not kids:
        return [[node]]
    out: list[list[int]] = []
    for k in kids:
        for p in _paths_from(tree, k):
            out.append([node] + p)
    return out


def _principal_path(tree: RootedAlkylTree) -> list[int]:
    paths = _paths_from(tree, tree.root)

    def key(p: list[int]) -> tuple:
        path_set = set(p)
        n_branch = len(tree.atoms - path_set)
        locs = _raw_branch_locants(tree, p)
        return (len(p), n_branch, tuple(-x for x in locs))

    return max(paths, key=key)


def _raw_branch_locants(tree: RootedAlkylTree, path: list[int]) -> list[int]:
    path_set = set(path)
    pos = {a: i + 1 for i, a in enumerate(path)}
    locs: list[int] = []
    for a in path:
        for c in tree.children.get(a) or ():
            if c not in path_set:
                locs.append(pos[a])
    return sorted(locs)


def _branch_size(tree: RootedAlkylTree, start: int) -> int:
    n, stack = 0, [start]
    while stack:
        cur = stack.pop()
        n += 1
        stack.extend(tree.children.get(cur) or ())
    return n


def _is_linear_branch(tree: RootedAlkylTree, start: int) -> bool:
    cur = start
    while True:
        kids = tree.children.get(cur) or ()
        if len(kids) > 1:
            return False
        if not kids:
            return True
        cur = kids[0]


def _collect_branches(
    tree: RootedAlkylTree, path: list[int],
) -> list[tuple[int, int]] | None:
    """(locant, n_carbons) list, or None if a branch is unsupported."""
    path_set = set(path)
    pos = {a: i + 1 for i, a in enumerate(path)}
    branches: list[tuple[int, int]] = []
    for a in path:
        for c in tree.children.get(a) or ():
            if c in path_set:
                continue
            size = _branch_size(tree, c)
            if size not in _BRANCH_EN or not _is_linear_branch(tree, c):
                return None
            branches.append((pos[a], size))
    return branches


def _format_en(chain_n: int, branches: list[tuple[int, int]]) -> str:
    stem = f"{_CHAIN_EN[chain_n]}yl"
    if not branches:
        return stem
    return f"{_prefix_en(branches)}{stem}"


def _format_zh(chain_n: int, branches: list[tuple[int, int]]) -> str:
    stem = f"{_CHAIN_ZH[chain_n]}基"
    if not branches:
        return stem
    return f"{_prefix_zh(branches)}{stem}"


def _prefix_en(branches: list[tuple[int, int]]) -> str:
    by_size: dict[int, list[int]] = {}
    for loc, size in branches:
        by_size.setdefault(size, []).append(loc)
    parts: list[str] = []
    for size in sorted(by_size):
        locs = sorted(by_size[size])
        name = _BRANCH_EN[size]
        if len(locs) == 1:
            parts.append(f"{locs[0]}-{name}")
        else:
            mult = _MULT_EN.get(len(locs), f"{len(locs)}")
            parts.append(f"{','.join(map(str, locs))}-{mult}{name}")
    return "-".join(parts)


def _prefix_zh(branches: list[tuple[int, int]]) -> str:
    by_size: dict[int, list[int]] = {}
    for loc, size in branches:
        by_size.setdefault(size, []).append(loc)
    parts: list[str] = []
    for size in sorted(by_size):
        locs = sorted(by_size[size])
        name = _BRANCH_ZH[size]
        if len(locs) == 1:
            parts.append(f"{locs[0]}-{name}")
        else:
            mult = _MULT_ZH.get(len(locs), f"{len(locs)}")
            parts.append(f"{','.join(map(str, locs))}-{mult}{name}")
    return "-".join(parts)


def name_rooted_alkyl(tree: RootedAlkylTree) -> tuple[str, str, bool] | None:
    """Return (en, zh, requires_parentheses) systematic -yl dual form."""
    path = _principal_path(tree)
    chain_n = len(path)
    if chain_n not in _CHAIN_EN:
        return None
    branches = _collect_branches(tree, path)
    if branches is None:
        return None
    if len(tree.atoms) != chain_n + sum(s for _, s in branches):
        return None
    return _format_en(chain_n, branches), _format_zh(chain_n, branches), bool(branches)
