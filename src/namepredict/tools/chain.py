"""链游走工具：同元素开链邻居枚举与最长开链路径搜索。"""

from __future__ import annotations

from rdkit.Chem import Mol


def _element_neighbors(mol: Mol, idx: int, z: int = 6, banned: set[int] = frozenset()) -> list[int]:
    """链游走用的开链同元素邻居（z = 原子序数），排除 banned 禁走原子。"""
    atom = mol.GetAtomWithIdx(idx)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == z and not n.GetIsAromatic() and not n.IsInRing()
        and n.GetIdx() not in banned
    ]


def _carbon_neighbors(mol: Mol, idx: int, banned: set[int] = frozenset()) -> list[int]:
    """链游走用的开链碳邻居，排除 banned 禁走碳。"""
    return _element_neighbors(mol, idx, 6, banned)


def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int],
              banned: set[int] = frozenset(), z: int = 6) -> list[int]:
    """深度优先寻找从 node 出发的最长开链路径（限同元素 z）。

    path 为进入时的前缀（末位即 node）。改用显式栈迭代：长链不再消耗
    Python 递归深度，也不再逐帧复制路径（原为 O(L²)）。邻居访问顺序与
    「先到的最长优先」判据同递归版，返回值逐一对应。
    """
    prefix = list(path)
    if not prefix or prefix[-1] != node:  # 前缀未含 node 时补上
        prefix.append(node)
    prefix_set = set(prefix)
    on_path = set(prefix)
    blocked = set(forbid)
    nbrs: dict[int, list[int]] = {}
    parent: dict[int, int] = {}
    depth0 = len(prefix)
    # 帧 = [节点, 邻居游标, 子树最优终点, 最优长度, 本节点深度]
    stack: list[list] = [[node, 0, node, depth0, depth0]]
    while stack:
        frame = stack[-1]
        u, i, end, best, depth = frame
        ns = nbrs.get(u)
        if ns is None:  # 每个节点的邻居表只取一次
            ns = nbrs[u] = _element_neighbors(mol, u, z, banned)
        if i < len(ns):
            frame[1] = i + 1
            nb = ns[i]
            if nb in on_path or nb in blocked:
                continue
            on_path.add(nb)
            parent[nb] = u
            stack.append([nb, 0, nb, depth + 1, depth + 1])
            continue
        stack.pop()
        if u not in prefix_set:  # 前缀内节点始终算已走过
            on_path.discard(u)
        if stack and best > stack[-1][3]:  # 严格大于：与递归版同序取先到者
            stack[-1][2], stack[-1][3] = end, best
    chain: list[int] = []
    while end != node:  # 沿父指针回溯出最优路径尾部
        chain.append(end)
        end = parent[end]
    chain.reverse()
    return prefix + chain


def _longest_from(mol: Mol, start: int, banned: set[int] = frozenset(), z: int = 6) -> list[int]:
    """从 start 出发避开 banned 的最长开链路径（限同元素 z）。"""
    return _dfs_path(mol, start, [start], set(), banned, z)
