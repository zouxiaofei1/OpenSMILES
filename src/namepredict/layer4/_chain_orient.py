"""L4 链位次与键端点位次的辅助工具函数。"""
from __future__ import annotations

from namepredict.tools.re import alkyl_alpha_key


def _edge_locants(chain: list[int], edge) -> tuple[int, int] | None:
    """返回边两端原子在链中的位次对；端点不在链中则 None。"""
    if not edge or edge[0] not in chain or edge[1] not in chain:
        return None
    return chain.index(edge[0]) + 1, chain.index(edge[1]) + 1

def _edge_min_locant(chain: list[int], edge) -> int | None:
    """返回边两端位次中较小的一个；无位次则 None。"""
    ends = _edge_locants(chain, edge)
    return min(ends) if ends else None

def _bond_min_locs(chain: list[int], bonds) -> tuple[int, ...] | None:
    """返回全部键较小端点位次的排序元组；有键无位次则 None。"""
    if not bonds:
        return None
    locs = [_edge_min_locant(chain, b) for b in bonds]
    if any(x is None for x in locs):
        return None
    return tuple(sorted(int(x) for x in locs))

def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    """返回 (基团字母序键, 位次) 排序对，用于字母序平局。"""
    return sorted(
        (alkyl_alpha_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents if s["attach_idx"] in chain
    )
