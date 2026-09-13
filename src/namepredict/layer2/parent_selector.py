"""排序、终态化并返回母体候选（P-44 选择入口）。"""
from __future__ import annotations


def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    """按评分降序排列候选母体。"""
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)


def _p45_2_prefix_count(info: dict, parent: dict) -> int:
    """P-45.2.1 键：以前缀引用的取代基团数目（=owned_atoms 边界外的 claim 个数，由 L3 iter_claims 枚举，无 cache 依赖）。"""
    mol = info["mol"]
    from namepredict.layer3.claimable_block import iter_claims
    return len(iter_claims(mol, parent.get("owned_atoms") or frozenset()))


def _reorder_p45_2(info: dict, cands: list[dict], *, tied: bool = False) -> list[dict]:
    """P-45.2 流水线：按前缀取代基团数目最多（P-45.2.1）稳定重排打平候选；tied 时只返回并列最大组。
    P-45.2.2/2.3 的位次需 L4 编号后才可得，只先接 P-45.2.1。"""
    if len(cands) <= 1:
        return cands
    keyed = sorted(((_p45_2_prefix_count(info, c), i, c) for i, c in enumerate(cands)),
                   key=lambda t: (-t[0], t[1]))
    if not tied:
        return [c for _, _, c in keyed]
    top = keyed[0][0]
    return [c for k, _, c in keyed if k == top]


def _finalize_ranked(info: dict, cands: list[dict]) -> list[dict]:
    """补齐契约/词干/编号并固化 owned_atoms。"""
    from namepredict.layer2.kind_registry import pack_parent_stem
    from namepredict.layer2.parent_candidate import with_principal_group_contract
    from namepredict.layer2.parent_ownership import finalize_parent_ownership
    mol = info.get("mol")
    return [
        finalize_parent_ownership(
            pack_parent_stem(with_principal_group_contract(c), mol), mol,
        )
        for c in _rank_candidates(info, cands)
    ]

def select_parent(info: dict) -> list[dict]:
    """P-44 评分降序->P-45.2-> P-45.2.1 并列最优"""
    from namepredict.layer2.candidates import _collect_candidates

    cands = _finalize_ranked(info, _collect_candidates(info))
    return _reorder_p45_2(info, cands, tied=True)
