"""排序、终态化并返回母体候选（P-44 选择入口）。"""
from __future__ import annotations


def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    """按评分降序排列候选母体。"""
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)


def _p45_2_prefix_count(info: dict, parent: dict) -> int:
    """P-45.2.1 键：以前缀引用的取代基团数目。

    母体 arm 选择在骨架/评分打平（如 N 上两条等长臂）时，IUPAC P-45.2.1
    要求优选以前缀引用的取代基团数目最多的母体结构。数目等于母体
    owned_atoms 边界之外的取代基组分（claim）个数，由 L3 claimable_block
    的 iter_claims 直接枚举（无需递归命名，无 name_mode/cache 依赖）。
    """
    mol = info["mol"]
    from namepredict.layer3.claimable_block import iter_claims
    return len(iter_claims(mol, parent.get("owned_atoms") or frozenset()))


def _reorder_p45_2(info: dict, cands: list[dict]) -> list[dict]:
    """P-45.2 流水线：按前缀取代基团数目最多（P-45.2.1）重排候选。

    稳定排序：同数目候选保持既有评分顺序；仅当候选数同、评分打平
    时（P-45.2.1 无法分出的骨架余隙）由 P-45.2.2/P-45.2.3 兜底——
    此处 P-45.2.2/2.3 的位次集合需 L4 编号后才可得，故只先接 P-45.2.1，
    完全打平的候选维持其余额顺序交由下游评分/覆盖门控处理。
    """
    if len(cands) <= 1:
        return cands
    return sorted(cands, key=lambda c: _p45_2_prefix_count(info, c), reverse=True)


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

def select_parent(info: dict, *, all_candidates: bool = False) -> dict | list[dict] | None:
    """排序后的母体候选，每个都以不可变 owned_atoms 完成最终确定。

    先按 P-44 简化评分降序，再经 P-45.2 流水线按前缀取代基团数目（P-45.2.1）
    重排打平候选，取首位为选定母体。
    """
    from namepredict.layer2.candidates import _collect_candidates

    cands = _finalize_ranked(info, _collect_candidates(info))
    cands = _reorder_p45_2(info, cands)
    if all_candidates:
        return cands
    return next(iter(cands), None)
