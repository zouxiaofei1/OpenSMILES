"""name_as_substituent: cut submol → free-name pipeline → P-29 -yl form.

No host (benzamide gate / n_block extract) wiring — pure cut→pipeline→yl.
"""
from __future__ import annotations

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.submol_build import build_cut_submol
from namepredict.layer3.yl_form import yl_form


def _locant_from_result(result, attach_new: int) -> int | None:
    chain = (result.meta or {}).get("parent_chain") or []
    if attach_new in chain:
        return chain.index(attach_new) + 1
    return None


def _locant_via_hetero(mol, attach_new: int, chain: list[int]) -> int | None:
    """When attach is a heteroatom (O/S/N) not in chain, use its chain-neighbor."""
    a = mol.GetAtomWithIdx(attach_new)
    if a.GetAtomicNum() not in (8, 16, 7):
        return None
    c_nbrs = [n.GetIdx() for n in a.GetNeighbors() if n.GetIdx() in chain]
    return chain.index(c_nbrs[0]) + 1 if len(c_nbrs) == 1 else None


def _yl_from_sub(
    sub, *, depth: int, name_mode: str = "general", cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    from namepredict.namer import _cache_put, _name_mol
    from rdkit import Chem

    result = None
    smiles = None
    if cache is not None:
        # 递归子结构命名与主分子共享缓存：命中则跳过整条 L2 递归。
        smiles = Chem.MolToSmiles(sub.mol)
        result = cache.get(smiles)
    if result is None:
        result = _name_mol(sub.mol, depth=depth, name_mode=name_mode, cache=cache)
        if cache is not None and result.success and result.en:
            _cache_put(cache, smiles, result)
    if not result.success or not result.en:
        return None
    loc = _locant_from_result(result, sub.attach_new)
    if loc is None:
        loc = _locant_via_hetero(sub.mol, sub.attach_new, (result.meta or {}).get("parent_chain") or [])
    return None if loc is None else yl_form(result.en, result.zh, loc)


def name_as_substituent(
    mol, attach_old: int, atoms, *, depth: int = 0, max_depth: int = 4, name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    """Cut atoms at attach_old, free-name the submol, emit -yl dual names."""
    if depth >= max_depth:
        return None
    sub = build_cut_submol(mol, frozenset(atoms), attach_old)
    return None if sub is None else _yl_from_sub(sub, depth=depth + 1, name_mode=name_mode, cache=cache)
