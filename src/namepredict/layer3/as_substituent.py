"""name_as_substituent: cut submol → free-name pipeline → P-29 -yl form.

No host (benzamide gate / n_block extract) wiring — pure cut→pipeline→yl.
"""
from __future__ import annotations

from rdkit.Chem import CanonicalRankAtoms

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.submol_build import build_cut_submol
from namepredict.layer5.free_to_yl import free_to_yl as yl_form


def _arene_yl_from_sub(sub, result) -> tuple[str, str, bool] | None:
    """Structured arene → yl for benzene-parent submol (P-29.6.2 phenyl).

    The free-name pipeline names a substituted benzene as e.g. 'chlorobenzene'
    (or 'phenol'/'aniline'/'benzonitrile' when the ring carries an OH/NH2/CN
    leaf — those are leaves of the phenyl substituent, not the free parent).
    The generic -k-yl fallback would produce 'chlorobenzen-1-yl'.  For any
    benzene ring we instead renumber with attach=1 (lowest locant set) and
    swap the stem to 'phenyl', letting name_ph_ring collect leaves on the cut
    submol (halo/Me/alkoxy/nitro/OH/NH2/CN/nested Ph all verified).
    """
    from namepredict.layer3.ring_namer import name_ph_ring

    mol = sub.mol
    attach_new = sub.attach_new
    ph = next(
        (set(r) for r in mol.GetRingInfo().AtomRings()
         if attach_new in r and len(r) == 6
         and all(mol.GetAtomWithIdx(i).GetIsAromatic()
                 and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in r)),
        None,
    )
    if ph is None:
        return None
    en, zh, _ = name_ph_ring(mol, ph, attach_new, -1, depth=1)
    if not en:
        return None
    return en, zh, en != "phenyl"


def _locant_from_result(mol, result, attach_new: int) -> int | None:
    """Locant of attach on the parent_chain (canonical ranks; cache-safe across cuts)."""
    chain = (result.meta or {}).get("parent_chain") or []
    if not chain:
        return None
    ranks = CanonicalRankAtoms(mol)
    ar = ranks[attach_new]
    if ar in chain:
        return chain.index(ar) + 1
    return None


def _locant_via_hetero(mol, attach_new: int, result) -> int | None:
    """When attach is a heteroatom (O/S/N) not in chain, use its chain-neighbor."""
    a = mol.GetAtomWithIdx(attach_new)
    if a.GetAtomicNum() not in (8, 16, 7):
        return None
    chain = (result.meta or {}).get("parent_chain") or []
    if not chain:
        return None
    ranks = CanonicalRankAtoms(mol)
    c_nbrs = [ranks[n.GetIdx()] for n in a.GetNeighbors() if ranks[n.GetIdx()] in chain]
    return chain.index(c_nbrs[0]) + 1 if len(c_nbrs) == 1 else None


def _yl_from_sub(
    sub, *, depth: int, name_mode: str = "general", cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    from namepredict.namer import _cache_put, _canonical_result, _name_mol
    from rdkit import Chem

    result = None
    smiles = None
    if cache is not None:
        # 递归子结构命名与主分子共享缓存：命中则跳过整条 L2 递归。
        smiles = Chem.MolToSmiles(sub.mol)
        result = cache.get(smiles)
    if result is None:
        # 用 canonical SMILES 重解析再命名：命名输入与缓存 key 严格一一对应，
        # 消除不同 cut 上下文（环断点/手性方向）泄漏进命名的差异。
        named_mol = sub.mol
        if smiles:
            canonical = Chem.MolFromSmiles(smiles)
            if canonical is not None:
                named_mol = canonical
        result = _name_mol(named_mol, depth=depth, name_mode=name_mode, cache=cache)
        if cache is not None and result.success and result.en:
            # 缓存条目的 parent_chain 必须是 canonical rank，跨 cut 复用才安全。
            result = _canonical_result(named_mol, result)
            _cache_put(cache, smiles, result)
    if not result.success or not result.en:
        return None
    arene = _arene_yl_from_sub(sub, result)
    if arene is not None:
        return arene
    loc = _locant_from_result(sub.mol, result, sub.attach_new)
    if loc is None:
        loc = _locant_via_hetero(sub.mol, sub.attach_new, result)
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
