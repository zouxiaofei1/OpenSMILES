"""L1 整分子磷酸母体族（kind=phosphate）谓词。
"""
from __future__ import annotations

from collections import deque

from rdkit.Chem import BondType, Mol

from namepredict.constants import (
    Br, C, Cl, F, I, N, O, P, S,
)

def _heavy(mol: Mol, a) -> list:
    """非氢邻居索引。"""
    return [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]

def detect_phosphate_whole(mol: Mol) -> dict | None:
    """检测磷酸整分子母体族。返回 {n_oh, n_om, n_arms, p_idx}；否则 None。"""
    p_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == P]
    if len(p_idxs) != 1:
        return None
    p_idx = p_idxs[0]
    p = mol.GetAtomWithIdx(p_idx)
    if p.GetTotalNumHs() != 0 or p.GetFormalCharge() != 0:
        return None
    nei = _heavy(mol, p)
    if len(nei) != 4 or any(mol.GetAtomWithIdx(i).GetAtomicNum() != O for i in nei):
        return None
    single_o: list[int] = []
    dbl_o: list[int] = []
    for i in nei:
        bt = mol.GetBondBetweenAtoms(p_idx, i).GetBondType()
        if bt == BondType.DOUBLE:
            dbl_o.append(i)
        elif bt == BondType.SINGLE:
            single_o.append(i)
        else:
            return None
    if len(dbl_o) != 1 or len(single_o) != 3:
        return None
    da = mol.GetAtomWithIdx(dbl_o[0])
    if da.GetFormalCharge() != 0 or da.GetTotalNumHs() != 0 or _heavy(mol, da) != [p_idx]:
        return None
    core = {p_idx, *nei}
    n_oh = n_om = n_arms = 0
    for o_idx in single_o:
        a = mol.GetAtomWithIdx(o_idx)
        heavy = _heavy(mol, a)
        if a.GetFormalCharge() == 0 and a.GetTotalNumHs() >= 1 and heavy == [p_idx]:
            n_oh += 1
            continue
        if a.GetFormalCharge() == -1 and a.GetTotalNumHs() == 0 and heavy == [p_idx]:
            n_om += 1
            continue
        # O–R：中性无 H、除 P 外另连 1 个重原子
        if a.GetFormalCharge() != 0 or a.GetTotalNumHs() != 0 or p_idx not in heavy:
            return None
        others = [j for j in heavy if j != p_idx]
        if len(others) != 1:
            return None
        start = others[0]
        # 从 start 取不穿过 core 的连通组分
        comp: set[int] = set()
        dq: deque[int] = deque([start])
        comp.add(start)
        while dq:
            i = dq.popleft()
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                if nb.GetAtomicNum() == 1:
                    continue
                j = nb.GetIdx()
                if j in core or j in comp:
                    continue
                comp.add(j)
                dq.append(j)
        if not comp or mol.GetAtomWithIdx(start).GetAtomicNum() != C:
            return None
        # 组分只贴 1 个 core 原子（桥 O）；不能连到 P 或其它 O
        attaches = [j for i in comp for j in _heavy(mol, mol.GetAtomWithIdx(i)) if j in core]
        if not attaches or len(set(attaches)) != 1 or attaches[0] != o_idx:
            return None
        n_arms += 1
    if n_oh + n_om + n_arms != 3:
        return None
    # 整分子纯度：重原子 = core ∪ 臂
    all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}
    arm_all = set()
    for o_idx in single_o:
        a = mol.GetAtomWithIdx(o_idx)
        heavy = _heavy(mol, a)
        if a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 0 and len(heavy) == 2:
            start = [j for j in heavy if j != p_idx][0]
            dq: deque[int] = deque([start])
            seen: set[int] = {start}
            while dq:
                i = dq.popleft()
                arm_all.add(i)
                for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                    if nb.GetAtomicNum() == 1 or nb.GetIdx() in core or nb.GetIdx() in seen:
                        continue
                    seen.add(nb.GetIdx())
                    dq.append(nb.GetIdx())
    if all_heavy != (core | arm_all):
        return None
    return {"n_oh": n_oh, "n_om": n_om, "n_arms": n_arms, "p_idx": p_idx}
