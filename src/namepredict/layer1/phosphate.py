"""L1 磷酸官能团（kind=phosphate）检测：分子中每个 P(=O)(O)₃ 中心产出一条条目。
"""
from __future__ import annotations

from collections import deque

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, O, P


def _heavy(mol: Mol, a) -> list:
    """非氢邻居索引。"""
    return [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]


def _arm_component(mol: Mol, start: int, core: set[int]) -> set[int]:
    """从臂根 start 取不穿过 core 的连通重原子组分。"""
    comp: set[int] = {start}
    dq: deque[int] = deque([start])
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
    return comp


def _one_phosphate(mol: Mol, p_idx: int) -> dict | None:
    """判定单个 P 是否为磷酸中心，是则返回 {p_idx,n_oh,n_om,n_arms}，否则 None。"""
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
    arm_all: set[int] = set()
    for o_idx in single_o:
        a = mol.GetAtomWithIdx(o_idx)
        heavy = _heavy(mol, a)
        if a.GetFormalCharge() == 0 and a.GetTotalNumHs() >= 1 and heavy == [p_idx]:
            n_oh += 1
            continue
        if a.GetFormalCharge() == -1 and a.GetTotalNumHs() == 0 and heavy == [p_idx]:
            n_om += 1
            continue
        if a.GetFormalCharge() != 0 or a.GetTotalNumHs() != 0 or p_idx not in heavy:  # O–R：中性无 H、除 P 外另连 1 个重原子
            return None
        others = [j for j in heavy if j != p_idx]
        if len(others) != 1 or mol.GetAtomWithIdx(others[0]).GetAtomicNum() != C:
            return None
        comp = _arm_component(mol, others[0], core)
        attaches = [j for i in comp for j in _heavy(mol, mol.GetAtomWithIdx(i)) if j in core]  # 组分只贴 1 个 core 原子（桥 O）；不能连到 P 或其它 O
        if not attaches or len(set(attaches)) != 1 or attaches[0] != o_idx:
            return None
        arm_all |= comp
        n_arms += 1
    if n_oh + n_om + n_arms != 3:
        return None
    all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}  # 整分子纯度：重原子 = core ∪ 臂（排除臂间成环、P–O–P 焦磷酸等）
    if all_heavy != (core | arm_all):
        return None
    return {"p_idx": p_idx, "n_oh": n_oh, "n_om": n_om, "n_arms": n_arms}


def phosphate_entries(mol: Mol) -> list[dict]:
    """分子中全部磷酸中心（P(=O)(O)₃）的条目列表。"""
    result = [e for a in mol.GetAtoms() if a.GetAtomicNum() == P  and (e := _one_phosphate(mol, a.GetIdx())) is not None]
    # print(result)
    return result
