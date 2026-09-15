"""L1 官能团分析器：枚举分子中各类官能团条目并汇总为分析结果 dict。
酰卤 R–C(=O)–X（X=F/Cl/Br/I，P-65.5）检测内联于本模块。
"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol
from collections import deque
from namepredict.constants import (
    C, H,
)
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer1.fg_local_smarts import match_local_fg

def _heavy(a) -> list:
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


def _phosphate_entry(mol: Mol, core: tuple[int, ...]) -> dict | None:
    """校验磷酸候选的非局部条件并组装条目（臂单点回接 + 整分子纯度）。"""
    p_idx, _o_dbl, *o_sgl = core
    core_set = set(core)
    n_oh = n_om = n_arms = 0
    arm_all: set[int] = set()
    for o_idx in o_sgl:
        o = mol.GetAtomWithIdx(o_idx)
        if o.GetFormalCharge() == 0 and o.GetTotalNumHs() >= 1:  # 羟基氧（SMARTS 已保证只连 P）
            n_oh += 1
            continue
        if o.GetFormalCharge() == -1:  # 羧酸盐式阴离子氧
            n_om += 1
            continue
        others = [j for j in _heavy(o) if j != p_idx]  # O–R：臂根须为碳（SMARTS 已保证）
        if len(others) != 1:
            return None
        comp = _arm_component(mol, others[0], core_set)
        attaches = [j for i in comp for j in _heavy(mol.GetAtomWithIdx(i)) if j in core_set]  # 组分只贴 1 个 core 原子（桥 O）；不能连到 P 或其它 O
        if not attaches or len(set(attaches)) != 1 or attaches[0] != o_idx:
            return None
        arm_all |= comp
        n_arms += 1
    all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}  # 整分子纯度：重原子 = core ∪ 臂（排除臂间成环、P–O–P 焦磷酸等）
    if all_heavy != (core_set | arm_all):
        return None
    return {"p_idx": p_idx, "n_oh": n_oh, "n_om": n_om, "n_arms": n_arms}


def phosphate_entries(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> list[dict]:
    """分子中全部磷酸中心（P(=O)(O)₃）的条目列表。"""
    if matches is None:
        matches = match_local_fg(mol).get("phosphate", [])
    out = [e for m in matches if (e := _phosphate_entry(mol, m)) is not None]
    return sorted(out, key=lambda e: e["p_idx"])
def _surr_idx(atom) -> list[int]:
    """周边原子：中心全部重原子邻居；碳中心不在环内时排除环内邻居。"""
    ring_excl = atom.GetAtomicNum() == C and not atom.IsInRing()
    return [n.GetIdx() for n in atom.GetNeighbors()
            if n.GetAtomicNum() != H and not (ring_excl and n.IsInRing())]

def _fg_entry(atom) -> dict:
    """组装官能团条目 dict（中心原子 + 重原子周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": _surr_idx(atom)}

def _is_cc_double(bond) -> bool:
    """判断键是否为 C=C 双键（排除芳香键）。"""
    if bond.GetBondType() != BondType.DOUBLE or bond.GetIsAromatic():
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == C and b.GetAtomicNum() == C

def _is_cc_triple(bond) -> bool:
    """判断键是否为 C≡C 三键。"""
    if bond.GetBondType() != BondType.TRIPLE:
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == C and b.GetAtomicNum() == C

def _bond_entry(bond) -> dict:
    """将 C=C/C≡C 键组装为条目 dict（两碳索引有序）。"""
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    return {"c1": min(a, b), "c2": max(a, b)}

def _filter_bond_entries(mol: Mol, pred, entry_fn) -> list[dict]:
    return [entry_fn(bond) for bond in mol.GetBonds() if pred(bond)]

def _ring_meta(mol: Mol) -> dict:
    """汇总环事实：环条目、环系与数量统计。"""
    from namepredict.layer1.ring_systems import build_ring_systems, sssr_rings
    rings = [{"atom_ids": r} for r in sssr_rings(mol)]
    systems = build_ring_systems(mol)
    return {
        "rings": rings, "n_rings": len(rings), "has_ring": bool(rings),
        "ring_systems": systems, "n_ring_systems": len(systems),
    }

def _radical_entry(atom) -> dict:
    """组装自由基位点条目 dict（哑原子所连重原子为中心，无周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": []}


_SUPPRESSIBLE = frozenset({"acid", "ester", "acyl_halide", "amide", "nitrile", "aldehyde"})  # 可被更高优先级 FG 整体压制的组合 FG：组合羰基 + 腈
_PRESENCE_SKIP = frozenset({"phosphate"})  # 磷酸不参与存在性判定（p41 与酯同为 9，纳入会改写压制结果）
_LEAF_DEMOTED = ("acid", "nitrile")  # 降级为"前缀叶"的组合 FG：整组碳排除出主链（P-61.1.3 carboxy/cyano）。其余组合 FG（酯/酰胺/醛/酰卤）降级为"氧代"——羰基碳留在链内，仅 O 作 oxo/formyl 前缀，由 L3 锚定叶识别。


def _arbitrate_parts(parts: dict) -> tuple[dict, frozenset[str]]:
    """P-41 仲裁：更高优先级 FG 使组合 FG 退出，叶型标 demoted。"""
    p41 = {sp.fg: sp.p41 for sp in FG_SPECS if sp.p41}
    present = {fg for fg in p41 if parts.get(fg) and fg not in _PRESENCE_SKIP}
    out = dict(parts)
    demoted: set[str] = set()
    for fg in _SUPPRESSIBLE:
        if out[fg] and any(p41[h] < p41[fg] for h in present if h != fg):
            if fg in _LEAF_DEMOTED:
                demoted |= {f"{fg}:{i}" for i in range(len(out[fg]))}
            else:
                out[fg] = []
    return out, frozenset(demoted)


# 局部环境 SMARTS 命中且无非局部判据的 FG 键；顺序沿用历史输出
_LOCAL_ENTRY_FGS = ("acid", "alcohol", "ester", "amide", "ketone", "amine", "thiol",
                    "nitrile", "acyl_halide")


def _local_entries(mol: Mol, hits: dict) -> dict:
    """按局部环境命中结果组装各 FG 的条目列表（中心原子升序）。"""
    return {fg: [_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get(fg, [])]
            for fg in _LOCAL_ENTRY_FGS}


def _detect_parts(mol: Mol) -> dict:
    """检测（未仲裁）分子中各类官能团条目。"""
    hits = match_local_fg(mol)
    acyl = [_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("acyl", [])]  # 须先于 aldehyde/radical
    heads = frozenset(e["center_idx"] for e in acyl)
    out = _local_entries(mol, hits)
    result = {**out,
        "radical": [_radical_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("radical", []) if t[0] not in heads],
        "acyl": acyl,
        "aldehyde": [e for e in (_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("aldehyde", [])) if e["center_idx"] not in heads],
        "phosphate": phosphate_entries(mol, hits.get("phosphate", []))}
    print(result)
    return result

def _collect_fgs(mol: Mol) -> dict:
    """聚合官能团条目并构建带类型清单（FG 唯一出口）。"""
    from namepredict.layer1.functional_group_inventory import build_inventory

    parts, demoted = _arbitrate_parts(_detect_parts(mol))
    return {"double_bonds": _filter_bond_entries(mol, _is_cc_double, _bond_entry),
            "triple_bonds": _filter_bond_entries(mol, _is_cc_triple, _bond_entry),
            "fg_inventory": build_inventory(parts, mol, demoted)}

def _info(mol: Mol, carbons: list[int]) -> dict:
    """组装分子分析结果 dict（碳信息 + 官能团 + 环事实）。"""
    base = {"mol": mol, "carbon_ids": carbons, "n_carbons": len(carbons)}
    return {**base, **_collect_fgs(mol), **_ring_meta(mol)}

def analyze(mol: Mol) -> dict:
    """分析分子并返回完整的官能团与结构信息 dict。"""
    result = _info(mol, [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == C])
    return result
