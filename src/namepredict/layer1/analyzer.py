"""L1 官能团分析器：枚举分子中各类官能团条目并汇总为分析结果 dict。
含氧酸中心（P/S 非碳骨架成键）的非局部判据与 oxo_kind 归一本模块。
羰基原语 _double_bonded_o_idxs/_alkoxy_c_of 供 L2 主基团表达式复用。
"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol
from collections import deque
from namepredict.constants import (
    C, H, N, O, S,
)
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer1.fg_local_smarts import match_local_fg

def _heavy(a) -> list:
    """非氢邻居索引。"""
    return [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]


def _double_bonded_o_idxs(carbon) -> list[int]:
    """返回碳上羰基双键氧的索引列表。"""
    mol = carbon.GetOwningMol()
    out: list[int] = []
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != O:
            continue
        b = mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
        if b is not None and b.GetBondType() == BondType.DOUBLE:
            out.append(n.GetIdx())
    return out


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    """返回氧上除羰基碳外的烷氧基碳索引。"""
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


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


_OXO_Z_ANCHORED = frozenset({"phosphate", "phosphonate", "sulfate"})  # 中心原子自任母体的 oxo_kind（碳骨架退为取代基）
_OXO_KIND_P = {(1, False): "phosphate", (1, True): "phosphonate"}  # P：(双键氧数, 有无直连碳臂) → kind
_OXO_KIND_S_ARM = {"halo": "sulfonyl_chloride", "n": "sulfonamide",   # S 带碳臂时：另一臂角色 → kind
                   "acid": "sulfonic", "o_arm": "sulfonate"}


def _oxo_roles(mol: Mol, z_idx: int) -> dict:
    """中心原子的邻居角色计数：双键氧/羟基氧/阴离子氧/O-臂/碳/卤素/氮/硫。"""
    roles = {"oxo": 0, "oh": 0, "om": 0, "o_arm": 0, "c": 0, "hal": 0, "n": 0, "s": 0}
    for n in mol.GetAtomWithIdx(z_idx).GetNeighbors():
        z, i = n.GetAtomicNum(), n.GetIdx()
        if z in (1, 0):
            continue
        if z == O:
            if mol.GetBondBetweenAtoms(z_idx, i).GetBondType() == BondType.DOUBLE:
                roles["oxo"] += 1
            elif n.GetFormalCharge() == -1:
                roles["om"] += 1
            elif n.GetTotalNumHs() >= 1:
                roles["oh"] += 1
            else:
                roles["o_arm"] += 1
        elif z == C:
            roles["c"] += 1
        elif z in (9, 17, 35, 53):
            roles["hal"] += 1
        elif z == N:
            roles["n"] += 1
        elif z == S:
            roles["s"] += 1
    return roles


def _oxo_kind(mol: Mol, z_idx: int) -> str | None:
    """按中心元素、双键氧数与臂角色归一 oxo_kind（P/S 同一张表驱动）。"""
    z = mol.GetAtomWithIdx(z_idx).GetAtomicNum()
    if z not in (15, 16):
        return None
    r = _oxo_roles(mol, z_idx)
    if z == 15:
        return _OXO_KIND_P.get((r["oxo"], r["c"] > 0))
    if r["oxo"] != 2:
        return None
    if not r["c"]:  # 无直连碳：中心自任母体，两臂须一酸式一 O-R（硫酸氢酯/硫酸酯）
        return "sulfate" if r["oh"] + r["om"] == 1 and r["o_arm"] == 1 else None
    arm = ("halo" if r["hal"] else "n" if r["n"] else "acid" if r["oh"] + r["om"] else "o_arm")
    return _OXO_KIND_S_ARM.get(arm)


def _oxo_arm_roots(mol: Mol, z_idx: int) -> tuple[list[int], list[int], list[int], list[int], list[int]]:
    """中心原子的 (直连碳臂, O-臂根碳, 羟基氧, 阴离子氧, 酸式氧) 索引。"""
    c_roots: list[int] = []
    o_roots: list[int] = []
    oh: list[int] = []
    om: list[int] = []
    acid: list[int] = []
    for n in mol.GetAtomWithIdx(z_idx).GetNeighbors():
        if n.GetAtomicNum() == 6:
            c_roots.append(n.GetIdx())
        elif n.GetAtomicNum() == O and mol.GetBondBetweenAtoms(z_idx, n.GetIdx()).GetBondType() != BondType.DOUBLE:
            if n.GetFormalCharge() == -1:
                om.append(n.GetIdx())
                acid.append(n.GetIdx())
            elif n.GetTotalNumHs() >= 1:
                oh.append(n.GetIdx())
                acid.append(n.GetIdx())
            else:
                o_roots.extend(j for j in _heavy(n) if j != z_idx)
    return c_roots, o_roots, oh, om, acid


def _arm_single_attach(mol: Mol, roots: list[int], core: set[int]) -> set[int] | None:
    """臂须单点回接且各臂互不相连：返回臂原子并集，违反返回 None。"""
    arm: set[int] = set()
    for root in roots:
        comp = _arm_component(mol, root, core)
        if comp & arm:
            return None
        arm |= comp
    return arm


def _oxoacid_entry(mol: Mol, core: tuple[int, ...]) -> dict | None:
    """校验含氧酸候选的非局部条件并归一 oxo_kind、锚点与负载。"""
    z_idx = core[0]
    kind = _oxo_kind(mol, z_idx)
    if kind is None:
        return None
    core_set = set(core)
    c_roots, o_roots, oh, om, _acid = _oxo_arm_roots(mol, z_idx)
    if kind in _OXO_Z_ANCHORED:  # 中心为母体：臂单点回接 + 整分子纯度（排除臂间成环、焦磷酸等）
        arm = _arm_single_attach(mol, c_roots + o_roots, core_set)
        if arm is None:
            return None
        all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}
        if all_heavy != (core_set | arm):
            return None
        return {"oxo_z": z_idx, "center_idx": z_idx, "oxo_kind": kind,
                "n_oh": len(oh), "n_om": len(om)}
    if len(c_roots) != 1:  # 碳锚定：直连碳臂唯一，该碳即母体碳骨架的入口
        return None
    return {"oxo_z": z_idx, "center_idx": c_roots[0], "oxo_kind": kind, "n_oh": len(oh), "n_om": len(om),
            "surr_idx": list(om)}  # surr_idx 只收阴离子氧：L2 据此补 anion 标志


# oxo_kind → P-41 类别键：磺酰胺属酰胺类 11，其余落 oxoacid 类 9
_OXO_CLASS_BY_KIND = {"sulfonamide": "sulfonamide"}


def oxoacid_entries(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> list[dict]:
    """分子中全部含氧酸中心（P/S 非碳骨架成键）的条目列表。"""
    if matches is None:
        matches = match_local_fg(mol).get("oxoacid", [])
    by_z: dict[int, dict] = {}
    for m in matches:
        e = _oxoacid_entry(mol, m)
        if e is not None:
            by_z[e["oxo_z"]] = e
    return [by_z[k] for k in sorted(by_z)]


def oxoacid_lists(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> dict[str, list[dict]]:
    """含氧酸条目按 P-41 类别键归位（同一检测器，类别由 oxo_kind 表决定）。"""
    out: dict[str, list[dict]] = {}
    for e in oxoacid_entries(mol, matches):
        out.setdefault(_OXO_CLASS_BY_KIND.get(e["oxo_kind"], "oxoacid"), []).append(e)
    return out


def phosphate_entries(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> list[dict]:
    """仅磷酸 kind（P(=O)(O⁄O⁻⁄OR)₃）的条目视图；膦酸与硫含氧酸不在内。"""
    return [e for e in oxoacid_entries(mol, matches) if e["oxo_kind"] == "phosphate"]
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
_PRESENCE_SKIP = frozenset({"oxoacid", "sulfonamide"})  # 含氧酸不参与存在性判定（纳入会改写压制结果）
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
        **oxoacid_lists(mol, hits.get("oxoacid", []))}
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
